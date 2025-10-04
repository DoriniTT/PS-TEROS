"""High-level pythonic scatter-gather workflow tying slab generation and relaxation."""

from __future__ import annotations

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, dynamic, namespace, task
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

MappingStrAny = t.Mapping[str, t.Any]
MappingStrStr = t.Mapping[str, str]


@task
def generate_slab_structures(
    bulk_structure: orm.StructureData | object,
    miller_indices: list[int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int | None = None,
) -> t.Annotated[dict[str, orm.StructureData], namespace(slabs=dynamic(orm.StructureData))]:
    """Generate slab terminations from the provided structure."""
    adaptor = AseAtomsAdaptor()

    if isinstance(bulk_structure, orm.StructureData):
        ase_structure = bulk_structure.get_ase()
    else:
        ase_structure = bulk_structure

    pymatgen_structure = adaptor.get_structure(ase_structure)

    if primitive:
        analyzer = SpacegroupAnalyzer(pymatgen_structure)
        pymatgen_structure = analyzer.get_primitive_standard_structure()

    generator = SlabGenerator(
        pymatgen_structure,
        tuple(miller_indices),
        min_slab_thickness,
        min_vacuum_thickness,
        center_slab=center_slab,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
        lll_reduce=lll_reduce,
    )

    slab_nodes: dict[str, orm.StructureData] = {}
    for index, slab in enumerate(generator.get_slabs(symmetrize=symmetrize)):
        orthogonal_slab = slab.get_orthogonal_c_slab()
        ase_slab = adaptor.get_atoms(orthogonal_slab)
        slab_nodes[f"slab_{index:02d}"] = orm.StructureData(ase=ase_slab)

    return {'slabs': slab_nodes}


@task
def extract_total_energy(energies: orm.Dict | dict) -> orm.Float:
    """Return the total energy stored in a VASP ``misc`` dictionary."""
    energy_dict = energies.get_dict() if isinstance(energies, orm.Dict) else dict(energies)
    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in VASP outputs. Available keys: {available}')


@task.graph
def relax_slabs_pythonic(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code_label: str,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[
    dict[str, object],
    namespace(
        relaxed_structures=dynamic(orm.StructureData),
        energies=dynamic(orm.Float),
    ),
]:
    """Relax each slab structure in parallel using the pythonic scatter pattern."""
    code = orm.load_code(code_label)
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task = task(vasp_wc)

    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}

    for label, structure in slabs.items():
        inputs: dict[str, t.Any] = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options),
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
        }
        if kpoints_spacing is not None:
            inputs['kpoints_spacing'] = kpoints_spacing

        relaxation = relax_task(**inputs)
        relaxed[label] = relaxation.structure
        energies_ns[label] = extract_total_energy(energies=relaxation.misc).result

    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
    }


@task.graph
def slab_relaxation_pythonic(
    bulk_structure: orm.StructureData,
    *,
    miller_indices: tuple[int, int, int] | list[int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    code_label: str,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[
    dict[str, object],
    namespace(
        generated_slabs=dynamic(orm.StructureData),
        relaxed_slabs=dynamic(orm.StructureData),
        slab_energies=dynamic(orm.Float),
    ),
]:
    """Build the pythonic scatter-gather slab workflow."""
    if isinstance(miller_indices, tuple):
        miller_indices = list(miller_indices)

    slab_namespace = generate_slab_structures(
        bulk_structure=bulk_structure,
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        symmetrize=symmetrize,
        primitive=primitive,
        in_unit_planes=in_unit_planes,
        max_normal_search=max_normal_search,
    ).slabs

    relaxation_outputs = relax_slabs_pythonic(
        slabs=slab_namespace,
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
    )

    return {
        'generated_slabs': slab_namespace,
        'relaxed_slabs': relaxation_outputs.relaxed_structures,
        'slab_energies': relaxation_outputs.energies,
    }


def build_pythonic_workgraph(**kwargs) -> WorkGraph:
    """Convenience wrapper returning a built workgraph ready for ``run``."""
    if 'miller_indices' in kwargs and isinstance(kwargs['miller_indices'], tuple):
        kwargs = dict(kwargs)
        kwargs['miller_indices'] = list(kwargs['miller_indices'])
    return slab_relaxation_pythonic.build(**kwargs)
