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


@task.calcfunction
def generate_slab_structures(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    min_slab_thickness: orm.Float,
    min_vacuum_thickness: orm.Float,
    lll_reduce: orm.Bool,
    center_slab: orm.Bool,
    symmetrize: orm.Bool,
    primitive: orm.Bool,
) -> t.Annotated[dict, namespace(slabs=dynamic(orm.StructureData))]:
    """Generate slab terminations from the provided structure."""
    adaptor = AseAtomsAdaptor()

    ase_structure = bulk_structure.get_ase()
    pymatgen_structure = adaptor.get_structure(ase_structure)

    if primitive.value:
        analyzer = SpacegroupAnalyzer(pymatgen_structure)
        pymatgen_structure = analyzer.get_primitive_standard_structure()

    generator = SlabGenerator(
        pymatgen_structure,
        tuple(miller_indices.get_list()),
        min_slab_thickness.value,
        min_vacuum_thickness.value,
        center_slab=center_slab.value,
        in_unit_planes=False,
        max_normal_search=None,
        lll_reduce=lll_reduce.value,
    )

    slab_nodes: dict[str, orm.StructureData] = {}
    for index, slab in enumerate(generator.get_slabs(symmetrize=symmetrize.value)):
        orthogonal_slab = slab.get_orthogonal_c_slab()
        ase_slab = adaptor.get_atoms(orthogonal_slab)
        slab_nodes[f"slab_{index:02d}"] = orm.StructureData(ase=ase_slab)

    return {'slabs': slab_nodes}


@task.calcfunction
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """Return the total energy stored in a VASP ``misc`` dictionary."""
    energy_dict = energies.get_dict()
    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ('energy_extrapolated', 'energy_no_entropy', 'energy'):
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    available = ', '.join(sorted(energy_dict.keys()))
    raise ValueError(f'Unable to find total energy in VASP outputs. Available keys: {available}')


@task.calcfunction
def generate_mock_scalars(count: orm.Int) -> t.Annotated[
    dict,
    namespace(values=dynamic(orm.Float)),
]:
    """Return a dynamic namespace of scalar values for scatter-gather testing."""
    values = {
        f'value_{index:02d}': orm.Float(float(index))
        for index in range(count.value)
    }
    return {'values': values}


@task.calcfunction
def shift_scalar(value: orm.Float, delta: orm.Float) -> orm.Float:
    """Apply a constant offset to a scalar value."""
    return orm.Float(value.value + delta.value)


@task.graph
def relax_slabs_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(relaxed_structures=dynamic(orm.StructureData), energies=dynamic(orm.Float))]:
    """
    Scatter phase: relax each slab structure in parallel.
    
    This task graph iterates over the input slabs dictionary and creates
    independent relaxation tasks that run in parallel. The loop wrapping
    is necessary because slabs is a future output that isn't available
    at graph construction time.
    """
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task_cls = task(vasp_wc)
    
    relaxed: dict[str, orm.StructureData] = {}
    energies_ns: dict[str, orm.Float] = {}

    # Scatter: create independent tasks for each slab (runs in parallel)
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
        
        relaxation = relax_task_cls(**inputs)
        relaxed[label] = relaxation.structure
        energies_ns[label] = extract_total_energy(energies=relaxation.misc).result

    # Gather: return collected results
    return {
        'relaxed_structures': relaxed,
        'energies': energies_ns,
    }


@task.graph
def scatter_shift_scalars(
    data: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    delta: float = 0.5,
) -> t.Annotated[dict, namespace(shifted=dynamic(orm.Float))]:
    """Scatter phase: applies shift_scalar to each value in parallel."""
    shifted = {}
    delta_node = orm.Float(delta)
    for key, value in data.items():
        shifted[key] = shift_scalar(value=value, delta=delta_node).result
    return {'shifted': shifted}


@task.graph
def slab_relaxation_scatter_gather(
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
    dict,
    namespace(
        generated_slabs=dynamic(orm.StructureData),
        relaxed_slabs=dynamic(orm.StructureData),
        slab_energies=dynamic(orm.Float),
    ),
]:
    """
    Full scatter-gather workflow for slab relaxation.
    
    This implements the scatter-gather pattern:
    1. Generate slabs from bulk structure
    2. Scatter: relax each slab in parallel
    3. Gather: collect relaxed structures and energies
    """
    if isinstance(miller_indices, tuple):
        miller_indices = list(miller_indices)

    # Load the code outside the nested graph
    code = orm.load_code(code_label)

    # Generate input slabs - convert parameters to AiiDA types
    slab_namespace = generate_slab_structures(
        bulk_structure=bulk_structure,
        miller_indices=orm.List(list=miller_indices),
        min_slab_thickness=orm.Float(min_slab_thickness),
        min_vacuum_thickness=orm.Float(min_vacuum_thickness),
        lll_reduce=orm.Bool(lll_reduce),
        center_slab=orm.Bool(center_slab),
        symmetrize=orm.Bool(symmetrize),
        primitive=orm.Bool(primitive),
    ).slabs

    # Scatter-gather: relax slabs in parallel and collect results
    relaxation_outputs = relax_slabs_scatter(
        slabs=slab_namespace,
        code=code,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        kpoints_spacing=kpoints_spacing,
        clean_workdir=clean_workdir,
    )

    # Return all outputs
    return {
        'generated_slabs': slab_namespace,
        'relaxed_slabs': relaxation_outputs.relaxed_structures,
        'slab_energies': relaxation_outputs.energies,
    }


@task.graph
def mock_scatter_gather(
    count: int = 3,
    delta: float = 0.5,
) -> t.Annotated[
    dict,
    namespace(
        source_values=dynamic(orm.Float),
        shifted_values=dynamic(orm.Float),
    ),
]:
    """
    Mock scatter-gather workflow for testing.
    
    Demonstrates the scatter-gather pattern with simple scalar operations:
    1. Generate scalar values
    2. Scatter: shift each value in parallel
    3. Gather: collect shifted values
    """
    # Generate inputs - direct access to .values namespace (no .result)
    data = generate_mock_scalars(count=orm.Int(count)).values
    
    # Scatter-gather: shift values in parallel
    shifted = scatter_shift_scalars(data=data, delta=delta).shifted
    
    # Return outputs
    return {
        'source_values': data,
        'shifted_values': shifted,
    }


def build_pythonic_workgraph(**kwargs) -> WorkGraph:
    """Convenience wrapper returning a built workgraph ready for ``run``."""
    if 'miller_indices' in kwargs and isinstance(kwargs['miller_indices'], tuple):
        kwargs = dict(kwargs)
        kwargs['miller_indices'] = list(kwargs['miller_indices'])
    return slab_relaxation_scatter_gather.build(**kwargs)


def build_mock_scatter_gather_workgraph(
    *,
    count: int = 3,
    delta: float = 0.5,
) -> WorkGraph:
    """Return a lightweight ``WorkGraph`` using scatter-gather over scalars."""
    return mock_scatter_gather.build(count=count, delta=delta)

