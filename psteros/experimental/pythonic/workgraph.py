"""High-level pythonic scatter-gather workflow tying slab generation and relaxation."""

from __future__ import annotations

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, dynamic, namespace, task
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from test_modules.pythonic.aiat_ternary import (
    compute_surface_energies_scatter,
    create_mock_bulk_energy,
    create_mock_formation_enthalpy,
    create_mock_reference_energies,
)

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
    compute_thermodynamics: bool = False,
    bulk_energy: orm.Float | None = None,
    reference_energies: orm.Dict | None = None,
    formation_enthalpy: orm.Float | None = None,
    sampling: int = 100,
) -> t.Annotated[
    dict,
    namespace(
        generated_slabs=dynamic(orm.StructureData),
        relaxed_slabs=dynamic(orm.StructureData),
        slab_energies=dynamic(orm.Float),
        surface_energies=dynamic(orm.Dict),
    ),
]:
    """
    Full scatter-gather workflow for slab relaxation with optional thermodynamics.
    
    This implements the scatter-gather pattern:
    1. Generate slabs from bulk structure
    2. Scatter: relax each slab in parallel
    3. Gather: collect relaxed structures and energies
    4. (Optional) Compute ab initio atomistic thermodynamics
    
    Args:
        bulk_structure: Bulk structure
        miller_indices: Miller indices for surface
        min_slab_thickness: Minimum slab thickness (Å)
        min_vacuum_thickness: Minimum vacuum thickness (Å)
        code_label: AiiDA code label
        potential_family: Pseudopotential family
        potential_mapping: Element to potential mapping
        parameters: VASP/DFT parameters
        options: Computation resources options
        kpoints_spacing: K-points spacing (Å⁻¹)
        lll_reduce: Apply LLL reduction
        center_slab: Center slab in cell
        symmetrize: Generate symmetrically distinct terminations
        primitive: Use primitive cell
        in_unit_planes: Thickness in unit planes
        max_normal_search: Max normal search
        clean_workdir: Clean remote directories
        compute_thermodynamics: Whether to compute surface energies
        bulk_energy: Bulk total energy (required if compute_thermodynamics=True)
        reference_energies: Reference energies dict (required if compute_thermodynamics=True)
        formation_enthalpy: Formation enthalpy (required if compute_thermodynamics=True)
        sampling: Grid resolution for chemical potential
        
    Returns:
        Dictionary with:
        - generated_slabs: Original slab structures
        - relaxed_slabs: Relaxed structures
        - slab_energies: Total energies
        - surface_energies: Surface energy data (if compute_thermodynamics=True)
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

    # Prepare outputs
    outputs = {
        'generated_slabs': slab_namespace,
        'relaxed_slabs': relaxation_outputs.relaxed_structures,
        'slab_energies': relaxation_outputs.energies,
    }

    # Optional: Compute ab initio atomistic thermodynamics
    if compute_thermodynamics:
        if bulk_energy is None or reference_energies is None or formation_enthalpy is None:
            raise ValueError(
                'compute_thermodynamics=True requires bulk_energy, '
                'reference_energies, and formation_enthalpy'
            )
        
        surface_energy_data = compute_surface_energies_scatter(
            slabs=relaxation_outputs.relaxed_structures,
            energies=relaxation_outputs.energies,
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            reference_energies=reference_energies,
            formation_enthalpy=formation_enthalpy,
            sampling=sampling,
        ).surface_energies
        
        outputs['surface_energies'] = surface_energy_data

    return outputs


@task.calcfunction
def create_mock_thermo_data(value: orm.Float) -> orm.Dict:
    """Create mock thermodynamics data from a shifted value."""
    return orm.Dict(dict={
        'phi': float(value.value) * 1.5,
        'gamma_at_reference': float(value.value) * 1.5,
        'mock_thermodynamics': True,
    })


@task.graph
def compute_mock_thermodynamics(
    data: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
) -> t.Annotated[dict, namespace(thermo=dynamic(orm.Dict))]:
    """Scatter-gather: compute mock thermodynamics for each value."""
    thermo = {}
    for key, value in data.items():
        thermo[key] = create_mock_thermo_data(value=value).result
    return {'thermo': thermo}


@task.graph
def mock_scatter_gather(
    count: int = 3,
    delta: float = 0.5,
    with_thermodynamics: bool = False,
) -> t.Annotated[
    dict,
    namespace(
        source_values=dynamic(orm.Float),
        shifted_values=dynamic(orm.Float),
        thermo_results=dynamic(orm.Dict),
    ),
]:
    """
    Mock scatter-gather workflow for testing.
    
    Demonstrates the scatter-gather pattern with simple scalar operations:
    1. Generate scalar values
    2. Scatter: shift each value in parallel
    3. Gather: collect shifted values
    4. (Optional) Mock thermodynamics computation
    """
    # Generate inputs - direct access to .values namespace (no .result)
    data = generate_mock_scalars(count=orm.Int(count)).values
    
    # Scatter-gather: shift values in parallel
    shifted = scatter_shift_scalars(data=data, delta=delta).shifted
    
    # Prepare outputs
    outputs = {
        'source_values': data,
        'shifted_values': shifted,
    }
    
    # Optional mock thermodynamics
    if with_thermodynamics:
        thermo = compute_mock_thermodynamics(data=shifted).thermo
        outputs['thermo_results'] = thermo
    
    return outputs


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
    with_thermodynamics: bool = False,
) -> WorkGraph:
    """Return a lightweight ``WorkGraph`` using scatter-gather over scalars."""
    return mock_scatter_gather.build(
        count=count,
        delta=delta,
        with_thermodynamics=with_thermodynamics,
    )


