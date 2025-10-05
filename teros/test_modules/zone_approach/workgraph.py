"""Workflow wiring for the slab relaxation Map-zone example."""

from __future__ import annotations

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import Map, WorkGraph, dynamic, namespace, task
from aiida_workgraph.engine import task_manager as _task_manager
from aiida_workgraph.utils import get_nested_dict

from collections.abc import Mapping
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

MappingStrAny = t.Mapping[str, t.Any]
MappingStrStr = t.Mapping[str, str]

# Temporary shim: aiida-workgraph 1.0.0b3 drops the ``source`` kwarg when the
# map zone iterates over a Task output namespace. Patch ``TaskManager`` so the
# ``source`` payload is recovered directly from the socket when missing.
_ORIG_EXECUTE_MAP_TASK = _task_manager.TaskManager.execute_map_task


def _patched_execute_map_task(self, task, kwargs):
    kwargs = dict(kwargs)
    source = kwargs.get('source')
    if not source:
        source = self.get_socket_value(task.inputs.source)
    if not source:
        links = getattr(task.inputs.source, '_links', [])
        if links:
            link = links[0]
            results = self.ctx._task_results.get(link.from_node.name)
            if results is not None:
                scoped_name = link.from_socket._scoped_name
                source = get_nested_dict(results, scoped_name, default=None)
                if source is None:
                    source = get_nested_dict(results, f'result.{scoped_name}', default=None)
                if source is None:
                    source = get_nested_dict(results, f'_outputs.{scoped_name}', default=None)
                if source is None and isinstance(results, Mapping):
                    result_namespace = results.get('result')
                    if isinstance(result_namespace, Mapping):
                        source = result_namespace.get(scoped_name)
    if isinstance(source, Mapping):
        source = dict(source)
    kwargs['source'] = source or {}
    return _ORIG_EXECUTE_MAP_TASK(self, task, kwargs)


if not getattr(_task_manager.TaskManager.execute_map_task, '_zone_patch', False):
    _task_manager.TaskManager.execute_map_task = _patched_execute_map_task
    _task_manager.TaskManager.execute_map_task._zone_patch = True


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
def relax_single_slab(
    structure: orm.StructureData,
    *,
    code_label: str,
    potential_family: str,
    potential_mapping: MappingStrStr,
    parameters: MappingStrAny,
    options: MappingStrAny,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[
    dict[str, object],
    namespace(relaxed_structure=orm.StructureData, energy=orm.Float),
]:
    """Relax a single slab and expose the relaxed structure and total energy."""
    code = orm.load_code(code_label)
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    relax_task = task(vasp_wc)

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
    energy = extract_total_energy(energies=relaxation.misc).result

    return {
        'relaxed_structure': relaxation.structure,
        'energy': energy,
    }


@task
def generate_mock_scalars(count: int = 3) -> t.Annotated[
    dict[str, orm.Float],
    namespace(values=dynamic(orm.Float)),
]:
    """Return a dynamic namespace of scalar values for Map-zone testing."""
    values = {
        f'value_{index:02d}': orm.Float(float(index))
        for index in range(count)
    }
    return {'values': values}


@task
def shift_scalar(value: orm.Float, delta: float = 1.0) -> orm.Float:
    """Apply a constant offset to a scalar value."""
    return orm.Float(value.value + float(delta))


def build_zone_workgraph(
    *,
    name: str = 'slab_relax_map_zone',
    bulk_structure: orm.StructureData,
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
) -> WorkGraph:
    """Return a ``WorkGraph`` that maps slab relaxations via the Map zone."""
    if isinstance(miller_indices, tuple):
        miller_indices = list(miller_indices)

    outputs_spec = namespace(
        generated_slabs=dynamic(orm.StructureData),
        relaxed_slabs=dynamic(orm.StructureData),
        slab_energies=dynamic(orm.Float),
    )

    with WorkGraph(name=name, outputs=outputs_spec) as wg:
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

        with Map(slab_namespace) as map_zone:
            relaxation = relax_single_slab(
                structure=map_zone.item.value,
                code_label=code_label,
                potential_family=potential_family,
                potential_mapping=potential_mapping,
                parameters=parameters,
                options=options,
                kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
            )

            map_zone.gather(
                {
                    'relaxed_slabs': relaxation.relaxed_structure,
                    'slab_energies': relaxation.energy,
                }
            )

        wg.outputs.generated_slabs = slab_namespace
        wg.outputs.relaxed_slabs = map_zone.outputs.relaxed_slabs
        wg.outputs.slab_energies = map_zone.outputs.slab_energies

    return wg


def build_mock_zone_workgraph(
    *,
    name: str = 'mock_map_zone',
    count: int = 3,
    delta: float = 0.5,
) -> WorkGraph:
    """Return a lightweight ``WorkGraph`` that maps over generated scalars."""
    outputs_spec = namespace(
        source_values=dynamic(orm.Float),
        shifted_values=dynamic(orm.Float),
    )

    with WorkGraph(name=name, outputs=outputs_spec) as wg:
        source_namespace = generate_mock_scalars(count=count).values

        with Map(source_namespace) as map_zone:
            shifted = shift_scalar(value=map_zone.item.value, delta=delta)
            map_zone.gather({'shifted_values': shifted})

        wg.outputs.source_values = source_namespace
        wg.outputs.shifted_values = map_zone.outputs.shifted_values

    return wg
