"""
PS-TEROS Core WorkGraph using Map zone for parallel slab relaxations.

This version uses the context manager paradigm (with WorkGraph() as wg:)
similar to the legacy active_map_zone pattern.
"""

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, Zone, task, spec
from teros.modules.hf import calculate_formation_enthalpy
from teros.modules.slabs import get_slabs
from typing import Any


@task
def load_structure(filepath: str):
    """Load structure from file using ASE."""
    from ase.io import read
    atoms = read(filepath)
    return orm.StructureData(ase=atoms)


@task
def extract_total_energy(energies):
    """Extract total energy from VASP energies output."""
    if isinstance(energies, orm.Dict):
        energy_dict = energies.get_dict()
    else:
        energy_dict = energies

    if 'total_energies' in energy_dict:
        energy_dict = energy_dict['total_energies']

    for key in ['energy_extrapolated', 'energy_no_entropy', 'energy']:
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    raise ValueError(f"Could not find energy. Available keys: {list(energy_dict.keys())}")


@task.graph(outputs=['relaxed_slabs', 'slab_energies'])
def relax_all_slabs(slabs, code, slab_parameters, slab_options, slab_kpoints_spacing,
                    potential_family, slab_potential_mapping, clean_workdir):
    """Graph task to relax all slabs in parallel using Zone."""
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Two-pass approach for parallel slab relaxations
    relaxed_slabs = {}
    slab_energies_dict = {}
    vasp_tasks = {}

    # Extract slabs from dynamic namespace
    if hasattr(slabs.slabs, "_get_keys"):
        slab_items = [(key, getattr(slabs.slabs, key)) for key in slabs.slabs._get_keys()]
    else:
        slab_items = list(slabs.slabs.items()) if isinstance(slabs.slabs, dict) else []

    # PASS 1: Create all VASP tasks (independent, run in parallel)
    for slab_id, slab_structure in slab_items:
        vasp_tasks[slab_id] = VaspTask(
            structure=slab_structure,
            code=code,
            parameters={'incar': slab_parameters},
            options=slab_options,
            kpoints_spacing=slab_kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=slab_potential_mapping,
            clean_workdir=clean_workdir,
        )

    # PASS 2: Extract energies and collect results
    for slab_id in vasp_tasks:
        vasp_task = vasp_tasks[slab_id]
        energy = extract_total_energy(energies=vasp_task.misc)

        relaxed_slabs[slab_id] = vasp_task.structure
        slab_energies_dict[slab_id] = energy.result

    return {
        'relaxed_slabs': relaxed_slabs,
        'slab_energies': slab_energies_dict,
    }


def build_core_workgraph_with_map(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    potential_family: str = 'PBE',
    bulk_potential_mapping: dict = None,
    metal_potential_mapping: dict = None,
    nonmetal_potential_mapping: dict = None,
    oxygen_potential_mapping: dict = None,
    kpoints_spacing: float = 0.3,
    bulk_parameters: dict = None,
    bulk_options: dict = None,
    metal_parameters: dict = None,
    metal_options: dict = None,
    nonmetal_parameters: dict = None,
    nonmetal_options: dict = None,
    oxygen_parameters: dict = None,
    oxygen_options: dict = None,
    clean_workdir: bool = True,
    slab_parameters: dict = None,
    slab_options: dict = None,
    slab_potential_mapping: dict = None,
    slab_kpoints_spacing: float = None,
    lll_reduce: bool = False,
    center_slab: bool = True,
    symmetrize: bool = False,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
    relax_slabs: bool = False,
    name: str = 'FormationEnthalpy',
):
    """
    Build WorkGraph using context manager paradigm with Map zone for slab relaxations.

    This follows the pattern from legacy/v2/teros/core/workgraph.py which uses
    with WorkGraph() as wg: and with active_map_zone() (now Map).
    """
    # Get VASP workchain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')

    # Load code
    code = orm.load_code(code_label)

    # Build WorkGraph using context manager
    with WorkGraph(name=name) as wg:
        # ===== BULK RELAXATION =====
        bulk_struct_task = wg.add_task(
            load_structure,
            name="load_bulk_structure",
            filepath=f"{structures_dir}/{bulk_name}"
        )

        bulk_vasp_task = wg.add_task(
            VaspWorkChain,
            name="bulk_relaxation",
            structure=bulk_struct_task.outputs.result,
            code=code,
            parameters={'incar': bulk_parameters or {}},
            options=bulk_options or {},
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=bulk_potential_mapping or {},
            clean_workdir=clean_workdir,
        )

        bulk_energy_task = wg.add_task(
            extract_total_energy,
            name="extract_bulk_energy",
            energies=bulk_vasp_task.outputs.misc
        )

        # ===== REFERENCE RELAXATIONS (in parallel) =====
        # Metal
        metal_struct_task = wg.add_task(
            load_structure,
            name="load_metal_structure",
            filepath=f"{structures_dir}/{metal_name}"
        )

        metal_vasp_task = wg.add_task(
            VaspWorkChain,
            name="metal_relaxation",
            structure=metal_struct_task.outputs.result,
            code=code,
            parameters={'incar': metal_parameters or {}},
            options=metal_options or {},
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=metal_potential_mapping or {},
            clean_workdir=clean_workdir,
        )

        metal_energy_task = wg.add_task(
            extract_total_energy,
            name="extract_metal_energy",
            energies=metal_vasp_task.outputs.misc
        )

        # Nonmetal
        nonmetal_struct_task = wg.add_task(
            load_structure,
            name="load_nonmetal_structure",
            filepath=f"{structures_dir}/{nonmetal_name}"
        )

        nonmetal_vasp_task = wg.add_task(
            VaspWorkChain,
            name="nonmetal_relaxation",
            structure=nonmetal_struct_task.outputs.result,
            code=code,
            parameters={'incar': nonmetal_parameters or {}},
            options=nonmetal_options or {},
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=nonmetal_potential_mapping or {},
            clean_workdir=clean_workdir,
        )

        nonmetal_energy_task = wg.add_task(
            extract_total_energy,
            name="extract_nonmetal_energy",
            energies=nonmetal_vasp_task.outputs.misc
        )

        # Oxygen
        oxygen_struct_task = wg.add_task(
            load_structure,
            name="load_oxygen_structure",
            filepath=f"{structures_dir}/{oxygen_name}"
        )

        oxygen_vasp_task = wg.add_task(
            VaspWorkChain,
            name="oxygen_relaxation",
            structure=oxygen_struct_task.outputs.result,
            code=code,
            parameters={'incar': oxygen_parameters or {}},
            options=oxygen_options or {},
            kpoints_spacing=kpoints_spacing,
            potential_family=potential_family,
            potential_mapping=oxygen_potential_mapping or {},
            clean_workdir=clean_workdir,
        )

        oxygen_energy_task = wg.add_task(
            extract_total_energy,
            name="extract_oxygen_energy",
            energies=oxygen_vasp_task.outputs.misc
        )

        # ===== FORMATION ENTHALPY =====
        formation_task = wg.add_task(
            calculate_formation_enthalpy,
            name="formation_enthalpy",
            bulk_structure=bulk_vasp_task.outputs.structure,
            bulk_energy=bulk_energy_task.outputs.result,
            metal_structure=metal_vasp_task.outputs.structure,
            metal_energy=metal_energy_task.outputs.result,
            nonmetal_structure=nonmetal_vasp_task.outputs.structure,
            nonmetal_energy=nonmetal_energy_task.outputs.result,
            oxygen_structure=oxygen_vasp_task.outputs.structure,
            oxygen_energy=oxygen_energy_task.outputs.result,
        )

        # ===== SLAB GENERATION =====
        slabs_task = wg.add_task(
            get_slabs,
            name="generate_slabs",
            relaxed_structure=bulk_vasp_task.outputs.structure,
            miller_indices=miller_indices,
            min_slab_thickness=min_slab_thickness,
            min_vacuum_thickness=min_vacuum_thickness,
            lll_reduce=lll_reduce,
            center_slab=center_slab,
            symmetrize=symmetrize,
            primitive=primitive,
            in_unit_planes=in_unit_planes,
            max_normal_search=max_normal_search,
        )

        # ===== SLAB RELAXATIONS =====
        if relax_slabs:
            # Use slab-specific parameters or fall back to bulk
            slab_params = slab_parameters if slab_parameters is not None else bulk_parameters
            slab_opts = slab_options if slab_options is not None else bulk_options
            slab_pot_map = slab_potential_mapping if slab_potential_mapping is not None else bulk_potential_mapping
            slab_kpts = slab_kpoints_spacing if slab_kpoints_spacing is not None else kpoints_spacing

            # Use graph task for parallel slab relaxations
            slab_relax_task = wg.add_task(
                relax_all_slabs,
                name="relax_all_slabs",
                slabs=slabs_task.outputs,
                code=code,
                slab_parameters=slab_params or {},
                slab_options=slab_opts or {},
                slab_kpoints_spacing=slab_kpts,
                potential_family=potential_family,
                slab_potential_mapping=slab_pot_map or {},
                clean_workdir=clean_workdir,
            )

            # Set outputs
            wg.outputs.relaxed_slabs = slab_relax_task.outputs.relaxed_slabs
            wg.outputs.slab_energies = slab_relax_task.outputs.slab_energies

        # Set workflow outputs
        wg.outputs.bulk_energy = bulk_energy_task.outputs.result
        wg.outputs.bulk_structure = bulk_vasp_task.outputs.structure
        wg.outputs.metal_energy = metal_energy_task.outputs.result
        wg.outputs.metal_structure = metal_vasp_task.outputs.structure
        wg.outputs.nonmetal_energy = nonmetal_energy_task.outputs.result
        wg.outputs.nonmetal_structure = nonmetal_vasp_task.outputs.structure
        wg.outputs.oxygen_energy = oxygen_energy_task.outputs.result
        wg.outputs.oxygen_structure = oxygen_vasp_task.outputs.structure
        wg.outputs.formation_enthalpy = formation_task.outputs.result
        wg.outputs.slab_structures = slabs_task.outputs.slabs

    return wg
