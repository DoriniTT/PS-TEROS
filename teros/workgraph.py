"""
PS-TEROS Core WorkGraph

This module contains the core workflow for PS-TEROS calculations.
One main WorkGraph with individual tasks running in parallel.
"""

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task
from ase.io import read
from typing import Annotated
from aiida_workgraph import spec


@task
def load_structure(filepath: str) -> orm.StructureData:
    """
    Load structure from a file using ASE.

    Args:
        filepath: Path to structure file

    Returns:
        StructureData node
    """
    atoms = read(filepath)
    return orm.StructureData(ase=atoms)


@task
def extract_total_energy(energies: orm.Dict) -> orm.Float:
    """
    Extract total energy from VASP energies output.

    Args:
        energies: Dictionary containing energy outputs from VASP

    Returns:
        Total energy as Float
    """
    energy_dict = energies.get_dict()

    # Try different keys in order of preference
    for key in ['energy_extrapolated', 'energy_no_entropy', 'energy']:
        if key in energy_dict:
            return orm.Float(energy_dict[key])

    raise ValueError(f"Could not find energy in energies dict. Available keys: {list(energy_dict.keys())}")


@task.graph(outputs=[
    'bulk_energy', 'metal_energy', 'nonmetal_energy', 'oxygen_energy',
    'bulk_structure', 'metal_structure', 'nonmetal_structure', 'oxygen_structure'
])
def formation_workgraph(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict,
    metal_potential_mapping: dict,
    nonmetal_potential_mapping: dict,
    oxygen_potential_mapping: dict,
    kpoints_spacing: float,
    bulk_parameters: dict,
    bulk_options: dict,
    metal_parameters: dict,
    metal_options: dict,
    nonmetal_parameters: dict,
    nonmetal_options: dict,
    oxygen_parameters: dict,
    oxygen_options: dict,
    clean_workdir: bool,
):
    """
    Core WorkGraph for formation enthalpy calculations of ternary oxides.

    This workflow relaxes the bulk compound and all reference elements in parallel.
    The workflow is general and works for any ternary oxide system (e.g., Ag3PO4, Fe2WO6, etc.).
    Each reference (metal, nonmetal, oxygen) can have its own specific calculation parameters.

    Args:
        structures_dir: Directory containing all structure files
        bulk_name: Filename of bulk structure (e.g., 'ag3po4.cif')
        metal_name: Filename of metal reference structure (e.g., 'Ag.cif')
        nonmetal_name: Filename of nonmetal reference structure (e.g., 'P.cif')
        oxygen_name: Filename of oxygen reference structure (e.g., 'O2.cif')
        code_label: Label of the VASP code in AiiDA
        potential_family: Name of the potential family
        bulk_potential_mapping: Mapping for bulk (e.g., {'Ag': 'Ag', 'P': 'P', 'O': 'O'})
        metal_potential_mapping: Mapping for metal (e.g., {'Ag': 'Ag'})
        nonmetal_potential_mapping: Mapping for nonmetal (e.g., {'P': 'P'})
        oxygen_potential_mapping: Mapping for oxygen (e.g., {'O': 'O'})
        kpoints_spacing: K-points spacing in A^-1 * 2pi
        bulk_parameters: VASP parameters for bulk
        bulk_options: Scheduler options for bulk
        metal_parameters: VASP parameters for metal reference
        metal_options: Scheduler options for metal reference
        nonmetal_parameters: VASP parameters for nonmetal reference
        nonmetal_options: Scheduler options for nonmetal reference
        oxygen_parameters: VASP parameters for oxygen reference
        oxygen_options: Scheduler options for oxygen reference
        clean_workdir: Whether to clean work directory

    Returns:
        Dictionary with energies and structures for all systems
    """
    # Load the code
    code = orm.load_code(code_label)

    # Get VASP workchain and wrap it as a task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # ===== BULK RELAXATION =====
    bulk_struct = load_structure(filepath=f"{structures_dir}/{bulk_name}")

    bulk_vasp = VaspTask(
        structure=bulk_struct.result,
        code=code,
        parameters={'incar': bulk_parameters},
        options=bulk_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=bulk_potential_mapping,
        clean_workdir=clean_workdir,
    )

    bulk_energy = extract_total_energy(energies=bulk_vasp.energies)

    # ===== METAL RELAXATION =====
    metal_struct = load_structure(filepath=f"{structures_dir}/{metal_name}")

    metal_vasp = VaspTask(
        structure=metal_struct.result,
        code=code,
        parameters={'incar': metal_parameters},
        options=metal_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=metal_potential_mapping,
        clean_workdir=clean_workdir,
    )

    metal_energy = extract_total_energy(energies=metal_vasp.energies)

    # ===== NONMETAL RELAXATION =====
    nonmetal_struct = load_structure(filepath=f"{structures_dir}/{nonmetal_name}")

    nonmetal_vasp = VaspTask(
        structure=nonmetal_struct.result,
        code=code,
        parameters={'incar': nonmetal_parameters},
        options=nonmetal_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=nonmetal_potential_mapping,
        clean_workdir=clean_workdir,
    )

    nonmetal_energy = extract_total_energy(energies=nonmetal_vasp.energies)

    # ===== OXYGEN RELAXATION =====
    oxygen_struct = load_structure(filepath=f"{structures_dir}/{oxygen_name}")

    oxygen_vasp = VaspTask(
        structure=oxygen_struct.result,
        code=code,
        parameters={'incar': oxygen_parameters},
        options=oxygen_options,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=oxygen_potential_mapping,
        clean_workdir=clean_workdir,
    )

    oxygen_energy = extract_total_energy(energies=oxygen_vasp.energies)

    # Return all outputs
    return {
        'bulk_energy': bulk_energy.result,
        'bulk_structure': bulk_vasp.structure,
        'metal_energy': metal_energy.result,
        'metal_structure': metal_vasp.structure,
        'nonmetal_energy': nonmetal_energy.result,
        'nonmetal_structure': nonmetal_vasp.structure,
        'oxygen_energy': oxygen_energy.result,
        'oxygen_structure': oxygen_vasp.structure,
    }


def build_formation_workgraph(
    structures_dir: str,
    bulk_name: str,
    metal_name: str,
    nonmetal_name: str,
    oxygen_name: str,
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
    name: str = 'FormationEnthalpy',
):
    """
    Build a formation enthalpy WorkGraph for any ternary oxide.

    This is a convenience wrapper that builds and returns a WorkGraph
    ready to calculate formation enthalpy.

    Returns:
        WorkGraph instance ready to be submitted
    """
    # Build the workgraph
    # Note: parameters will be wrapped with {'incar': ...} inside the graph
    wg = formation_workgraph.build(
        structures_dir=structures_dir,
        bulk_name=bulk_name,
        metal_name=metal_name,
        nonmetal_name=nonmetal_name,
        oxygen_name=oxygen_name,
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping=bulk_potential_mapping or {},
        metal_potential_mapping=metal_potential_mapping or {},
        nonmetal_potential_mapping=nonmetal_potential_mapping or {},
        oxygen_potential_mapping=oxygen_potential_mapping or {},
        kpoints_spacing=kpoints_spacing,
        bulk_parameters=bulk_parameters or {},
        bulk_options=bulk_options or {},
        metal_parameters=metal_parameters or {},
        metal_options=metal_options or {},
        nonmetal_parameters=nonmetal_parameters or {},
        nonmetal_options=nonmetal_options or {},
        oxygen_parameters=oxygen_parameters or {},
        oxygen_options=oxygen_options or {},
        clean_workdir=clean_workdir,
    )

    # Set the name
    wg.name = name

    return wg


if __name__ == '__main__':
    """
    Simple test/example of the WorkGraph.
    For full examples, see examples/formation/
    """
    print("PS-TEROS Core WorkGraph Module")
    print("=" * 50)
    print("\nThis creates a WorkGraph using @task.graph decorator:")
    print("  - load_structure tasks (using @task)")
    print("  - VASP relaxation tasks (run in parallel)")
    print("  - extract_total_energy tasks (using @task)")
    print("\nFor examples, see:")
    print("  - examples/formation/formation.py (formation enthalpy)")
