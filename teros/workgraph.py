"""
PS-TEROS Core WorkGraph

This module contains the core workflow for PS-TEROS calculations.
The core_workgraph function is the main entry point that orchestrates
all calculations using helper functions from the modules package.
"""

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task
from teros.modules.helper_functions import (
    get_structure_from_file,
    prepare_vasp_inputs,
)


@task.graph(outputs=["relaxed_structure", "misc", "energies"])
def core_workgraph(
    structure_path: str = None,
    structure: orm.StructureData = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    kpoints_spacing: float = 0.3,
    potential_family: str = 'PBE.54',
    potential_mapping: dict = None,
    parameters: dict = None,
    options: dict = None,
    clean_workdir: bool = True,
):
    """
    Core WorkGraph for PS-TEROS calculations.

    This is the main workflow function that orchestrates all calculations.
    It is designed to be simple, clean, and focused on the workflow logic.

    The workflow:
    1. Load structure from file (if needed)
    2. Prepare VASP inputs using helper functions
    3. Run VASP relaxation using the vasp.vasp workchain
    4. Return relaxed structure and outputs

    Args:
        structure_path: Path to structure file (CIF format)
        structure: Pre-loaded structure (alternative to structure_path)
        code_label: Label of the VASP code in AiiDA (default: VASP-VTST-6.4.3@bohr)
        kpoints_spacing: K-points spacing in A^-1 * 2pi (default: 0.3)
        potential_family: Name of the potential family (default: PBE.54)
        potential_mapping: Mapping of elements to potentials (auto-generated if None)
        parameters: VASP input parameters (uses defaults if None)
        options: Scheduler options (uses defaults if None)
        clean_workdir: Whether to clean work directory after completion (default: True)

    Returns:
        Dictionary with:
            - relaxed_structure: The relaxed structure
            - misc: Miscellaneous outputs from VASP
            - energies: Energy outputs from VASP
    """
    from aiida_workgraph import task as task_decorator

    # Step 1: Load structure from file if needed
    if structure is None and structure_path is not None:
        struct_output = get_structure_from_file(filepath=orm.Str(structure_path))
        struct_input = struct_output.result
    elif structure is not None:
        struct_input = structure
    else:
        raise ValueError("Either structure_path or structure must be provided")

    # Step 2: Prepare all VASP inputs using helper function
    # This is plain Python that runs during graph construction
    inputs = prepare_vasp_inputs(
        code_label=code_label,
        structure=structure,
        structure_path=structure_path,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        clean_workdir=clean_workdir,
    )

    # Step 3: Run VASP relaxation
    # Load and wrap the VASP workchain as a task
    VaspWorkChain = WorkflowFactory('vasp.vasp')
    VaspTask = task_decorator(VaspWorkChain)

    # Execute the VASP relaxation task
    vasp_node = VaspTask(
        structure=struct_input,
        code=inputs['code'],
        parameters=inputs['params_dict'],
        options=inputs['options_dict'],
        kpoints_spacing=inputs['kpoints_spacing'],
        potential_family=inputs['potential_family'],
        potential_mapping=inputs['pot_mapping_dict'],
        clean_workdir=inputs['clean_workdir'],
    )

    # Step 4: Return the outputs we want to expose
    return {
        'relaxed_structure': vasp_node.structure,
        'misc': vasp_node.misc,
        'energies': vasp_node.energies,
    }


# Convenience function to build the core workgraph
def build_workgraph(
    structure_path: str = None,
    structure: orm.StructureData = None,
    code_label: str = 'VASP-VTST-6.4.3@bohr',
    kpoints_spacing: float = 0.3,
    potential_family: str = 'PBE.54',
    potential_mapping: dict = None,
    parameters: dict = None,
    options: dict = None,
    clean_workdir: bool = True,
    name: str = 'PS-TEROS',
):
    """
    Build a PS-TEROS WorkGraph instance.

    This is a convenience wrapper around core_workgraph that builds
    and returns a WorkGraph ready to be submitted.

    Args:
        structure_path: Path to structure file (CIF format)
        structure: Pre-loaded structure (alternative to structure_path)
        code_label: Label of the VASP code in AiiDA
        kpoints_spacing: K-points spacing in A^-1 * 2pi
        potential_family: Name of the potential family
        potential_mapping: Mapping of elements to potentials
        parameters: VASP input parameters (if None, use defaults)
        options: Scheduler options (if None, use defaults)
        clean_workdir: Whether to clean work directory after completion
        name: Name of the WorkGraph (default: PS-TEROS)

    Returns:
        WorkGraph instance ready to be submitted

    Example:
        >>> wg = build_workgraph(
        ...     structure_path='/path/to/structure.cif',
        ...     code_label='VASP-6.4.3@cluster',
        ...     potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
        ...     name='MyRelaxation'
        ... )
        >>> wg.submit()
    """
    # Build the workgraph using the core_workgraph factory
    wg = core_workgraph.build(
        structure_path=structure_path,
        structure=structure,
        code_label=code_label,
        kpoints_spacing=kpoints_spacing,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=parameters,
        options=options,
        clean_workdir=clean_workdir,
    )

    # Set the name
    wg.name = name

    return wg


if __name__ == '__main__':
    """
    Simple test/example of the core WorkGraph.
    For full examples, see examples/relaxation/relaxation.py
    """
    print("PS-TEROS Core WorkGraph Module")
    print("=" * 50)
    print("\nMain function: core_workgraph")
    print("  - This is a @task.graph decorated function")
    print("  - Use core_workgraph.build() to create a WorkGraph\n")
    print("Convenience function: build_workgraph")
    print("  - Wrapper around core_workgraph.build()")
    print("  - Adds custom naming support\n")
    print("For examples, see examples/relaxation/relaxation.py")
