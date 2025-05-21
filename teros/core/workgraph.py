"""
Generic DFT workflow for structure relaxations using AiiDA Workgraphs.

This module defines a workflow that can be used with multiple DFT codes
(CP2K, Quantum ESPRESSO, VASP, etc.) by changing the imported builders.
It focuses on defining tasks and executing the workflow,
while all code-specific configuration is imported from builder modules.
"""
from aiida_workgraph import WorkGraph, task, active_map_zone
from teros.functions.slabs import get_slabs
from teros.functions.thermodynamics.formation import calculate_formation_enthalpy
from teros.functions.thermodynamics.binary import calculate_surface_energy_binary
from teros.functions.thermodynamics.ternary import calculate_surface_energy_ternar
from aiida.orm import Dict, Int, Float, Bool, List
from aiida import load_profile
load_profile()

# Module docstring to explain the workflow
"""
Teros: Modular workflow for surface thermodynamics of oxides using AiiDA Workgraphs.

This module automates bulk and slab relaxations, reference calculations, formation enthalpy,
and surface Gibbs free energy (γ) for binary and ternary oxides. It is compatible with
Quantum ESPRESSO, CP2K, and VASP, and uses AiiDA WorkGraph for workflow construction.

Main steps:
  1. Bulk relaxation.
  2. Slab generation and relaxation.
  3. Reference relaxations (O2, metals).
  4. Formation enthalpy calculation.
  5. Surface thermodynamics (γ) for each slab.
  6. WorkGraph assembly for automated execution.

Usage:
    wg = create_teros_workgraph(dft_workchain=MyDFTWorkChain, builder_bulk=..., 
                             builder_slab=..., reference_builders=..., code="VASP")
    wg.submit()
"""

@task.graph_builder(outputs=[{"name": "result", "from": "ctx.result"}])
def create_teros_workgraph(dft_workchain, builder_bulk, builder_slab, reference_builders,
                          workgraph_name="teros_surface_workflow", code="VASP",
                          manual_slabs=None,
                          sampling: Int = None,
                          # Slab generation parameters
                          miller_indices: List = None,
                          min_slab_thickness: Float = None,
                          min_vacuum_thickness: Float = None,
                          lll_reduce: Bool = None,
                          center_slab: Bool = None,
                          symmetrize: Bool = None,
                          primitive: Bool = None,
                          in_unit_planes: Bool = None,
                          max_normal_search: Int = None):
    """
    Factory function to create a TEROS WorkGraph for surface thermodynamics calculations.

    This function assembles an AiiDA WorkGraph that performs:
    1. Bulk structure relaxation.
    2. Slab generation (if not manually provided) and relaxation for specified Miller indices.
    3. Relaxation of reference structures (e.g., elemental phases, O2 molecule).
    4. Calculation of formation enthalpy for the bulk material.
    5. Calculation of surface energies for the generated slabs, potentially as a function
       of chemical potentials for ternary systems.

    :param dft_workchain: The AiiDA DFT workchain class to be used for all relaxations
                          (e.g., ``WorkflowFactory('vasp.vasp')`` or ``WorkflowFactory('quantumespresso.pw.relax')``).
    :type dft_workchain: aiida.plugins.WorkflowFactory
    :param builder_bulk: A pre-configured AiiDA ``ProcessBuilder`` instance for the bulk relaxation.
                         This builder should have the structure, DFT parameters, and resource options set.
    :type builder_bulk: aiida.engine.ProcessBuilder
    :param builder_slab: A pre-configured AiiDA ``ProcessBuilder`` instance for slab relaxations.
                         This builder should have DFT parameters and resource options set. The structure
                         for each slab will be injected by the workgraph.
    :type builder_slab: aiida.engine.ProcessBuilder
    :param reference_builders: A dictionary mapping reference material identifiers (e.g., 'O2', 'Ag')
                               to their pre-configured AiiDA ``ProcessBuilder`` instances for relaxation.
                               Example: ``{'O2': o2_builder, 'Ag': ag_builder}``
    :type reference_builders: dict[str, aiida.engine.ProcessBuilder]
    :param workgraph_name: Name for the generated WorkGraph.
    :type workgraph_name: str, optional
    :param code: DFT code identifier, used to handle code-specific output parsing.
                 Supported codes: "VASP", "QUANTUM_ESPRESSO", "CP2K".
    :type code: str, optional
    :param manual_slabs: A dictionary mapping slab identifiers (str) to AiiDA ``StructureData`` nodes
                         for manually created/provided slab structures. If provided, automated slab
                         generation via ``get_slabs`` is skipped. Default is None (automatic generation).
    :type manual_slabs: dict[str, aiida.orm.StructureData], optional
    :param sampling: An AiiDA ``Int`` node specifying the number of sampling points for the
                     chemical potential grid in ternary surface phase diagrams.
                     Only used for ternary compounds. Default is None.
    :type sampling: aiida.orm.Int, optional
    :param miller_indices: An AiiDA ``List`` of Miller indices for which slabs should be generated
                           (e.g., ``List(list=[[1,0,0], [1,1,0]])``). Used if ``manual_slabs`` is not provided.
    :type miller_indices: aiida.orm.List, optional
    :param min_slab_thickness: An AiiDA ``Float`` for the minimum slab thickness in Angstroms for slab generation.
    :type min_slab_thickness: aiida.orm.Float, optional
    :param min_vacuum_thickness: An AiiDA ``Float`` for the minimum vacuum thickness in Angstroms for slab generation.
    :type min_vacuum_thickness: aiida.orm.Float, optional
    :param lll_reduce: An AiiDA ``Bool`` indicating whether to perform LLL reduction on the cell for slab generation.
    :type lll_reduce: aiida.orm.Bool, optional
    :param center_slab: An AiiDA ``Bool`` indicating whether to center the slab in the simulation cell.
    :type center_slab: aiida.orm.Bool, optional
    :param symmetrize: An AiiDA ``Bool`` indicating whether to generate symmetrically distinct slab terminations.
    :type symmetrize: aiida.orm.Bool, optional
    :param primitive: An AiiDA ``Bool`` indicating whether to use the primitive cell for slab generation.
    :type primitive: aiida.orm.Bool, optional
    :param in_unit_planes: An AiiDA ``Bool`` for restricting slab generation to unit planes.
    :type in_unit_planes: aiida.orm.Bool, optional
    :param max_normal_search: An AiiDA ``Int`` for the maximum normal search distance in Pymatgen's slab generation.
                              Default is None, using Pymatgen's default.
    :type max_normal_search: aiida.orm.Int, optional

    :raises ValueError: If ``reference_builders`` is not provided or empty.
    :raises ValueError: If the bulk structure cannot be extracted from ``builder_bulk``.

    :return: The configured AiiDA WorkGraph instance.
    :rtype: aiida_workgraph.WorkGraph
    """
    # Block 1: WorkGraph Context Setup
    # ---------------------------------
    # Create a new WorkGraph context with the specified workflow name.
    # All tasks and dependencies will be defined within this context.
    with WorkGraph(workgraph_name) as wg:
        # Block 2: Bulk Relaxation Task
        # -----------------------------
        bulk_task = wg.add_task(dft_workchain, name="bulk_relaxation")
        bulk_task.set_from_builder(builder_bulk)

        # Extract the bulk structure from the builder for later use (e.g., to determine elements).
        bulk_structure = builder_bulk.structure if hasattr(builder_bulk, 'structure') else None
        # For CP2K, the structure may be nested under the 'cp2k' attribute.
        if bulk_structure is None and hasattr(builder_bulk, 'cp2k') and 'structure' in builder_bulk.cp2k:
            bulk_structure = builder_bulk.cp2k.structure

        # If the structure could not be extracted, raise an error.
        if bulk_structure is None:
            raise ValueError("Could not extract structure from bulk builder")

        # Get the list of elements present in the bulk structure for later logic.
        elements = list(set(bulk_structure.get_ase().get_chemical_symbols()))

        # Determine the correct output node for the relaxed structure depending on the DFT code.
        if code == 'VASP':
            output_struct_node = bulk_task.outputs.structure
        else:
            output_struct_node = bulk_task.outputs.output_structure

        if manual_slabs:
            # Block 3 (Alternative): Use manually provided slab structures
            # ----------------------------------------------------------
            # When manual slabs are provided, skip the get_slabs function
            slab_structures = manual_slabs
            print(f"Using {len(manual_slabs)} manually provided slab structures")
        else:
            # Block 3: Slab Generation Task
            # -----------------------------
            # Add a task to generate slab structures from the relaxed bulk structure.
            generate_slabs_task = wg.add_task(
                get_slabs,
                name="generate_slabs",
                relaxed_structure=output_struct_node,
                # Pass slab generation parameters
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
            
            # Use the generated slabs for the rest of the workflow
            slab_structures = generate_slabs_task.outputs.structures

        # Block 4: Slab Relaxation Tasks (Active Map Zone)
        # ------------------------------------------------
        with active_map_zone(slab_structures) as map_zone:
            # Depending on the DFT code, set up the slab relaxation task with the correct input keys.
            if code == 'QUANTUM_ESPRESSO':
                slab_task = map_zone.add_task(
                    dft_workchain,
                    name="slab_relaxation",
                    structure=map_zone.item,
                )
            elif code == 'CP2K':
                slab_task = map_zone.add_task(
                    dft_workchain,
                    name="slab_relaxation",
                    cp2k={'structure': map_zone.item}
                )
            elif code == 'VASP':
                slab_task = map_zone.add_task(
                    dft_workchain,
                    name="slab_relaxation",
                    structure=map_zone.item,
                )
            # Set the builder for the slab relaxation task (common for all codes).
            slab_task.set_from_builder(builder_slab)

        # Block 5: Reference Calculations for Formation Enthalpy
        # -----------------------------------------------------
        reference_tasks = {}

        # Ensure that reference builders are provided.
        if not reference_builders:
            raise ValueError("Reference builders must be provided. They are required for formation enthalpy calculations.")

        # For each reference element, add a relaxation task and store it for later use.
        for element, builder in reference_builders.items():
            task_name = f"{element.lower()}_relaxation"
            ref_task = wg.add_task(dft_workchain, name=task_name)
            ref_task.set_from_builder(builder)
            reference_tasks[element] = ref_task

        # Block 6: Prepare Inputs for Formation Enthalpy Calculation
        # ---------------------------------------------------------
        enthalpy_kwargs = {}

        # For VASP, the output keys are 'structure' and 'misc'.
        if code == 'VASP':
            bulk_struct_out = bulk_task.outputs.structure
            bulk_params_out = bulk_task.outputs.misc

            for element, task in reference_tasks.items():
                struct_key = f"{element.lower()}_structure"
                param_key = f"{element.lower()}_parameters"
                enthalpy_kwargs[struct_key] = task.outputs.structure
                enthalpy_kwargs[param_key] = task.outputs.misc
        # For other codes, use 'output_structure' and 'output_parameters'.
        else:
            bulk_struct_out = bulk_task.outputs.output_structure
            bulk_params_out = bulk_task.outputs.output_parameters

            for element, task in reference_tasks.items():
                struct_key = f"{element.lower()}_structure"
                param_key = f"{element.lower()}_parameters"
                enthalpy_kwargs[struct_key] = task.outputs.output_structure
                enthalpy_kwargs[param_key] = task.outputs.output_parameters

        # Block 7: Formation Enthalpy Calculation Task
        # --------------------------------------------
        enthalpy_task = wg.add_task(
            calculate_formation_enthalpy,
            name="formation_enthalpy",
            bulk_structure=bulk_struct_out,
            bulk_parameters=bulk_params_out,
            code=code,
            **enthalpy_kwargs
        )

        # Block 8: Prepare Inputs for Surface Thermodynamics Task
        # ------------------------------------------------------
        if code == 'VASP':
            bulk_struct_out = bulk_task.outputs.structure
            bulk_params_out = bulk_task.outputs.misc
            slab_structures_out = slab_task.outputs.structure
            slab_parameters_out = slab_task.outputs.misc
        else:
            bulk_struct_out = bulk_task.outputs.output_structure
            bulk_params_out = bulk_task.outputs.output_parameters
            slab_structures_out = slab_task.outputs.output_structure
            slab_parameters_out = slab_task.outputs.output_parameters

        # Block 9: Determine Compound Type (Binary or Ternary)
        # ----------------------------------------------------
        bulk_atoms = bulk_structure.get_ase()
        elements = set(bulk_atoms.get_chemical_symbols())
        num_elements = len(elements)
        is_binary = num_elements == 2
        is_ternary = num_elements == 3

        # Block 10: Surface Thermodynamics Calculation Task
        # -------------------------------------------------
        if is_binary:
            surface_thermo_task = wg.add_task(
                calculate_surface_energy_binary,
                name="surface_thermodynamics",
                bulk_structure=bulk_struct_out,
                bulk_parameters=bulk_params_out,
                slab_structures=slab_structures_out,
                slab_parameters=slab_parameters_out,
                formation_enthalpy=enthalpy_task.outputs.formation_enthalpy,
                code=code
            )
        elif is_ternary:
            surface_thermo_task = wg.add_task(
                calculate_surface_energy_ternar,
                name="surface_thermodynamics",
                bulk_structure=bulk_struct_out,
                bulk_parameters=bulk_params_out,
                sampling=sampling,
                slab_structures=slab_structures_out,
                slab_parameters=slab_parameters_out,
                formation_enthalpy=enthalpy_task.outputs.formation_enthalpy,
                code=code
            )

    return wg

# Example of how to use the module:
# from aiida.plugins import WorkflowFactory
# 
# # Configure builders for bulk and slab calculations
# builder_bulk = ...  # Configure your bulk builder
# builder_slab = ...  # Configure your slab builder
# reference_builders = {
#     'O': o2_builder,      # O2 molecule reference
#     'Metal1': metal1_builder,  # First metal reference
#     'Metal2': metal2_builder,  # Second metal reference (for ternary compounds)
# }
# 
# # Get the workgraph with a specific DFT workchain
# dft_workchain = WorkflowFactory('quantumespresso.pw.relax')
# workgraph = create_teros_workgraph(
#     dft_workchain=dft_workchain, 
#     builder_bulk=builder_bulk,
#     builder_slab=builder_slab,
#     reference_builders=reference_builders, 
#     code="QUANTUM_ESPRESSO"
# )
# 
# # Submit the workgraph
# from aiida.engine import submit
# results = submit(workgraph)