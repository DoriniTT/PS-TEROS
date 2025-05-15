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
from teros.functions.thermodynamics.ternary import calculate_surface_energy_ternary

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
                          workgraph_name="teros_surface_workflow", code="VASP"):
    """
    Factory function to create a generic WorkGraph for bulk and slab relaxations.

    Args:
        dft_workchain: AiiDA workchain class for DFT relaxations
        builder_bulk: Preconfigured builder for bulk relaxation
        builder_slab: Preconfigured builder for slab relaxation
        reference_builders: Dictionary of reference builders for formation enthalpy
                        Format: {element: builder_function}
        workgraph_name: Name for the WorkGraph
        code: DFT code identifier ("QUANTUM_ESPRESSO", "CP2K", "VASP")

    Returns:
        WorkGraph: Configured workflow
    """
    # Block 1: WorkGraph Context Setup
    # ---------------------------------
    # Create a new WorkGraph context with the specified workflow name.
    # All tasks and dependencies will be defined within this context.
    with WorkGraph(workgraph_name) as wg:
        # Block 2: Bulk Relaxation Task
        # -----------------------------
        # Add a task to relax the bulk structure using the provided DFT workchain.
        # The builder_bulk contains all necessary inputs for the bulk relaxation.
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
        # VASP uses 'structure', others use 'output_structure'.
        if code == 'VASP':
            output_struct_node = bulk_task.outputs.structure
        else:
            output_struct_node = bulk_task.outputs.output_structure

        # Block 3: Slab Generation Task
        # -----------------------------
        # Add a task to generate slab structures from the relaxed bulk structure.
        # This uses the get_slabs function.
        generate_slabs_task = wg.add_task(
            get_slabs,
            name="generate_slabs",
            relaxed_structure=output_struct_node
        )

        # Block 4: Slab Relaxation Tasks (Active Map Zone)
        # ------------------------------------------------
        # For each generated slab, add a relaxation task.
        # The active_map_zone context allows mapping over all generated slabs.
        with active_map_zone(generate_slabs_task.outputs.structures) as map_zone:
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
        # For each reference system (e.g., O2, metals), add a relaxation task.
        # These are needed to compute the formation enthalpy of the bulk.
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
        # Collect the outputs from the bulk and reference relaxation tasks to use as inputs
        # for the formation enthalpy calculation.
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
        # Add a task to calculate the formation enthalpy of the bulk structure.
        # This uses the outputs from the bulk and reference relaxation tasks.
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
        # Collect the outputs from the slab relaxation tasks for use in the surface thermodynamics calculation.
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
        # Use the number of unique elements in the bulk structure to determine if the compound is binary or ternary.
        bulk_atoms = bulk_structure.get_ase()
        elements = set(bulk_atoms.get_chemical_symbols())
        num_elements = len(elements)
        is_binary = num_elements == 2
        is_ternary = num_elements == 3

        # Block 10: Surface Thermodynamics Calculation Task
        # -------------------------------------------------
        # Add a task to calculate the surface thermodynamics (Gibbs free energy per area).
        # Use the binary or ternary implementation depending on the compound type.
        if is_binary:
            # For binary compounds (e.g., Ag2O), use the binary implementation.
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
            # For ternary and more complex compounds, use the ternary implementation.
            surface_thermo_task = wg.add_task(
                calculate_surface_energy_ternary,
                name="surface_thermodynamics",
                bulk_structure=bulk_struct_out,
                bulk_parameters=bulk_params_out,
                slab_structures=slab_structures_out,
                slab_parameters=slab_parameters_out,
                formation_enthalpy=enthalpy_task.outputs.formation_enthalpy,
                code=code
            )

    # Return the fully constructed WorkGraph object.
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