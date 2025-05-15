from aiida.orm import Dict
from aiida_workgraph import task
from collections import Counter
from aiida.orm import Dict

@task.calcfunction(outputs=[{"name": "formation_enthalpy"}])
def calculate_formation_enthalpy(bulk_structure, bulk_parameters, code=None, **reference_systems):
    """
    Calculate the enthalpy of formation for any oxide material.

    Args:
        bulk_structure: AiiDA StructureData of the bulk oxide
        bulk_parameters: Dict with energy data for the oxide
        code: DFT code identifier ("QUANTUM_ESPRESSO", "CP2K", "VASP")
        **reference_systems: Dictionary of reference systems with structure and parameters
                            Format: {element_name_structure: structure, element_name_parameters: parameters}

    Returns:
        Dict with enthalpy of formation in eV
    """
    # Step 1: Extract the atomic structure and count elements in the bulk
    bulk_atoms = bulk_structure.get_ase()
    element_counts = Counter(bulk_atoms.get_chemical_symbols())
    elements = list(element_counts.keys())

    # Step 2: Determine the number of formula units in the bulk cell
    formula_units = 1
    if all(count % 2 == 0 for count in element_counts.values()):
        formula_units = 2
    elif all(count % 3 == 0 for count in element_counts.values()):
        formula_units = 3
    elif all(count % 4 == 0 for count in element_counts.values()):
        formula_units = 4

    # Step 3: Extract the total energy of the bulk structure
    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        e_bulk = bulk_parameters.get_dict()['energy']
    else:  # VASP and others
        e_bulk = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']

    # Step 4: Extract reference energies for each element
    element_energy_per_atom = {}
    for element in elements:
        # Construct the expected keys for the structure and parameters in the reference_systems dict.
        struct_key = f"{element.lower()}_structure"
        param_key = f"{element.lower()}_parameters"

        # Special handling for oxygen: often referenced as O2 molecule.
        if struct_key not in reference_systems or param_key not in reference_systems:
            if element == 'O':
                struct_key = "o2_structure"
                param_key = "o2_parameters"
            # Try uppercase keys as a fallback.
            if struct_key not in reference_systems:
                struct_key = f"{element.upper()}_structure"
            if param_key not in reference_systems:
                param_key = f"{element.upper()}_parameters"

        # Check if both structure and parameters are available for this element.
        if struct_key in reference_systems and param_key in reference_systems:
            structure = reference_systems[struct_key]
            parameters = reference_systems[param_key]
            # Count how many atoms of this element are present in the reference structure.
            atom_symbols = structure.get_ase().get_chemical_symbols()
            count = atom_symbols.count(element)
            if count > 0:
                # Extract the total energy of the reference system.
                if code in ['QUANTUM_ESPRESSO', 'CP2K']:
                    energy = parameters.get_dict()['energy']
                else:
                    energy = parameters.get_dict()['total_energies']['energy_extrapolated']
                # Calculate the energy per atom for this element.
                element_energy_per_atom[element] = energy / count
            else:
                raise ValueError(f"Reference structure for {element} does not contain any {element} atoms")
        else:
            # If reference data is missing, raise an error and show available keys for debugging.
            available_keys = ', '.join(reference_systems.keys())
            raise ValueError(f"Missing reference data for element {element}. Looking for keys '{struct_key}' and '{param_key}'.\nAvailable keys: {available_keys}")

    # Step 5: Calculate the formation energy of the bulk
    formation_energy = e_bulk
    for element, count in element_counts.items():
        formation_energy -= count * element_energy_per_atom[element]

    # Step 6: Normalize the formation energy per formula unit
    formation_energy_per_fu = formation_energy / formula_units

    # Step 7: Convert the formation energy to kJ/mol
    formation_energy_kjmol = formation_energy_per_fu * 96.485

    # Step 8: Calculate the formation energy per atom
    total_atoms = sum(element_counts.values())
    formation_energy_per_atom = formation_energy_per_fu / (total_atoms / formula_units)

    # Step 9: Prepare the results dictionary
    results = {
        'formation_enthalpy_ev': formation_energy_per_fu,
        'formation_enthalpy_kjmol': formation_energy_kjmol,
        'formation_enthalpy_ev_per_atom': formation_energy_per_atom,
        'bulk_energy': e_bulk,
        'formula_units': formula_units,
        'elements': elements,
        'element_counts': {element: int(count) for element, count in element_counts.items()},
    }
    for element, energy in element_energy_per_atom.items():
        results[f'{element.lower()}_energy_per_atom'] = energy

    # Step 10: Return the results as an AiiDA Dict node
    return {'formation_enthalpy': Dict(dict=results)}