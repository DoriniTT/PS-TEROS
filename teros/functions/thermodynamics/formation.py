from aiida.orm import Dict
from aiida_workgraph import task
from collections import Counter
from aiida.orm import Dict

@task.calcfunction(outputs=[{"name": "formation_enthalpy"}])
def calculate_formation_enthalpy(bulk_structure, bulk_parameters, code=None, **reference_systems):
    """
    Calculate the enthalpy of formation for a material, typically an oxide.

    The formation enthalpy is calculated with respect to the provided elemental references
    (or O2 for oxygen). The function handles parsing energies from different DFT codes
    and normalizes the formation enthalpy per formula unit and per atom.

    :param bulk_structure: AiiDA ``StructureData`` node of the bulk material.
    :type bulk_structure: aiida.orm.StructureData
    :param bulk_parameters: AiiDA ``Dict`` node containing the output parameters (including total energy)
                            of the DFT calculation for the bulk material.
    :type bulk_parameters: aiida.orm.Dict
    :param code: DFT code identifier, used to correctly parse energy from `bulk_parameters`
                 and reference parameters. Supported: "VASP", "QUANTUM_ESPRESSO", "CP2K".
    :type code: str, optional
    :param reference_systems: Keyword arguments where each key-value pair represents a reference system.
                              For each element in the bulk material (e.g., 'Ag', 'P', 'O'),
                              two keyword arguments are expected:
                              - ``<element_lower>_structure``: AiiDA ``StructureData`` for the reference (e.g., ``ag_structure``).
                              - ``<element_lower>_parameters``: AiiDA ``Dict`` with output parameters for the reference (e.g., ``ag_parameters``).
                              For oxygen, the keys are typically ``o2_structure`` and ``o2_parameters``.
                              The function also checks for uppercase element names as a fallback (e.g., ``AG_structure``).
    :type reference_systems: dict[str, aiida.orm.StructureData or aiida.orm.Dict]

    :raises ValueError: If reference data (structure or parameters) is missing for an element
                        in the bulk material.
    :raises ValueError: If a reference structure for an element does not contain that element.

    :return: A dictionary where the key 'formation_enthalpy' maps to an AiiDA ``Dict`` node.
             This Dict node contains the following key-value pairs:
             - ``formation_enthalpy_ev`` (float): Formation enthalpy in eV per formula unit.
             - ``formation_enthalpy_kjmol`` (float): Formation enthalpy in kJ/mol per formula unit.
             - ``formation_enthalpy_ev_per_atom`` (float): Formation enthalpy in eV per atom.
             - ``bulk_energy`` (float): Total energy of the bulk material in eV.
             - ``formula_units`` (int): Number of formula units in the provided bulk cell.
             - ``elements`` (list[str]): List of unique chemical symbols in the bulk material.
             - ``element_counts`` (dict[str, int]): Count of each element in the bulk cell.
             - ``<element_lower>_energy_per_atom`` (float): Energy per atom for each reference element (e.g., ``ag_energy_per_atom``).
    :rtype: dict[str, aiida.orm.Dict]
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