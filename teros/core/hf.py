"""
Formation Enthalpy Calculation Module

This module provides functions to calculate the enthalpy of formation for materials,
particularly ternary oxides, from DFT relaxation results.
"""

from aiida import orm
from aiida_workgraph import task
from collections import Counter
from typing import Union

from .constants import EV_TO_KJ_PER_MOL


@task.calcfunction
def calculate_formation_enthalpy(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    metal_structure: orm.StructureData,
    metal_energy: orm.Float,
    nonmetal_structure: orm.StructureData,
    nonmetal_energy: orm.Float,
    oxygen_structure: orm.StructureData,
    oxygen_energy: orm.Float,
) -> orm.Dict:
    """
    Calculate the enthalpy of formation for binary or ternary oxides.

    This function computes the formation enthalpy with respect to elemental references
    (metal, nonmetal/metal2, and O2 for oxygen). The formation enthalpy is normalized per
    formula unit and per atom, and converted to both eV and kJ/mol.

    For ternary oxides (M-N-O):
        ΔH_f = E_bulk - (n_M * E_M/atom + n_N * E_N/atom + n_O * E_O2/atom)
    
    For binary oxides (M-O):
        ΔH_f = E_bulk - (n_M * E_M/atom + n_O * E_O2/atom)
        (nonmetal reference is ignored even if provided)

    Args:
        bulk_structure: Relaxed structure of the bulk oxide
        bulk_energy: Total energy of the bulk system (eV)
        metal_structure: Relaxed structure of the metal reference
        metal_energy: Total energy of the metal reference (eV)
        nonmetal_structure: Relaxed structure of the nonmetal reference (or dummy for binary)
        nonmetal_energy: Total energy of the nonmetal reference (eV) (or dummy for binary)
        oxygen_structure: Relaxed structure of O2 molecule
        oxygen_energy: Total energy of O2 molecule (eV)

    Returns:
        Dict node containing:
            - oxide_type (str): 'binary' or 'ternary'
            - formation_enthalpy_ev (float): Formation enthalpy in eV per formula unit
            - formation_enthalpy_kjmol (float): Formation enthalpy in kJ/mol per formula unit
            - formation_enthalpy_ev_per_atom (float): Formation enthalpy in eV per atom
            - bulk_energy (float): Total energy of bulk (eV)
            - formula_units (int): Number of formula units in the bulk cell
            - elements (list): List of unique elements in the bulk
            - element_counts (dict): Count of each element in the bulk cell
            - metal_symbol (str): Chemical symbol of the metal
            - nonmetal_symbol (str): Chemical symbol of the nonmetal (or None for binary)
            - metal_energy_per_atom (float): Energy per atom of metal reference (eV)
            - nonmetal_energy_per_atom (float): Energy per atom of nonmetal reference (eV) (or None)
            - oxygen_energy_per_atom (float): Energy per atom from O2 (eV)

    Raises:
        ValueError: If element identification fails or structures are inconsistent
    """
    # Extract bulk structure information
    bulk_atoms = bulk_structure.get_ase()
    element_counts = Counter(bulk_atoms.get_chemical_symbols())
    elements = sorted(element_counts.keys())  # Sort for reproducibility

    # Identify metal and nonmetal elements
    metal_atoms = metal_structure.get_ase()
    nonmetal_atoms = nonmetal_structure.get_ase()
    oxygen_atoms = oxygen_structure.get_ase()

    metal_symbol = metal_atoms.get_chemical_symbols()[0]
    nonmetal_symbol = nonmetal_atoms.get_chemical_symbols()[0]
    oxygen_symbol = 'O'

    # Verify oxygen structure
    if oxygen_symbol not in oxygen_atoms.get_chemical_symbols():
        raise ValueError("Oxygen reference structure does not contain O atoms")

    # Determine oxide type based on bulk elements
    if len(elements) == 2 and oxygen_symbol in elements:
        # Binary oxide: M-O
        oxide_type = 'binary'
        expected_elements = sorted([metal_symbol, oxygen_symbol])
        
        if elements != expected_elements:
            raise ValueError(
                f"Binary oxide: Bulk elements {elements} do not match expected elements {expected_elements}. "
                f"Metal: {metal_symbol}, Oxygen: {oxygen_symbol}"
            )
        
        # For binary, nonmetal is not used
        use_nonmetal = False
        
    elif len(elements) == 3 and oxygen_symbol in elements:
        # Ternary oxide: M-N-O
        oxide_type = 'ternary'
        expected_elements = sorted([metal_symbol, nonmetal_symbol, oxygen_symbol])
        
        if elements != expected_elements:
            raise ValueError(
                f"Ternary oxide: Bulk elements {elements} do not match expected elements {expected_elements}. "
                f"Metal: {metal_symbol}, Nonmetal: {nonmetal_symbol}, Oxygen: {oxygen_symbol}"
            )
        
        use_nonmetal = True
        
    else:
        raise ValueError(
            f"Unexpected bulk composition: {elements}. "
            f"Expected either binary (M-O) or ternary (M-N-O) oxide."
        )

    # Determine number of formula units in the bulk cell
    # Use greatest common divisor of all element counts
    from math import gcd
    from functools import reduce

    def gcd_list(numbers):
        return reduce(gcd, numbers)

    formula_units = gcd_list(list(element_counts.values()))

    # Calculate energy per atom for each reference
    metal_count = len([s for s in metal_atoms.get_chemical_symbols() if s == metal_symbol])
    oxygen_count = len([s for s in oxygen_atoms.get_chemical_symbols() if s == oxygen_symbol])

    if metal_count == 0:
        raise ValueError(f"Metal reference structure does not contain {metal_symbol} atoms")
    if oxygen_count == 0:
        raise ValueError("Oxygen reference structure does not contain O atoms")

    metal_energy_per_atom = metal_energy.value / metal_count
    oxygen_energy_per_atom = oxygen_energy.value / oxygen_count

    # Calculate formation energy
    e_bulk = bulk_energy.value
    formation_energy = e_bulk
    formation_energy -= element_counts[metal_symbol] * metal_energy_per_atom
    formation_energy -= element_counts[oxygen_symbol] * oxygen_energy_per_atom

    # For ternary, subtract nonmetal contribution
    nonmetal_energy_per_atom_val = None
    if use_nonmetal:
        nonmetal_count = len([s for s in nonmetal_atoms.get_chemical_symbols() if s == nonmetal_symbol])
        if nonmetal_count == 0:
            raise ValueError(f"Nonmetal reference structure does not contain {nonmetal_symbol} atoms")
        
        nonmetal_energy_per_atom_val = nonmetal_energy.value / nonmetal_count
        formation_energy -= element_counts[nonmetal_symbol] * nonmetal_energy_per_atom_val

    # Normalize per formula unit
    formation_energy_per_fu = formation_energy / formula_units

    # Convert to kJ/mol
    formation_energy_kjmol = formation_energy_per_fu * EV_TO_KJ_PER_MOL

    # Calculate per atom
    total_atoms = sum(element_counts.values())
    atoms_per_formula_unit = total_atoms / formula_units
    formation_energy_per_atom = formation_energy_per_fu / atoms_per_formula_unit

    # Prepare results
    results = {
        'oxide_type': oxide_type,
        'formation_enthalpy_ev': formation_energy_per_fu,
        'formation_enthalpy_kjmol': formation_energy_kjmol,
        'formation_enthalpy_ev_per_atom': formation_energy_per_atom,
        'bulk_energy': e_bulk,
        'formula_units': int(formula_units),
        'elements': elements,
        'element_counts': {element: int(count) for element, count in element_counts.items()},
        'metal_symbol': metal_symbol,
        'nonmetal_symbol': nonmetal_symbol if use_nonmetal else None,
        'metal_energy_per_atom': metal_energy_per_atom,
        'nonmetal_energy_per_atom': nonmetal_energy_per_atom_val if use_nonmetal else None,
        'oxygen_energy_per_atom': oxygen_energy_per_atom,
    }

    return orm.Dict(dict=results)
