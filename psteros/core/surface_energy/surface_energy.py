"""
Surface Energy Calculations for Metals and Intermetallics.

This module implements the simple surface energy formula:
γ = (E_slab - N·E_bulk_atom) / (2A)

This formula applies to:
- Elemental metals (Au, Ag, Cu, Pt, etc.)
- Stoichiometric and symmetric intermetallic surfaces (PdIn, AuCu, NiAl, etc.)

For stoichiometric slabs where the surface composition matches the bulk
stoichiometry and both surfaces are equivalent (symmetric slab), no chemical
potential dependencies are needed. The simple formula above is sufficient.

Note: For non-stoichiometric or asymmetric intermetallic surfaces, chemical
potential-dependent formulations are required (to be implemented separately).
"""

from __future__ import annotations

import typing as t
from collections import Counter

import numpy as np
from aiida import orm
from aiida_workgraph import task, dynamic

# Conversion factor: 1 eV/Ų = 16.0218 J/m²
EV_PER_A2_TO_J_PER_M2 = 16.0218


def identify_compound_type(structure: orm.StructureData) -> dict:
    """
    Identify whether a structure is an elemental metal or intermetallic.

    Args:
        structure: The structure to analyze

    Returns:
        Dictionary containing:
        - compound_type: 'elemental' or 'intermetallic'
        - elements: List of element symbols
        - composition: Dict mapping element to count
        - formula: Chemical formula string (e.g., 'Au', 'PdIn', 'Au3Cu')
        - stoichiometry: Dict mapping element to stoichiometric ratio
    """
    ase_atoms = structure.get_ase()
    symbols = ase_atoms.get_chemical_symbols()
    composition = dict(Counter(symbols))
    elements = sorted(composition.keys())

    # Determine compound type
    if len(elements) == 1:
        compound_type = 'elemental'
        formula = elements[0]
    else:
        compound_type = 'intermetallic'
        # Build formula string (e.g., 'PdIn' or 'Au3Cu')
        # Find GCD to reduce stoichiometry
        from math import gcd
        from functools import reduce
        counts = list(composition.values())
        common_divisor = reduce(gcd, counts)
        reduced_counts = {el: composition[el] // common_divisor for el in elements}

        # Build formula
        formula_parts = []
        for el in elements:
            count = reduced_counts[el]
            if count == 1:
                formula_parts.append(el)
            else:
                formula_parts.append(f'{el}{count}')
        formula = ''.join(formula_parts)

    # Calculate stoichiometric ratios
    total_atoms = sum(composition.values())
    stoichiometry = {el: composition[el] / total_atoms for el in elements}

    return {
        'compound_type': compound_type,
        'elements': elements,
        'composition': composition,
        'formula': formula,
        'stoichiometry': stoichiometry,
    }

@task.calcfunction
def calculate_metal_surface_energy(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
) -> orm.Dict:
    """
    Calculate surface energy for a metal or stoichiometric intermetallic slab.

    Uses the formula: γ = (E_slab - N_slab·E_bulk/N_bulk) / (2A)

    This formula is valid for:
    - Elemental metals (Au, Ag, Cu, Pt, etc.)
    - Stoichiometric and symmetric intermetallic surfaces (PdIn, AuCu, NiAl, etc.)

    Args:
        bulk_structure: Bulk crystal structure
        bulk_energy: Total energy of bulk structure (eV)
        slab_structure: Slab structure
        slab_energy: Total energy of slab (eV)

    Returns:
        Dictionary containing:
        - gamma_eV_A2: Surface energy in eV/Ų
        - gamma_J_m2: Surface energy in J/m²
        - area_A2: Surface area in Ų
        - E_slab_eV: Slab total energy
        - E_bulk_per_atom_eV: Bulk energy per atom
        - N_slab: Number of atoms in slab
        - N_bulk: Number of atoms in bulk
        - compound_type: 'elemental' or 'intermetallic'
        - formula: Chemical formula (e.g., 'Au', 'PdIn', 'Au3Cu')
        - elements: List of element symbols
        - composition: Dict with element counts in bulk
        - slab_composition: Dict with element counts in slab
        - is_stoichiometric: Whether slab has same stoichiometry as bulk
    """
    # Get bulk info and composition
    bulk_ase = bulk_structure.get_ase()
    n_bulk = len(bulk_ase)
    E_bulk_per_atom = bulk_energy.value / n_bulk
    bulk_comp_info = identify_compound_type(bulk_structure)

    # Get slab info and composition
    slab_ase = slab_structure.get_ase()
    n_slab = len(slab_ase)
    slab_symbols = slab_ase.get_chemical_symbols()
    slab_composition = dict(Counter(slab_symbols))

    # Check if slab is stoichiometric (same composition ratio as bulk)
    bulk_stoich = bulk_comp_info['stoichiometry']
    slab_total = sum(slab_composition.values())
    slab_stoich = {el: slab_composition.get(el, 0) / slab_total for el in bulk_comp_info['elements']}

    # Compare stoichiometries (with small tolerance for numerical precision)
    is_stoichiometric = all(
        abs(bulk_stoich.get(el, 0) - slab_stoich.get(el, 0)) < 0.01
        for el in set(list(bulk_stoich.keys()) + list(slab_stoich.keys()))
    )

    # Calculate surface area from cell vectors (a × b)
    cell = slab_ase.get_cell()
    a_vec = cell[0]
    b_vec = cell[1]
    cross = np.cross(a_vec, b_vec)
    area_A2 = float(np.linalg.norm(cross))

    # Calculate surface energy: γ = (E_slab - N·E_bulk/atom) / (2A)
    gamma_eV_A2 = (slab_energy.value - n_slab * E_bulk_per_atom) / (2 * area_A2)
    gamma_J_m2 = gamma_eV_A2 * EV_PER_A2_TO_J_PER_M2

    return orm.Dict(dict={
        'gamma_eV_A2': float(gamma_eV_A2),
        'gamma_J_m2': float(gamma_J_m2),
        'area_A2': float(area_A2),
        'E_slab_eV': float(slab_energy.value),
        'E_bulk_per_atom_eV': float(E_bulk_per_atom),
        'N_slab': int(n_slab),
        'N_bulk': int(n_bulk),
        'compound_type': bulk_comp_info['compound_type'],
        'formula': bulk_comp_info['formula'],
        'elements': bulk_comp_info['elements'],
        'composition': bulk_comp_info['composition'],
        'slab_composition': slab_composition,
        'is_stoichiometric': is_stoichiometric,
    })


@task.calcfunction
def calculate_bulk_energy_per_atom(
    bulk_energy: orm.Float,
    bulk_structure: orm.StructureData,
) -> orm.Float:
    """
    Calculate bulk energy per atom.
    
    Args:
        bulk_energy: Total energy of bulk structure
        bulk_structure: Bulk crystal structure
        
    Returns:
        Energy per atom as Float
    """
    n_atoms = len(bulk_structure.get_ase())
    return orm.Float(bulk_energy.value / n_atoms)


@task.graph
def compute_metal_surface_energies_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
) -> t.Annotated[dict, dict]:
    """
    Scatter-gather pattern for computing surface energies for multiple slabs.
    
    Args:
        slabs: Dictionary of slab structures keyed by termination ID
        energies: Dictionary of slab energies keyed by termination ID
        bulk_structure: Bulk crystal structure
        bulk_energy: Total energy of bulk structure
        
    Returns:
        Dictionary with 'surface_energies' containing results for each slab
    """
    surface_results = {}
    
    # Compute surface energy for each slab in parallel
    for key, slab_structure in slabs.items():
        slab_energy = energies[key]
        
        surface_data = calculate_metal_surface_energy(
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energy,
        ).result
        
        surface_results[key] = surface_data
    
    return {'surface_energies': surface_results}
