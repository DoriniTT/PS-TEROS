"""
Surface Energy Calculations for Metals.

This module implements the simple surface energy formula for elemental metals:
γ = (E_slab - N·E_bulk_atom) / (2A)

No chemical potential dependencies since metals are elemental systems.
"""

from __future__ import annotations

import typing as t
from collections import Counter

import numpy as np
from aiida import orm
from aiida_workgraph import task, dynamic


# Conversion factor: 1 eV/Ų = 16.0218 J/m²
EV_PER_A2_TO_J_PER_M2 = 16.0218


@task.calcfunction
def calculate_metal_surface_energy(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
) -> orm.Dict:
    """
    Calculate surface energy for a metal slab.
    
    Uses the formula: γ = (E_slab - N_slab·E_bulk/N_bulk) / (2A)
    
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
        - element: Element symbol
    """
    # Get bulk info
    bulk_ase = bulk_structure.get_ase()
    n_bulk = len(bulk_ase)
    E_bulk_per_atom = bulk_energy.value / n_bulk
    
    # Get slab info
    slab_ase = slab_structure.get_ase()
    n_slab = len(slab_ase)
    
    # Calculate surface area from cell vectors (a × b)
    cell = slab_ase.get_cell()
    a_vec = cell[0]
    b_vec = cell[1]
    cross = np.cross(a_vec, b_vec)
    area_A2 = float(np.linalg.norm(cross))
    
    # Calculate surface energy: γ = (E_slab - N·E_bulk/atom) / (2A)
    gamma_eV_A2 = (slab_energy.value - n_slab * E_bulk_per_atom) / (2 * area_A2)
    gamma_J_m2 = gamma_eV_A2 * EV_PER_A2_TO_J_PER_M2
    
    # Get element (for single-element metals)
    symbols = set(bulk_ase.get_chemical_symbols())
    element = list(symbols)[0] if len(symbols) == 1 else ','.join(sorted(symbols))
    
    return orm.Dict(dict={
        'gamma_eV_A2': float(gamma_eV_A2),
        'gamma_J_m2': float(gamma_J_m2),
        'area_A2': float(area_A2),
        'E_slab_eV': float(slab_energy.value),
        'E_bulk_per_atom_eV': float(E_bulk_per_atom),
        'N_slab': int(n_slab),
        'N_bulk': int(n_bulk),
        'element': element,
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
