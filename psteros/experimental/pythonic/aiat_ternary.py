"""Ab initio atomistic thermodynamics for ternary oxide surfaces.

This module implements surface energy calculations as a function of chemical potential
for ternary oxide slabs (M-N-O systems), following the ab initio atomistic thermodynamics
framework.

Updated for current aiida-workgraph using scatter-gather pattern.
"""

from __future__ import annotations

import typing as t
from collections import Counter
from functools import reduce
from math import gcd

import numpy as np
from aiida import orm
from aiida_workgraph import dynamic, namespace, task


@task.calcfunction
def calculate_surface_energy_ternary(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Float,
    sampling: orm.Int,
) -> orm.Dict:
    """
    Compute γ(Δμ_M, Δμ_O) surface energy for a single ternary oxide slab.
    
    This calcfunction computes surface energy as a function of chemical potential
    deviations (Δμ) for a ternary oxide system M-N-O, where:
    - M: first metal (independent variable)
    - N: second metal (reference)
    - O: oxygen
    
    The surface energy γ is calculated using:
    γ(Δμ_M, Δμ_O) = φ - Γ_M·Δμ_M - Γ_O·Δμ_O
    
    where:
    - φ: reference surface energy at bulk equilibrium
    - Γ_M, Γ_O: surface excess of M and O relative to N
    - Δμ_M, Δμ_O: chemical potential deviations from reference state
    
    Args:
        bulk_structure: Bulk structure containing M, N, O
        bulk_energy: Total energy of bulk structure (eV)
        slab_structure: Slab structure
        slab_energy: Total energy of slab (eV)
        reference_energies: Dict with '{element}_energy_per_atom' keys
        formation_enthalpy: Formation enthalpy of bulk (eV)
        sampling: Number of grid points for Δμ sampling
        
    Returns:
        Dictionary containing:
        - phi: Reference surface energy
        - Gamma_M_vs_Nref: Surface excess of M
        - Gamma_O_vs_Nref: Surface excess of O
        - gamma_values_grid: Dict of γ values at grid points
        - gamma_values_fixed_muM_zero: γ values with Δμ_M = 0
        - gamma_at_reference: γ at reference conditions
        - area_A2: Slab surface area (Ų)
        - element_M_independent: Element M symbol
        - element_N_reference: Element N symbol
        - bulk_stoichiometry_MxNyOz: Reduced stoichiometry
        - slab_atom_counts: Atom counts in slab
        - reference_energies_per_atom: Reference energies used
        - E_slab_eV: Slab energy
        - E_bulk_fu_eV: Bulk energy per formula unit
        
    Raises:
        ValueError: If bulk is not a ternary oxide or missing reference data
    """
    grid_points = sampling.value
    if grid_points <= 0:
        raise ValueError('Sampling must be a positive integer')
    
    # Extract bulk composition
    bulk_ase = bulk_structure.get_ase()
    bulk_counts = Counter(bulk_ase.get_chemical_symbols())
    
    if 'O' not in bulk_counts:
        raise ValueError('The bulk structure contains no oxygen; expected a ternary oxide.')
    
    metal_elements = sorted(element for element in bulk_counts if element != 'O')
    if len(metal_elements) != 2:
        raise ValueError(
            f'Expected exactly two distinct metal species; found: {metal_elements}'
        )
    
    element_M = metal_elements[0]
    element_N_ref = metal_elements[1]
    element_O = 'O'
    
    # Get reduced stoichiometry
    stoichiometric_counts = [
        bulk_counts[element_M],
        bulk_counts[element_N_ref],
        bulk_counts[element_O]
    ]
    common_divisor = reduce(gcd, stoichiometric_counts)
    x_M = bulk_counts[element_M] // common_divisor
    y_N = bulk_counts[element_N_ref] // common_divisor
    z_O = bulk_counts[element_O] // common_divisor
    
    # Bulk energy per formula unit
    formula_units_in_bulk = bulk_counts[element_N_ref] / y_N
    bulk_energy_per_fu = bulk_energy.value / formula_units_in_bulk
    
    # Extract reference energies
    ref_data = reference_energies.get_dict()
    required_keys = [
        f'{element_M.lower()}_energy_per_atom',
        f'{element_N_ref.lower()}_energy_per_atom',
        f'{element_O.lower()}_energy_per_atom',
    ]
    for key in required_keys:
        if key not in ref_data:
            raise ValueError(f'Missing reference energy "{key}" in reference_energies')
    
    ref_energies = {
        element_M: float(ref_data[f'{element_M.lower()}_energy_per_atom']),
        element_N_ref: float(ref_data[f'{element_N_ref.lower()}_energy_per_atom']),
        element_O: float(ref_data[f'{element_O.lower()}_energy_per_atom']),
    }
    
    # Define chemical potential ranges
    delta_h = formation_enthalpy.value
    delta_mu_M_min = min(0.0, delta_h / bulk_counts[element_M])
    delta_mu_O_min = min(0.0, delta_h / bulk_counts[element_O])
    
    delta_mu_M_range = np.linspace(delta_mu_M_min, 0.0, grid_points)
    delta_mu_O_range = np.linspace(delta_mu_O_min, 0.0, grid_points)
    
    # Process slab
    slab_ase = slab_structure.get_ase()
    slab_counts = Counter(slab_ase.get_chemical_symbols())
    area = float(np.linalg.norm(np.cross(slab_ase.cell[0], slab_ase.cell[1])))
    
    N_M_slab = slab_counts.get(element_M, 0)
    N_N_slab = slab_counts.get(element_N_ref, 0)
    N_O_slab = slab_counts.get(element_O, 0)
    
    # Calculate surface excesses (Γ) relative to N
    excess_M = N_M_slab - (x_M * N_N_slab / y_N)
    excess_O = N_O_slab - (z_O * N_N_slab / y_N)
    gamma_M = excess_M / (2.0 * area)  # Factor of 2 for two surfaces
    gamma_O = excess_O / (2.0 * area)
    
    # Calculate φ (reference surface energy)
    bulk_equivalent_energy = (N_N_slab / y_N) * bulk_energy_per_fu
    term1 = (slab_energy.value - bulk_equivalent_energy) / (2.0 * area)
    term2 = -gamma_M * ref_energies[element_M]
    term3 = -gamma_O * ref_energies[element_O]
    phi = term1 + term2 + term3
    
    # Compute γ(Δμ_M, Δμ_O) grid
    gamma_grid: dict[str, float] = {}
    for delta_mu_M in delta_mu_M_range:
        for delta_mu_O in delta_mu_O_range:
            gamma = phi - gamma_M * float(delta_mu_M) - gamma_O * float(delta_mu_O)
            gamma_grid[f'muM_{delta_mu_M:.4f}_muO_{delta_mu_O:.4f}'] = float(gamma)
    
    # Compute γ with Δμ_M = 0 (special case)
    gamma_muM_zero: dict[str, float] = {}
    for delta_mu_O in delta_mu_O_range:
        gamma = phi - gamma_O * float(delta_mu_O)
        gamma_muM_zero[f'muM_0.0000_muO_{delta_mu_O:.4f}'] = float(gamma)
    
    return orm.Dict(
        dict={
            'phi': float(phi),
            'Gamma_M_vs_Nref': float(gamma_M),
            'Gamma_O_vs_Nref': float(gamma_O),
            'gamma_values_grid': gamma_grid,
            'gamma_values_fixed_muM_zero': gamma_muM_zero,
            'gamma_at_reference': float(phi),
            'area_A2': float(area),
            'element_M_independent': element_M,
            'element_N_reference': element_N_ref,
            'bulk_stoichiometry_MxNyOz': {
                f'x_{element_M}': int(x_M),
                f'x_{element_N_ref}': int(y_N),
                f'x_{element_O}': int(z_O),
            },
            'slab_atom_counts': {
                f'N_{element_M}': int(N_M_slab),
                f'N_{element_N_ref}': int(N_N_slab),
                f'N_{element_O}': int(N_O_slab),
            },
            'reference_energies_per_atom': {k: float(v) for k, v in ref_energies.items()},
            'E_slab_eV': float(slab_energy.value),
            'E_bulk_fu_eV': float(bulk_energy_per_fu),
            'formation_enthalpy_eV': float(delta_h),
            'delta_mu_M_range': [float(x) for x in delta_mu_M_range],
            'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
        }
    )


@task.graph
def compute_surface_energies_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Float,
    sampling: int = 100,
) -> t.Annotated[dict, namespace(surface_energies=dynamic(orm.Dict))]:
    """
    Scatter-gather pattern for computing surface energies.
    
    Computes surface energy γ(Δμ_M, Δμ_O) for each slab in parallel.
    
    Args:
        slabs: Dynamic namespace of slab structures
        energies: Dynamic namespace of slab energies
        bulk_structure: Bulk structure
        bulk_energy: Bulk total energy
        reference_energies: Reference energies for elements
        formation_enthalpy: Formation enthalpy of bulk
        sampling: Grid resolution for chemical potential
        
    Returns:
        Dictionary with 'surface_energies' namespace containing
        surface energy data for each slab
    """
    surface_results = {}
    sampling_node = orm.Int(sampling)
    
    # Scatter: compute surface energy for each slab in parallel
    for key, slab_structure in slabs.items():
        slab_energy = energies[key]
        
        surface_data = calculate_surface_energy_ternary(
            bulk_structure=bulk_structure,
            bulk_energy=bulk_energy,
            slab_structure=slab_structure,
            slab_energy=slab_energy,
            reference_energies=reference_energies,
            formation_enthalpy=formation_enthalpy,
            sampling=sampling_node,
        ).result
        
        surface_results[key] = surface_data
    
    # Gather: return collected results
    return {'surface_energies': surface_results}


@task.calcfunction
def create_mock_reference_energies() -> orm.Dict:
    """
    Create mock reference energies for testing.
    
    Returns typical reference energies for Ag-P-O system.
    """
    return orm.Dict(dict={
        'ag_energy_per_atom': -2.83,  # eV/atom (bulk Ag metal)
        'p_energy_per_atom': -5.35,   # eV/atom (white P)
        'o_energy_per_atom': -4.95,   # eV/atom (O2 molecule / 2)
    })


@task.calcfunction
def create_mock_bulk_energy(bulk_structure: orm.StructureData) -> orm.Float:
    """
    Create mock bulk energy for testing.
    
    Returns a fake energy proportional to number of atoms.
    """
    n_atoms = len(bulk_structure.get_ase())
    # Typical energy per atom ~ -5 eV
    return orm.Float(-5.0 * n_atoms)


@task.calcfunction
def create_mock_formation_enthalpy() -> orm.Float:
    """
    Create mock formation enthalpy for testing.
    
    Returns typical formation enthalpy for ternary oxide.
    """
    return orm.Float(-2.5)  # eV per formula unit
