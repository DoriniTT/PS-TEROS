"""Ab initio atomistic thermodynamics for oxide surfaces.

This module implements surface energy calculations as a function of chemical potential
for both binary and ternary oxide slabs, following the ab initio atomistic thermodynamics
framework.
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
def identify_oxide_type(bulk_structure: orm.StructureData) -> orm.Str:
    """
    Identify whether a bulk structure is a binary or ternary oxide.
    
    Args:
        bulk_structure: Bulk structure to analyze
        
    Returns:
        Str node with value 'binary' or 'ternary'
        
    Raises:
        ValueError: If structure is not a binary or ternary oxide
    """
    bulk_ase = bulk_structure.get_ase()
    bulk_counts = Counter(bulk_ase.get_chemical_symbols())
    
    if 'O' not in bulk_counts:
        raise ValueError('Structure contains no oxygen; not an oxide.')
    
    metal_elements = sorted(element for element in bulk_counts if element != 'O')
    
    if len(metal_elements) == 1:
        oxide_type = 'binary'
    elif len(metal_elements) == 2:
        oxide_type = 'ternary'
    else:
        raise ValueError(
            f'Found {len(metal_elements)} non-oxygen elements: {metal_elements}. '
            'Only binary (1 metal) and ternary (2 metals) oxides are supported.'
        )
    
    return orm.Str(oxide_type)


@task.calcfunction
def calculate_surface_energy_ternary(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Dict,
    sampling: orm.Int,
) -> orm.Dict:
    """
    Compute γ(Δμ_M, Δμ_O) surface energy for a single ternary oxide slab.
    
    This calcfunction computes surface energy as a function of chemical potential
    deviations (Δμ) for a ternary oxide system M-N-O, where:
    - M: first metal (element A, independent variable)
    - N: second metal (element B, reference)
    - O: oxygen
    
    The surface energy γ is calculated using:
    γ(Δμ_M, Δμ_O) = φ - Γ_M·Δμ_M - Γ_O·Δμ_O
    
    where:
    - φ: reference surface energy at bulk equilibrium
    - Γ_M, Γ_O: surface excess of M and O relative to N
    - Δμ_M, Δμ_O: chemical potential deviations from reference state
    
    This function also computes the alternative formulation with B as independent
    variable and returns both results in the output dictionary.
    
    Args:
        bulk_structure: Bulk structure containing M, N, O
        bulk_energy: Total energy of bulk structure (eV)
        slab_structure: Slab structure
        slab_energy: Total energy of slab (eV)
        reference_energies: Dict from formation_enthalpy calculation with '*_energy_per_atom' keys
        formation_enthalpy: Dict from formation_enthalpy calculation with 'formation_enthalpy_ev' key
        sampling: Number of grid points for Δμ sampling
        
    Returns:
        Dictionary containing surface energy data for both A-based and B-based formulations
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
    
    # Extract reference energies and formation enthalpy
    ref_data = reference_energies.get_dict()
    formation_data = formation_enthalpy.get_dict()
    
    ref_energies = {
        element_M: float(ref_data['metal_energy_per_atom']),
        element_N_ref: float(ref_data['nonmetal_energy_per_atom']),
        element_O: float(ref_data['oxygen_energy_per_atom']),
    }
    
    delta_h = float(formation_data['formation_enthalpy_ev'])
    
    # Calculate surface area
    slab_ase = slab_structure.get_ase()
    cell = slab_ase.get_cell()
    a_vec = cell[0]
    b_vec = cell[1]
    cross = np.cross(a_vec, b_vec)
    area = float(np.linalg.norm(cross))
    
    # Slab atom counts
    slab_counts = Counter(slab_ase.get_chemical_symbols())
    N_M_slab = slab_counts.get(element_M, 0)
    N_N_slab = slab_counts.get(element_N_ref, 0)
    N_O_slab = slab_counts.get(element_O, 0)
    
    # Surface excess relative to N
    gamma_M = (N_M_slab - (x_M / y_N) * N_N_slab) / (2 * area)
    gamma_O = (N_O_slab - (z_O / y_N) * N_N_slab) / (2 * area)
    
    # Reference surface energy φ
    # Formula: φ = [E_slab - N_M·E_M - N_N·E_N - N_O·E_O - (N_N/y_N)·ΔH_f] / (2A)
    # This corresponds to γ at Δμ_M = Δμ_O = 0 (equilibrium with pure elements)
    phi = (
        slab_energy.value
        - N_M_slab * ref_energies[element_M]
        - N_N_slab * ref_energies[element_N_ref]
        - N_O_slab * ref_energies[element_O]
        - (N_N_slab / y_N) * delta_h
    ) / (2 * area)
    
    # Generate grid for chemical potentials
    # Physical convention: Δμ = μ - μ_ref where μ_ref = E_element
    # - Element-rich limit: Δμ = 0 (equilibrium with pure element)
    # - Decomposition limit: Δμ = ΔH_f / stoich (stability boundary)
    # Since ΔH_f < 0 for stable oxides, Δμ ranges from negative to 0
    delta_mu_M_range = np.linspace(delta_h / x_M, 0, grid_points)
    delta_mu_O_range = np.linspace(delta_h / z_O, 0, grid_points)
    
    # Compute γ(Δμ_M, Δμ_O) on a 2D grid
    # Store as 2D array where gamma_grid[i][j] = γ(Δμ_M[i], Δμ_O[j])
    gamma_grid_2d = []
    for delta_mu_M in delta_mu_M_range:
        gamma_row = []
        for delta_mu_O in delta_mu_O_range:
            gamma = phi - gamma_M * float(delta_mu_M) - gamma_O * float(delta_mu_O)
            gamma_row.append(float(gamma))
        gamma_grid_2d.append(gamma_row)
    
    # Compute γ with Δμ_M = 0 (1D slice for oxygen-rich/poor conditions)
    gamma_at_muM_zero = []
    for delta_mu_O in delta_mu_O_range:
        gamma = phi - gamma_O * float(delta_mu_O)
        gamma_at_muM_zero.append(float(gamma))
    
    # Compute γ with Δμ_O = 0 (1D slice for metal-rich/poor conditions)
    gamma_at_muO_zero = []
    for delta_mu_M in delta_mu_M_range:
        gamma = phi - gamma_M * float(delta_mu_M)
        gamma_at_muO_zero.append(float(gamma))
    
    # ========================================================================
    # ALTERNATIVE FORMULATION: B as independent variable (N becomes reference)
    # ========================================================================
    # In the standard formulation above, we have:
    #   - M (element A) as independent variable
    #   - N (element B) as reference (eliminated via bulk equilibrium)
    #
    # Now we compute the alternative formulation with:
    #   - N (element B) as independent variable
    #   - M (element A) as reference (eliminated via bulk equilibrium)
    #
    # Following the derivation: γ(Δμ_B, Δμ_O) = φ_B - Γ_B·Δμ_B - Γ_O·Δμ_O
    
    # Calculate stoichiometric reference for B-based formulation
    # Number of N (B) atoms needed to maintain bulk stoichiometry with N_M (A) atoms:
    # N_M / x_M = N_N_stoich / y_N  =>  N_N_stoich = N_M * y_N / x_M
    N_N_stoich = (N_M_slab * y_N) / x_M
    N_O_stoich_B = (N_M_slab * z_O) / x_M
    
    # Surface excesses for B-based formulation (per unit area)
    gamma_N = (N_N_stoich - N_N_slab) / (2 * area)  # Γ_B
    gamma_O_B = (N_O_stoich_B - N_O_slab) / (2 * area)  # Γ_O (with respect to A as reference)
    
    # Reference surface energy for B-based formulation (at Δμ_B=0, Δμ_O=0)
    # φ_B = [E_slab - N_M·E_A - N_N·E_B - N_O·(E_O2/2) - (N_M/x_M)·ΔH_f] / (2A)
    phi_B = (
        slab_energy.value 
        - N_M_slab * ref_energies[element_M]
        - N_N_slab * ref_energies[element_N_ref]
        - N_O_slab * ref_energies['O']
        - (N_M_slab / x_M) * delta_h
    ) / (2 * area)
    
    # Chemical potential ranges for B (element N)
    # Maximum: Δμ_N = 0 (N-rich limit, equilibrium with element N reference)
    delta_mu_N_max = 0.0

    # Minimum: from bulk stability constraint
    # At most restrictive condition (Δμ_M = 0, Δμ_O = 0):
    # y_N·Δμ_N >= ΔH_f  =>  Δμ_N >= ΔH_f/y_N
    # Since ΔH_f is negative for stable compounds, this gives a negative minimum
    delta_mu_N_min = delta_h / y_N

    # Range for N: from N-poor (ΔH_f/y_N) to N-rich (0)
    delta_mu_N_range = np.linspace(delta_mu_N_min, delta_mu_N_max, grid_points)
    
    # Compute 2D surface energy grid: γ(Δμ_N, Δμ_O)
    gamma_grid_2d_B = []
    for delta_mu_N in delta_mu_N_range:
        gamma_row_B = []
        for delta_mu_O in delta_mu_O_range:
            gamma_B = phi_B - gamma_N * float(delta_mu_N) - gamma_O_B * float(delta_mu_O)
            gamma_row_B.append(float(gamma_B))
        gamma_grid_2d_B.append(gamma_row_B)
    
    # Compute γ with Δμ_N = 0 (1D slice for N-rich conditions)
    gamma_at_muN_zero = []
    for delta_mu_O in delta_mu_O_range:
        gamma_B = phi_B - gamma_O_B * float(delta_mu_O)
        gamma_at_muN_zero.append(float(gamma_B))
    
    # Compute γ with Δμ_O = 0 (1D slice for O-rich conditions)
    gamma_at_muO_zero_B = []
    for delta_mu_N in delta_mu_N_range:
        gamma_B = phi_B - gamma_N * float(delta_mu_N)
        gamma_at_muO_zero_B.append(float(gamma_B))
    
    return orm.Dict(
        dict={
            # ===== Primary formulation (A-based, for consistency with binary) =====
            # Use 'primary' key for the default formulation
            'primary': {
                'phi': float(phi),
                'Gamma_M': float(gamma_M),
                'Gamma_O': float(gamma_O),
                'delta_mu_M_range': [float(x) for x in delta_mu_M_range],
                'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
                'gamma_grid': gamma_grid_2d,
                'gamma_at_muM_zero': gamma_at_muM_zero,
                'gamma_at_muO_zero': gamma_at_muO_zero,
                'gamma_at_reference': float(phi),
                'element_M_independent': element_M,
                'element_N_reference': element_N_ref,
            },

            # ===== A-based formulation (original) =====
            'A_based': {
                'phi': float(phi),
                'Gamma_A': float(gamma_M),
                'Gamma_O': float(gamma_O),
                'delta_mu_A_range': [float(x) for x in delta_mu_M_range],
                'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
                'gamma_grid': gamma_grid_2d,
                'gamma_at_muA_zero': gamma_at_muM_zero,
                'gamma_at_muO_zero': gamma_at_muO_zero,
                'gamma_at_reference': float(phi),
                'element_A_independent': element_M,
                'element_B_reference': element_N_ref,
            },

            # ===== B-based formulation (alternative) =====
            'B_based': {
                'phi': float(phi_B),
                'Gamma_B': float(gamma_N),
                'Gamma_O': float(gamma_O_B),
                'delta_mu_B_range': [float(x) for x in delta_mu_N_range],
                'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
                'gamma_grid': gamma_grid_2d_B,
                'gamma_at_muB_zero': gamma_at_muN_zero,
                'gamma_at_muO_zero': gamma_at_muO_zero_B,
                'gamma_at_reference': float(phi_B),
                'element_B_independent': element_N_ref,
                'element_A_reference': element_M,
            },

            # ===== Common system information =====
            'oxide_type': 'ternary',
            'area_A2': float(area),
            'bulk_stoichiometry': {
                f'x_{element_M}': int(x_M),
                f'y_{element_N_ref}': int(y_N),
                'z_O': int(z_O),
            },
            'slab_atom_counts': {
                f'N_{element_M}': int(N_M_slab),
                f'N_{element_N_ref}': int(N_N_slab),
                'N_O': int(N_O_slab),
            },
            'reference_energies_per_atom': {k: float(v) for k, v in ref_energies.items()},
            'E_slab_eV': float(slab_energy.value),
            'E_bulk_per_fu_eV': float(bulk_energy_per_fu),
            'formation_enthalpy_eV': float(delta_h),

            # ===== Legacy keys for backward compatibility =====
            'bulk_stoichiometry_AxByOz': {
                f'x_{element_M}': int(x_M),
                f'y_{element_N_ref}': int(y_N),
                f'z_O': int(z_O),
            },
            'E_bulk_fu_eV': float(bulk_energy_per_fu),
        }
    )


@task.calcfunction
def calculate_surface_energy_binary(
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    slab_structure: orm.StructureData,
    slab_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Dict,
    sampling: orm.Int,
) -> orm.Dict:
    """
    Compute γ(Δμ_O) surface energy for a single binary oxide slab.

    Following Reuter & Scheffler convention (PRB 65, 035406):
    - Δμ_O = μ_O - (1/2)E_O2 (deviation from O2 reference)
    - O-rich limit: Δμ_O = 0 (equilibrium with O2 gas)
    - O-poor limit: Δμ_O = ΔH_f / y (oxide decomposition into metal)

    For a binary oxide M_x O_y, the surface energy is calculated as:
    γ(Δμ_O) = φ - Γ_O·Δμ_O

    where:
    - φ: reference surface energy at Δμ_O = 0 (O-rich limit)
    - Γ_O: surface oxygen excess per unit area
    - Δμ_O: oxygen chemical potential deviation from O2 reference

    Args:
        bulk_structure: Bulk structure containing M, O
        bulk_energy: Total energy of bulk structure (eV)
        slab_structure: Slab structure
        slab_energy: Total energy of slab (eV)
        reference_energies: Dict with 'metal_energy_per_atom', 'oxygen_energy_per_atom'
        formation_enthalpy: Dict with 'formation_enthalpy_ev'
        sampling: Number of grid points for Δμ_O sampling

    Returns:
        Dictionary containing surface energy data
    """
    grid_points = sampling.value
    if grid_points <= 0:
        raise ValueError('Sampling must be a positive integer')
    
    # Extract bulk composition
    bulk_ase = bulk_structure.get_ase()
    bulk_counts = Counter(bulk_ase.get_chemical_symbols())
    
    if 'O' not in bulk_counts:
        raise ValueError('The bulk structure contains no oxygen; expected a binary oxide.')
    
    metal_elements = [element for element in bulk_counts if element != 'O']
    if len(metal_elements) != 1:
        raise ValueError(
            f'Expected exactly one metal species; found: {metal_elements}'
        )
    
    element_M = metal_elements[0]
    element_O = 'O'
    
    # Get stoichiometry (x and y in M_x O_y)
    x = bulk_counts[element_M]
    y = bulk_counts[element_O]
    
    # Get reduced stoichiometry
    common_divisor = gcd(x, y)
    x_reduced = x // common_divisor
    y_reduced = y // common_divisor
    
    # Extract reference energies and formation enthalpy
    ref_data = reference_energies.get_dict()
    formation_data = formation_enthalpy.get_dict()
    
    E_M_ref = float(ref_data['metal_energy_per_atom'])
    E_O_ref = float(ref_data['oxygen_energy_per_atom'])
    delta_h = float(formation_data['formation_enthalpy_ev'])
    
    # Calculate surface area
    slab_ase = slab_structure.get_ase()
    cell = slab_ase.get_cell()
    a_vec = cell[0]
    b_vec = cell[1]
    cross = np.cross(a_vec, b_vec)
    area = float(np.linalg.norm(cross))
    
    # Slab atom counts
    slab_counts = Counter(slab_ase.get_chemical_symbols())
    N_M_slab = slab_counts.get(element_M, 0)
    N_O_slab = slab_counts.get(element_O, 0)
    
    # Stoichiometric imbalance: Δ_O = expected_O - actual_O
    # Expected oxygen based on metal count and bulk stoichiometry
    expected_O = (y / x) * N_M_slab
    stoichiometric_imbalance = expected_O - N_O_slab

    # ========== CORRECTED: Δμ_O bounds (Reuter & Scheffler convention) ==========
    # Δμ_O = μ_O - (1/2)E_O2, referenced to O2 molecule at T=0K
    #
    # O-poor limit: oxide decomposes into metal + O2
    #   Δμ_O_min = ΔH_f / y (negative for stable oxide)
    delta_mu_O_min = delta_h / y_reduced
    #
    # O-rich limit: equilibrium with O2 gas (reference state)
    #   Δμ_O_max = 0
    delta_mu_O_max = 0.0

    # Generate chemical potential range (O-poor to O-rich)
    delta_mu_O_range = np.linspace(delta_mu_O_min, delta_mu_O_max, grid_points)

    # ========== CORRECTED: Reference surface energy φ at Δμ_O = 0 (O-rich) ==========
    # φ = (1/2A) × [E_slab - N_M×(E_bulk/x) - N_O×E_O_ref]
    # where E_O_ref = (1/2)E_O2 is the oxygen reference energy per atom
    phi = (
        slab_energy.value
        - N_M_slab * (bulk_energy.value / x)
        - N_O_slab * E_O_ref
    ) / (2 * area)

    # ========== CORRECTED: Surface oxygen excess (Reuter & Scheffler convention) ==========
    # Γ_O = (N_O - (y/x)×N_M) / (2A) = -stoichiometric_imbalance / (2A)
    # Positive Γ_O means O-rich surface, negative means O-poor surface
    Gamma_O = -stoichiometric_imbalance / (2 * area)

    # ========== CORRECTED: Compute γ(Δμ_O) ==========
    # γ(Δμ_O) = φ - Γ_O × Δμ_O
    gamma_array = []
    for delta_mu_O in delta_mu_O_range:
        gamma = phi - Gamma_O * float(delta_mu_O)
        gamma_array.append(float(gamma))
    
    # Special values
    gamma_O_poor = gamma_array[0]   # At Δμ_O_min (O-poor/decomposition limit)
    gamma_O_rich = gamma_array[-1]  # At Δμ_O_max = 0 (O-rich/O2 equilibrium)
    
    # Calculate bulk energy per formula unit for consistency with ternary
    formula_units_in_bulk = bulk_counts[element_M] / x_reduced
    bulk_energy_per_fu = bulk_energy.value / formula_units_in_bulk

    return orm.Dict(
        dict={
            # ===== Primary calculation (M as reference) =====
            'primary': {
                'phi': float(phi),
                'Gamma_O': float(Gamma_O),
                'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
                'gamma_array': gamma_array,
                'gamma_O_poor': float(gamma_O_poor),
                'gamma_O_rich': float(gamma_O_rich),
                'gamma_at_reference': float(phi),
                'element_M': element_M,
            },

            # ===== Common system information =====
            'oxide_type': 'binary',
            'area_A2': float(area),
            'bulk_stoichiometry': {
                f'x_{element_M}': int(x_reduced),
                'y_O': int(y_reduced),
            },
            'slab_atom_counts': {
                f'N_{element_M}': int(N_M_slab),
                'N_O': int(N_O_slab),
            },
            'reference_energies_per_atom': {
                element_M: float(E_M_ref),
                'O': float(E_O_ref),
            },
            'E_slab_eV': float(slab_energy.value),
            'E_bulk_per_fu_eV': float(bulk_energy_per_fu),
            'formation_enthalpy_eV': float(delta_h),

            # ===== Binary-specific data =====
            'stoichiometric_imbalance': float(stoichiometric_imbalance),
            'delta_mu_O_min': float(delta_mu_O_min),  # O-poor limit (= ΔH_f / y)
            'delta_mu_O_max': float(delta_mu_O_max),  # O-rich limit (= 0)

            # ===== Legacy keys for backward compatibility =====
            'phi': float(phi),
            'Gamma_O': float(Gamma_O),
            'delta_mu_O_range': [float(x) for x in delta_mu_O_range],
            'gamma_array': gamma_array,
            'gamma_O_poor': float(gamma_O_poor),
            'gamma_O_rich': float(gamma_O_rich),
            'gamma_at_reference': float(phi),
            'element_M': element_M,
            'E_bulk_eV': float(bulk_energy.value),
            'bulk_stoichiometry_MxOy': {
                f'x_{element_M}': int(x_reduced),
                f'y_O': int(y_reduced),
            },
        }
    )


@task.graph
def compute_surface_energies_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)],
    bulk_structure: orm.StructureData,
    bulk_energy: orm.Float,
    reference_energies: orm.Dict,
    formation_enthalpy: orm.Dict,
    oxide_type: orm.Str,
    sampling: int = 100,
) -> t.Annotated[dict, namespace(surface_energies=dynamic(orm.Dict))]:
    """
    Scatter-gather pattern for computing surface energies.
    
    Automatically selects the appropriate calculation method (binary or ternary)
    based on oxide_type and computes surface energy for each slab in parallel.
    
    Args:
        slabs: Dynamic namespace of slab structures
        energies: Dynamic namespace of slab energies
        bulk_structure: Bulk structure
        bulk_energy: Bulk total energy
        reference_energies: Dict from formation_enthalpy calculation
        formation_enthalpy: Dict from formation_enthalpy calculation
        oxide_type: Type of oxide ('binary' or 'ternary') as Str node
        sampling: Grid resolution for chemical potential
        
    Returns:
        Dictionary with 'surface_energies' namespace containing
        surface energy data for each slab
    """
    surface_results = {}
    sampling_node = orm.Int(sampling)
    
    # Extract string value from Str node
    oxide_type_str = oxide_type.value
    
    # Select the appropriate calculation function based on oxide type
    if oxide_type_str == 'ternary':
        calc_func = calculate_surface_energy_ternary
    elif oxide_type_str == 'binary':
        calc_func = calculate_surface_energy_binary
    else:
        raise ValueError(f'Unknown oxide_type: {oxide_type_str}. Must be "binary" or "ternary".')
    
    # Scatter: compute surface energy for each slab in parallel
    for key, slab_structure in slabs.items():
        slab_energy = energies[key]
        
        surface_data = calc_func(
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
