from ase.io import read
from aiida.orm import Dict, load_node, Int
from aiida_workgraph import task
from collections import Counter
from math import gcd
from functools import reduce
import numpy as np


@task.calcfunction()
def calculate_surface_energy_ternar(bulk_structure, bulk_parameters, sampling=None, 
                                  formation_enthalpy=None, code=None, **kwargs):
    """
    Calculate surface Gibbs free energies (γ) for slabs of a ternary oxide material.

    ================================================================================================
    THEORETICAL FORMULATION
    ================================================================================================
    
    For a ternary oxide M_x N_y O_z (where M and N are metals), the surface energy is:
    
        γ(Δμ_M, Δμ_O) = φ - Γ_M × Δμ_M - Γ_O × Δμ_O
    
    Where:
    ------
    1. φ (surface energy at reference chemical potentials):
       φ = (1/2A) × [E_slab - (N_N,slab/y) × E_bulk,fu] - Γ_M × E_M,ref - Γ_O × E_O,ref
    
    2. Γ_M (surface excess of metal M relative to metal N):
       Γ_M = (1/2A) × (N_M,slab - x × N_N,slab/y)
    
    3. Γ_O (surface excess of oxygen relative to metal N):
       Γ_O = (1/2A) × (N_O,slab - z × N_N,slab/y)
    
    Variables:
    ----------
    - E_slab: Total energy of the slab
    - E_bulk,fu: Energy per formula unit of the bulk
    - N_i,slab: Number of atoms of element i in the slab
    - x, y, z: Stoichiometric coefficients in M_x N_y O_z
    - A: Surface area of one side of the slab (factor of 2 accounts for both surfaces)
    - E_i,ref: Reference energy per atom for element i
    - Δμ_i: Chemical potential of element i relative to its reference state
    
    Stability constraint (from bulk formation):
    ------------------------------------------
    x × Δμ_M + z × Δμ_O ≤ ΔH_f(M_x N_y O_z)
    
    ================================================================================================
    
    :param bulk_structure: AiiDA StructureData of the bulk ternary material
    :param bulk_parameters: AiiDA Dict with bulk DFT calculation results
    :param sampling: Number of points for chemical potential grid (default: 100)
    :param formation_enthalpy: AiiDA Dict with formation enthalpy data
    :param code: DFT code identifier ("QUANTUM_ESPRESSO", "CP2K", "VASP")
    :param kwargs: Should contain slab_structures and slab_parameters dicts
    
    :return: Dictionary mapping slab IDs to surface energy results
    """
    
    # ============================================================================================
    # SECTION 1: EXTRACT INPUT DATA
    # ============================================================================================
    
    slab_structures = kwargs.get('slab_structures', {})
    slab_parameters = kwargs.get('slab_parameters', {})
    
    # Set default sampling if not provided
    if sampling is None:
        sampling = Int(100)
    grid_points = sampling.value
    
    # ============================================================================================
    # SECTION 2: ANALYZE BULK COMPOSITION AND STOICHIOMETRY
    # ============================================================================================
    
    # Get bulk composition
    bulk_ase = bulk_structure.get_ase()
    bulk_atom_counts = Counter(bulk_ase.get_chemical_symbols())
    
    # Identify elements (expecting M, N, O for ternary oxide)
    all_elements = sorted(list(bulk_atom_counts.keys()))
    if 'O' not in all_elements:
        raise ValueError("No oxygen found. This workflow is for oxide materials.")
    
    # Separate metals from oxygen
    metal_elements = sorted([el for el in all_elements if el != 'O'])
    if len(metal_elements) != 2:
        raise ValueError(f"Expected 2 metal elements for ternary oxide, found {len(metal_elements)}: {metal_elements}")
    
    # Define element roles
    element_M = metal_elements[0]      # First metal (independent chemical potential)
    element_N_ref = metal_elements[1]  # Second metal (reference for surface excess)
    element_O = 'O'                    # Oxygen
    
    # Calculate stoichiometric coefficients (x, y, z in M_x N_y O_z)
    element_list = [element_M, element_N_ref, element_O]
    counts_for_gcd = [bulk_atom_counts[el] for el in element_list]
    common_divisor = reduce(gcd, counts_for_gcd)
    
    x_M = bulk_atom_counts[element_M] // common_divisor
    y_N = bulk_atom_counts[element_N_ref] // common_divisor
    z_O = bulk_atom_counts[element_O] // common_divisor
    
    print(f"\n{'='*80}")
    print(f"BULK COMPOSITION: {element_M}_{x_M} {element_N_ref}_{y_N} {element_O}_{z_O}")
    print(f"{'='*80}")
    
    # ============================================================================================
    # SECTION 3: EXTRACT ENERGIES FROM DFT CALCULATIONS
    # ============================================================================================
    
    # Extract bulk energy
    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        E_bulk_total = bulk_parameters.get_dict()['energy']
    else:  # VASP
        E_bulk_total = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']
    
    # Calculate energy per formula unit
    formula_units_in_bulk = bulk_atom_counts[element_N_ref] / y_N
    E_bulk_per_fu = E_bulk_total / formula_units_in_bulk
    
    # Extract reference energies
    formation_data = formation_enthalpy.get_dict() if formation_enthalpy else {}
    reference_energies = {}
    for element in element_list:
        energy_key = f'{element.lower()}_energy_per_atom'
        if energy_key not in formation_data:
            raise ValueError(f"Reference energy for {element} not found")
        reference_energies[element] = formation_data[energy_key]
    
    E_M_ref = reference_energies[element_M]
    E_O_ref = reference_energies[element_O]
    
    # Get formation enthalpy if available
    formation_enthalpy_value = formation_data.get('formation_enthalpy_ev')
    
    print(f"\nENERGY REFERENCES:")
    print(f"  E_bulk (per f.u.) = {E_bulk_per_fu:.6f} eV")
    print(f"  E_{element_M} (ref) = {E_M_ref:.6f} eV/atom")
    print(f"  E_{element_O} (ref) = {E_O_ref:.6f} eV/atom")
    if formation_enthalpy_value:
        print(f"  ΔH_f = {formation_enthalpy_value:.6f} eV/f.u.")
    
    # ============================================================================================
    # SECTION 4: DEFINE CHEMICAL POTENTIAL RANGES
    # ============================================================================================
    
    # Default bounds
    delta_mu_M_min = -5.0  # eV
    delta_mu_O_min = -5.0  # eV
    
    # Apply stability constraints if formation enthalpy is available
    if formation_enthalpy_value is not None:
        # From stability: x·Δμ_M + z·Δμ_O ≤ ΔH_f
        # When Δμ_O = 0: Δμ_M ≤ ΔH_f/x
        # When Δμ_M = 0: Δμ_O ≤ ΔH_f/z
        
        if x_M > 0:
            delta_mu_M_min = min(0.0, formation_enthalpy_value / x_M)
        if z_O > 0:
            delta_mu_O_min = min(0.0, formation_enthalpy_value / z_O)
    
    # Ensure bounds don't exceed 0
    delta_mu_M_min = min(delta_mu_M_min, 0.0)
    delta_mu_O_min = min(delta_mu_O_min, 0.0)
    
    # Create chemical potential grids
    delta_mu_M_range = np.linspace(delta_mu_M_min, 0.0, grid_points)
    delta_mu_O_range = np.linspace(delta_mu_O_min, 0.0, grid_points)
    
    print(f"\nCHEMICAL POTENTIAL RANGES:")
    print(f"  Δμ_{element_M}: [{delta_mu_M_min:.3f}, 0.000] eV")
    print(f"  Δμ_{element_O}: [{delta_mu_O_min:.3f}, 0.000] eV")
    
    # ============================================================================================
    # SECTION 5: PROCESS EACH SLAB
    # ============================================================================================
    
    results = {}
    
    for i, ((_, slab_data), (_, structure)) in enumerate(zip(slab_parameters.items(), slab_structures.items())):
        slab_id = f's_{i}'
        
        print(f"\n{'-'*80}")
        print(f"PROCESSING SLAB: {slab_id}")
        print(f"{'-'*80}")
        
        # ----------------------------------------------------------------------------------------
        # 5.1: Extract slab properties
        # ----------------------------------------------------------------------------------------
        
        slab_ase = structure.get_ase()
        slab_atom_counts = Counter(slab_ase.get_chemical_symbols())
        
        # Calculate surface area (cross product of lattice vectors a and b)
        area = np.linalg.norm(np.cross(slab_ase.cell[0], slab_ase.cell[1]))
        
        # Extract slab energy
        if code in ['QUANTUM_ESPRESSO', 'CP2K']:
            E_slab = slab_data.get_dict()['energy']
        else:  # VASP
            E_slab = slab_data.get_dict()['total_energies']['energy_extrapolated']
        
        if E_slab is None:
            print(f"  WARNING: No energy found for {slab_id}, skipping...")
            continue
        
        # Get atom counts in slab
        N_M_slab = slab_atom_counts.get(element_M, 0)
        N_N_slab = slab_atom_counts.get(element_N_ref, 0)
        N_O_slab = slab_atom_counts.get(element_O, 0)
        
        print(f"\n  Slab composition: {element_M}_{N_M_slab} {element_N_ref}_{N_N_slab} {element_O}_{N_O_slab}")
        print(f"  Surface area: {area:.3f} Å²")
        print(f"  E_slab = {E_slab:.6f} eV")
        
        # ----------------------------------------------------------------------------------------
        # 5.2: Calculate surface excess (Gamma) terms
        # ----------------------------------------------------------------------------------------
        
        print(f"\n  SURFACE EXCESS CALCULATIONS:")
        print(f"  " + "="*60)
        
        # Γ_M = (1/2A) × (N_M,slab - x × N_N,slab/y)
        excess_M = N_M_slab - (x_M * N_N_slab / y_N)
        Gamma_M = excess_M / (2 * area)
        
        print(f"  Γ_{element_M} = (1/2A) × (N_{element_M},slab - x × N_{element_N_ref},slab/y)")
        print(f"       = (1/{2*area:.3f}) × ({N_M_slab} - {x_M} × {N_N_slab}/{y_N})")
        print(f"       = (1/{2*area:.3f}) × {excess_M:.3f}")
        print(f"       = {Gamma_M:.6f} atoms/Å²")
        
        # Γ_O = (1/2A) × (N_O,slab - z × N_N,slab/y)
        excess_O = N_O_slab - (z_O * N_N_slab / y_N)
        Gamma_O = excess_O / (2 * area)
        
        print(f"\n  Γ_{element_O} = (1/2A) × (N_{element_O},slab - z × N_{element_N_ref},slab/y)")
        print(f"      = (1/{2*area:.3f}) × ({N_O_slab} - {z_O} × {N_N_slab}/{y_N})")
        print(f"      = (1/{2*area:.3f}) × {excess_O:.3f}")
        print(f"      = {Gamma_O:.6f} atoms/Å²")
        
        # ----------------------------------------------------------------------------------------
        # 5.3: Calculate φ (surface energy at reference chemical potentials)
        # ----------------------------------------------------------------------------------------
        
        print(f"\n  PHI (φ) CALCULATION:")
        print(f"  " + "="*60)
        
        # φ = (1/2A) × [E_slab - (N_N,slab/y) × E_bulk,fu] - Γ_M × E_M,ref - Γ_O × E_O,ref
        
        # Term 1: Energy difference contribution
        bulk_equivalent_energy = (N_N_slab / y_N) * E_bulk_per_fu
        energy_diff = E_slab - bulk_equivalent_energy
        term1 = energy_diff / (2 * area)
        
        # Term 2: Metal reference contribution
        term2 = -Gamma_M * E_M_ref
        
        # Term 3: Oxygen reference contribution
        term3 = -Gamma_O * E_O_ref
        
        # Total φ
        phi = term1 + term2 + term3
        
        print(f"  φ = (1/2A) × [E_slab - (N_{element_N_ref},slab/y) × E_bulk,fu] - Γ_{element_M} × E_{element_M},ref - Γ_{element_O} × E_{element_O},ref")
        print(f"\n  Breaking down the calculation:")
        print(f"    Term 1 = (1/{2*area:.3f}) × [{E_slab:.6f} - ({N_N_slab}/{y_N}) × {E_bulk_per_fu:.6f}]")
        print(f"           = (1/{2*area:.3f}) × {energy_diff:.6f}")
        print(f"           = {term1:.6f} eV/Å²")
        print(f"\n    Term 2 = -{Gamma_M:.6f} × {E_M_ref:.6f}")
        print(f"           = {term2:.6f} eV/Å²")
        print(f"\n    Term 3 = -{Gamma_O:.6f} × {E_O_ref:.6f}")
        print(f"           = {term3:.6f} eV/Å²")
        print(f"\n  φ = {term1:.6f} + {term2:.6f} + {term3:.6f} = {phi:.6f} eV/Å²")
        
        # ----------------------------------------------------------------------------------------
        # 5.4: Calculate surface energy over chemical potential grid
        # ----------------------------------------------------------------------------------------
        
        print(f"\n  SURFACE ENERGY EQUATION:")
        print(f"  γ(Δμ_{element_M}, Δμ_{element_O}) = {phi:.6f} - {Gamma_M:.6f} × Δμ_{element_M} - {Gamma_O:.6f} × Δμ_{element_O}")
        
        # Calculate γ for full 2D grid
        gamma_values_grid = {}
        for delta_mu_M in delta_mu_M_range:
            for delta_mu_O in delta_mu_O_range:
                gamma = phi - Gamma_M * delta_mu_M - Gamma_O * delta_mu_O
                key = f"muM_{delta_mu_M:.4f}_muO_{delta_mu_O:.4f}"
                gamma_values_grid[key] = float(gamma)
        
        # Calculate γ for fixed Δμ_M = 0 (common analysis case)
        gamma_values_fixed_muM_zero = {}
        for delta_mu_O in delta_mu_O_range:
            gamma = phi - Gamma_O * delta_mu_O  # Δμ_M = 0
            key = f"muM_0.0000_muO_{delta_mu_O:.4f}"
            gamma_values_fixed_muM_zero[key] = float(gamma)
        
        # Example values
        gamma_at_zero = phi  # Both chemical potentials at zero
        print(f"\n  Example: γ(0, 0) = {gamma_at_zero:.6f} eV/Å²")
        
        # ----------------------------------------------------------------------------------------
        # 5.5: Store results
        # ----------------------------------------------------------------------------------------
        
        results[slab_id] = Dict(dict={
            # Primary results
            'phi': float(phi),
            'Gamma_M_vs_Nref': float(Gamma_M),
            'Gamma_O_vs_Nref': float(Gamma_O),
            'gamma_values_grid': gamma_values_grid,
            'gamma_values_fixed_muM_zero': gamma_values_fixed_muM_zero,
            
            # Structural information
            'area_A2': float(area),
            'element_M_independent': element_M,
            'element_N_reference': element_N_ref,
            'bulk_stoichiometry_MxNyOz': {
                f"x_{element_M}": int(x_M),
                f"x_{element_N_ref}": int(y_N),
                f"x_{element_O}": int(z_O)
            },
            'slab_atom_counts': {
                f"N_{element_M}": int(N_M_slab),
                f"N_{element_N_ref}": int(N_N_slab),
                f"N_{element_O}": int(N_O_slab)
            },
            
            # Energy data
            'reference_energies_per_atom': {k: float(v) for k, v in reference_energies.items()},
            'E_slab_eV': float(E_slab),
            'E_bulk_fu_eV': float(E_bulk_per_fu),
        })
        
        if formation_enthalpy_value is not None:
            results[slab_id].dict['formation_enthalpy_MxNyOz_eV'] = float(formation_enthalpy_value)
    
    print(f"\n{'='*80}")
    print(f"CALCULATION COMPLETE")
    print(f"{'='*80}\n")
    
    return results