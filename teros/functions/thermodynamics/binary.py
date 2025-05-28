from aiida.orm import Dict
from aiida_workgraph import task
from collections import Counter
import numpy as np


@task.calcfunction()
def calculate_surface_energy_binary(bulk_structure, bulk_parameters, formation_enthalpy=None, 
                                   code=None, **kwargs):
    """
    Calculate surface Gibbs free energies (γ) for slabs of a binary oxide material.

    ================================================================================================
    THEORETICAL FORMULATION
    ================================================================================================
    
    For a binary oxide M_x O_y (where M is a metal), the surface energy is calculated as:
    
    1. Oxygen-poor limit (μ_O = μ_O^min):
       ----------------------------------------
       γ_O-poor = (1/2A) × [E_slab - N_M,slab × (E_bulk/N_M,bulk) + Δ_O × μ_O^min]
       
       where:
       - Δ_O = (y/x) × N_M,slab - N_O,slab  (oxygen deficiency/excess)
       - μ_O^min = (E_bulk - x × E_M,ref) / y  (decomposition limit)
    
    2. Oxygen-rich limit:
       ----------------------------------------
       γ_O-rich = γ_O-poor - (1/2A) × Δ_O × (ΔH_f / y)
       
       where:
       - ΔH_f is the formation enthalpy of M_x O_y
    
    Key Terms:
    ----------
    - E_slab: Total energy of the slab
    - E_bulk: Total energy of the bulk unit cell
    - N_i,slab: Number of atoms of element i in the slab
    - N_M,bulk: Number of metal atoms in the bulk unit cell (= x × formula units)
    - x, y: Stoichiometric coefficients in M_x O_y
    - A: Surface area of one side of the slab (factor of 2 accounts for both surfaces)
    - E_M,ref: Reference energy per metal atom (from pure metal)
    - Δ_O: Stoichiometric imbalance (positive means O deficient, negative means O rich)
    
    Chemical Potential Bounds:
    --------------------------
    The oxygen chemical potential range is:
    - μ_O^min: Set by decomposition to metal + O2
    - μ_O^max: min(0, μ_O^min + ΔH_f/y)
    
    Physical Interpretation:
    ------------------------
    - γ_O-poor: Surface energy when oxygen leaves the surface (reducing conditions)
    - γ_O-rich: Surface energy in equilibrium with O2 gas (oxidizing conditions)
    - The term Δ_O × (ΔH_f/y) accounts for the energy change when the surface 
      stoichiometry differs from the bulk
    
    ================================================================================================
    
    :param bulk_structure: AiiDA StructureData of the bulk binary oxide
    :param bulk_parameters: AiiDA Dict with bulk DFT calculation results
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
    
    # ============================================================================================
    # SECTION 2: ANALYZE BULK COMPOSITION AND STOICHIOMETRY
    # ============================================================================================
    
    # Get bulk composition
    bulk_atoms = bulk_structure.get_ase()
    elements_count = Counter(bulk_atoms.get_chemical_symbols())
    elements = list(elements_count.keys())
    
    # Validate binary oxide
    if len(elements) != 2:
        raise ValueError(f"Found {len(elements)} elements: {elements}. This workflow requires a binary compound.")
    
    if 'O' not in elements:
        raise ValueError(f"No oxygen atoms found in structure: {elements}. This workflow is designed for oxide materials.")
    
    # Identify metal and oxygen
    metal = [element for element in elements if element != 'O'][0]
    
    # Get stoichiometric coefficients (x and y in M_x O_y)
    x = elements_count[metal]  # Total metal atoms in bulk cell
    y = elements_count['O']    # Total oxygen atoms in bulk cell
    
    print(f"\n{'='*80}")
    print(f"BULK COMPOSITION: {metal}_{x} O_{y}")
    print(f"{'='*80}")
    
    # ============================================================================================
    # SECTION 3: EXTRACT ENERGIES FROM DFT CALCULATIONS
    # ============================================================================================
    
    # Extract bulk energy
    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        E_bulk = bulk_parameters.get_dict()['energy']
    else:  # VASP
        E_bulk = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']
    
    # Extract reference energies and formation data
    formation_dict = formation_enthalpy.get_dict()
    element_energy_per_atom = {
        element: formation_dict.get(f'{element.lower()}_energy_per_atom', 0.0)
        for element in elements
    }
    
    E_metal_ref = element_energy_per_atom[metal]
    E_O_ref = element_energy_per_atom['O']  # Per O atom (typically E_O2/2)
    E_O2 = E_O_ref * 2  # O2 molecule energy
    
    # Formation energy is given per atom in the dict, need total for formula unit
    formation_energy = formation_dict.get('formation_enthalpy_ev', 0.0) * (x + y)
    
    print(f"\nENERGY REFERENCES:")
    print(f"  E_bulk = {E_bulk:.6f} eV")
    print(f"  E_{metal} (ref) = {E_metal_ref:.6f} eV/atom")
    print(f"  E_O (ref) = {E_O_ref:.6f} eV/atom")
    print(f"  ΔH_f = {formation_energy:.6f} eV (total for {metal}_{x}O_{y})")
    
    # ============================================================================================
    # SECTION 4: CALCULATE OXYGEN CHEMICAL POTENTIAL BOUNDS
    # ============================================================================================
    
    print(f"\nOXYGEN CHEMICAL POTENTIAL BOUNDS:")
    print(f"{'='*60}")
    
    # Lower bound: decomposition limit
    # At equilibrium: E_bulk = x × E_metal + y × μ_O
    mu_O_min = (E_bulk - x * E_metal_ref) / y
    
    print(f"\nLower bound (O-poor, decomposition limit):")
    print(f"  μ_O^min = (E_bulk - x × E_{metal}) / y")
    print(f"         = ({E_bulk:.6f} - {x} × {E_metal_ref:.6f}) / {y}")
    print(f"         = {mu_O_min:.6f} eV")
    
    # Upper bound from formation energy
    mu_O_from_formation = mu_O_min + formation_energy / y
    mu_O_max = min(0.0, mu_O_from_formation)
    
    print(f"\nUpper bound (O-rich limit):")
    print(f"  From formation energy: μ_O ≤ μ_O^min + ΔH_f/y")
    print(f"                        μ_O ≤ {mu_O_min:.6f} + {formation_energy:.6f}/{y}")
    print(f"                        μ_O ≤ {mu_O_from_formation:.6f} eV")
    print(f"  μ_O^max = min(0, {mu_O_from_formation:.6f}) = {mu_O_max:.6f} eV")
    
    # Create range for analysis
    delta_mu_O_range = np.linspace(mu_O_min, mu_O_max, 10)
    
    # ============================================================================================
    # SECTION 5: PROCESS EACH SLAB
    # ============================================================================================
    
    results = {}
    
    for i, ((_, misc), (_, structure)) in enumerate(zip(slab_parameters.items(), slab_structures.items())):
        label = f's_{str(i)}'
        
        print(f"\n{'-'*80}")
        print(f"PROCESSING SLAB: {label}")
        print(f"{'-'*80}")
        
        # ----------------------------------------------------------------------------------------
        # 5.1: Extract slab properties
        # ----------------------------------------------------------------------------------------
        
        slab_atoms = structure.get_ase()
        slab_elements_count = Counter(slab_atoms.get_chemical_symbols())
        
        # Calculate surface area
        area = slab_atoms.get_volume() / slab_atoms.get_cell()[2][2]
        
        # Extract slab energy
        if code in ['QUANTUM_ESPRESSO', 'CP2K']:
            E_slab = misc.get_dict()['energy']
        else:  # VASP
            E_slab = misc.get_dict()['total_energies']['energy_extrapolated']
        
        if E_slab is None:
            continue
        
        # Get atom counts
        N_M_slab = slab_elements_count.get(metal, 0)
        N_O_slab = slab_elements_count.get('O', 0)
        N_M_bulk = x  # Total metal atoms in bulk cell
        
        print(f"\nSlab composition: {metal}_{N_M_slab} O_{N_O_slab}")
        print(f"Surface area: {area:.3f} Å²")
        print(f"E_slab = {E_slab:.6f} eV")
        
        # ----------------------------------------------------------------------------------------
        # 5.2: Calculate stoichiometric imbalance
        # ----------------------------------------------------------------------------------------
        
        print(f"\nSTOICHIOMETRIC ANALYSIS:")
        print(f"{'='*60}")
        
        # Expected oxygen based on metal count and bulk stoichiometry
        expected_O = (y/x) * N_M_slab
        
        # Stoichiometric imbalance (positive = O deficient, negative = O rich)
        stoichiometric_imbalance = expected_O - N_O_slab
        
        print(f"Expected O atoms (from bulk ratio): (y/x) × N_M = ({y}/{x}) × {N_M_slab} = {expected_O:.2f}")
        print(f"Actual O atoms in slab: {N_O_slab}")
        print(f"Stoichiometric imbalance Δ_O = {expected_O:.2f} - {N_O_slab} = {stoichiometric_imbalance:.2f}")
        if stoichiometric_imbalance > 0:
            print(f"  → Slab is OXYGEN DEFICIENT by {stoichiometric_imbalance:.2f} atoms")
        elif stoichiometric_imbalance < 0:
            print(f"  → Slab is OXYGEN RICH by {-stoichiometric_imbalance:.2f} atoms")
        else:
            print(f"  → Slab is STOICHIOMETRIC")
        
        # ----------------------------------------------------------------------------------------
        # 5.3: Calculate surface energy at O-poor limit
        # ----------------------------------------------------------------------------------------
        
        print(f"\nO-POOR LIMIT CALCULATION:")
        print(f"{'='*60}")
        
        # Term 1: Slab energy
        term1 = E_slab
        print(f"Term 1 (slab energy): {term1:.6f} eV")
        
        # Term 2: Bulk reference energy
        term2 = N_M_slab * (E_bulk / N_M_bulk)
        print(f"Term 2 (bulk reference): N_M × (E_bulk/N_M_bulk)")
        print(f"                       = {N_M_slab} × ({E_bulk:.6f}/{N_M_bulk})")
        print(f"                       = {term2:.6f} eV")
        
        # Term 3: Chemical potential contribution
        term3 = stoichiometric_imbalance * mu_O_min
        print(f"Term 3 (chem potential): Δ_O × μ_O^min")
        print(f"                       = {stoichiometric_imbalance:.2f} × {mu_O_min:.6f}")
        print(f"                       = {term3:.6f} eV")
        
        # Total surface energy at O-poor limit
        gamma_O_poor = (term1 - term2 + term3) / (2 * area)
        
        print(f"\nγ_O-poor = (Term1 - Term2 + Term3) / (2 × A)")
        print(f"        = ({term1:.6f} - {term2:.6f} + {term3:.6f}) / (2 × {area:.3f})")
        print(f"        = {gamma_O_poor:.6f} eV/Å²")
        
        # ----------------------------------------------------------------------------------------
        # 5.4: Calculate surface energy at O-rich limit
        # ----------------------------------------------------------------------------------------
        
        print(f"\nO-RICH LIMIT CALCULATION:")
        print(f"{'='*60}")
        
        # Correction term for O-rich conditions
        correction = (stoichiometric_imbalance * formation_energy / y) / (2 * area)
        
        print(f"Correction = Δ_O × (ΔH_f/y) / (2 × A)")
        print(f"          = {stoichiometric_imbalance:.2f} × ({formation_energy:.6f}/{y}) / (2 × {area:.3f})")
        print(f"          = {correction:.6f} eV/Å²")
        
        gamma_O_rich = gamma_O_poor - correction
        
        print(f"\nγ_O-rich = γ_O-poor - correction")
        print(f"        = {gamma_O_poor:.6f} - {correction:.6f}")
        print(f"        = {gamma_O_rich:.6f} eV/Å²")
        
        # ----------------------------------------------------------------------------------------
        # 5.5: Summary and physical interpretation
        # ----------------------------------------------------------------------------------------
        
        print(f"\nSURFACE ENERGY SUMMARY:")
        print(f"{'='*60}")
        print(f"γ_O-poor = {gamma_O_poor:.6f} eV/Å² (reducing conditions)")
        print(f"γ_O-rich = {gamma_O_rich:.6f} eV/Å² (oxidizing conditions)")
        print(f"Difference = {gamma_O_rich - gamma_O_poor:.6f} eV/Å²")
        
        if gamma_O_poor < gamma_O_rich:
            print(f"→ This surface is more stable under O-poor conditions")
        else:
            print(f"→ This surface is more stable under O-rich conditions")
        
        # ----------------------------------------------------------------------------------------
        # 5.6: Store results
        # ----------------------------------------------------------------------------------------
        
        results_dict = {
            'slab_energy': float(E_slab),
            'bulk_energy': float(E_bulk),
            'area': float(area),
            'mu_O_min': float(mu_O_min),
            'mu_O_max': float(mu_O_max),
            'gamma_O_poor': float(gamma_O_poor),
            'gamma_O_rich': float(gamma_O_rich),
            'stoichiometric_imbalance': float(stoichiometric_imbalance),
            'element_counts': {element: int(count) for element, count in slab_elements_count.items()},
            'reference_energies': {element: float(energy) for element, energy in element_energy_per_atom.items()},
            'formation_enthalpy_ev': float(formation_energy),
        }
        
        results[label] = Dict(dict=results_dict)
    
    print(f"\n{'='*80}")
    print(f"CALCULATION COMPLETE")
    print(f"{'='*80}\n")
    
    return results