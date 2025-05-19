from ase.io import read
from aiida.orm import Dict, load_node, Int
from aiida_workgraph import task
from collections import Counter
from math import gcd
from functools import reduce
import numpy as np
# from aiida.orm import Dict # Duplicate import

@task.calcfunction()
def calculate_surface_energy_ternar(bulk_structure, bulk_parameters, sampling=None, formation_enthalpy=None, code=None, **kwargs):
    """
    Calculate the surface Gibbs free energy per unit area (γ) for each slab in a ternary oxide system.
    This version aligns with the formulation where Gamma terms are defined using direct stoichiometric
    multipliers of N_N_ref_slab, and phi's bulk energy subtraction is generalized.

    γ = φ - Γ_NM * Δμ_M - Γ_NO * Δμ_O
    φ = (1/2A) * [E_slab - (N_N_ref_slab / x_N_ref_coeff) * E_bulk_fu]
        - Γ_NM * E_M_bulk_atom
        - Γ_NO * E_O_atom_ref
    Γ_NM = (1/2A) * (N_M_slab - x_M_coeff * N_N_ref_slab)
    Γ_NO = (1/2A) * (N_O_slab - x_O_coeff * N_N_ref_slab)

    where M (el_M) is the first non-oxygen element, N_ref (el_N_ref) is the second non-oxygen element
    (reference for Γ terms) from the bulk formula M_xM_coeff N_ref_xN_ref_coeff O_xO_coeff.
    E_O_atom_ref is typically (1/2 * E_O2_molecule).
    """
    print("--- Starting calculate_surface_energy_ternar (Updated Formulation) ---") # LOG
    slab_structures = kwargs.get('slab_structures', {})
    slab_parameters = kwargs.get('slab_parameters', {})

    # --- Bulk Properties ---
    bulk_ase = bulk_structure.get_ase()
    bulk_elements_counter_raw = Counter(bulk_ase.get_chemical_symbols())
    print(f"[LOG] Raw bulk elements count: {bulk_elements_counter_raw}") # LOG

    unique_elements_bulk_formula = sorted(list(bulk_elements_counter_raw.keys()))
    el_O = 'O'
    if el_O not in unique_elements_bulk_formula:
        raise ValueError("No oxygen atoms found in bulk structure. This workflow is designed for oxide materials.")

    metal_elements = sorted([el for el in unique_elements_bulk_formula if el != el_O])
    if len(metal_elements) != 2:
        raise ValueError(
            f"Expected 2 metal elements for a ternary oxide, found {len(metal_elements)}: {metal_elements}. "
            f"Elements in bulk: {unique_elements_bulk_formula}"
        )

    el_M = metal_elements[0]
    el_N_ref = metal_elements[1]
    print(f"[LOG] Identified elements: el_M = {el_M}, el_N_ref (reference) = {el_N_ref}, el_O = {el_O}") # LOG

    counts_for_gcd = [bulk_elements_counter_raw[el] for el in [el_M, el_N_ref, el_O] if el in bulk_elements_counter_raw]
    if not counts_for_gcd or len(counts_for_gcd) < 3: # Ensure all three elements are present for GCD
        raise ValueError(f"Bulk structure does not contain all expected elements M, N_ref, O for GCD calculation. Counts: {bulk_elements_counter_raw}")
    common_divisor = reduce(gcd, counts_for_gcd)
    print(f"[LOG] GCD for bulk formula: {common_divisor}") # LOG

    bulk_stoich = {
        el_M: bulk_elements_counter_raw[el_M] // common_divisor,
        el_N_ref: bulk_elements_counter_raw[el_N_ref] // common_divisor,
        el_O: bulk_elements_counter_raw[el_O] // common_divisor,
    }
    x_M_coeff = bulk_stoich[el_M]
    x_N_ref_coeff = bulk_stoich[el_N_ref]
    x_O_coeff = bulk_stoich[el_O]
    print(f"[LOG] Bulk stoichiometric coefficients: x_{el_M}={x_M_coeff}, x_{el_N_ref}={x_N_ref_coeff}, x_{el_O}={x_O_coeff}") # LOG

    if x_N_ref_coeff == 0: # Denominator for bulk energy term and potentially old Gamma
        raise ValueError(f"Stoichiometric coefficient for reference element {el_N_ref} (x_N_ref_coeff) is zero.")

    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        E_bulk_cell = bulk_parameters.get_dict()['energy']
    else: # VASP and others
        E_bulk_cell = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']
    print(f"[LOG] Total energy of bulk simulation cell (E_bulk_cell): {E_bulk_cell} eV") # LOG

    if x_N_ref_coeff > 0 :
        num_fu_in_bulk_cell = bulk_elements_counter_raw[el_N_ref] / x_N_ref_coeff
    elif x_M_coeff > 0:
         num_fu_in_bulk_cell = bulk_elements_counter_raw[el_M] / x_M_coeff
    elif x_O_coeff > 0:
        num_fu_in_bulk_cell = bulk_elements_counter_raw[el_O] / x_O_coeff
    else:
        raise ValueError("Cannot determine number of formula units in bulk cell due to zero stoichiometric coefficients.")

    print(f"[LOG] Number of formula units in bulk cell: {num_fu_in_bulk_cell}") # LOG
    E_bulk_fu = E_bulk_cell / num_fu_in_bulk_cell
    print(f"[LOG] Energy per formula unit of bulk (E_bulk_fu): {E_bulk_fu} eV") # LOG

    formation_data_dict = formation_enthalpy.get_dict() if formation_enthalpy else {}
    element_ref_energies = {}
    for el_symbol in [el_M, el_N_ref, el_O]: # N_ref needed if it's part of chem pots, but not for O2 ref.
        energy = formation_data_dict.get(f'{el_symbol.lower()}_energy_per_atom')
        if energy is None:
            raise ValueError(f"Reference energy for element {el_symbol} not found in formation_enthalpy node.")
        element_ref_energies[el_symbol] = energy
    print(f"[LOG] Elemental reference energies (per atom): {element_ref_energies}") # LOG

    E_M_bulk_atom = element_ref_energies[el_M]
    E_O_atom_ref = element_ref_energies[el_O] # This is E_O_atom_ref, often 1/2 * E_O2_molecule

    results = {}
    first_slab_logged = False

    for i, ((slab_key_param, slab_misc_data), (slab_key_struct, slab_structure_data)) in enumerate(zip(slab_parameters.items(), slab_structures.items())):
        label = f's_{str(i)}'
        print(f"\n--- Processing Slab: {label} (param key: {slab_key_param}, struct key: {slab_key_struct}) ---") # LOG

        slab_ase = slab_structure_data.get_ase()
        slab_atom_counts = Counter(slab_ase.get_chemical_symbols())
        area = np.linalg.norm(np.cross(slab_ase.cell[0], slab_ase.cell[1]))

        if not first_slab_logged:
            print(f"[LOG][{label}] Slab atoms count: {slab_atom_counts}") # LOG
            print(f"[LOG][{label}] Calculated surface area (A): {area} Å^2 (Using 2*A: {2*area} Å^2)") # LOG

        if code in ['QUANTUM_ESPRESSO', 'CP2K']:
            E_slab = slab_misc_data.get_dict()['energy']
        else: # VASP and others
            E_slab = slab_misc_data.get_dict()['total_energies']['energy_extrapolated']

        if E_slab is None:
            results[label] = Dict(dict={'error': f'E_slab for {label} is None'})
            if not first_slab_logged: print(f"[LOG][{label}] E_slab is None. Skipping.") # LOG
            continue
        if not first_slab_logged: print(f"[LOG][{label}] Slab energy (E_slab): {E_slab} eV") # LOG

        N_M_slab = slab_atom_counts.get(el_M, 0)
        N_N_ref_slab = slab_atom_counts.get(el_N_ref, 0)
        N_O_slab = slab_atom_counts.get(el_O, 0)
        if not first_slab_logged:
            print(f"[LOG][{label}] N_M_slab ({el_M}): {N_M_slab}") # LOG
            print(f"[LOG][{label}] N_N_ref_slab ({el_N_ref}): {N_N_ref_slab}") # LOG
            print(f"[LOG][{label}] N_O_slab ({el_O}): {N_O_slab}") # LOG

        # --- Calculate Gamma terms (surface excess concentrations) ---
        term_in_Gamma_NM = (N_M_slab - x_M_coeff * N_N_ref_slab)
        Gamma_NM = term_in_Gamma_NM / (2 * area)

        term_in_Gamma_NO = (N_O_slab - x_O_coeff * N_N_ref_slab)
        Gamma_NO = term_in_Gamma_NO / (2 * area)

        if not first_slab_logged:
            print(f"[LOG][{label}] Term in Γ_NM (N_M_slab - x_M_coeff * N_N_ref_slab): {term_in_Gamma_NM}") # LOG
            print(f"[LOG][{label}] Γ_NM ({el_M} wrt {el_N_ref}): {Gamma_NM} atoms/Å^2") # LOG
            print(f"[LOG][{label}] Term in Γ_NO (N_O_slab - x_O_coeff * N_N_ref_slab): {term_in_Gamma_NO}") # LOG
            print(f"[LOG][{label}] Γ_NO ({el_O} wrt {el_N_ref}): {Gamma_NO} atoms/Å^2") # LOG

        # --- Calculate phi term ---
        if x_N_ref_coeff == 0: # Should have been caught earlier, but as a safeguard for division
             raise ValueError(f"x_N_ref_coeff for {el_N_ref} is zero, division by zero in phi calculation.")
        equivalent_bulk_energy_for_phi = (N_N_ref_slab / x_N_ref_coeff) * E_bulk_fu
        phi_part1_numerator = E_slab - equivalent_bulk_energy_for_phi
        phi_part1 = phi_part1_numerator / (2 * area)

        phi_part2 = - Gamma_NM * E_M_bulk_atom
        phi_part3 = - Gamma_NO * E_O_atom_ref
        phi_value = phi_part1 + phi_part2 + phi_part3

        if not first_slab_logged:
            print(f"[LOG][{label}] For φ calculation:") # LOG
            print(f"[LOG][{label}]   (N_N_ref_slab / x_{el_N_ref}_coeff) * E_bulk_fu = {equivalent_bulk_energy_for_phi}") # LOG
            print(f"[LOG][{label}]   E_slab - ((N_N_ref_slab / x_{el_N_ref}_coeff) * E_bulk_fu) (numerator for part1) = {phi_part1_numerator}") # LOG
            print(f"[LOG][{label}]   φ Part 1 (energy term): {phi_part1} eV/Å^2") # LOG
            print(f"[LOG][{label}]   Γ_NM * E_M_bulk_atom = {Gamma_NM * E_M_bulk_atom}") # LOG
            print(f"[LOG][{label}]   φ Part 2 (chem pot M term): {phi_part2} eV/Å^2") # LOG
            print(f"[LOG][{label}]   Γ_NO * E_O_atom_ref = {Gamma_NO * E_O_atom_ref}") # LOG
            print(f"[LOG][{label}]   φ Part 3 (chem pot O term): {phi_part3} eV/Å^2") # LOG
            print(f"[LOG][{label}] Total φ value: {phi_value} eV/Å^2") # LOG

        # --- Chemical potential ranges and stability conditions ---
        formation_energy_MxNyOz = formation_data_dict.get('formation_enthalpy_ev')
        if not first_slab_logged and formation_energy_MxNyOz is not None:
             print(f"[LOG][{label}] Formation enthalpy of bulk (ΔH_f(M_{x_M_coeff}{el_N_ref}_{x_N_ref_coeff}{el_O}_{x_O_coeff})): {formation_energy_MxNyOz} eV") # LOG

        # Initialize with a default practical lower bound
        default_min_mu = -5.0
        delta_mu_M_min = default_min_mu
        delta_mu_O_min = default_min_mu

        if formation_energy_MxNyOz is not None:
            # Calculate the stability limit for delta_mu_M (when delta_mu_O = 0)
            if x_M_coeff > 0:
                limit_M_from_FE = formation_energy_MxNyOz / x_M_coeff
                # The actual lower bound is this limit, but not exceeding 0.0
                delta_mu_M_min = min(0.0, limit_M_from_FE)

            # Calculate the stability limit for delta_mu_O (when delta_mu_M = 0)
            if x_O_coeff > 0:
                limit_O_from_FE = formation_energy_MxNyOz / x_O_coeff
                # The actual lower bound is this limit, but not exceeding 0.0
                delta_mu_O_min = min(0.0, limit_O_from_FE)
        
        # Ensure that even if default values were used (e.g., FE not provided),
        # or if calculated limits were positive (e.g. positive FE),
        # the final min values for chemical potentials do not exceed 0.0.
        delta_mu_M_min = min(delta_mu_M_min, 0.0)
        delta_mu_O_min = min(delta_mu_O_min, 0.0)

        if not first_slab_logged:
            print(f"[LOG][{label}] Δμ_{el_M} range for plots: [{delta_mu_M_min}, 0.0] eV") # LOG
            print(f"[LOG][{label}] Δμ_{el_O} range for plots: [{delta_mu_O_min}, 0.0] eV") # LOG

        if sampling is None:
            sampling = Int(100)
        delta_mu_M_range = np.linspace(delta_mu_M_min, 0.0, sampling.value) # User defined sampling
        delta_mu_O_range = np.linspace(delta_mu_O_min, 0.0, sampling.value) # User defined sampling

        gamma_values_grid = {}
        for delta_mu_M_val in delta_mu_M_range:
            for delta_mu_O_val in delta_mu_O_range:
                gamma = phi_value - Gamma_NM * delta_mu_M_val - Gamma_NO * delta_mu_O_val
                key = f"muM_{delta_mu_M_val:.4f}_muO_{delta_mu_O_val:.4f}"
                gamma_values_grid[key] = float(gamma)

        # --- Calculate gamma for fixed delta_mu_M = 0 ---
        gamma_values_fixed_muM_zero = {}
        delta_mu_M_fixed_val = 0.0
        for delta_mu_O_val in delta_mu_O_range:
            gamma = phi_value - Gamma_NM * delta_mu_M_fixed_val - Gamma_NO * delta_mu_O_val
            key = f"muM_0.0000_muO_{delta_mu_O_val:.4f}" # muM is fixed at 0.0
            gamma_values_fixed_muM_zero[key] = float(gamma)

        if not first_slab_logged:
            print(f"[LOG][{label}] Example γ (at Δμ_{el_M}=0, Δμ_{el_O}=0, if stable): {gamma_values_grid.get('muM_0.0000_muO_0.0000', 'N/A (not in stable region or not calculated)')}") # LOG
            first_slab_logged = True

        n_i_original_calc = {}
        for el, bulk_c in bulk_stoich.items():
            if bulk_c > 0:
                n_i_original_calc[el] = slab_atom_counts.get(el, 0) / bulk_c
            else:
                n_i_original_calc[el] = 0.0
        n_values_original = [v for v in n_i_original_calc.values() if v > 0]
        formula_units_original_def = min(n_values_original) if n_values_original else 1.0
        excess_atoms_original_def = {
            el: slab_atom_counts.get(el, 0) - formula_units_original_def * bulk_stoich.get(el,0)
            for el in [el_M, el_N_ref, el_O]
        }

        current_results_dict = {
            'phi': float(phi_value),
            'Gamma_M_vs_Nref': float(Gamma_NM),
            'Gamma_O_vs_Nref': float(Gamma_NO),
            'gamma_values_grid': gamma_values_grid,
            'gamma_values_fixed_muM_zero': gamma_values_fixed_muM_zero,
            'area_A2': float(area),
            'element_M_independent': el_M,
            'element_N_reference': el_N_ref,
            'bulk_stoichiometry_MxNyOz': {
                f"x_{el_M}": int(x_M_coeff),
                f"x_{el_N_ref}": int(x_N_ref_coeff),
                f"x_{el_O}": int(x_O_coeff)
            },
            'slab_atom_counts': {
                f"N_{el_M}": int(N_M_slab),
                f"N_{el_N_ref}": int(N_N_ref_slab),
                f"N_{el_O}": int(N_O_slab)
            },
            'reference_energies_per_atom': {k:float(v) for k,v in element_ref_energies.items()},
            'E_slab_eV': float(E_slab),
            'E_bulk_fu_eV': float(E_bulk_fu),
            'formula_units_original_definition': float(formula_units_original_def),
            'excess_atoms_original_definition': {k:float(v) for k,v in excess_atoms_original_def.items()},
        }
        if formation_energy_MxNyOz is not None:
            current_results_dict['formation_enthalpy_MxNyOz_eV'] = float(formation_energy_MxNyOz)

        results[label] = Dict(dict=current_results_dict)

    return results

if __name__ == "__main__":
    from aiida.orm import StructureData, Dict
    from ase import Atoms

    mock_bulk_structure = StructureData(ase=read('/home/thiagotd/git/aiida_teros/examples/structures/bulk/Ag6O8P2_optimized.cif'))

    mock_bulk_parameters = Dict(dict={
        'total_energies': {
            'energy_extrapolated': -83.74683825
        }
    })

    mock_formation_enthalpy_data = {
        'formation_enthalpy_ev': -8.800103665,
        'ag_energy_per_atom': -2.717776815,
        'p_energy_per_atom': -5.200724015,
        'o_energy_per_atom': -4.92981525
    }
    mock_formation_enthalpy = Dict(dict=mock_formation_enthalpy_data)

    mock_code = 'VASP'

    mock_slab1_structure = load_node(23123)
    mock_slab1_parameters = Dict(dict={
        'total_energies': {
            'energy_extrapolated': -246.86343315
        }
    })

    mock_slab_structures = {'slab1': mock_slab1_structure}
    mock_slab_parameters = {'slab1_params': mock_slab1_parameters}

    print("Testing calculate_surface_energy_ternar...")
    try:
        results_dict = calculate_surface_energy_ternar(
            bulk_structure=mock_bulk_structure,
            bulk_parameters=mock_bulk_parameters,
            formation_enthalpy=mock_formation_enthalpy,
            code=mock_code,
            slab_structures=mock_slab_structures,
            slab_parameters=mock_slab_parameters
        )

        for slab_label, result_data_node in results_dict.items():
            print(f"\nResults for {slab_label}:")
            result_data = result_data_node.get_dict()
            for key, value in result_data.items():
                if isinstance(value, dict):
                    print(f"  {key}:")
                    for sub_key, sub_value in value.items():
                        print(f"    {sub_key}: {sub_value}")
                else:
                    print(f"  {key}: {value}")
    
    except ValueError as e:
        print(f"An error occurred during the test: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")