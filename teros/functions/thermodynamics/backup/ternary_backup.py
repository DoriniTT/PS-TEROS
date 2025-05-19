from aiida.orm import Dict
from aiida_workgraph import task
from collections import Counter
from math import gcd
from functools import reduce
import numpy as np
from aiida.orm import Dict

@task.calcfunction()
def calculate_surface_energy_ternary(bulk_structure, bulk_parameters, formation_enthalpy=None, code=None, **kwargs):
    """
    Calculate the surface Gibbs free energy per unit area (γ) for each slab in a ternary oxide system.

    Args:
        bulk_structure: AiiDA StructureData of the bulk
        bulk_parameters: Dict with energy data
        formation_enthalpy: Dict with formation enthalpy data
        code: DFT code identifier ("QUANTUM_ESPRESSO", "CP2K", "VASP")
        **kwargs: Contains slab_structures and slab_parameters

    Returns:
        Dict with γ for each slab (in eV/Å²)
    """
    slab_structures = kwargs.get('slab_structures', {})
    slab_parameters = kwargs.get('slab_parameters', {})

    bulk_atoms = bulk_structure.get_ase()
    elements_count_raw = Counter(bulk_atoms.get_chemical_symbols())

    formula_gcd = 1
    elements_count = {el: cnt // formula_gcd for el, cnt in elements_count_raw.items()}
    elements = list(elements_count.keys())

    oxygen_idx = elements.index('O') if 'O' in elements else None
    metal_indices = [i for i, el in enumerate(elements) if el != 'O']

    if oxygen_idx is None:
        raise ValueError("No oxygen atoms found in structure. This workflow is designed for oxide materials.")
    if len(elements) < 3:
        raise ValueError(f"Less than 3 elements found: {elements}. This workflow is designed for ternary oxides.")

    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        E_bulk = bulk_parameters.get_dict()['energy']
    else:
        E_bulk = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']

    formation_dict = formation_enthalpy.get_dict()
    formation_energy = formation_dict.get('formation_enthalpy_ev')
    element_energy_per_atom = {
        element: formation_dict.get(f'{element.lower()}_energy_per_atom', 0.0)
        for element in elements
    }
    E_O2 = element_energy_per_atom.get('O', 0.0) * 2

    results = {}

    for i, ((_, misc), (_, structure)) in enumerate(zip(slab_parameters.items(), slab_structures.items())):
        label = f's_{str(i)}'
        slab_atoms = structure.get_ase()
        slab_elements_count = Counter(slab_atoms.get_chemical_symbols())
        area = slab_atoms.get_volume() / slab_atoms.get_cell()[2][2]

        if code in ['QUANTUM_ESPRESSO', 'CP2K']:
            E_slab = misc.get_dict()['energy']
        else:
            E_slab = misc.get_dict()['total_energies']['energy_extrapolated']
        if E_slab is None:
            continue

        n_i = {}
        for element, bulk_count in elements_count.items():
            if bulk_count > 0:
                n_i[element] = slab_elements_count.get(element, 0) / bulk_count
            else:
                n_i[element] = 0.0

        n_values = [v for v in n_i.values() if v > 0]
        if n_values and max(n_values) - min(n_values) < 1e-3:
            formula_units = n_values[0]
            stoichiometric_slab = True
        else:
            formula_units = min(n_values) if n_values else 1
            stoichiometric_slab = False

        excess_atoms = {
            element: slab_elements_count.get(element, 0) - formula_units * bulk_count
            for element, bulk_count in elements_count.items()
        }
        excess_contribution = sum(
            excess * element_energy_per_atom.get(element, 0)
            for element, excess in excess_atoms.items()
        )

        theta = (E_slab - formula_units * E_bulk - excess_contribution) / (2 * area)

        primary_metal = next((el for el in elements if el != 'O'), None)
        if primary_metal is None:
            raise ValueError("No primary metal found in elements.")

        counts = list(elements_count.values())
        formula_gcd = reduce(gcd, counts)
        x = elements_count[primary_metal] // formula_gcd
        z = elements_count['O'] // formula_gcd
        
        if formation_energy is not None:
            delta_mu_M_min = formation_energy / x if x > 0 else -5.0
            delta_mu_O_min = formation_energy / z if z > 0 else -2.0
        else:
            delta_mu_M_min = -5.0
            delta_mu_O_min = -2.0
        
        delta_mu_M_range = np.linspace(delta_mu_M_min, 0.0, 10)
        delta_mu_O_range = np.linspace(delta_mu_O_min, 0.0, 10)
        
        gamma_values_grid = {}
        for delta_mu_M in delta_mu_M_range:
            for delta_mu_O in delta_mu_O_range:
                if x*delta_mu_M + z*delta_mu_O <= formation_energy or formation_energy is None:
                    gamma = theta
                    if primary_metal in excess_atoms:
                        gamma -= (excess_atoms[primary_metal] * delta_mu_M) / (2 * area)
                    if 'O' in excess_atoms:
                        gamma -= (excess_atoms['O'] * delta_mu_O) / (2 * area)
                    key = f"muM={float(delta_mu_M):.4f},muO={float(delta_mu_O):.4f}"
                    gamma_values_grid[key] = float(gamma)

        gamma_values_fixed_metals = {}
        for delta_mu_O in delta_mu_O_range:
            gamma = theta
            if 'O' in excess_atoms:
                gamma -= (excess_atoms['O'] * delta_mu_O) / (2 * area)
            gamma_values_fixed_metals[float(delta_mu_O)] = float(gamma)

        results_dict = {
            'theta': float(theta),
            'gamma_values_fixed_metals': gamma_values_fixed_metals,
            'gamma_values_grid': gamma_values_grid,
            'area': float(area),
            'formula_units': float(formula_units),
            'element_counts': {element: int(count) for element, count in slab_elements_count.items()},
            'excess_atoms': {element: float(excess) for element, excess in excess_atoms.items()},
            'reference_energies': {element: float(energy) for element, energy in element_energy_per_atom.items()},
            'stoichiometric_slab': stoichiometric_slab,
            'n_i_per_element': n_i,
        }
        if formation_energy is not None:
            results_dict['formation_enthalpy_ev'] = float(formation_energy)
        results[label] = Dict(dict=results_dict)

    return results