from aiida.orm import Dict
from aiida_workgraph import task
from collections import Counter
import numpy as np
from aiida.orm import Dict

@task.calcfunction()
def calculate_surface_energy_binary(bulk_structure, bulk_parameters, formation_enthalpy=None, code=None, **kwargs):
    """
    Calculate surface Gibbs free energies (γ) for slabs of a binary oxide material.

    This function determines γ under oxygen-poor and oxygen-rich conditions,
    which correspond to the lower and upper bounds of the allowed oxygen chemical potential.

    :param bulk_structure: AiiDA ``StructureData`` node of the bulk material.
    :type bulk_structure: aiida.orm.StructureData
    :param bulk_parameters: AiiDA ``Dict`` node with energy data from the bulk DFT calculation.
    :type bulk_parameters: aiida.orm.Dict
    :param formation_enthalpy: AiiDA ``Dict`` node containing formation enthalpy data for the bulk material,
                               as returned by ``calculate_formation_enthalpy``. This includes reference energies.
    :type formation_enthalpy: aiida.orm.Dict, optional
    :param code: DFT code identifier ("QUANTUM_ESPRESSO", "CP2K", "VASP"), used for energy parsing.
    :type code: str, optional
    :param kwargs: Should contain:
                   - ``slab_structures`` (dict[str, aiida.orm.StructureData]): A dictionary mapping slab identifiers
                     to their AiiDA ``StructureData`` nodes.
                   - ``slab_parameters`` (dict[str, aiida.orm.Dict]): A dictionary mapping slab identifiers
                     to AiiDA ``Dict`` nodes containing energy data from slab DFT calculations.
    :type kwargs: dict

    :raises ValueError: If the material is not a binary compound or does not contain oxygen.

    :return: A dictionary where keys are slab identifiers (e.g., "s_0", "s_1") and values are
             AiiDA ``Dict`` nodes. Each inner ``Dict`` contains:
             - ``slab_energy`` (float): Total energy of the slab.
             - ``bulk_energy`` (float): Total energy of the bulk material.
             - ``area`` (float): Surface area of the slab in Å².
             - ``mu_O_min`` (float): Minimum oxygen chemical potential (O-poor limit).
             - ``mu_O_max`` (float): Maximum oxygen chemical potential (O-rich limit).
             - ``gamma_O_poor`` (float): Surface energy at the O-poor limit (eV/Å²).
             - ``gamma_O_rich`` (float): Surface energy at the O-rich limit (eV/Å²).
             - ``element_counts`` (dict[str, int]): Elemental composition of the slab.
             - ``reference_energies`` (dict[str, float]): Reference energies per atom for elements.
             - ``formation_enthalpy_ev`` (float): Formation enthalpy of the bulk material.
    :rtype: dict[str, aiida.orm.Dict]
    """
    # Block 1: Extract and Validate Input Parameters
    # ---------------------------------------------
    slab_structures = kwargs.get('slab_structures', {})
    slab_parameters = kwargs.get('slab_parameters', {})

    bulk_atoms = bulk_structure.get_ase()
    elements_count = Counter(bulk_atoms.get_chemical_symbols())
    elements = list(elements_count.keys())
    
    if len(elements) != 2:
        raise ValueError(f"Found {len(elements)} elements: {elements}. This workflow requires a binary compound.")
    
    if 'O' not in elements:
        raise ValueError(f"No oxygen atoms found in structure: {elements}. This workflow is designed for oxide materials.")
    
    metal = [element for element in elements if element != 'O'][0]
    
    # Block 2: Extract Energies and Reference Data
    # -------------------------------------------
    if code in ['QUANTUM_ESPRESSO', 'CP2K']:
        E_bulk = bulk_parameters.get_dict()['energy']
    else:  # VASP
        E_bulk = bulk_parameters.get_dict()['total_energies']['energy_extrapolated']
    
    formation_dict = formation_enthalpy.get_dict()
    element_energy_per_atom = {
        element: formation_dict.get(f'{element.lower()}_energy_per_atom', 0.0)
        for element in elements
    }
    E_O2 = element_energy_per_atom.get('O', 0.0) * 2

    # Block 3: Calculate Oxygen Chemical Potential Bounds
    # --------------------------------------------------
    x = elements_count[metal]
    y = elements_count['O']
    formation_energy = formation_dict.get('formation_enthalpy_ev') * (x + y)

    mu_O_min = (E_bulk - x * element_energy_per_atom[metal]) / y
    mu_O_from_formation = mu_O_min + formation_energy / y
    mu_O_max = min(0.0, mu_O_from_formation)
    
    delta_mu_O_range = np.linspace(mu_O_min, mu_O_max, 10)
    
    # Block 4: Process Each Slab Structure and Calculate Properties
    # ------------------------------------------------------------
    results = {}
    
    for i, ((_, misc), (_, structure)) in enumerate(zip(slab_parameters.items(), slab_structures.items())):        
        label = f's_{str(i)}'
        slab_atoms = structure.get_ase()
        slab_elements_count = Counter(slab_atoms.get_chemical_symbols())
        
        area = slab_atoms.get_volume() / slab_atoms.get_cell()[2][2]
        
        if code in ['QUANTUM_ESPRESSO', 'CP2K']:
            E_slab = misc.get_dict()['energy']
        else:  # VASP
            E_slab = misc.get_dict()['total_energies']['energy_extrapolated']
        
        if E_slab is None:
            continue
        
        N_M_slab = slab_elements_count.get(metal, 0)
        N_O_slab = slab_elements_count.get('O', 0)
        
        N_M_bulk = x
        
        mu_O_poor = mu_O_min
        mu_M_bulk = element_energy_per_atom[metal]
        
        term1 = E_slab
        term2 = N_M_slab * (E_bulk / N_M_bulk)
        term3 = ((y/x) * N_M_slab - N_O_slab) * mu_O_min
        gamma_O_poor = (term1 - term2 + term3) / (2 * area)
        
        #stoichiometric_imbalance = (y/x * N_M_slab - N_O_slab)
        #gamma_O_rich = gamma_O_poor + (1 / (2 * area)) * stoichiometric_imbalance * formation_energy

        stoichiometric_imbalance = (y/x * N_M_slab - N_O_slab)
        gamma_O_rich = gamma_O_poor - (1 / (2 * area)) * stoichiometric_imbalance * formation_energy / y

        results_dict = {
            'slab_energy': float(E_slab),
            'bulk_energy': float(E_bulk),
            'area': float(area),
            'mu_O_min': float(mu_O_min),
            'mu_O_max': float(mu_O_max),
            'gamma_O_poor': float(gamma_O_poor),
            'gamma_O_rich': float(gamma_O_rich),
            'element_counts': {element: int(count) for element, count in slab_elements_count.items()},
            'reference_energies': {element: float(energy) for element, energy in element_energy_per_atom.items()},
        }
        
        results_dict['formation_enthalpy_ev'] = float(formation_energy)
        
        results[label] = Dict(dict=results_dict)

    return results