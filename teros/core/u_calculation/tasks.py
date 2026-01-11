"""AiiDA calcfunction tasks for Hubbard U calculation.

These tasks extract d-electron occupations from VASP output and calculate
the Hubbard U parameter using linear response theory.

Expected VASP misc output structure (when LORBIT=11 and LDAU=True, LDAUTYPE=3):
    misc = {
        'site_magnetization': {
            'sphere_0': {'s': 0.1, 'p': 0.2, 'd': 5.5, 'f': 0.0, 'total': 5.8},
            'sphere_1': {'s': 0.1, 'p': 0.2, 'd': 5.4, 'f': 0.0, 'total': 5.7},
            ...
        },
        # OR (alternative format)
        'charge_per_site': [
            {'s': 0.1, 'p': 0.2, 'd': 5.5, 'f': 0.0},
            {'s': 0.1, 'p': 0.2, 'd': 5.4, 'f': 0.0},
            ...
        ],
        'symbols': ['Fe', 'Fe', 'O', 'O', ...],  # Atom types per site
        ...
    }
"""

import typing as t

from aiida import orm
from aiida_workgraph import task

from .utils import linear_regression


@task.calcfunction
def extract_d_electron_occupation(
    misc: orm.Dict,
    target_species: orm.Str,
    structure: orm.StructureData,
) -> orm.Dict:
    """
    Extract d-electron occupation for target species from VASP misc output.

    Requires LORBIT=11 in VASP INCAR to produce orbital-resolved occupations.
    The occupation is summed over all atoms of the target species.

    Args:
        misc: VASP misc output Dict
        target_species: Element symbol to extract occupations for (e.g., 'Fe')
        structure: StructureData to identify atom types

    Returns:
        Dict with:
            - total_d_occupation: Sum of d-occupations for target species
            - per_atom_d_occupation: List of d-occupation per atom
            - atom_indices: 0-based indices of target atoms
            - atom_count: Number of target atoms

    Raises:
        ValueError: If d-occupation data not found in misc output
    """
    misc_dict = misc.get_dict()
    species = target_species.value

    # Get atom symbols from structure
    ase_struct = structure.get_ase()
    symbols = ase_struct.get_chemical_symbols()

    # Find indices of target species
    target_indices = [i for i, s in enumerate(symbols) if s == species]

    if not target_indices:
        raise ValueError(
            f"Target species '{species}' not found in structure. "
            f"Available: {sorted(set(symbols))}"
        )

    # Try to extract d-occupations from different possible formats
    per_atom_d_occ = []

    # Format 1: site_magnetization with sphere_N keys
    if 'site_magnetization' in misc_dict:
        site_mag = misc_dict['site_magnetization']
        for idx in target_indices:
            sphere_key = f'sphere_{idx}'
            if sphere_key in site_mag:
                sphere_data = site_mag[sphere_key]
                if 'd' in sphere_data:
                    per_atom_d_occ.append(float(sphere_data['d']))
                else:
                    raise ValueError(
                        f"'d' occupation not found in {sphere_key}. "
                        f"Available keys: {list(sphere_data.keys())}. "
                        f"Ensure LORBIT=11 is set in INCAR."
                    )
            else:
                raise ValueError(
                    f"{sphere_key} not found in site_magnetization. "
                    f"Available: {list(site_mag.keys())}"
                )

    # Format 2: charge_per_site as list of dicts
    elif 'charge_per_site' in misc_dict:
        charges = misc_dict['charge_per_site']
        for idx in target_indices:
            if idx < len(charges):
                site_charges = charges[idx]
                if 'd' in site_charges:
                    per_atom_d_occ.append(float(site_charges['d']))
                else:
                    raise ValueError(
                        f"'d' occupation not found for site {idx}. "
                        f"Ensure LORBIT=11 is set in INCAR."
                    )
            else:
                raise ValueError(
                    f"Site index {idx} out of range for charge_per_site "
                    f"(length {len(charges)})"
                )

    # Format 3: projectors/local_charges (some VASP versions)
    elif 'projectors' in misc_dict:
        projectors = misc_dict['projectors']
        for idx in target_indices:
            if idx < len(projectors):
                proj = projectors[idx]
                if 'd' in proj:
                    per_atom_d_occ.append(float(proj['d']))
                elif 'd_total' in proj:
                    per_atom_d_occ.append(float(proj['d_total']))
                else:
                    raise ValueError(
                        f"'d' occupation not found in projectors for site {idx}. "
                        f"Available: {list(proj.keys())}"
                    )
            else:
                raise ValueError(
                    f"Site index {idx} out of range for projectors"
                )

    else:
        available_keys = list(misc_dict.keys())
        raise ValueError(
            f"Could not find orbital occupation data in misc output. "
            f"Expected 'site_magnetization', 'charge_per_site', or 'projectors'. "
            f"Available keys: {available_keys}. "
            f"Ensure LORBIT=11 is set in INCAR and parser is configured correctly."
        )

    total_d_occ = sum(per_atom_d_occ)

    return orm.Dict(dict={
        'total_d_occupation': total_d_occ,
        'per_atom_d_occupation': per_atom_d_occ,
        'atom_indices': target_indices,
        'atom_count': len(target_indices),
        'target_species': species,
    })


@task.calcfunction
def calculate_occupation_response(
    ground_state_occupation: orm.Dict,
    nscf_occupation: orm.Dict,
    scf_occupation: orm.Dict,
    potential_value: orm.Float,
) -> orm.Dict:
    """
    Calculate response functions (chi, chi_0) for a single potential value.

    The linear response theory relates occupation change to applied potential:
        ΔN_NSCF = χ₀ * V  (non-self-consistent response)
        ΔN_SCF = χ * V    (self-consistent response)

    Args:
        ground_state_occupation: D-occupation from ground state (no potential)
        nscf_occupation: D-occupation from non-SCF response calculation
        scf_occupation: D-occupation from SCF response calculation
        potential_value: Applied potential V (eV)

    Returns:
        Dict with:
            - potential: Applied potential V (eV)
            - ground_state_d: Ground state d-occupation
            - nscf_d: Non-SCF response d-occupation
            - scf_d: SCF response d-occupation
            - delta_n_nscf: N_nscf - N_ground
            - delta_n_scf: N_scf - N_ground
    """
    gs_dict = ground_state_occupation.get_dict()
    nscf_dict = nscf_occupation.get_dict()
    scf_dict = scf_occupation.get_dict()
    V = potential_value.value

    gs_d = gs_dict['total_d_occupation']
    nscf_d = nscf_dict['total_d_occupation']
    scf_d = scf_dict['total_d_occupation']

    delta_n_nscf = nscf_d - gs_d
    delta_n_scf = scf_d - gs_d

    return orm.Dict(dict={
        'potential': V,
        'ground_state_d': gs_d,
        'nscf_d': nscf_d,
        'scf_d': scf_d,
        'delta_n_nscf': delta_n_nscf,
        'delta_n_scf': delta_n_scf,
    })


@task.calcfunction
def calculate_hubbard_u_single_point(
    response: orm.Dict,
) -> orm.Dict:
    """
    Calculate Hubbard U from a single potential value (quick estimate).

    Uses U = 1/χ - 1/χ₀ where:
        χ = ΔN_SCF / V
        χ₀ = ΔN_NSCF / V

    Note: This is less accurate than linear regression with multiple potentials.

    Args:
        response: Response dict from calculate_occupation_response

    Returns:
        Dict with U value and intermediate quantities

    Raises:
        ValueError: If potential is zero (cannot calculate chi)
    """
    resp_dict = response.get_dict()
    V = resp_dict['potential']

    if abs(V) < 1e-10:
        raise ValueError(
            "Cannot calculate U from V=0 response. "
            "Use a non-zero potential or multiple potentials with linear regression."
        )

    delta_n_nscf = resp_dict['delta_n_nscf']
    delta_n_scf = resp_dict['delta_n_scf']

    chi_0 = delta_n_nscf / V
    chi = delta_n_scf / V

    # Avoid division by zero
    if abs(chi) < 1e-10 or abs(chi_0) < 1e-10:
        raise ValueError(
            f"Response is too small to calculate U. "
            f"chi={chi:.6f}, chi_0={chi_0:.6f}. "
            f"Try a larger potential value."
        )

    U = (1.0 / chi) - (1.0 / chi_0)

    return orm.Dict(dict={
        'U': U,
        'chi': chi,
        'chi_0': chi_0,
        'potential': V,
        'delta_n_scf': delta_n_scf,
        'delta_n_nscf': delta_n_nscf,
    })


@task.calcfunction
def calculate_hubbard_u_linear_regression(
    responses: orm.List,
) -> orm.Dict:
    """
    Calculate Hubbard U from multiple potential values using linear regression.

    Performs linear fits:
        ΔN_NSCF = χ₀ * V + c₀
        ΔN_SCF = χ * V + c

    Then calculates: U = 1/χ - 1/χ₀

    This is more accurate than single-point calculation as it averages
    over multiple perturbation strengths.

    Args:
        responses: List of response dicts from calculate_occupation_response
            Each dict should have: potential, delta_n_nscf, delta_n_scf

    Returns:
        Dict with:
            - U: Final Hubbard U value (eV)
            - chi_slope: SCF response slope
            - chi_0_slope: NSCF response slope
            - chi_intercept: SCF intercept
            - chi_0_intercept: NSCF intercept
            - chi_r2: SCF fit R²
            - chi_0_r2: NSCF fit R²
            - potential_values: List of potentials used
            - delta_n_scf_values: List of SCF occupation changes
            - delta_n_nscf_values: List of NSCF occupation changes

    Raises:
        ValueError: If too few data points or regression fails
    """
    responses_list = responses.get_list()

    if len(responses_list) < 2:
        raise ValueError(
            f"Need at least 2 potential values for linear regression, "
            f"got {len(responses_list)}"
        )

    # Extract data arrays
    potentials = []
    delta_n_nscf_vals = []
    delta_n_scf_vals = []

    for resp in responses_list:
        potentials.append(resp['potential'])
        delta_n_nscf_vals.append(resp['delta_n_nscf'])
        delta_n_scf_vals.append(resp['delta_n_scf'])

    # Linear regression for chi_0 (NSCF response)
    chi_0_slope, chi_0_intercept, chi_0_r2 = linear_regression(
        potentials, delta_n_nscf_vals
    )

    # Linear regression for chi (SCF response)
    chi_slope, chi_intercept, chi_r2 = linear_regression(
        potentials, delta_n_scf_vals
    )

    # Calculate U
    if abs(chi_slope) < 1e-10 or abs(chi_0_slope) < 1e-10:
        raise ValueError(
            f"Response slopes are too small to calculate U. "
            f"chi_slope={chi_slope:.6f}, chi_0_slope={chi_0_slope:.6f}. "
            f"This may indicate insufficient perturbation or numerical issues."
        )

    U = (1.0 / chi_slope) - (1.0 / chi_0_slope)

    return orm.Dict(dict={
        'U': U,
        'chi_slope': chi_slope,
        'chi_0_slope': chi_0_slope,
        'chi_intercept': chi_intercept,
        'chi_0_intercept': chi_0_intercept,
        'chi_r2': chi_r2,
        'chi_0_r2': chi_0_r2,
        'potential_values': potentials,
        'delta_n_scf_values': delta_n_scf_vals,
        'delta_n_nscf_values': delta_n_nscf_vals,
        'n_points': len(potentials),
    })


@task.calcfunction
def gather_responses(
    **kwargs: orm.Dict,
) -> orm.List:
    """
    Gather response Dicts into a List for linear regression.

    Args:
        **kwargs: Response dicts keyed by label (e.g., V_0, V_1, ...)

    Returns:
        List of response dicts sorted by potential value
    """
    responses = []

    for key, response_node in kwargs.items():
        if isinstance(response_node, orm.Dict):
            responses.append(response_node.get_dict())
        else:
            responses.append(dict(response_node))

    # Sort by potential value for consistent ordering
    responses.sort(key=lambda r: r.get('potential', 0.0))

    return orm.List(list=responses)
