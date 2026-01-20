"""
Surface energy calculator for hydroxylated surfaces.

Implements equations 4-10 from Section S3 for calculating surface Gibbs free energies.
"""

from collections import Counter
import numpy as np
from aiida_workgraph import task
from aiida.orm import Float, Int


def _get_formula_dict(structure):
    """
    Extract chemical formula as dictionary from structure.

    Args:
        structure: AiiDA StructureData

    Returns:
        dict: Element symbols to counts, e.g. {'Ag': 12, 'P': 4, 'O': 16}
    """
    # Get chemical symbols
    symbols = [site.kind_name for site in structure.sites]

    # Count occurrences
    formula_dict = dict(Counter(symbols))

    return formula_dict


def analyze_composition(slab_structure, bulk_structure):
    """
    Determine stoichiometry parameters (n, x, y) from slab composition.

    DEPRECATED: This function uses parametric formula Ag₃ₙPₙO₄ₙ₊ₓ₋₂ᵧH₂ₓ which
    only works for stoichiometric surfaces. Use analyze_composition_general()
    instead for non-stoichiometric surfaces.

    Algorithm:
        1. Parse bulk formula (e.g., Ag3PO4 → {'Ag': 3, 'P': 1, 'O': 4})
        2. Parse slab formula (e.g., Ag12P4O18H2)
        3. Determine n from element ratios
        4. Calculate x from hydrogen count: x = n_h / 2
        5. Calculate y from oxygen deficit: y = (n_o_deficit + x) / 2

    Args:
        slab_structure: AiiDA StructureData for modified surface
        bulk_structure: AiiDA StructureData for bulk unit cell

    Returns:
        dict: {
            'n': int,              # Bulk formula units in slab
            'x': int,              # OH groups added
            'y': int,              # Net O atoms removed
            'n_h': int,            # Total H atoms
            'n_o_deficit': int,    # O deficit vs n×bulk
            'formulas': {
                'bulk': str,
                'slab': str
            }
        }

    Raises:
        ValueError: If composition analysis fails or gives non-integer results
    """
    # Parse formulas
    bulk_formula = _get_formula_dict(bulk_structure)
    slab_formula = _get_formula_dict(slab_structure)

    # Determine n (number of bulk units in slab)
    # Use elements present in bulk (exclude H which isn't in bulk)
    ratios = []
    for element, bulk_count in bulk_formula.items():
        if element in slab_formula:
            ratio = slab_formula[element] / bulk_count
            ratios.append((element, ratio))

    # All ratios should be equal for stoichiometric slab
    # Use most common ratio as n
    n_values = [r[1] for r in ratios]
    n = round(sum(n_values) / len(n_values))  # Average and round

    # Validate n is positive integer
    if n <= 0:
        raise ValueError(f"Invalid n={n}, must be positive")

    # Calculate hydrogen addition
    n_h = slab_formula.get("H", 0)

    # Validate n_h is even (required for H_{2x} formula)
    if n_h % 2 != 0:
        raise ValueError(
            f"Number of H atoms must be even for H_{{2x}} formula, got n_h={n_h}"
        )

    # Calculate oxygen deficit
    expected_o = n * bulk_formula["O"]
    actual_o = slab_formula["O"]
    n_o_deficit = expected_o - actual_o

    # Derive x and y
    # From modified formula: Ag_{3n}P_nO_{4n+x-2y}H_{2x}
    # n_h = 2x → x = n_h / 2
    # n_o_deficit = -(x - 2y) → y = (n_o_deficit + x) / 2
    x = n_h // 2
    y = (n_o_deficit + x) // 2

    # Validate results
    if (n_o_deficit + x) % 2 != 0:
        raise ValueError(
            f"O balance gives non-integer y: n_o_deficit={n_o_deficit}, x={x}"
        )

    if (n_o_deficit + x) != 2 * y:
        raise ValueError(
            f"O balance inconsistent: n_o_deficit={n_o_deficit}, x={x}, y={y}"
        )

    # Format output
    bulk_str = "".join(f"{elem}{count}" for elem, count in sorted(bulk_formula.items()))
    slab_str = "".join(f"{elem}{count}" for elem, count in sorted(slab_formula.items()))

    return {
        "n": n,
        "x": x,
        "y": y,
        "n_h": n_h,
        "n_o_deficit": n_o_deficit,
        "formulas": {
            "bulk": bulk_str,
            "slab": slab_str,
        },
    }


def analyze_composition_general(slab_structure, pristine_structure):
    """
    Determine composition changes relative to pristine slab (general approach).

    This approach does NOT assume any specific parametric formula like Ag₃ₙPₙO₄ₙ₊ₓ₋₂ᵧH₂ₓ.
    It works for ANY surface composition, stoichiometric or not.

    Formation reaction: pristine + Δn_H2O H₂O → modified + Δn_O2 O₂

    Where:
        - Δn_H2O = (H_modified - H_pristine) / 2 (H₂O molecules added)
        - Δn_O2 = (O_pristine - O_modified + Δn_H2O) / 2 (O₂ molecules released)

    Args:
        slab_structure: AiiDA StructureData for modified surface
        pristine_structure: AiiDA StructureData for unmodified slab

    Returns:
        dict: {
            'n_h2o': float,        # H₂O molecules added (can be negative)
            'n_o2': float,         # O₂ molecules released (can be negative)
            'delta_h': int,        # Change in H atoms
            'delta_o': int,        # Change in O atoms
            'formulas': {
                'pristine': str,
                'modified': str
            }
        }

    Raises:
        ValueError: If composition analysis fails
    """
    # Parse formulas
    pristine_formula = _get_formula_dict(pristine_structure)
    slab_formula = _get_formula_dict(slab_structure)

    # Calculate changes in H and O
    h_pristine = pristine_formula.get("H", 0)
    h_slab = slab_formula.get("H", 0)
    delta_h = h_slab - h_pristine

    o_pristine = pristine_formula["O"]
    o_slab = slab_formula["O"]
    delta_o = o_slab - o_pristine

    # Validate H change is even (required for H₂O)
    if delta_h % 2 != 0:
        raise ValueError(f"H change must be even for H₂O, got delta_H={delta_h}")

    # Calculate H₂O and O₂ exchange
    n_h2o = delta_h / 2.0
    n_o2 = (o_pristine - o_slab + n_h2o) / 2.0

    # Format output
    pristine_str = "".join(
        f"{elem}{count}" for elem, count in sorted(pristine_formula.items())
    )
    slab_str = "".join(f"{elem}{count}" for elem, count in sorted(slab_formula.items()))

    return {
        "n_h2o": n_h2o,
        "n_o2": n_o2,
        "delta_h": delta_h,
        "delta_o": delta_o,
        "formulas": {
            "pristine": pristine_str,
            "modified": slab_str,
        },
    }


# Internal helper functions (not @calcfunction) for use within other functions
def _calc_delta_g_general_reaction1_internal(
    e_modified, e_pristine, n_h2o, n_o2, delta_mu_h2o, delta_mu_o2
):
    """Internal helper: Calculate ΔG for Reaction 1 (H2O/O2). Uses plain floats."""
    delta_g = e_modified - e_pristine - n_h2o * delta_mu_h2o + n_o2 * delta_mu_o2
    return delta_g


def _calc_delta_g_general_reaction2_internal(
    e_modified, e_pristine, n_h2o, n_o2, delta_mu_h2, delta_mu_h2o
):
    """Internal helper: Calculate ΔG for Reaction 2 (H2/H2O). Uses plain floats."""
    delta_g = (
        e_modified
        - e_pristine
        - n_h2o * delta_mu_h2
        + (n_o2 - 0.5 * n_h2o) * delta_mu_h2o
    )
    return delta_g


def _calc_delta_g_general_reaction3_internal(
    e_modified, e_pristine, n_h2o, n_o2, delta_mu_h2, delta_mu_o2
):
    """Internal helper: Calculate ΔG for Reaction 3 (H2/O2). Uses plain floats."""
    delta_g = (
        e_modified
        - e_pristine
        - n_h2o * delta_mu_h2
        + (n_o2 - 0.5 * n_h2o) * delta_mu_o2
    )
    return delta_g


@task.calcfunction
def calc_delta_g_general_reaction1(
    E_modified, E_pristine, n_h2o, n_o2, delta_mu_h2o, delta_mu_o2
):
    """
    Calculate formation energy using general approach - Reaction 1: H2O/O2 reservoirs.

    Formation reaction: pristine + n_H2O H₂O → modified + n_O2 O₂

    ΔG = E_modified - E_pristine - n_H2O·μ_H2O + n_O2·μ_O2

    Args:
        E_modified: Float - Total energy of modified slab (eV)
        E_pristine: Float - Total energy of pristine slab (eV)
        n_h2o: Float - Number of H₂O molecules added
        n_o2: Float - Number of O₂ molecules released
        delta_mu_h2o: Float - Δμ(H2O) from JANAF (eV/molecule)
        delta_mu_o2: Float - Δμ(O2) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV
    """
    delta_g = (
        E_modified.value
        - E_pristine.value
        - n_h2o.value * delta_mu_h2o.value
        + n_o2.value * delta_mu_o2.value
    )

    return Float(delta_g)


@task.calcfunction
def calc_delta_g_general_reaction2(
    E_modified, E_pristine, n_h2o, n_o2, delta_mu_h2, delta_mu_h2o
):
    """
    Calculate formation energy using general approach - Reaction 2: H2/H2O reservoirs.

    Formation reaction: pristine + n_H2O H₂O → modified + n_O2 O₂
    But H₂O = H₂ + (1/2)O₂, so:
    pristine + n_H2O H₂ → modified + (n_O2 - n_H2O) O₂ + n_H2O H₂O

    Simplifying: pristine + n_H2O H₂ → modified + (n_O2 - 0.5·n_H2O) O₂

    Args:
        E_modified: Float - Total energy of modified slab (eV)
        E_pristine: Float - Total energy of pristine slab (eV)
        n_h2o: Float - H₂O equivalent (from H balance)
        n_o2: Float - O₂ released (from O balance with H₂O)
        delta_mu_h2: Float - Δμ(H2) from JANAF (eV/molecule)
        delta_mu_h2o: Float - Δμ(H2O) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV
    """
    # Convert reaction to H2 reservoir
    # pristine + n_h2o·H2 → modified + (n_o2 - 0.5·n_h2o)·H2O
    delta_g = (
        E_modified.value
        - E_pristine.value
        - n_h2o.value * delta_mu_h2.value
        + (n_o2.value - 0.5 * n_h2o.value) * delta_mu_h2o.value
    )

    return Float(delta_g)


@task.calcfunction
def calc_delta_g_general_reaction3(
    E_modified, E_pristine, n_h2o, n_o2, delta_mu_h2, delta_mu_o2
):
    """
    Calculate formation energy using general approach - Reaction 3: H2/O2 reservoirs.

    Formation reaction: pristine + n_H2O H₂O → modified + n_O2 O₂
    But H₂O = H₂ + (1/2)O₂, so:
    pristine + n_H2O H₂ + 0.5·n_H2O O₂ → modified + n_O2 O₂

    Simplifying: pristine + n_H2O H₂ → modified + (n_O2 - 0.5·n_H2O) O₂

    Args:
        E_modified: Float - Total energy of modified slab (eV)
        E_pristine: Float - Total energy of pristine slab (eV)
        n_h2o: Float - H₂O equivalent (from H balance)
        n_o2: Float - O₂ released (from O balance with H₂O)
        delta_mu_h2: Float - Δμ(H2) from JANAF (eV/molecule)
        delta_mu_o2: Float - Δμ(O2) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV
    """
    # pristine + n_h2o·H2 → modified + (n_o2 - 0.5·n_h2o)·O2
    delta_g = (
        E_modified.value
        - E_pristine.value
        - n_h2o.value * delta_mu_h2.value
        + (n_o2.value - 0.5 * n_h2o.value) * delta_mu_o2.value
    )

    return Float(delta_g)


@task.calcfunction
def calc_delta_g_reaction1(E_slab, E_bulk, n, x, y, delta_mu_h2o, delta_mu_o2):
    """
    Calculate formation energy using Reaction 1: H2O/O2 reservoirs.

    Reaction: n Ag3PO4 + x H2O ↔ surface + y O2
    Use case: Water-rich oxidizing environment

    Equation 4: ΔG = E + 2y·μ_O - n·E_bulk - x·μ_H2O

    Args:
        E_slab: Float - Total energy of modified slab (eV)
        E_bulk: Float - Bulk energy per formula unit (eV/f.u.)
        n: Int - Number of bulk formula units in slab
        x: Int - Number of OH groups added
        y: Int - Net oxygen atoms removed (can be negative)
        delta_mu_h2o: Float - Δμ(H2O) from JANAF (eV/molecule)
        delta_mu_o2: Float - Δμ(O2) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV

    Raises:
        ValueError: If validation fails (n <= 0 or E_bulk >= 0)
    """
    # Validate inputs
    if n.value <= 0:
        raise ValueError(f"n must be positive, got {n.value}")

    if E_bulk.value >= 0:
        raise ValueError(
            f"E_bulk should be negative (binding energy), got {E_bulk.value} eV. "
            f"Check that bulk energy is per formula unit, not total energy."
        )

    # Calculate atomic oxygen chemical potential from O2
    # μ(O) = (1/2)·μ(O2) = (1/2)·Δμ(O2)
    mu_o = 0.5 * delta_mu_o2.value

    # Equation 4: ΔG = E + 2y·μ_O - n·E_bulk - x·μ_H2O
    delta_g = (
        E_slab.value
        + 2 * y.value * mu_o
        - n.value * E_bulk.value
        - x.value * delta_mu_h2o.value
    )

    return Float(delta_g)


@task.calcfunction
def calc_delta_g_reaction2(E_slab, E_bulk, n, x, y, delta_mu_h2, delta_mu_h2o):
    """
    Calculate formation energy using Reaction 2: H2/H2O reservoirs.

    Reaction: n Ag3PO4 + x H2 ↔ surface + (2y-x) H2O
    Use case: Reducing environment with water

    Equation 5: ΔG = E + (2y-x)·μ_H2O - n·E_bulk - x·μ_H2

    Args:
        E_slab: Float - Total energy of modified slab (eV)
        E_bulk: Float - Bulk energy per formula unit (eV/f.u.)
        n: Int - Number of bulk formula units in slab
        x: Int - Number of OH groups added
        y: Int - Net oxygen atoms removed
        delta_mu_h2: Float - Δμ(H2) from JANAF (eV/molecule)
        delta_mu_h2o: Float - Δμ(H2O) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV
    """
    # Validate inputs
    if n.value <= 0:
        raise ValueError(f"n must be positive, got {n.value}")

    if E_bulk.value >= 0:
        raise ValueError(f"E_bulk should be negative, got {E_bulk.value}")

    # Equation 5
    delta_g = (
        E_slab.value
        + (2 * y.value - x.value) * delta_mu_h2o.value
        - n.value * E_bulk.value
        - x.value * delta_mu_h2.value
    )

    return Float(delta_g)


@task.calcfunction
def calc_delta_g_reaction3(E_slab, E_bulk, n, x, y, delta_mu_h2, delta_mu_o2):
    """
    Calculate formation energy using Reaction 3: H2/O2 reservoirs.

    Reaction: n Ag3PO4 + x H2 ↔ surface + (y-0.5x) O2
    Use case: Reducing environment

    Equation 6: ΔG = E + (2y-x)·μ_O - n·E_bulk - x·μ_H2

    Args:
        E_slab: Float - Total energy of modified slab (eV)
        E_bulk: Float - Bulk energy per formula unit (eV/f.u.)
        n: Int - Number of bulk formula units in slab
        x: Int - Number of OH groups added
        y: Int - Net oxygen atoms removed
        delta_mu_h2: Float - Δμ(H2) from JANAF (eV/molecule)
        delta_mu_o2: Float - Δμ(O2) from JANAF (eV/molecule)

    Returns:
        Float: Formation energy ΔG in eV
    """
    # Validate inputs
    if n.value <= 0:
        raise ValueError(f"n must be positive, got {n.value}")

    if E_bulk.value >= 0:
        raise ValueError(f"E_bulk should be negative, got {E_bulk.value}")

    # Calculate atomic oxygen chemical potential
    mu_o = 0.5 * delta_mu_o2.value

    # Equation 6
    delta_g = (
        E_slab.value
        + (2 * y.value - x.value) * mu_o
        - n.value * E_bulk.value
        - x.value * delta_mu_h2.value
    )

    return Float(delta_g)


def select_reaction_function(which_reaction):
    """
    Return the appropriate formation energy calculation function.

    Args:
        which_reaction: int (1, 2, or 3)
            1: H2O/O2 reservoirs (oxidizing, water-rich)
            2: H2/H2O reservoirs (reducing, water-present)
            3: H2/O2 reservoirs (reducing, dry)

    Returns:
        callable: The corresponding calc_delta_g function

    Raises:
        ValueError: If which_reaction not in {1, 2, 3}
    """
    reactions = {
        1: calc_delta_g_reaction1,
        2: calc_delta_g_reaction2,
        3: calc_delta_g_reaction3,
    }

    if which_reaction not in reactions:
        raise ValueError(f"which_reaction must be 1, 2, or 3, got {which_reaction}")

    return reactions[which_reaction]


def get_surface_area(structure):
    """
    Calculate surface area from structure cell parameters.

    Assumes surface is perpendicular to c-axis (standard slab convention).
    Area is calculated as magnitude of cross product of in-plane lattice vectors.

    Args:
        structure: AiiDA StructureData

    Returns:
        Float: Surface area in m² (converted from ų)
    """
    # Get cell vectors (3x3 array in ų)
    cell = np.array(structure.cell)

    # In-plane vectors (a and b)
    a = cell[0]
    b = cell[1]

    # Cross product gives area vector
    cross = np.cross(a, b)
    area_angstrom2 = np.linalg.norm(cross)

    # Convert ų to m²: 1 ų = 1e-20 m²
    area_m2 = area_angstrom2 * 1e-20

    return Float(area_m2)


# Unit conversion constant
EV_TO_J = 1.602176634e-19  # eV to Joules


@task.calcfunction
def calc_gamma_s(delta_g, area):
    """
    Calculate average surface energy of both surfaces.

    Equation 9: γ_s = ΔG / (2A)

    Where:
        γ_s = average energy of top and bottom surfaces
        ΔG = formation energy (eV)
        A = surface area (m²)
        Factor of 2 accounts for two surfaces in slab

    Args:
        delta_g: Float - Formation energy (eV)
        area: Float - Surface area (m²)

    Returns:
        Float: Average surface energy γ_s in J/m²

    Raises:
        ValueError: If area <= 0
    """
    # Validate
    if area.value <= 0:
        raise ValueError(f"Surface area must be positive, got {area.value} m²")

    # Convert energy from eV to J
    delta_g_J = delta_g.value * EV_TO_J

    # Equation 9: γ_s = ΔG / (2A)
    gamma_s = delta_g_J / (2 * area.value)

    return Float(gamma_s)


@task.calcfunction
def calc_gamma(gamma_s_modified, gamma_0_pristine):
    """
    Calculate corrected surface energy (removes pristine contribution).

    Equation 10: γ = 2γ_s - γ_0

    Physical meaning:
        γ_s = (γ_top + γ_bottom) / 2

        For single-sided modification:
            γ_top = modified surface energy (what we want)
            γ_bottom = γ_0 (pristine, unchanged)

        Solving for γ_top:
            γ_top = 2γ_s - γ_0

    Args:
        gamma_s_modified: Float - Average surface energy of modified structure (J/m²)
        gamma_0_pristine: Float - Pristine surface energy reference (J/m²)

    Returns:
        Float: Corrected surface energy γ in J/m²
    """
    # Equation 10
    gamma = 2 * gamma_s_modified.value - gamma_0_pristine.value

    return Float(gamma)


def calculate_surface_energies_general(
    structures_dict,
    energies_dict,
    temperature,
    pressures,
    which_reaction,
):
    """
    Calculate surface energies using general composition-based approach.

    This approach does NOT require parametric formula Ag₃ₙPₙO₄ₙ₊ₓ₋₂ᵧH₂ₓ.
    Works for ANY surface composition, including non-stoichiometric surfaces.

    Uses reference slab as reference state:
        Formation reaction: reference + Δn_H2O H₂O → modified + Δn_O2 O₂

    Reference structure selection:
        - Prefers structure with H=0 (true pristine)
        - Falls back to minimum-H structure if no H=0 exists
        - This handles 'combine' mode where all structures have some H

    Main steps:
        1. Identify reference structure (minimum H count)
        2. Calculate composition changes relative to reference
        3. Calculate formation energies using chemical potential corrections
        4. Calculate surface energies

    Args:
        structures_dict: dict - {name: StructureData}
        energies_dict: dict - {name: Float (eV)}
        temperature: float - Temperature in Kelvin
        pressures: dict - {'H2O': float, 'O2': float, 'H2': float} in bar
        which_reaction: int - Formation reaction to use (1, 2, or 3)

    Returns:
        dict: {
            'surface_energies': {name: gamma (J/m²)},
            'formation_energies': {name: DeltaG (eV)},
            'composition_changes': {name: {'n_h2o': float, 'n_o2': float, ...}},
            'reference_data': {
                'temperature': float (K),
                'pressures': dict,
                'reaction_used': int,
                'mu_corrections': dict,
                'pristine_name': str,
                'E_pristine': float (eV),
            }
        }

    Raises:
        ValueError: If no pristine structure found or if required pressure missing
    """
    from teros.core.surface_hydroxylation.thermodynamics import JanafDatabase

    # Initialize JANAF database
    janaf_db = JanafDatabase()

    # Get chemical potential corrections
    mu_corrections = {}
    for species in ["H2O", "H2", "O2"]:
        if species in pressures:
            mu_corrections[species] = janaf_db.get_mu_correction(
                species=species, T=temperature, P=pressures[species]
            )

    # Step 1: Identify reference structure (minimum H count)
    # For combine mode where all structures have H, use minimum-H structure as reference
    min_h = float("inf")
    pristine_name = None

    for name, structure in structures_dict.items():
        formula_dict = _get_formula_dict(structure)
        n_h = formula_dict.get("H", 0)

        if n_h < min_h:
            min_h = n_h
            pristine_name = name

    if pristine_name is None:
        raise ValueError("No structures found in structures_dict")

    # Log reference selection
    if min_h == 0:
        print(f"Using pristine structure (H=0): {pristine_name}")
    else:
        print(f"Using minimum-H structure (H={min_h}) as reference: {pristine_name}")
        print(f"  (No true pristine H=0 found - typical for 'combine' mode)")

    pristine_structure = structures_dict[pristine_name]
    pristine_energy = energies_dict[pristine_name]

    # Step 2: Analyze all structures relative to pristine
    composition_changes = {}
    for name, structure in structures_dict.items():
        comp = analyze_composition_general(structure, pristine_structure)
        composition_changes[name] = comp

    # Step 3: Calculate formation energies
    formation_energies = {}

    # Select reaction functions (use internal helpers to avoid @calcfunction nesting)
    if which_reaction == 1:
        required_species = ["H2O", "O2"]
    elif which_reaction == 2:
        required_species = ["H2", "H2O"]
    elif which_reaction == 3:
        required_species = ["H2", "O2"]
    else:
        raise ValueError(f"which_reaction must be 1, 2, or 3, got {which_reaction}")

    # Check required species are present
    for species in required_species:
        if species not in mu_corrections:
            raise ValueError(
                f"Reaction {which_reaction} requires {species} pressure, but it was not provided"
            )

    # Extract plain float values for calculations
    # energies_dict values can be Float nodes or plain floats
    def get_value(val):
        return val.value if hasattr(val, "value") else val

    pristine_energy_val = get_value(pristine_energy)

    # Calculate ΔG for each structure using internal helper functions
    for name, comp in composition_changes.items():
        e_modified = get_value(energies_dict[name])

        if which_reaction == 1:
            delta_g = _calc_delta_g_general_reaction1_internal(
                e_modified=e_modified,
                e_pristine=pristine_energy_val,
                n_h2o=comp["n_h2o"],
                n_o2=comp["n_o2"],
                delta_mu_h2o=mu_corrections["H2O"],
                delta_mu_o2=mu_corrections["O2"],
            )
        elif which_reaction == 2:
            delta_g = _calc_delta_g_general_reaction2_internal(
                e_modified=e_modified,
                e_pristine=pristine_energy_val,
                n_h2o=comp["n_h2o"],
                n_o2=comp["n_o2"],
                delta_mu_h2=mu_corrections["H2"],
                delta_mu_h2o=mu_corrections["H2O"],
            )
        else:  # reaction 3
            delta_g = _calc_delta_g_general_reaction3_internal(
                e_modified=e_modified,
                e_pristine=pristine_energy_val,
                n_h2o=comp["n_h2o"],
                n_o2=comp["n_o2"],
                delta_mu_h2=mu_corrections["H2"],
                delta_mu_o2=mu_corrections["O2"],
            )

        formation_energies[name] = delta_g

    # Step 4: Calculate surface energies
    surface_energies = {}
    for name, structure in structures_dict.items():
        area = get_surface_area(structure)
        delta_g_value = formation_energies[name]

        # γ = ΔG / (2A)
        # Note: No pristine correction needed - ΔG is already relative to pristine!
        # Calculate directly to avoid calling @calcfunction from within @task
        # Conversion: eV to J/m²
        J_PER_EV = 1.602176634e-19  # J/eV
        gamma_j_m2 = (delta_g_value * J_PER_EV) / (2 * area)  # J/m²
        surface_energies[name] = gamma_j_m2

    return {
        "surface_energies": surface_energies,
        "formation_energies": formation_energies,
        "composition_changes": composition_changes,
        "reference_data": {
            "temperature": temperature,
            "pressures": pressures,
            "reaction_used": which_reaction,
            "mu_corrections": mu_corrections,
            "reference_name": pristine_name,
            "reference_H_count": int(min_h),
            "E_reference": pristine_energy_val,
        },
    }


def calculate_surface_energies(
    structures_dict,
    energies_dict,
    bulk_structure,
    bulk_energy,
    temperature,
    pressures,
    which_reaction,
):
    """
    Calculate surface energies for all structures in workflow.

    DEPRECATED: This function uses parametric formula Ag₃ₙPₙO₄ₙ₊ₓ₋₂ᵧH₂ₓ which
    only works for stoichiometric surfaces. Use calculate_surface_energies_general()
    instead for non-stoichiometric surfaces.

    Main orchestration function that:
        1. Gets chemical potential corrections from JANAF database
        2. Identifies pristine structure (x=0, y=0) for γ_0 reference
        3. Calculates formation energies for all structures
        4. Calculates surface energies with pristine correction

    Args:
        structures_dict: dict - {name: StructureData}
        energies_dict: dict - {name: Float (eV)}
        bulk_structure: StructureData - Bulk unit cell
        bulk_energy: Float - Bulk energy per formula unit (eV/f.u.)
        temperature: float - Temperature in Kelvin
        pressures: dict - {'H2O': float, 'O2': float, 'H2': float} in bar
        which_reaction: int - Formation reaction to use (1, 2, or 3)

    Returns:
        dict: {
            'surface_energies': {name: gamma (J/m²)},
            'formation_energies': {name: DeltaG (eV)},
            'stoichiometry': {name: {'n': int, 'x': int, 'y': int, ...}},
            'reference_data': {
                'temperature': float (K),
                'pressures': dict,
                'reaction_used': int,
                'mu_corrections': dict,
                'E_bulk': float (eV/f.u.),
                'gamma_0': float (J/m²),
            }
        }

    Raises:
        ValueError: If no pristine structure found or if required pressure missing
    """
    # DEBUG: Print to confirm this function is being called
    print("=" * 80)
    print(
        "DEBUG: calculate_surface_energies (OLD) called - redirecting to general version"
    )
    print("=" * 80)

    # Try using general approach instead
    return calculate_surface_energies_general(
        structures_dict=structures_dict,
        energies_dict=energies_dict,
        temperature=temperature,
        pressures=pressures,
        which_reaction=which_reaction,
    )
    from teros.core.surface_hydroxylation.thermodynamics import JanafDatabase

    # Initialize JANAF database
    janaf_db = JanafDatabase()

    # Get chemical potential corrections
    mu_corrections = {}
    for species in ["H2O", "H2", "O2"]:
        if species in pressures:
            mu_corrections[species] = janaf_db.get_mu_correction(
                species=species, T=temperature, P=pressures[species]
            )

    # Select reaction function
    calc_delta_g = select_reaction_function(which_reaction)

    # Prepare outputs
    surface_energies = {}
    formation_energies = {}
    stoichiometry_data = {}

    # Step 1: Identify pristine structure (no H atoms) for reference
    # For non-stoichiometric slabs, we can't use x=0,y=0 because pristine
    # itself may have oxygen excess/deficit vs bulk stoichiometry
    pristine_name = None
    for name, structure in structures_dict.items():
        # Get H count directly without full analysis
        formula_dict = _get_formula_dict(structure)
        n_h = formula_dict.get("H", 0)

        if n_h == 0:
            pristine_name = name
            break

    if pristine_name is None:
        raise ValueError(
            "No pristine structure found (n_H=0). "
            "Set include_empty=True in surface_params or ensure unmodified slab is included."
        )

    # Step 2: Analyze all structures using pristine as oxygen reference
    pristine_structure = structures_dict[pristine_name]

    for name, structure in structures_dict.items():
        # Analyze composition using pristine as reference for non-stoichiometric slabs
        comp = analyze_composition(structure, bulk_structure, pristine_structure)
        stoichiometry_data[name] = comp

    # Calculate γ_0 (pristine surface energy)
    pristine_comp = stoichiometry_data[pristine_name]
    pristine_area = get_surface_area(pristine_structure)

    # For pristine: x=0, y=0, so formation energy simplifies
    # Prepare arguments based on reaction
    if which_reaction == 1:
        pristine_delta_g = calc_delta_g(
            E_slab=energies_dict[pristine_name],
            E_bulk=bulk_energy,
            n=Int(pristine_comp["n"]),
            x=Int(0),
            y=Int(0),
            delta_mu_h2o=Float(mu_corrections.get("H2O", 0.0)),
            delta_mu_o2=Float(mu_corrections.get("O2", 0.0)),
        )
    elif which_reaction == 2:
        pristine_delta_g = calc_delta_g(
            E_slab=energies_dict[pristine_name],
            E_bulk=bulk_energy,
            n=Int(pristine_comp["n"]),
            x=Int(0),
            y=Int(0),
            delta_mu_h2=Float(mu_corrections.get("H2", 0.0)),
            delta_mu_h2o=Float(mu_corrections.get("H2O", 0.0)),
        )
    else:  # reaction 3
        pristine_delta_g = calc_delta_g(
            E_slab=energies_dict[pristine_name],
            E_bulk=bulk_energy,
            n=Int(pristine_comp["n"]),
            x=Int(0),
            y=Int(0),
            delta_mu_h2=Float(mu_corrections.get("H2", 0.0)),
            delta_mu_o2=Float(mu_corrections.get("O2", 0.0)),
        )

    pristine_gamma_s = calc_gamma_s(pristine_delta_g, pristine_area)
    gamma_0 = pristine_gamma_s  # For pristine, γ_s = γ_0 (no correction needed)

    # Process all structures
    for name, structure in structures_dict.items():
        comp = stoichiometry_data[name]
        area = get_surface_area(structure)

        # Calculate formation energy based on reaction
        if which_reaction == 1:
            delta_g = calc_delta_g(
                E_slab=energies_dict[name],
                E_bulk=bulk_energy,
                n=Int(comp["n"]),
                x=Int(comp["x"]),
                y=Int(comp["y"]),
                delta_mu_h2o=Float(mu_corrections["H2O"]),
                delta_mu_o2=Float(mu_corrections["O2"]),
            )
        elif which_reaction == 2:
            delta_g = calc_delta_g(
                E_slab=energies_dict[name],
                E_bulk=bulk_energy,
                n=Int(comp["n"]),
                x=Int(comp["x"]),
                y=Int(comp["y"]),
                delta_mu_h2=Float(mu_corrections["H2"]),
                delta_mu_h2o=Float(mu_corrections["H2O"]),
            )
        else:  # reaction 3
            delta_g = calc_delta_g(
                E_slab=energies_dict[name],
                E_bulk=bulk_energy,
                n=Int(comp["n"]),
                x=Int(comp["x"]),
                y=Int(comp["y"]),
                delta_mu_h2=Float(mu_corrections["H2"]),
                delta_mu_o2=Float(mu_corrections["O2"]),
            )

        formation_energies[name] = delta_g.value

        # Calculate surface energies
        gamma_s = calc_gamma_s(delta_g, area)
        gamma = calc_gamma(gamma_s, gamma_0)

        surface_energies[name] = gamma.value

    # Return complete results
    return {
        "surface_energies": surface_energies,
        "formation_energies": formation_energies,
        "stoichiometry": stoichiometry_data,
        "reference_data": {
            "temperature": temperature,
            "pressures": pressures,
            "reaction_used": which_reaction,
            "mu_corrections": {k: v for k, v in mu_corrections.items()},
            "E_bulk": bulk_energy.value,
            "gamma_0": gamma_0.value,
        },
    }
