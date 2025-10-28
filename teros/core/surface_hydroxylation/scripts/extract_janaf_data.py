#!/usr/bin/env python
"""
Extract JANAF thermochemical data for H2O, H2, O2.

Reads manually copied JANAF table data and outputs thermodynamics_data.json.
"""

import json
from pathlib import Path

# Unit conversion constants
KJ_TO_EV = 0.010364269  # kJ/mol to eV/molecule
J_TO_KJ = 0.001  # J to kJ

# JANAF table URLs for reference
JANAF_URLS = {
    'H2O': 'https://janaf.nist.gov/tables/H-064.html',
    'H2': 'https://janaf.nist.gov/tables/H-050.html',
    'O2': 'https://janaf.nist.gov/tables/O-029.html',
}

def calculate_delta_mu(H_kJ: float, S_J_K: float, T: float) -> float:
    """
    Calculate chemical potential correction.

    Δμ⁰(T) = [H(T) - H(0K)] - T·S(T)

    Args:
        H_kJ: Enthalpy H(T) - H(0K) in kJ/mol
        S_J_K: Entropy S(T) in J/(mol·K)
        T: Temperature in K

    Returns:
        Δμ⁰ in eV/molecule
    """
    S_kJ_K = S_J_K * J_TO_KJ  # Convert to kJ/(mol·K)
    delta_mu_kJ = H_kJ - T * S_kJ_K
    delta_mu_eV = delta_mu_kJ * KJ_TO_EV
    return delta_mu_eV


def process_molecule_data(raw_data: list[tuple]) -> dict:
    """
    Process raw JANAF data for one molecule.

    Args:
        raw_data: List of (T, H, S) tuples

    Returns:
        Dict mapping temperature to {H, S, delta_mu}
    """
    processed = {}

    for T, H, S in raw_data:
        delta_mu = calculate_delta_mu(H, S, T)
        processed[str(int(T))] = {
            'H': H,
            'S': S,
            'delta_mu': round(delta_mu, 6)
        }

    return processed


def main():
    """Extract and process JANAF data."""
    print("JANAF Data Extraction Script")
    print("=" * 50)

    # Test with H2O at 298 K
    # From JANAF: H(298) - H(0) = 9.904 kJ/mol, S(298) = 188.834 J/(mol·K)
    test_data = [(298, 9.904, 188.834)]
    result = process_molecule_data(test_data)

    print(f"\nTest calculation for H2O at 298 K:")
    print(f"  H = {result['298']['H']} kJ/mol")
    print(f"  S = {result['298']['S']} J/(mol·K)")
    print(f"  Δμ⁰ = {result['298']['delta_mu']} eV")
    print(f"  (Expected Δμ⁰ ≈ -0.564 eV)")

if __name__ == '__main__':
    main()
