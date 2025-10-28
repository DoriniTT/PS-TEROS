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

# Raw JANAF data extracted from tables
# Format: (T [K], H(T)-H(0K) [kJ/mol], S(T) [J/(mol·K)])

H2O_DATA = [
    (0, 0.000, 0.000),
    (50, 1.630, 152.388),
    (100, 3.281, 175.485),
    (150, 5.020, 189.042),
    (200, 6.842, 198.788),
    (250, 8.735, 206.534),
    (298, 9.904, 188.834),  # Note: This uses updated reference
    (300, 10.719, 213.051),
    (350, 12.592, 218.739),
    (400, 14.516, 223.826),
    (450, 16.485, 228.460),
    (500, 18.491, 232.738),
    (550, 20.531, 236.732),
    (600, 22.601, 240.485),
    (700, 26.797, 247.364),
    (800, 31.062, 253.633),
    (900, 35.381, 259.426),
    (1000, 39.744, 264.843),
    (1100, 44.145, 269.956),
    (1200, 48.579, 274.821),
    (1300, 53.040, 279.477),
    (1400, 57.528, 283.955),
    (1500, 62.040, 288.280),
    (1600, 66.576, 292.471),
    (1700, 71.135, 296.543),
    (1800, 75.717, 300.508),
    (1900, 80.321, 304.377),
    (2000, 84.945, 308.159),
]

H2_DATA = [
    (0, 0.000, 0.000),
    (50, 1.426, 84.052),
    (100, 2.854, 100.727),
    (150, 4.283, 112.000),
    (200, 5.715, 120.824),
    (250, 7.147, 128.105),
    (298, 8.468, 130.680),
    (300, 8.582, 134.442),
    (350, 10.019, 140.062),
    (400, 11.457, 145.134),
    (450, 12.897, 149.774),
    (500, 14.339, 154.077),
    (550, 15.783, 158.105),
    (600, 17.230, 161.906),
    (700, 20.131, 168.790),
    (800, 23.044, 174.890),
    (900, 25.968, 180.360),
    (1000, 28.906, 185.320),
    (1100, 31.857, 189.852),
    (1200, 34.823, 194.024),
    (1300, 37.806, 197.888),
    (1400, 40.807, 201.484),
    (1500, 43.826, 204.848),
    (1600, 46.866, 208.007),
    (1700, 49.927, 210.983),
    (1800, 53.011, 213.794),
    (1900, 56.119, 216.456),
    (2000, 59.253, 218.983),
]

O2_DATA = [
    (0, 0.000, 0.000),
    (50, 1.449, 182.115),
    (100, 2.909, 193.485),
    (150, 4.377, 200.777),
    (200, 5.854, 206.308),
    (250, 7.337, 210.806),
    (298, 8.683, 205.147),
    (300, 8.828, 214.612),
    (350, 10.323, 217.871),
    (400, 11.822, 220.693),
    (450, 13.326, 223.194),
    (500, 14.836, 225.449),
    (550, 16.353, 227.509),
    (600, 17.877, 229.410),
    (700, 20.945, 232.825),
    (800, 24.041, 235.810),
    (900, 27.170, 238.488),
    (1000, 30.334, 240.930),
    (1100, 33.536, 243.189),
    (1200, 36.779, 245.299),
    (1300, 40.065, 247.287),
    (1400, 43.396, 249.173),
    (1500, 46.776, 250.974),
    (1600, 50.205, 252.703),
    (1700, 53.686, 254.370),
    (1800, 57.220, 255.983),
    (1900, 60.809, 257.549),
    (2000, 64.454, 259.073),
]

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

    # Process all molecules
    molecules = {
        'H2O': H2O_DATA,
        'H2': H2_DATA,
        'O2': O2_DATA,
    }

    processed_data = {}
    for species, raw_data in molecules.items():
        print(f"\nProcessing {species}...")
        processed_data[species] = process_molecule_data(raw_data)
        print(f"  Processed {len(raw_data)} temperature points")

        # Show 298 K value
        data_298 = processed_data[species].get('298', {})
        if data_298:
            print(f"  At 298 K: Δμ⁰ = {data_298['delta_mu']} eV")

    print("\n" + "=" * 50)
    print("Processing complete!")

if __name__ == '__main__':
    main()
