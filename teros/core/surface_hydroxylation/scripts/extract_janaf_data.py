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

def main():
    """Extract and process JANAF data."""
    print("JANAF Data Extraction Script")
    print("=" * 50)

if __name__ == '__main__':
    main()
