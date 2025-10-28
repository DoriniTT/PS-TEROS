"""
JANAF thermodynamics database for chemical potential corrections.

Provides temperature-dependent chemical potential corrections for gas-phase
species (H2O, H2, O2) from NIST-JANAF thermochemical tables.
"""

import json
from pathlib import Path
from typing import Optional


# Unit conversion constants
KJ_TO_EV = 0.010364269  # kJ/mol to eV/molecule
KB_EV = 8.617333e-5     # Boltzmann constant in eV/K

# Supported species
SPECIES = ['H2O', 'H2', 'O2']


class JanafDatabase:
    """
    Thermodynamic database from NIST-JANAF tables.

    Provides chemical potential corrections for gas-phase species
    used in surface energy calculations.
    """

    def __init__(self):
        """
        Load thermodynamics_data.json from module directory.

        Raises:
            FileNotFoundError: If JSON data file not found
            ValueError: If JSON schema is invalid
        """
        self._load_data()

    def _load_data(self):
        """Load and validate JSON data file."""
        json_path = Path(__file__).parent / "thermodynamics_data.json"

        if not json_path.exists():
            raise FileNotFoundError(
                f"Thermodynamics data file not found: {json_path}\n"
                f"Expected location: teros/core/surface_hydroxylation/thermodynamics_data.json"
            )

        with open(json_path) as f:
            self._data = json.load(f)

        # Validate schema
        required_keys = ['metadata', 'molecules']
        missing = [k for k in required_keys if k not in self._data]
        if missing:
            raise ValueError(f"Missing required keys in JSON: {missing}")

    @property
    def available_species(self) -> list:
        """Get list of available molecular species."""
        return list(self._data['molecules'].keys())
