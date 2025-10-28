"""
JANAF thermodynamics database for chemical potential corrections.

Provides temperature-dependent chemical potential corrections for gas-phase
species (H2O, H2, O2) from NIST-JANAF thermochemical tables.
"""

import json
import math
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

    def get_mu_correction(
        self,
        species: str,
        T: float,
        P: float = 1.0
    ) -> float:
        """
        Get chemical potential correction Δμ⁰(T, P).

        Args:
            species: Molecular species ('H2O', 'H2', 'O2')
            T: Temperature in Kelvin (must be in discrete points: 0, 50, 100, ..., 2000)
            P: Pressure in bar (default 1.0 = standard pressure)

        Returns:
            Chemical potential correction in eV/molecule

        Raises:
            ValueError: If species invalid or T not in database

        Example:
            >>> db = JanafDatabase()
            >>> db.get_mu_correction('H2O', T=298)
            -0.4806
        """
        # Validate species
        self._validate_species(species)

        # Validate temperature
        T_str = str(int(T))
        if T_str not in self._data['molecules'][species]['data']:
            available_temps = sorted([int(t) for t in self._data['molecules'][species]['data'].keys()])
            raise ValueError(
                f"Temperature {T} K not found in database for {species}. "
                f"Available: {min(available_temps)}-{max(available_temps)} K in 50 K steps."
            )

        # Lookup pre-calculated delta_mu
        delta_mu = self._data['molecules'][species]['data'][T_str]['delta_mu']

        # Add pressure correction if P != 1.0
        if P != 1.0:
            pressure_correction = KB_EV * T * math.log(P)
            delta_mu += pressure_correction

        return delta_mu

    def _validate_species(self, species: str) -> None:
        """Validate species is supported."""
        if species not in SPECIES:
            raise ValueError(
                f"Species '{species}' not found. Available: {SPECIES}"
            )
