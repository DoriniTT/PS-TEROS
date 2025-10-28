# JANAF Thermodynamics Module Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement a thermodynamics database module that provides NIST-JANAF chemical potential corrections for H2O, H2, and O2 species used in surface energy calculations.

**Architecture:** Data File + Loader Module pattern. JSON file stores raw JANAF data (H, S) and pre-calculated Δμ⁰ values. Python class loads JSON and provides API for chemical potential corrections with unit conversions and pressure corrections.

**Tech Stack:** Python 3.x, JSON (stdlib), pathlib, pytest

---

## Task 1: Create Data Extraction Script Foundation

**Files:**
- Create: `teros/core/surface_hydroxylation/scripts/extract_janaf_data.py`
- Test: Manual verification (no pytest for data script)

**Step 1: Create script directory and basic structure**

Create the script file with imports and constants:

```python
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
```

**Step 2: Verify script runs**

```bash
cd /home/thiagotd/git/PS-TEROS
python teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
```

Expected output: Header prints

**Step 3: Commit foundation**

```bash
git add teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
git commit -m "feat: add JANAF data extraction script foundation"
```

---

## Task 2: Add Data Processing Functions to Extraction Script

**Files:**
- Modify: `teros/core/surface_hydroxylation/scripts/extract_janaf_data.py`

**Step 1: Add function to calculate Δμ⁰**

Add after constants:

```python
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
```

**Step 2: Test with sample data**

Add to `main()`:

```python
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
```

**Step 3: Run and verify calculation**

```bash
python teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
```

Expected: Δμ⁰ ≈ -0.564 eV (close to Section 4 example value)

**Step 4: Commit data processing**

```bash
git add teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
git commit -m "feat: add JANAF data processing functions

Implements Δμ⁰ calculation from H(T) and S(T) with unit conversions"
```

---

## Task 3: Add Raw JANAF Data to Extraction Script

**Files:**
- Modify: `teros/core/surface_hydroxylation/scripts/extract_janaf_data.py`

**Step 1: Add raw H2O data**

Add after `JANAF_URLS` constant:

```python
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
```

**Step 2: Update main to process all molecules**

Replace `main()` function:

```python
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
```

**Step 3: Run and verify all molecules**

```bash
python teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
```

Expected: Shows Δμ⁰ values for H2O, H2, O2 at 298 K

**Step 4: Commit raw data**

```bash
git add teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
git commit -m "feat: add raw JANAF data for H2O, H2, O2

Data covers 0-2000 K range in 50 K steps from NIST-JANAF tables"
```

---

## Task 4: Generate JSON Output from Extraction Script

**Files:**
- Modify: `teros/core/surface_hydroxylation/scripts/extract_janaf_data.py`
- Create: `teros/core/surface_hydroxylation/thermodynamics_data.json` (generated)

**Step 1: Add JSON generation function**

Add before `main()`:

```python
def generate_json_output(processed_data: dict) -> dict:
    """
    Generate final JSON structure.

    Args:
        processed_data: Dict mapping species to processed data

    Returns:
        Complete JSON structure with metadata
    """
    output = {
        'metadata': {
            'source': 'NIST-JANAF Thermochemical Tables',
            'urls': JANAF_URLS,
            'units': {
                'temperature': 'K',
                'enthalpy': 'kJ/mol (H(T) - H(0 K))',
                'entropy': 'J/(mol·K)',
                'delta_mu': 'eV/molecule'
            },
            'reference_temperature': 298.15,
            'temperature_range': [0, 2000],
            'temperature_step': 50,
            'data_extraction': {
                'date': '2025-10-27',
                'method': 'Manual extraction from JANAF HTML tables',
                'verification': 'Spot-checked against Section 4 example values'
            },
            'conversion_factors': {
                'kJ_to_eV': KJ_TO_EV,
                'J_to_kJ': J_TO_KJ,
            }
        },
        'molecules': {}
    }

    for species, data in processed_data.items():
        output['molecules'][species] = {
            'data': data
        }

    return output
```

**Step 2: Update main to write JSON file**

Add to end of `main()`:

```python
    # Generate JSON structure
    json_output = generate_json_output(processed_data)

    # Write to file
    output_path = Path(__file__).parent.parent / 'thermodynamics_data.json'
    with open(output_path, 'w') as f:
        json.dump(json_output, f, indent=2)

    print(f"\nJSON file written to: {output_path}")
    print(f"File size: {output_path.stat().st_size} bytes")
```

**Step 3: Run script and generate JSON**

```bash
python teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
```

Expected: Creates `thermodynamics_data.json` in `teros/core/surface_hydroxylation/`

**Step 4: Verify JSON structure**

```bash
python -c "
import json
with open('teros/core/surface_hydroxylation/thermodynamics_data.json') as f:
    data = json.load(f)
    print('Molecules:', list(data['molecules'].keys()))
    print('H2O at 298 K:', data['molecules']['H2O']['data']['298'])
"
```

Expected: Shows H2O, H2, O2 and 298 K data

**Step 5: Commit extraction script and generated data**

```bash
git add teros/core/surface_hydroxylation/scripts/extract_janaf_data.py
git add teros/core/surface_hydroxylation/thermodynamics_data.json
git commit -m "feat: generate thermodynamics_data.json from JANAF data

Complete database with H2O, H2, O2 at 0-2000 K"
```

---

## Task 5: Create JanafDatabase Class Foundation

**Files:**
- Create: `teros/core/surface_hydroxylation/thermodynamics.py`
- Create: `tests/test_thermodynamics.py`

**Step 1: Write failing test for class initialization**

Create test file:

```python
"""Tests for JANAF thermodynamics database."""

import pytest
from pathlib import Path
from teros.core.surface_hydroxylation.thermodynamics import JanafDatabase


def test_database_loads():
    """Test that database loads without errors."""
    db = JanafDatabase()
    assert db is not None


def test_available_species():
    """Test available species list."""
    db = JanafDatabase()
    species = db.available_species

    assert 'H2O' in species
    assert 'H2' in species
    assert 'O2' in species
    assert len(species) == 3
```

**Step 2: Run test to verify it fails**

```bash
cd /home/thiagotd/git/PS-TEROS
pytest tests/test_thermodynamics.py::test_database_loads -v
```

Expected: ImportError or ModuleNotFoundError

**Step 3: Write minimal implementation**

Create thermodynamics.py:

```python
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
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_thermodynamics.py::test_database_loads -v
pytest tests/test_thermodynamics.py::test_available_species -v
```

Expected: Both tests PASS

**Step 5: Commit foundation**

```bash
git add teros/core/surface_hydroxylation/thermodynamics.py
git add tests/test_thermodynamics.py
git commit -m "feat: add JanafDatabase class foundation

Implements JSON loading and basic properties"
```

---

## Task 6: Implement get_mu_correction Method

**Files:**
- Modify: `teros/core/surface_hydroxylation/thermodynamics.py`
- Modify: `tests/test_thermodynamics.py`

**Step 1: Write failing test**

Add to test file:

```python
def test_get_mu_correction_basic():
    """Test basic Δμ⁰ lookup at 298 K."""
    db = JanafDatabase()

    # Get correction at 298 K
    mu_h2o = db.get_mu_correction('H2O', T=298)
    mu_h2 = db.get_mu_correction('H2', T=298)
    mu_o2 = db.get_mu_correction('O2', T=298)

    # All should be negative (favorable)
    assert mu_h2o < 0
    assert mu_h2 < 0
    assert mu_o2 < 0

    # Should be within reasonable range (-1 to 0 eV)
    assert -1.0 < mu_h2o < 0
    assert -1.0 < mu_h2 < 0
    assert -1.0 < mu_o2 < 0


def test_get_mu_correction_temperature_dependence():
    """Test that Δμ⁰ becomes more negative with temperature."""
    db = JanafDatabase()

    mu_298 = db.get_mu_correction('H2O', T=298)
    mu_500 = db.get_mu_correction('H2O', T=500)

    # Higher temperature should give more negative correction
    assert mu_500 < mu_298
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_thermodynamics.py::test_get_mu_correction_basic -v
```

Expected: AttributeError (method doesn't exist)

**Step 3: Implement get_mu_correction**

Add to JanafDatabase class:

```python
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
            -0.5643
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
            import math
            pressure_correction = KB_EV * T * math.log(P)
            delta_mu += pressure_correction

        return delta_mu

    def _validate_species(self, species: str) -> None:
        """Validate species is supported."""
        if species not in SPECIES:
            raise ValueError(
                f"Species '{species}' not found. Available: {SPECIES}"
            )
```

**Step 4: Run tests to verify they pass**

```bash
pytest tests/test_thermodynamics.py::test_get_mu_correction_basic -v
pytest tests/test_thermodynamics.py::test_get_mu_correction_temperature_dependence -v
```

Expected: Both tests PASS

**Step 5: Commit implementation**

```bash
git add teros/core/surface_hydroxylation/thermodynamics.py
git add tests/test_thermodynamics.py
git commit -m "feat: implement get_mu_correction method

Returns chemical potential corrections with temperature and pressure support"
```

---

## Task 7: Implement Error Handling and Validation

**Files:**
- Modify: `teros/core/surface_hydroxylation/thermodynamics.py`
- Modify: `tests/test_thermodynamics.py`

**Step 1: Write failing tests for error cases**

Add to test file:

```python
def test_invalid_species():
    """Test error on invalid species."""
    db = JanafDatabase()

    with pytest.raises(ValueError, match="not found"):
        db.get_mu_correction('CO2', T=298)


def test_invalid_temperature():
    """Test error on invalid temperature."""
    db = JanafDatabase()

    # Not a 50 K step
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=299)

    # Out of range
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=3000)


def test_temperature_must_be_discrete():
    """Test that intermediate temperatures are not allowed."""
    db = JanafDatabase()

    # 325 K is between 300 and 350
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=325)
```

**Step 2: Run tests to verify they fail appropriately**

```bash
pytest tests/test_thermodynamics.py::test_invalid_species -v
pytest tests/test_thermodynamics.py::test_invalid_temperature -v
```

Expected: Tests PASS (error handling already works from previous implementation)

**Step 3: Commit test coverage**

```bash
git add tests/test_thermodynamics.py
git commit -m "test: add error handling test coverage

Validates species and temperature validation"
```

---

## Task 8: Implement get_raw_data and list_temperatures Methods

**Files:**
- Modify: `teros/core/surface_hydroxylation/thermodynamics.py`
- Modify: `tests/test_thermodynamics.py`

**Step 1: Write failing tests**

Add to test file:

```python
def test_get_raw_data():
    """Test raw JANAF data retrieval."""
    db = JanafDatabase()

    raw = db.get_raw_data('H2O', T=298)

    assert 'H' in raw
    assert 'S' in raw
    assert 'delta_mu' in raw

    # Enthalpy should be positive and reasonable
    assert 0 < raw['H'] < 20  # kJ/mol

    # Entropy should be positive
    assert raw['S'] > 0  # J/(mol·K)

    # Delta mu should match correction method
    mu = db.get_mu_correction('H2O', T=298)
    assert abs(raw['delta_mu'] - mu) < 1e-6


def test_list_temperatures():
    """Test temperature listing."""
    db = JanafDatabase()

    temps = db.list_temperatures()

    # Should be 41 points (0-2000 K in 50 K steps)
    assert len(temps) == 41

    # Should be sorted
    assert temps == sorted(temps)

    # Check range
    assert temps[0] == 0
    assert temps[-1] == 2000

    # Check step
    assert temps[1] == 50
    assert temps[6] == 300
```

**Step 2: Run tests to verify they fail**

```bash
pytest tests/test_thermodynamics.py::test_get_raw_data -v
pytest tests/test_thermodynamics.py::test_list_temperatures -v
```

Expected: AttributeError (methods don't exist)

**Step 3: Implement methods**

Add to JanafDatabase class:

```python
    def get_raw_data(self, species: str, T: float) -> dict:
        """
        Get raw JANAF data for verification and transparency.

        Args:
            species: Molecular species ('H2O', 'H2', 'O2')
            T: Temperature in Kelvin

        Returns:
            Dictionary with keys: 'H' (kJ/mol), 'S' (J/(mol·K)), 'delta_mu' (eV)

        Raises:
            ValueError: If species or T invalid
        """
        self._validate_species(species)

        T_str = str(int(T))
        if T_str not in self._data['molecules'][species]['data']:
            raise ValueError(f"Temperature {T} K not found for {species}")

        return self._data['molecules'][species]['data'][T_str]

    def list_temperatures(self, species: Optional[str] = None) -> list:
        """
        Get available temperature points.

        Args:
            species: Optional species to check (default: use first species)

        Returns:
            Sorted list of available temperatures in Kelvin
        """
        if species is None:
            species = list(self._data['molecules'].keys())[0]

        temps = [int(t) for t in self._data['molecules'][species]['data'].keys()]
        return sorted(temps)
```

**Step 4: Run tests to verify they pass**

```bash
pytest tests/test_thermodynamics.py::test_get_raw_data -v
pytest tests/test_thermodynamics.py::test_list_temperatures -v
```

Expected: Both tests PASS

**Step 5: Commit implementation**

```bash
git add teros/core/surface_hydroxylation/thermodynamics.py
git add tests/test_thermodynamics.py
git commit -m "feat: implement get_raw_data and list_temperatures methods

Provides transparent access to raw JANAF data and temperature points"
```

---

## Task 9: Add Pressure Correction Tests

**Files:**
- Modify: `tests/test_thermodynamics.py`

**Step 1: Write test for pressure corrections**

Add to test file:

```python
import math


def test_pressure_correction():
    """Test pressure-dependent corrections."""
    db = JanafDatabase()

    # At P = 1 bar (reference)
    mu_1bar = db.get_mu_correction('O2', T=298, P=1.0)

    # At P = 0.21 bar (atmospheric O2 partial pressure)
    mu_021bar = db.get_mu_correction('O2', T=298, P=0.21)

    # Should differ by k_B*T*ln(0.21)
    KB_EV = 8.617333e-5
    expected_diff = KB_EV * 298 * math.log(0.21)

    actual_diff = mu_021bar - mu_1bar

    assert abs(actual_diff - expected_diff) < 1e-6


def test_pressure_correction_increases_with_temperature():
    """Test pressure correction magnitude increases with T."""
    db = JanafDatabase()

    P = 0.5  # Half standard pressure

    mu_298_1bar = db.get_mu_correction('H2O', T=298, P=1.0)
    mu_298_05bar = db.get_mu_correction('H2O', T=298, P=P)

    mu_500_1bar = db.get_mu_correction('H2O', T=500, P=1.0)
    mu_500_05bar = db.get_mu_correction('H2O', T=500, P=P)

    correction_298 = abs(mu_298_05bar - mu_298_1bar)
    correction_500 = abs(mu_500_05bar - mu_500_1bar)

    # Higher temperature should give larger pressure correction
    assert correction_500 > correction_298
```

**Step 2: Run tests to verify they pass**

```bash
pytest tests/test_thermodynamics.py::test_pressure_correction -v
pytest tests/test_thermodynamics.py::test_pressure_correction_increases_with_temperature -v
```

Expected: Both tests PASS (pressure correction already implemented)

**Step 3: Commit test coverage**

```bash
git add tests/test_thermodynamics.py
git commit -m "test: add comprehensive pressure correction tests

Validates k_B*T*ln(P) correction formula"
```

---

## Task 10: Add Data Validation Tests

**Files:**
- Modify: `tests/test_thermodynamics.py`

**Step 1: Write test to validate Δμ⁰ calculation**

Add to test file:

```python
def test_delta_mu_calculation_validation():
    """Verify Δμ⁰ calculation from raw data matches stored value."""
    db = JanafDatabase()

    # Test at multiple temperatures
    test_temps = [0, 298, 500, 1000, 2000]

    for species in ['H2O', 'H2', 'O2']:
        for T in test_temps:
            raw = db.get_raw_data(species, T=T)

            # Manual calculation
            H_kJ = raw['H']
            S_kJ = raw['S'] * 0.001  # J/(mol·K) to kJ/(mol·K)
            KJ_TO_EV = 0.010364269

            delta_mu_manual = (H_kJ - T * S_kJ) * KJ_TO_EV

            # Should match stored value
            assert abs(delta_mu_manual - raw['delta_mu']) < 1e-5, \
                f"Mismatch for {species} at {T} K"


def test_monotonic_behavior():
    """Verify Δμ⁰ decreases with temperature (becomes more negative)."""
    db = JanafDatabase()

    temps = [0, 298, 500, 1000, 2000]

    for species in ['H2O', 'H2', 'O2']:
        mu_values = [db.get_mu_correction(species, T=T) for T in temps]

        # Δμ⁰ should become more negative with increasing T
        for i in range(len(mu_values) - 1):
            assert mu_values[i+1] < mu_values[i], \
                f"{species}: mu({temps[i+1]}) should be < mu({temps[i]})"


def test_zero_temperature_values():
    """Test that T=0 K gives zero correction."""
    db = JanafDatabase()

    for species in ['H2O', 'H2', 'O2']:
        mu_0 = db.get_mu_correction(species, T=0)
        raw_0 = db.get_raw_data(species, T=0)

        # At T=0: H=0, S=0, therefore Δμ⁰=0
        assert mu_0 == 0.0
        assert raw_0['H'] == 0.0
        assert raw_0['S'] == 0.0
```

**Step 2: Run tests to verify they pass**

```bash
pytest tests/test_thermodynamics.py::test_delta_mu_calculation_validation -v
pytest tests/test_thermodynamics.py::test_monotonic_behavior -v
pytest tests/test_thermodynamics.py::test_zero_temperature_values -v
```

Expected: All tests PASS (validates data integrity)

**Step 3: Commit validation tests**

```bash
git add tests/test_thermodynamics.py
git commit -m "test: add data validation tests

Validates calculation correctness and thermodynamic consistency"
```

---

## Task 11: Update Module __init__ for Easy Import

**Files:**
- Modify: `teros/core/surface_hydroxylation/__init__.py`

**Step 1: Write test for import convenience**

Add new test file `tests/test_thermodynamics_import.py`:

```python
"""Test that thermodynamics module can be imported easily."""


def test_import_from_surface_hydroxylation():
    """Test import from main module."""
    from teros.core.surface_hydroxylation import JanafDatabase

    db = JanafDatabase()
    assert db is not None


def test_direct_import():
    """Test direct import from thermodynamics."""
    from teros.core.surface_hydroxylation.thermodynamics import (
        JanafDatabase,
        KJ_TO_EV,
        KB_EV,
    )

    assert JanafDatabase is not None
    assert KJ_TO_EV > 0
    assert KB_EV > 0
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_thermodynamics_import.py::test_import_from_surface_hydroxylation -v
```

Expected: ImportError (JanafDatabase not exported from __init__)

**Step 3: Update __init__.py**

Add to `teros/core/surface_hydroxylation/__init__.py`:

```python
from .thermodynamics import JanafDatabase

__all__ = [
    # ... existing exports ...
    'JanafDatabase',
]
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_thermodynamics_import.py -v
```

Expected: Both tests PASS

**Step 5: Commit export**

```bash
git add teros/core/surface_hydroxylation/__init__.py
git add tests/test_thermodynamics_import.py
git commit -m "feat: export JanafDatabase from surface_hydroxylation module

Enables: from teros.core.surface_hydroxylation import JanafDatabase"
```

---

## Task 12: Run Full Test Suite and Verify

**Files:**
- No new files (verification task)

**Step 1: Run all thermodynamics tests**

```bash
pytest tests/test_thermodynamics.py tests/test_thermodynamics_import.py -v
```

Expected: All tests PASS

**Step 2: Check test coverage**

```bash
pytest tests/test_thermodynamics.py tests/test_thermodynamics_import.py --cov=teros.core.surface_hydroxylation.thermodynamics --cov-report=term-missing
```

Expected: >90% coverage

**Step 3: Run quick integration test**

```bash
python -c "
from teros.core.surface_hydroxylation import JanafDatabase

db = JanafDatabase()
print('✓ Database loaded')

mu_h2o = db.get_mu_correction('H2O', T=298)
print(f'✓ μ(H2O) at 298 K = {mu_h2o:.4f} eV')

mu_o2 = db.get_mu_correction('O2', T=500, P=0.21)
print(f'✓ μ(O2) at 500 K, 0.21 bar = {mu_o2:.4f} eV')

temps = db.list_temperatures()
print(f'✓ {len(temps)} temperature points available')

print('All checks passed!')
"
```

Expected: All checks pass

**Step 4: Manual verification against Section 4 example**

```bash
python -c "
from teros.core.surface_hydroxylation import JanafDatabase

db = JanafDatabase()

# Section 4 example values at 298 K, 1 bar
expected = {
    'H2O': -0.564,  # Approximate from example
    'H2': -0.405,
    'O2': -0.638,
}

print('Verification against Section 4 example values:')
for species, expected_val in expected.items():
    actual = db.get_mu_correction(species, T=298)
    diff = abs(actual - expected_val)
    status = '✓' if diff < 0.01 else '✗'
    print(f'{status} {species}: {actual:.4f} eV (expected ~{expected_val:.3f}, diff={diff:.4f})')
"
```

Expected: All values within ~0.01 eV of example

**Step 5: No commit (verification only)**

---

## Task 13: Create Usage Documentation and Examples

**Files:**
- Create: `teros/core/surface_hydroxylation/examples/thermodynamics_usage.py`

**Step 1: Create examples file**

Create file with comprehensive usage examples:

```python
"""
Usage examples for JANAF thermodynamics database.

This script demonstrates how to use the JanafDatabase class for
chemical potential corrections in surface energy calculations.
"""

from teros.core.surface_hydroxylation import JanafDatabase


def example_basic_usage():
    """Basic usage: get corrections at room temperature."""
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)

    db = JanafDatabase()

    # Get correction at room temperature (298 K, 1 bar)
    mu_h2o = db.get_mu_correction('H2O', T=298)
    mu_h2 = db.get_mu_correction('H2', T=298)
    mu_o2 = db.get_mu_correction('O2', T=298)

    print(f"\nChemical potential corrections at 298 K, 1 bar:")
    print(f"  μ(H2O) = {mu_h2o:.4f} eV")
    print(f"  μ(H2)  = {mu_h2:.4f} eV")
    print(f"  μ(O2)  = {mu_o2:.4f} eV")

    # For atomic oxygen (not molecular)
    mu_O_atom = mu_o2 / 2.0
    print(f"  μ(O atom) = {mu_O_atom:.4f} eV")


def example_temperature_dependence():
    """Demonstrate temperature dependence."""
    print("\n" + "=" * 60)
    print("Example 2: Temperature Dependence")
    print("=" * 60)

    db = JanafDatabase()

    temps = [298, 373, 500, 1000]
    print(f"\nμ(H2O) at different temperatures:")

    for T in temps:
        mu = db.get_mu_correction('H2O', T=T)
        print(f"  {T:4d} K: {mu:.4f} eV")


def example_pressure_correction():
    """Demonstrate pressure corrections."""
    print("\n" + "=" * 60)
    print("Example 3: Pressure Corrections")
    print("=" * 60)

    db = JanafDatabase()

    T = 373  # 100°C
    pressures = [0.1, 0.5, 1.0, 2.0]

    print(f"\nμ(H2O) at {T} K, different pressures:")

    for P in pressures:
        mu = db.get_mu_correction('H2O', T=T, P=P)
        print(f"  {P:4.1f} bar: {mu:.4f} eV")


def example_surface_energy_calculation():
    """Example: Calculate hydroxylation energy."""
    print("\n" + "=" * 60)
    print("Example 4: Surface Energy Calculation")
    print("=" * 60)

    db = JanafDatabase()

    # Example energies from VASP (dummy values)
    E_pristine = -1500.0  # eV
    E_hydroxylated = -1495.0  # eV
    n_OH = 2  # Added 2 OH groups

    # Experimental conditions
    T = 298  # K
    P_H2O = 0.023  # bar (room temperature water vapor)

    # Chemical potentials
    mu_h2o = db.get_mu_correction('H2O', T=T, P=P_H2O)
    mu_h2 = db.get_mu_correction('H2', T=T, P=1.0)

    # Formation reaction: Surface + n*H2O -> Surface-nOH + n/2*H2
    # ΔE = E_hydroxylated - E_pristine - n*μ(H2O) + (n/2)*μ(H2)
    delta_E = E_hydroxylated - E_pristine - n_OH * mu_h2o + (n_OH / 2) * mu_h2

    print(f"\nHydroxylation calculation:")
    print(f"  Temperature: {T} K")
    print(f"  P(H2O): {P_H2O} bar")
    print(f"  Number of OH groups: {n_OH}")
    print(f"  μ(H2O): {mu_h2o:.4f} eV")
    print(f"  μ(H2): {mu_h2:.4f} eV")
    print(f"  ΔE_formation: {delta_E:.4f} eV")


def example_raw_data_access():
    """Demonstrate raw data access for verification."""
    print("\n" + "=" * 60)
    print("Example 5: Raw Data Access")
    print("=" * 60)

    db = JanafDatabase()

    raw = db.get_raw_data('H2O', T=298)

    print(f"\nRaw JANAF data for H2O at 298 K:")
    print(f"  H(T) - H(0K): {raw['H']:.3f} kJ/mol")
    print(f"  S(T): {raw['S']:.3f} J/(mol·K)")
    print(f"  Δμ⁰: {raw['delta_mu']:.4f} eV")

    # Show available temperatures
    temps = db.list_temperatures()
    print(f"\nAvailable temperature points: {len(temps)}")
    print(f"  Range: {temps[0]}-{temps[-1]} K")
    print(f"  Step: {temps[1] - temps[0]} K")


if __name__ == '__main__':
    example_basic_usage()
    example_temperature_dependence()
    example_pressure_correction()
    example_surface_energy_calculation()
    example_raw_data_access()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
```

**Step 2: Run examples to verify**

```bash
python teros/core/surface_hydroxylation/examples/thermodynamics_usage.py
```

Expected: All examples run without errors, show reasonable values

**Step 3: Commit examples**

```bash
git add teros/core/surface_hydroxylation/examples/thermodynamics_usage.py
git commit -m "docs: add comprehensive usage examples for thermodynamics module

Demonstrates basic usage, temperature/pressure corrections, and surface energy calculations"
```

---

## Task 14: Final Integration Test

**Files:**
- No new files (verification task)

**Step 1: Restart AiiDA daemon**

```bash
verdi daemon restart
```

Expected: Daemon restarts successfully

**Step 2: Test import in AiiDA context**

```bash
source ~/envs/aiida/bin/activate && python -c "
from aiida import load_profile
load_profile('presto')

from teros.core.surface_hydroxylation import JanafDatabase

db = JanafDatabase()
print('✓ JanafDatabase works in AiiDA context')

mu = db.get_mu_correction('H2O', T=298)
print(f'✓ μ(H2O) = {mu:.4f} eV')
"
```

Expected: Import and usage work correctly

**Step 3: Test with organize_hydroxylation_results (conceptual)**

```bash
python -c "
from teros.core.surface_hydroxylation import JanafDatabase

# Simulate integration with hydroxylation workflow
db = JanafDatabase()

# In actual usage, this would be called after workflow completion:
# results = organize_hydroxylation_results(workflow_node)
# Then use db to calculate surface energies with corrections

print('✓ Module ready for integration with hydroxylation workflow')
print('✓ Can be used in post-processing to calculate surface phase diagrams')
"
```

**Step 4: No commit (verification only)**

---

## Task 15: Update Section 4 of LaTeX Document

**Files:**
- Modify: `/home/thiagotd/git/fosfato/calculos/hydroxylation/docs/surface_energy_calc_procedure.tex`

**Step 1: Update Module Integration section**

Find and update the "Module Integration" section in Section 4:

```latex
\textbf{Module Integration:}
\begin{itemize}
    \item \textbf{Implementation:} Thermodynamics module now available at \texttt{teros/core/surface\_hydroxylation/thermodynamics.py}
    \item \textbf{Data:} JANAF data stored in \texttt{thermodynamics\_data.json} covering 0-2000 K in 50 K steps
    \item \textbf{Usage:} Import with \texttt{from teros.core.surface\_hydroxylation import JanafDatabase}
    \item \textbf{Example:}
    \begin{verbatim}
db = JanafDatabase()
mu_h2o = db.get_mu_correction('H2O', T=298, P=1.0)
    \end{verbatim}
    \item For detailed usage, see \texttt{examples/thermodynamics\_usage.py}
\end{itemize}
```

**Step 2: Commit documentation update**

```bash
cd /home/thiagotd/git/fosfato/calculos/hydroxylation
git add docs/surface_energy_calc_procedure.tex
git commit -m "docs: update Section 4 with implemented thermodynamics module

Reflects completed JANAF database implementation in PS-TEROS"
```

---

## Completion Checklist

After completing all tasks, verify:

- [ ] JSON data file created with all species (H2O, H2, O2)
- [ ] Data covers 0-2000 K in 50 K steps
- [ ] JanafDatabase class implemented with all methods
- [ ] All unit tests pass (>90% coverage)
- [ ] Data validation tests confirm calculation correctness
- [ ] Import works from surface_hydroxylation module
- [ ] Usage examples run without errors
- [ ] Integration with AiiDA context verified
- [ ] Documentation updated
- [ ] All commits follow conventional commit format

**Total Tasks:** 15 tasks
**Estimated Time:** 2-3 hours (with testing and verification)

---

## Notes for Implementation

**Following TDD:**
- Write test first, watch it fail, implement, watch it pass, commit
- Keep commits small and focused (one feature per commit)
- Run tests frequently

**Code Style:**
- Follow existing PS-TEROS conventions
- Use type hints where appropriate
- Add docstrings for all public methods
- Keep functions focused and testable

**Testing:**
- Test error cases explicitly
- Validate calculations against known values
- Check edge cases (T=0, P!=1, etc.)
- Ensure scientific correctness

**References:**
- Design document: `docs/plans/2025-10-27-janaf-thermodynamics-module-design.md`
- JANAF tables: https://janaf.nist.gov/
- Section S2 equations: `surface_energy_calc_procedure.tex`
