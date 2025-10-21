# OOH Surface Modification Tool

## Overview

This tool modifies a LaMnO3 surface slab structure by adding an OOH radical to the most exposed oxygen atom coordinated with the topmost surface Mn atom.

## Usage

```bash
python add_ooh_to_surface.py
```

## Input

- `LaMnO3_100_A4_surface.cif`: Input surface slab structure

## Output

- `LaMnO3_100_A4_surface_OOH.cif`: Modified structure with OOH radical

## Methodology

1. Identify topmost Mn atom (highest z-coordinate)
2. Find O atoms coordinated to this Mn (within 2.5 Angstrom)
3. Select most exposed O (highest z-coordinate)
4. Add OOH radical:
   - O-O bond: 1.45 Angstrom
   - O-H bond: 0.97 Angstrom
   - Orientation: perpendicular to surface (+z direction)

## Dependencies

- pymatgen
- ASE
- numpy

## Testing

Run all tests:
```bash
pytest test_*.py -v
```

## Module Structure

The implementation is split into modular components:

- `structure_loader.py`: Load and validate CIF structures
- `surface_detector.py`: Identify topmost surface Mn atom
- `coordination_finder.py`: Find O atoms coordinated to Mn
- `oxygen_selector.py`: Select most exposed O atom
- `ooh_constructor.py`: Construct OOH radical geometry
- `structure_exporter.py`: Export modified structure to CIF
- `add_ooh_to_surface.py`: Main script orchestrating the workflow

Each module has corresponding tests in `test_*.py` files.

## Example Output

```
============================================================
OOH Surface Modification Script
============================================================

Loading structure from: LaMnO3_100_A4_surface.cif
✓ Structure loaded successfully
  Total atoms: 37
  Elements: O, Mn, La
  Mn atoms: 7
  O atoms: 22

Finding surface Mn atom...
✓ Surface Mn identified
  Index: 14
  Position (Cartesian): [2.3807, 1.3922, 15.2740]

Finding coordinated O atoms...
✓ Found 4 coordinated O atoms
  O1: index=35, distance=1.881 Angstrom
  O2: index=34, distance=1.910 Angstrom
  O3: index=36, distance=1.987 Angstrom
  O4: index=29, distance=2.044 Angstrom

Selecting most exposed O atom...
✓ Most exposed O selected
  Index: 35
  Position (Cartesian): [3.9085, 0.8299, 16.2157]
  Z-coordinate: 16.2157 Angstrom

Constructing OOH radical...
✓ OOH radical constructed
  Original atoms: 37
  Modified atoms: 39
  Added: 1 O + 1 H
  O-O bond length: 1.450 Angstrom
  O-H bond length: 0.970 Angstrom

Exporting modified structure to: LaMnO3_100_A4_surface_OOH.cif
✓ Structure exported successfully

============================================================
OOH modification completed successfully!
============================================================
```
