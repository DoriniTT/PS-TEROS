# Charge-Neutral Terminations for Semiconductor Surfaces

## Summary

`termination_mode='charge_neutral'` makes PS-TEROS generate only slabs that a
semiconductor or insulator can actually have: two equivalent faces and zero
net formal charge. The default mode (`'pymatgen'`) is unchanged.

## Why

Pymatgen's `SlabGenerator` returns every distinct cut. With
`symmetrize=True` it makes the two faces equivalent by deleting sites one by
one, without regard to valence. For an ionic semiconductor most of those
slabs carry a net formal charge. DFT then puts the surplus electrons in the
conduction band (or holes in the valence band) and the "surface" is a metal.

Ag3PO4(110) is a clear example. With `min_slab_size=10`,
`symmetrize=True` and Ag +1, P +5, O −2, none of the six slabs is neutral
and four split PO4 groups:

| Slab | Formal charge | PO4 intact |
| --- | --- | --- |
| Ag12P6O20 | +2 | no |
| Ag16P6O24 | −2 | yes |
| Ag14P4O20 | −6 | no |
| Ag16P6O20 | +6 | no |
| Ag14P4O16 | +2 | yes |
| Ag18P4O20 | −2 | no |

At 15 Å the list includes Ag20P6O24 (+2), the Ag-rich S0 slab. Its 8 extra
electrons per 2×2 cell sit in the conduction band in CP2K.

## Criterion

A slab is kept if, with nominal oxidation states:

1. its net formal charge is zero (the electron-counting rule), and
2. a symmetry operation maps one face onto the other (no dipole).

Stoichiometry is **not** required. A neutral, off-stoichiometric slab (an
extra Ag2O or P2O5 unit) is closed-shell and is handled by the
chemical-potential dependent surface thermodynamics already in PS-TEROS. A
stoichiometric slab is the special case with zero surface excess.

Formal charge is necessary, not sufficient. Confirm the band gap of each
candidate with a single-point calculation before production work.

## Algorithm

`psteros.core.terminations.find_charge_neutral_terminations`:

1. Assign oxidation states (given, or guessed by pymatgen) and check the bulk
   is neutral.
2. If no bulk operation reverses the (hkl) normal, stop: the direction is
   polar (e.g. zinc blende (111), wurtzite (0001)) and no slab of it has
   equivalent faces.
3. Build a thick oriented slab and split it into units: single atoms, or
   finite clusters joined by `unit_bonds` (P–O for PO4). Units are never
   split.
4. Group units into planes and take every contiguous window of planes.
5. Keep windows with a face-reversing symmetry operation.
6. If a window has formal charge Q ≠ 0, remove surface units in
   symmetry-related pairs until Q = 0:
   - layer by layer: the outermost planes completely, then part of the next
     plane, so no vacancy is left below a kept surface unit;
   - with the fewest units per face;
   - each removal pattern once per symmetry orbit.
7. Keep slabs within one bulk repeat above `min_slab_thickness` and remove
   duplicates: identical slabs, and slabs with the same face environment at
   a larger thickness.

For Ag3PO4(110) in the 1×1 cell this gives two terminations:

| Termination | Built from | Removed per face |
| --- | --- | --- |
| Ag18P6O24 | Ag20P6O24 (+2) | one Ag |
| Ag12P4O16 | Ag12P6O24 (−6) | one PO4 |

Ag18P6O24 has the same layer sequence as an independently relaxed
closed-shell slab (0.14 Å RMS displacement, consistent with relaxation).

## Usage

### In the core workflow

```python
wg = build_core_workgraph(
    ...,
    miller_indices=[1, 1, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    termination_mode='charge_neutral',
    oxidation_states={'Ag': 1, 'P': 5, 'O': -2},
    unit_bonds=[['P', 'O', 1.9]],      # keep PO4 whole
    termination_supercell=[1, 1],      # e.g. [2, 2] for quarter coverages
)
```

The slabs are returned as `term_0`, `term_1`, … like before, so relaxation,
thermodynamics and the other stages work unchanged. `symmetrize` and
`center_slab` are implied in this mode.

### Standalone (no AiiDA)

```python
from pymatgen.core import Structure
from psteros.core.terminations import find_charge_neutral_terminations

bulk = Structure.from_file('ag3po4.cif')
for termination in find_charge_neutral_terminations(
    bulk, (1, 1, 0), 10.0, 15.0,
    oxidation_states={'Ag': 1, 'P': 5, 'O': -2},
    unit_bonds={('P', 'O'): 1.9},
):
    print(termination.to_dict())
    termination.structure.to(filename=f'{termination.formula}.vasp')
```

`classify_slab(slab, bulk, oxidation_states)` reports formal charge, face
symmetry and stoichiometry of any slab, e.g. one supplied via
`input_slabs`.

## Options

| Option | Meaning |
| --- | --- |
| `oxidation_states` | Element → nominal oxidation state. Guessed if omitted. |
| `unit_bonds` | Bonds never broken (polyanions such as PO4, SO4, CO3). Do not list framework cation–anion bonds: they form an extended network and are rejected. |
| `supercell` / `termination_supercell` | In-plane repetition. Needed when a 1×1 cell cannot be neutralised, e.g. Ag3PO4(100) (window charges ±1) needs 1×2. |
| `surface_depth` | How deep the partially emptied plane may lie. Default: half the bulk repeat. |
| `max_removed_per_face` | Default: 4 per in-plane cell. |
| `max_variants_per_window` | Cap on removal patterns per window (default 20). Large supercells with partial coverages can have thousands; a warning is issued when the list is cut. |

## Limits

- One oxidation state per element. Reduced surfaces with mixed valence
  (e.g. SnO2 with surface Sn2+) are reported as charged under fixed Sn4+.
- Units are only removed, never added.
- Covalent semiconductors (Si, GaAs) need dangling-bond counting and surface
  reconstructions, which this mode does not attempt.
- Polar directions are reported, not handled. They need an asymmetric slab
  with a dipole correction or a passivated back face.
