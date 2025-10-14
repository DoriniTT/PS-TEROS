# Input Structures Directory

Place your pre-generated slab structure files here.

## Supported Formats

- CIF (`.cif`)
- POSCAR/CONTCAR (VASP format)
- XYZ (`.xyz`)
- Any format supported by ASE

## Naming Convention

Name your files sequentially for easy loading:
- `slab_term_0.cif`
- `slab_term_1.cif`
- `slab_term_2.cif`
- etc.

Or use descriptive names:
- `ag3po4_100_term_0.cif`
- `fe2wo6_110_oxygen_rich.cif`
- etc.

## Structure Requirements

Your slab structures should:
1. Have sufficient vacuum spacing (typically 10-15 Å)
2. Have the c-axis perpendicular to the surface
3. Be centered in the c direction (if using `center_slab` behavior)
4. Have appropriate cell parameters for DFT calculations

## Example

See the parent directory's `slabs_input_relax.py` for an example of how to load and use these structures.
