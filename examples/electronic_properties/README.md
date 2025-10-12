# Electronic Properties (DOS and Bands) Example

This example demonstrates how to compute **density of states (DOS)** and **band structure** for relaxed bulk structures in PS-TEROS.

## Overview

The electronic properties feature uses AiiDA-VASP's `vasp.v2.bands` workchain to automatically:
1. Perform an SCF calculation with charge density output
2. Compute band structure along high-symmetry k-paths (using SeeKpath)
3. Compute density of states with the tetrahedron method

## Files

- `bulk_dos_bands_ag2o.py` - Complete working example for Ag2O bulk
- `README.md` - This documentation file

## Features

### Electronic Properties Builder

The `get_electronic_properties_defaults()` builder provides material-agnostic parameters for DOS and band structure calculations:

```python
from teros.core.builders import get_electronic_properties_defaults

ep_defaults = get_electronic_properties_defaults(
    energy_cutoff=500,              # ENCUT (eV)
    electronic_convergence=1e-5,    # EDIFF
    ncore=4,                        # Cores per band
    ispin=2,                        # Spin-polarized
    kpoints_mesh_density=0.3,       # SCF k-mesh density
    band_kpoints_distance=0.2,      # Band path density
    dos_kpoints_distance=0.2,       # DOS k-mesh density
    nedos=2000,                     # DOS grid points
    band_mode="seekpath-aiida",     # Band path algorithm
)
```

### Workflow Integration

Enable electronic properties in any PS-TEROS workflow:

```python
from teros.core.workgraph import build_core_workgraph

wg = build_core_workgraph(
    structures_dir="structures",
    bulk_name="ag2o.cif",
    # ... other bulk parameters ...
    compute_electronic_properties_bulk=True,  # Enable DOS and bands
    bands_parameters=ep_defaults,             # Pass electronic properties params
    band_settings=ep_defaults['band_settings'],  # Pass band workflow settings
    bands_options=bulk_options,               # Scheduler options
)
```

## Key Parameters

### SCF Stage (Critical for Restart)
- `LWAVE=True` - Must be True (bands calculation needs WAVECAR)
- `LCHARG=True` - Must be True (DOS calculation needs CHGCAR)
- `LORBIT=11` - L(m)-decomposed charge densities
- `ISMEAR=0` - Gaussian smearing
- `SIGMA=0.05` - Smearing width (eV)

### Band Structure Stage
- `ISTART=1` - Read from WAVECAR
- `ICHARG=11` - Non-self-consistent calculation
- `ISMEAR=0` - Gaussian smearing
- `SIGMA=0.01` - Small smearing for smooth bands
- `LREAL=False` - Reciprocal space (more accurate)

### DOS Stage
- `ISTART=1` - Read from WAVECAR
- `ICHARG=11` - Non-self-consistent calculation
- `ISMEAR=-5` - Tetrahedron method with Blöchl corrections
- `NEDOS=2000` - Number of DOS grid points
- `LREAL=False` - Reciprocal space

### Band Workflow Settings
- `band_mode` - Algorithm for k-path: "seekpath-aiida", "pymatgen", "latimer-munro", "bradcrack"
- `symprec` - Symmetry precision for path determination (default: 1e-4)
- `band_kpoints_distance` - Spacing between k-points on band path (default: 0.2)
- `line_density` - Points per Angstrom along high-symmetry lines (default: 0.2)
- `run_dos` - Always True (compute DOS with bands)
- `dos_kpoints_distance` - K-mesh density for DOS (default: 0.2)

## Usage

### Prerequisites

1. AiiDA profile configured and daemon running:
```bash
verdi profile set-default psteros
verdi status
verdi daemon start
```

2. Structure files in `examples/structures/`:
   - `ag2o.cif` - Bulk Ag2O
   - `Ag.cif` - Metal reference
   - `O2.cif` - Oxygen molecule

3. VASP code configured in AiiDA:
```bash
verdi code list
# Should show: VASP-VTST-6.4.3@bohr
```

### Running the Example

```bash
# Activate environment
source ~/envs/psteros/bin/activate

# Clear cache (if you made code changes)
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
verdi daemon restart

# Run the example
cd /home/thiagotd/git/PS-TEROS/examples/electronic_properties
python bulk_dos_bands_ag2o.py
```

### Monitoring Progress

```bash
# Get the workflow PK from the script output, then:
verdi process show <PK>
verdi process report <PK>

# List all child processes
verdi process list -a -p <PK>
```

### Expected Runtime

- Bulk relaxation: ~5-10 minutes
- SCF calculation: ~3-5 minutes
- Band structure: ~2-3 minutes
- DOS calculation: ~2-3 minutes
- **Total: ~15-25 minutes** (depends on system size and resources)

## Accessing Results

### After Workflow Completes

```python
from aiida import load_node

# Load the workflow
wg = load_node(<PK>)

# Check if completed successfully
print(wg.is_finished_ok)  # Should be True

# Get relaxed bulk structure
relaxed_bulk = wg.outputs.bulk_structure
bulk_energy = wg.outputs.bulk_energy.value  # eV

# Get band structure (BandsData)
bands = wg.outputs.bulk_bands
print(f"Number of bands: {len(bands.get_bands())}")
print(f"Number of kpoints: {bands.get_kpoints().shape[0]}")

# Get DOS (ArrayData)
dos = wg.outputs.bulk_dos
dos_dict = dos.get_dict()
energies = dos_dict['energy']  # Energy grid (eV)
total_dos = dos_dict['dos']    # Total DOS (states/eV)

# Get additional outputs
misc = wg.outputs.bulk_electronic_properties_misc
```

### Plotting Results

#### Band Structure

```python
import matplotlib.pyplot as plt
from aiida.tools.visualization import plot_bands

# Simple plot
plot_bands(bands)
plt.savefig('ag2o_bands.png')

# Advanced: Extract data for custom plotting
bands_array = bands.get_bands()  # Shape: (num_kpoints, num_bands)
kpoints = bands.get_kpoints()

plt.figure(figsize=(8, 6))
for i in range(bands_array.shape[1]):
    plt.plot(range(len(kpoints)), bands_array[:, i], 'b-', linewidth=0.5)
plt.xlabel('K-point index')
plt.ylabel('Energy (eV)')
plt.title('Ag2O Band Structure')
plt.savefig('ag2o_bands_custom.png')
```

#### Density of States

```python
import matplotlib.pyplot as plt

dos = wg.outputs.bulk_dos
dos_dict = dos.get_dict()

energies = dos_dict['energy']
total_dos = dos_dict['dos']

plt.figure(figsize=(8, 6))
plt.plot(total_dos, energies, 'b-', linewidth=1.5)
plt.xlabel('DOS (states/eV)')
plt.ylabel('Energy (eV)')
plt.title('Ag2O Density of States')
plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)  # Fermi level
plt.grid(True, alpha=0.3)
plt.savefig('ag2o_dos.png')
```

### Exporting Structures

```python
# Export relaxed bulk to file
relaxed_bulk = wg.outputs.bulk_structure
atoms = relaxed_bulk.get_ase()
atoms.write('ag2o_bulk_relaxed.cif')
atoms.write('ag2o_bulk_relaxed.vasp', format='vasp', direct=True, vasp5=True)
```

## Customization

### Changing DOS Resolution

```python
ep_defaults = get_electronic_properties_defaults(
    nedos=3000,  # More points for finer DOS grid
)
```

### Denser K-Point Meshes

```python
ep_defaults = get_electronic_properties_defaults(
    kpoints_mesh_density=0.2,    # Denser SCF mesh
    dos_kpoints_distance=0.15,   # Denser DOS mesh
)
```

### Different Band Path Algorithm

```python
ep_defaults = get_electronic_properties_defaults(
    band_mode="pymatgen",  # Use pymatgen instead of seekpath
)
```

### Custom INCAR Parameters

```python
ep_defaults = get_electronic_properties_defaults()

# Modify specific parameters
ep_defaults['scf']['NELM'] = 200  # More SCF iterations
ep_defaults['bands']['SIGMA'] = 0.005  # Smaller smearing
ep_defaults['dos']['NEDOS'] = 5000  # Very fine DOS grid
```

## Troubleshooting

### Common Issues

**1. "BandsWorkChain not found"**
- Make sure aiida-vasp is installed: `pip show aiida-vasp`
- Check plugin entry points: `verdi plugin list aiida.workflows`

**2. "SCF did not converge"**
- Increase NELM: `ep_defaults['scf']['NELM'] = 200`
- Use denser k-mesh: `kpoints_mesh_density=0.2`
- Check structure quality (overlapping atoms, etc.)

**3. "Bands calculation failed"**
- Verify LWAVE=True and LCHARG=True in SCF stage
- Check that SCF completed successfully
- Look at detailed error: `verdi process report <PK>`

**4. "DOS looks noisy"**
- Increase NEDOS: `nedos=3000` or higher
- Use denser k-mesh: `dos_kpoints_distance=0.15`
- Verify ISMEAR=-5 (tetrahedron method)

### Debugging

```bash
# Check workflow status
verdi process show <PK>

# Get detailed error report
verdi process report <PK>

# Check specific calculation
verdi process list -a -p <PK>
verdi calcjob outputcat <CALCJOB_PK>

# View daemon log
verdi daemon logshow
```

## Implementation Details

### Workflow Structure

```
PS-TEROS Workflow
├── Bulk Relaxation (VaspWorkChain)
│   └── Outputs: relaxed structure, energy
│
├── Reference Calculations (Metal, O2)
│   └── Outputs: reference energies
│
├── Formation Enthalpy
│   └── Outputs: ΔH_f
│
└── Electronic Properties (vasp.v2.bands) ← NEW!
    ├── SCF (internal)
    │   └── LWAVE=True, LCHARG=True
    ├── Bands (internal)
    │   └── Reads WAVECAR, computes along k-path
    └── DOS (internal)
        └── Reads CHGCAR, tetrahedron method
```

### File Organization

The electronic properties feature consists of:

1. **Builder**: `teros/core/builders/electronic_properties_builder.py`
   - Material-agnostic parameter defaults
   - Returns structured dict for vasp.v2.bands

2. **Core Integration**: `teros/core/workgraph.py`
   - `core_workgraph()` - Adds electronic properties task
   - `build_core_workgraph()` - Manual task setup with namespace handling
   - `build_core_workgraph_with_map()` - Parameter forwarding

3. **Example**: `examples/electronic_properties/`
   - Working example script
   - Documentation

## Future Extensions

Potential enhancements for this feature:

1. **Slab Electronic Properties**
   - Add similar functionality for slab structures
   - Useful for surface state analysis

2. **Projected DOS**
   - Extract atom/orbital-projected DOS from LORBIT=11 output
   - Visualize contribution of each element

3. **Band Gap Analysis**
   - Automatic detection of direct/indirect gaps
   - Valence band maximum and conduction band minimum

4. **Fermi Surface**
   - Add support for Fermi surface calculations
   - Requires additional k-point sampling

5. **Hybrid Functionals**
   - Support for HSE06, PBE0 calculations
   - Requires `hybrid_reuse_wavecar=True` setting

## References

- [AiiDA-VASP Bands Workchain Documentation](https://aiida-vasp.readthedocs.io/en/latest/workchains/bands.html)
- [SeeKpath for k-path determination](https://seekpath.readthedocs.io/)
- [VASP Manual - Band Structure](https://www.vasp.at/wiki/index.php/Band-structure_calculation)
- [VASP Manual - Density of States](https://www.vasp.at/wiki/index.php/Density_of_states)

## Support

For issues or questions:
- Check the main PS-TEROS documentation
- Review AiiDA-VASP documentation
- Check `verdi process report <PK>` for detailed error messages
