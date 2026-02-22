# COHP Analysis with LOBSTER

This directory contains an example of using the COHP (Crystal Orbital Hamilton Population) brick for bonding analysis via LOBSTER.

## What is COHP?

COHP (Crystal Orbital Hamilton Population) is a bonding analysis method that projects plane-wave DFT results onto local atomic orbitals. It provides insight into:
- Chemical bonding between atoms
- Bonding vs antibonding character at different energies
- Integrated COHP (ICOHP) values quantifying bond strength

COHP analysis is performed by LOBSTER, a post-processing tool for VASP (and other codes).

## Requirements

### Software
- **LOBSTER binary**: Install to `~/.local/bin/lobster` or add to PATH
  - Download from: https://www.cohp.de/
  - Requires registration (free for academic use)
- **VASP**: For the SCF calculation
- **AiiDA** with PS-TEROS installed

### VASP Prerequisites
The SCF calculation MUST have:
- `ISYM = -1` (no symmetry reduction)
- `LWAVE = True` (write WAVECAR)
- Dense k-point mesh (recommended: 0.02-0.03 Å⁻¹)
- Retrieved files: WAVECAR, CONTCAR, POTCAR, INCAR, DOSCAR

## Example Workflow

The example script `run_cohp_sno2.py` demonstrates a two-stage workflow:

1. **SCF stage**: Static DFT calculation with LOBSTER prerequisites
2. **COHP stage**: LOBSTER bonding analysis

### Stage Configuration

```python
stages = [
    {
        'name': 'scf',
        'type': 'vasp',
        'incar': {
            'encut': 520,
            'ediff': 1e-6,
            'ismear': 0,
            'sigma': 0.05,
            'ibrion': -1,
            'nsw': 0,
            'isym': -1,      # Required for LOBSTER
            'lwave': True,   # Required for LOBSTER
            'prec': 'Accurate',
            'lreal': 'Auto',
        },
        'kpoints_spacing': 0.02,  # Dense mesh recommended
        'retrieve': ['WAVECAR', 'CONTCAR', 'POTCAR', 'INCAR', 'DOSCAR'],
    },
    {
        'name': 'cohp',
        'type': 'cohp',
        'cohp_from': 'scf',
        'cohp_start_energy': -15.0,   # eV below Fermi level
        'cohp_end_energy': 10.0,      # eV above Fermi level
        'basis_set': 'pbeVaspFit2015', # LOBSTER basis set
        'bond_cutoff': 3.0,           # Å, max bond distance
        'save_projection': True,
        'compute_coop': True,
    },
]
```

### COHP Stage Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `cohp_from` | Yes | — | Name of SCF stage with VASP outputs |
| `cohp_start_energy` | No | -15.0 | Energy window start (eV) |
| `cohp_end_energy` | No | 10.0 | Energy window end (eV) |
| `basis_set` | No | `'pbeVaspFit2015'` | LOBSTER basis set name |
| `bond_cutoff` | No | 3.0 | Maximum bond distance (Å) |
| `save_projection` | No | True | Save projection data |
| `compute_coop` | No | True | Compute COOP in addition to COHP |

### Basis Sets

Common LOBSTER basis sets:
- `pbeVaspFit2015` - PBE functional (VASP 5.x)
- `pbeVaspFit2021` - Updated PBE basis
- `LDAVaspFit2015` - LDA functional

## Running the Example

1. Install LOBSTER binary:
   ```bash
   # Download from https://www.cohp.de/
   mv lobster-x.x.x ~/.local/bin/lobster
   chmod +x ~/.local/bin/lobster
   ```

2. Activate AiiDA environment:
   ```bash
   source ~/envs/aiida/bin/activate
   verdi daemon restart
   ```

3. Submit workflow:
   ```bash
   python examples/lego/cohp/run_cohp_sno2.py
   ```

4. Monitor progress:
   ```bash
   verdi process show <PK>
   verdi process report <PK>
   ```

5. Get results:
   ```python
   from teros.core.lego import print_sequential_results
   print_sequential_results(result)
   ```

## Output Files

The COHP stage produces:

1. **cohp_data** (Dict): Parsed ICOHPLIST.lobster with bond information
   - List of bonds with atom pairs, distances, ICOHP values
   - Total bond count

2. **icohp** (SinglefileData): ICOHPLIST.lobster file
   - Integrated COHP values for each bond
   - Text format, human-readable

3. **cohpcar** (SinglefileData): COHPCAR.lobster file
   - Energy-resolved COHP data
   - Used for plotting COHP curves

4. **doscar** (SinglefileData): DOSCAR.lobster file
   - DOS from LOBSTER projection
   - Can differ from VASP DOS

## Interpreting Results

The `print_stage_results()` function displays:
- Total number of bonds analyzed
- Top 5 strongest bonding interactions (by |ICOHP|)

Example output:
```
[2] cohp (COHP)
    Bonds analyzed: 48
    Strongest bonding interactions (top 5 by |ICOHP|):
      1. Sn-O: d=2.089 Å, ICOHP=-4.2315 eV
      2. Sn-O: d=2.089 Å, ICOHP=-4.2315 eV
      3. Sn-O: d=2.089 Å, ICOHP=-4.1892 eV
      4. Sn-O: d=2.089 Å, ICOHP=-4.1892 eV
      5. O-O: d=2.654 Å, ICOHP=-0.5621 eV
```

### ICOHP Interpretation
- **Negative ICOHP**: Bonding interaction
- **Positive ICOHP**: Antibonding interaction
- **Magnitude**: Larger |ICOHP| = stronger interaction
- Typical ranges:
  - Strong covalent: -5 to -2 eV
  - Moderate covalent: -2 to -0.5 eV
  - Weak interaction: -0.5 to 0 eV

## Common Issues

### LOBSTER not found
```
FileNotFoundError: LOBSTER binary not found
```
**Solution**: Install LOBSTER to `~/.local/bin/lobster` or add to PATH

### WAVECAR not retrieved
```
FileNotFoundError: WAVECAR not found in retrieved FolderData
```
**Solution**: Add `'WAVECAR'` to the SCF stage's `retrieve` list

### Missing ISYM=-1
```
Warning: Stage 'scf' missing prerequisite: INCAR key 'isym' should be -1
```
**Solution**: Add `'isym': -1` to SCF stage INCAR

### LOBSTER crashes
Check:
- POTCAR is compatible with basis set
- k-point mesh is not too coarse
- WAVECAR is from a converged calculation

## Advanced Usage

### Custom lobsterin Settings

Add custom LOBSTER directives:
```python
{
    'name': 'cohp',
    'type': 'cohp',
    'cohp_from': 'scf',
    'custom_lines': [
        'userecommendedbasisfunctions',
        'skipdos',
        'skipCOOP',
    ],
}
```

### Extracting Specific Bonds

After the calculation completes, parse `cohp_data` for specific bonds:
```python
from teros.core.lego import get_stage_results

cohp_result = get_stage_results(result, 'cohp')
cohp_data = cohp_result['cohp_data']

# Find all Sn-O bonds
sn_o_bonds = [
    bond for bond in cohp_data['bonds']
    if (bond['atom1'].startswith('Sn') and bond['atom2'].startswith('O'))
    or (bond['atom1'].startswith('O') and bond['atom2'].startswith('Sn'))
]

for bond in sn_o_bonds:
    print(f"{bond['atom1']}-{bond['atom2']}: {bond['icohp']:.4f} eV")
```

## References

1. Dronskowski, R. & Bloechl, P. E. *Crystal orbital Hamilton populations (COHP): energy-resolved visualization of chemical bonding in solids based on density-functional calculations.* J. Phys. Chem. **97**, 8617–8624 (1993).

2. Deringer, V. L., Tchougréeff, A. L. & Dronskowski, R. *Crystal orbital Hamilton population (COHP) analysis as projected from plane-wave basis sets.* J. Phys. Chem. A **115**, 5461–5466 (2011).

3. Maintz, S. et al. *LOBSTER: A tool to extract chemical bonding from plane-wave based DFT.* J. Comput. Chem. **37**, 1030–1035 (2016).

4. LOBSTER website: https://www.cohp.de/
