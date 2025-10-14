# Complete PS-TEROS Workflow Example

This example demonstrates **ALL features** of the PS-TEROS `build_core_workgraph` function in a single comprehensive workflow.

## What This Example Tests

This workflow exercises every capability of the PS-TEROS core workgraph:

### 1. **Bulk and Reference Relaxations** ✓
- Relaxes the bulk structure (Ag₃PO₄)
- Relaxes all reference structures in parallel:
  - Metal reference (Ag)
  - Nonmetal reference (P)
  - Oxygen reference (O₂)
- Each with optimized VASP parameters for their material type

### 2. **Formation Enthalpy Calculation** ✓
- Calculates ΔH_f of Ag₃PO₄ using the relaxed energies
- Uses proper stoichiometry for ternary oxide

### 3. **Slab Generation** ✓
- Generates slab structures from the relaxed bulk
- Miller index: (1,0,0)
- Minimum slab thickness: 15 Å
- Minimum vacuum thickness: 15 Å
- Creates multiple terminations

### 4. **Slab Relaxation with Unrelaxed SCF** ✓
- Performs SCF calculation on unrelaxed slabs
- Performs full relaxation on slabs
- Both calculations run in parallel for all terminations

### 5. **Relaxation Energy Calculation** ✓
- Computes E_relaxed - E_unrelaxed for each slab
- Quantifies the energy gained by relaxation
- **Flag**: `compute_relaxation_energy=True` (default)

### 6. **Cleavage Energy Calculation** ✓
- Calculates cleavage energies for complementary slab pairs
- Important for understanding surface formation
- **Flag**: `compute_cleavage=True` (default)

### 7. **Surface Thermodynamics** ✓
- Identifies oxide type (binary/ternary, metal/nonmetal identity)
- Performs chemical potential sampling:
  - μ_O range: from metal-rich to oxygen-rich conditions
  - μ_P range: from P-rich to P-poor conditions
  - Grid sampling: 50×50 points
- Calculates surface energies γ(μ_O, μ_P) for each termination
- **Flag**: `compute_thermodynamics=True` (default)

## Material System

**Ag₃PO₄** - Silver phosphate (ternary oxide)
- **Bulk**: Ag₃PO₄ (semiconductor, photocatalyst)
- **Metal reference**: Ag (FCC metal)
- **Nonmetal reference**: P (black phosphorus)
- **Oxygen reference**: O₂ (molecule)
- **Surface**: (1,0,0) terminations

## Files

```
complete/
├── README.md              # This file
├── complete_example.py    # Main executable script
└── structures/            # Structure files
    ├── ag3po4.cif        # Bulk structure
    ├── Ag.cif            # Metal reference
    ├── P.cif             # Nonmetal reference
    └── O2.cif            # Oxygen reference
```

## Usage

### Prerequisites

1. **AiiDA profile**: Make sure 'psteros' profile is active
   ```bash
   verdi profile set-default psteros
   verdi status
   ```

2. **Daemon running**: Start if not already running
   ```bash
   verdi daemon start
   ```

3. **Code configured**: VASP code must be set up
   ```bash
   verdi code show VASP-VTST-6.4.3@bohr
   ```

### Running the Example

```bash
# Activate environment
source ~/envs/aiida/bin/activate

# Navigate to example directory
cd /home/thiagotd/git/worktree/PS-TEROS/feature-relax-energy/examples/complete

# Run the complete example
python complete_example.py
```

### Expected Output

The script will print detailed information about:
1. Structure files being used
2. Parameters for each calculation type
3. Calculation flags (all enabled)
4. WorkGraph structure
5. Submission confirmation
6. Monitoring instructions

A file `ag3po4_complete_workgraph.html` will be created showing the full workflow visualization.

## Monitoring the Workflow

### Check Status
```bash
verdi process show <PK>
verdi process report <PK>
```

### List All Processes
```bash
verdi process list
verdi process list -a -p1  # Detailed view
```

### Monitor Real-Time
```bash
watch -n 5 'verdi process list'
```

## Expected Outputs

Upon successful completion, the workflow produces:

### Energies
- `bulk_energy`: Total energy of relaxed Ag₃PO₄ (Float)
- `metal_energy`: Total energy of relaxed Ag (Float)
- `nonmetal_energy`: Total energy of relaxed P (Float)
- `oxygen_energy`: Total energy of relaxed O₂ (Float)

### Structures
- `bulk_structure`: Relaxed Ag₃PO₄ structure (StructureData)
- `metal_structure`: Relaxed Ag structure (StructureData)
- `nonmetal_structure`: Relaxed P structure (StructureData)
- `oxygen_structure`: Relaxed O₂ structure (StructureData)

### Formation Enthalpy
- `formation_enthalpy`: ΔH_f in eV/formula unit (Dict)
  - Contains: `delta_hf`, stoichiometry, and component energies

### Slab Data
- `slab_structures`: All generated (1,0,0) terminations (namespace)
- `unrelaxed_slab_energies`: SCF energies for each slab (namespace)
- `slab_energies`: Relaxed energies for each slab (namespace)
- `relaxed_slabs`: Relaxed slab structures (namespace)

### Relaxation Energies
- `relaxation_energies`: E_relax - E_unrelaxed for each slab (namespace)
  - Shows how much energy is gained by relaxation
  - Larger values indicate more surface reconstruction

### Cleavage Energies
- `cleavage_energies`: For complementary termination pairs (namespace)
  - E_cleavage = (E_slab_A + E_slab_B - E_bulk) / (2 × Area)
  - Measures energy cost to cleave the bulk

### Surface Thermodynamics
- `surface_energies`: γ(μ_O, μ_P) for each termination (namespace)
  - Chemical potential dependent surface energies
  - Includes stability regions
  - Format: arrays with shape (n_samples_O, n_samples_P)

## Calculation Flags

All flags can be controlled in the script:

```python
relax_slabs = True                      # Enable slab relaxation
compute_relaxation_energy = True        # Calculate relaxation energies (default: True)
compute_cleavage = True                 # Calculate cleavage energies (default: True)
compute_thermodynamics = True           # Calculate surface energies (default: True)
thermodynamics_sampling = 50            # Grid resolution for μ sampling
```

To disable any feature, set its flag to `False`.

## Computational Cost

This is a **full production calculation** that will:
- Run 4 bulk/reference relaxations in parallel (~1-4 hours each)
- Generate multiple slab terminations
- Run SCF + relaxation for each slab (~2-6 hours each)
- Perform thermodynamics calculations

**Total estimated time**: Several hours to a day, depending on:
- System size
- Number of slab terminations generated
- Cluster load and resources

## Comparison with Simple Example

| Feature | Simple (relaxation.py) | Complete (this example) |
|---------|----------------------|-------------------------|
| Bulk relaxation | ✓ | ✓ |
| Reference relaxations | ✗ | ✓ |
| Formation enthalpy | ✗ | ✓ |
| Slab generation | ✗ | ✓ |
| Slab relaxation | ✗ | ✓ |
| Relaxation energies | ✗ | ✓ |
| Cleavage energies | ✗ | ✓ |
| Surface thermodynamics | ✗ | ✓ |
| Estimated time | ~1 hour | Several hours |

## Tips

1. **Test with a subset first**: Reduce `thermodynamics_sampling` to 10 for faster testing
2. **Check intermediate results**: Monitor individual VASP calculations
3. **Restart capability**: If calculation fails, use `restart_from_node` parameter
4. **Visualization**: Open the HTML file to understand the workflow structure

## Troubleshooting

### Daemon not running
```bash
verdi daemon start
```

### Profile not set
```bash
verdi profile set-default psteros
```

### Python cache issues
```bash
find . -type d -name __pycache__ -exec rm -rf {} +
find . -name "*.pyc" -delete
verdi daemon restart
```

### Check calculation details
```bash
verdi process show <PK>
verdi calcjob outputcat <PK>  # For VASP calculations
```

## Further Information

- See `examples/relaxation/` for a minimal bulk-only example
- See `examples/formation/` for formation enthalpy without slabs
- See PS-TEROS documentation for detailed API reference

## Citation

If you use this workflow in your research, please cite:
[Citation information to be added]
