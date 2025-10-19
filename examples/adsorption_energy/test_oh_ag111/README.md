# Adsorption Energy Example: OH on Ag(111)

## Overview

This example demonstrates how to calculate adsorption energies using the PSTEROS adsorption energy module.

## Workflow

1. **Input**: Substrate+adsorbate structure (Ag slab with OH)
2. **Separation**: Automatically separates into:
   - Bare Ag substrate
   - Isolated OH molecule
   - Complete Ag+OH system
3. **Relaxation**: VASP relaxes all three structures in parallel
4. **Calculation**: Computes E_ads = E_complete - E_substrate - E_molecule

## Usage

### Prerequisites

1. Configure AiiDA VASP code:
   ```bash
   verdi code create core.code.installed
   ```

2. Set up pseudopotentials:
   ```bash
   verdi data core.upf uploadfamily --path /path/to/potentials --name PBE
   ```

### Run the workflow

```bash
source ~/envs/aiida/bin/activate
cd examples/adsorption_energy/test_oh_ag111
python run_adsorption_energy.py
```

### Monitor progress

```bash
verdi process list
verdi process show <PK>
verdi process report <PK>
```

## Expected Results

- **Substrate energy**: Total energy of bare Ag(111) slab
- **Molecule energy**: Total energy of isolated OH in same cell
- **Complete energy**: Total energy of Ag(111)+OH system
- **Adsorption energy**: E_ads (negative = favorable adsorption)

## Notes

- Adjust VASP parameters for production calculations
- Use appropriate k-point sampling for your system
- Check convergence with respect to slab thickness and vacuum
- This example uses minimal parameters for fast testing
