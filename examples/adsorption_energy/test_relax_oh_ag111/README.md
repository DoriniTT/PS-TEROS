# Test: Adsorption Energy with Relaxation

This test validates the new relaxation workflow for adsorption energy calculations.

## System

- **Substrate:** Ag(111) surface (4 layers, 2x2 supercell)
- **Adsorbate:** OH radical
- **Adsorption site:** Top site (directly above Ag atom)

## Workflow

**Phase 1: Relaxation**
- Relax complete Ag+OH system
- Plugin: `vasp.v2.relax`
- Settings: NSW=50, IBRION=2, ISIF=2

**Phase 2: Separation**
- Separate relaxed structure using connectivity analysis
- Components: Ag substrate, OH molecule, complete system

**Phase 3: SCF Calculations**
- Single-point calculations for all three components
- Plugin: `vasp.v2.vasp`
- Settings: NSW=0 (enforced), IBRION=-1 (enforced)

## Running

```bash
source ~/envs/aiida/bin/activate
python run_relax_adsorption.py
```

## Expected Results

- Relaxation should converge in < 50 steps
- Relaxed structure should show OH bond length ~0.97 Å
- Adsorption energy: E_ads ≈ -2.5 to -3.5 eV (exothermic)

## Validation

Check that:
1. Relaxation completed successfully
2. Structure separation identified OH correctly
3. All 3 SCF calculations completed
4. Adsorption energy is negative (favorable)

Access outputs:
```python
wg = load_node(PK)
relaxed = wg.outputs.relaxed_complete_structures['oh_ag']
E_ads = wg.outputs.adsorption_energies['oh_ag']
```
