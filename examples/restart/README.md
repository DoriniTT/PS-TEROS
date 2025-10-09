# PS-TEROS Restart Feature Documentation

## Overview

The PS-TEROS restart feature allows you to continue slab relaxation calculations from a previous run that failed, didn't converge, or needs different parameters. This is particularly useful when:

- **Calculations hit time/resource limits**: NSW wasn't enough, walltime exceeded
- **Convergence too strict**: EDIFFG was too tight for the system
- **Algorithm changes needed**: Switching IBRION, adjusting POTIM
- **Iterative refinement**: Progressively tightening convergence criteria
- **VASP errors**: Recoverable errors that can be fixed by continuing

## Key Features

✅ **Full Feature Parity**: Restart mode supports all features (thermodynamics, cleavage energies, etc.)  
✅ **Automatic Data Extraction**: Automatically loads structures and RemoteData from previous runs  
✅ **VASP Restart**: Uses WAVECAR and CONTCAR from previous calculations  
✅ **Iterative Restart**: Can chain restarts indefinitely  
✅ **Complete Provenance**: All outputs properly tracked in AiiDA database  
✅ **New RemoteData**: Creates new RemoteData for potential future restarts

## Quick Start

### 1. Run Initial Calculation

First, run a PS-TEROS calculation that may fail or not converge:

```bash
cd examples/restart
source ~/envs/psteros/bin/activate
python slabs_relax_ag2o_restart.py
```

Note the PK of the workgraph (e.g., 22223).

### 2. Check if Restart is Possible

```bash
verdi process show 22223
```

Look for `slab_remote` outputs:
```
Outputs             PK     Type
------------------  -----  -------------
slab_remote
    term_0          22310  RemoteData
    term_1          22311  RemoteData
```

If you see `slab_remote` outputs, you can restart from this run!

### 3. Restart the Calculation

Edit `restart_example.py` to set the previous run PK:

```python
PREVIOUS_RUN_PK = 22223  # Your previous run PK
```

Optionally modify parameters (e.g., increase NSW, adjust EDIFFG):

```python
slab_parameters = {
    "NSW": 200,        # Increase from 100
    "EDIFFG": -0.05,   # Relax from -0.1
    # ... other parameters
}
```

Run the restart:

```bash
python restart_example.py
```

### 4. Monitor Progress

```bash
verdi process show <NEW_PK>
verdi process report <NEW_PK>
```

## How It Works

### Under the Hood

When you provide `restart_from_node=<PK>`:

1. **Data Extraction**: 
   - Loads previous workgraph node
   - Extracts `slab_structures` (starting structures)
   - Extracts `slab_remote` (RemoteData with WAVECAR, CONTCAR, etc.)

2. **VASP Task Creation**:
   - Creates individual VaspWorkChain tasks (not using scatter-gather)
   - Each task gets the correct `restart_folder` input
   - VASP reads WAVECAR (wavefunctions) and CONTCAR (positions) from restart_folder

3. **Calculation Continuation**:
   - VASP continues from the last ionic step
   - Uses provided parameters (can be different from original)
   - Bulk/reference calculations are re-run (not restarted)

4. **Output Collection**:
   - All outputs collected via collector task
   - New RemoteData created for each slab
   - Can be used for another restart if needed!

### VASP Restart Behavior

When `restart_folder` is provided to VaspWorkChain:
- **WAVECAR**: Electronic wavefunctions from previous SCF
- **CONTCAR**: Last ionic positions as starting structure  
- **CHGCAR**: Charge density (if available)
- Continues optimization from last step instead of starting from scratch

This saves time and can help with convergence in difficult cases.

## Examples

### Example 1: Increase NSW (More Ionic Steps)

Original calculation didn't converge in 100 steps:

```python
# In restart_example.py
slab_parameters = {
    # ... same parameters as original ...
    "NSW": 300,  # Increased from 100
}

wg = build_core_workgraph(
    # ... same parameters ...
    restart_from_node=22223,
)
```

### Example 2: Relax Convergence Criteria

Original EDIFFG=-0.01 was too strict:

```python
slab_parameters = {
    # ... same parameters ...
    "EDIFFG": -0.05,  # Relaxed from -0.01
}

wg = build_core_workgraph(
    # ... same parameters ...
    restart_from_node=22223,
)
```

### Example 3: Change Algorithm

RMM-DIIS (IBRION=1) converging faster than CG (IBRION=2):

```python
slab_parameters = {
    # ... same parameters ...
    "IBRION": 1,   # Changed from 2
    "POTIM": 0.5,  # Adjust step size
}

wg = build_core_workgraph(
    # ... same parameters ...
    restart_from_node=22223,
)
```

### Example 4: Chain Multiple Restarts

```python
# First restart (increase NSW)
wg1 = build_core_workgraph(..., restart_from_node=22223)
wg1.submit()
# PK: 22844

# Wait for completion...

# Second restart (tighten convergence)
wg2 = build_core_workgraph(
    ...,
    slab_parameters={"EDIFFG": -0.03, "NSW": 400},
    restart_from_node=22844,  # Restart from first restart!
)
wg2.submit()
```

## API Reference

### `build_core_workgraph(..., restart_from_node=<PK>)`

**Parameters:**

- `restart_from_node` (int, optional): PK of previous PS-TEROS workgraph to restart from
  - Must have `slab_structures` and `slab_remote` outputs
  - Automatically extracts structures and RemoteData
  - Overrides `input_slabs` parameter

**Behavior with restart_from_node:**

1. Loads previous workgraph node
2. Extracts slab structures → uses as `input_slabs`
3. Extracts RemoteData → passes as `restart_folder` to each VASP task
4. Creates individual VASP tasks (not using relax_slabs_scatter)
5. Enables thermodynamics and cleavage calculations
6. Returns new RemoteData for potential future restarts

**Outputs (identical to normal mode):**

- `bulk_energy`, `bulk_structure`
- `metal_energy`, `metal_structure`
- `nonmetal_energy`, `nonmetal_structure`
- `oxygen_energy`, `oxygen_structure`
- `formation_enthalpy`
- `slab_structures` - Input slab structures
- `relaxed_slabs` - Relaxed structures from restart
- `slab_energies` - Energies from restart
- `slab_remote` - **NEW RemoteData** from restart (can restart again!)
- `surface_energies` - Surface energy calculations
- `cleavage_energies` - Cleavage energy calculations

### Helper Function: `extract_restart_folders_from_node(node_pk)`

Manually extract RemoteData from a previous run:

```python
from teros.core.slabs import extract_restart_folders_from_node

restart_folders = extract_restart_folders_from_node(22223)
# Returns: {'term_0': RemoteData(22310), 'term_1': RemoteData(22311)}
```

**Parameters:**
- `node_pk` (int): PK of PS-TEROS workgraph

**Returns:**
- `dict`: Mapping of slab labels to RemoteData nodes

**Raises:**
- `ValueError`: If node doesn't have `slab_remote` outputs

## Troubleshooting

### Error: "Node does not have 'slab_remote' outputs"

**Cause**: Previous run was with an older version of PS-TEROS that didn't expose RemoteData.

**Solution**: 
1. Re-run with updated PS-TEROS code, OR
2. Manually extract RemoteData from individual VASP calculations

### Error: "Could not load node <PK>"

**Cause**: Node doesn't exist or wrong PK.

**Solution**: 
1. Check PK with `verdi process list -a`
2. Verify it's a PS-TEROS workgraph with `verdi process show <PK>`

### Remote Files Cleaned

**Cause**: `clean_workdir=True` was used in previous run.

**Solution**:
- Remote files may have been deleted
- Cannot restart from cleaned directories
- Run new calculation from scratch

### Restart Still Not Converging

**Solutions**:
1. **Increase NSW**: Give more ionic steps
2. **Relax EDIFFG**: Use less strict convergence (e.g., -0.05 instead of -0.01)
3. **Change algorithm**: Try IBRION=1 (RMM-DIIS) instead of IBRION=2 (CG)
4. **Adjust POTIM**: Reduce step size if oscillating
5. **Check structure**: May have fundamental instability

## Best Practices

### 1. Keep Remote Directories

```python
wg = build_core_workgraph(
    # ...
    clean_workdir=False,  # Keep files for restart!
)
```

### 2. Progressive Refinement

Start with loose convergence, progressively tighten:

```python
# Run 1: Quick convergence
params1 = {"NSW": 100, "EDIFFG": -0.1}

# Run 2: Tighter (restart from Run 1)
params2 = {"NSW": 200, "EDIFFG": -0.05}

# Run 3: Final (restart from Run 2)
params3 = {"NSW": 300, "EDIFFG": -0.01}
```

### 3. Monitor Convergence

Before restarting, check why it failed:

```bash
# Check last ionic steps
verdi process report <PK>

# Check VASP outputs
verdi calcjob outputcat <VASP_PK> --path OUTCAR | tail -100
```

### 4. Document Parameters

Keep track of parameters used in each restart:

```python
wg.name = f"Ag2O_Restart_from_{PREVIOUS_PK}_NSW{NSW}_EDIFFG{EDIFFG}"
```

## Limitations

### Current Limitations

1. **Bulk/References Not Restarted**: Only slab calculations use restart functionality
   - Bulk, metal, oxygen calculations are re-run from scratch
   - This is usually fine as they converge quickly

2. **Same Slabs Required**: Must restart with same slab structures
   - Automatically handled by extracting from previous run
   - Cannot add/remove slabs in restart

3. **VASP Compatibility**: VASP must be able to read files from previous run
   - Different VASP versions may have compatibility issues
   - Check VASP documentation for WAVECAR/CONTCAR compatibility

### Not Limitations

These work perfectly:

✅ Thermodynamics calculations  
✅ Cleavage energy calculations  
✅ Full output provenance  
✅ Chain multiple restarts  
✅ Access all outputs

## Files in This Directory

- **`README.md`** (this file): Complete documentation
- **`restart_example.py`**: Minimal restart example
- **`slabs_relax_ag2o_restart.py`**: Full Ag2O example with restart capability
- **`TUTORIAL.md`**: Step-by-step tutorial with screenshots
- **`ADVANCED.md`**: Advanced usage patterns and tips

## Getting Help

### Check Status

```bash
verdi process show <PK>
verdi process report <PK>
verdi process list -a -p1
```

### Access Outputs

```python
from aiida import orm, load_profile
load_profile('psteros')

wg = orm.load_node(<PK>)

# Check available outputs
print(list(wg.outputs))

# Access specific outputs
relaxed_term_0 = wg.outputs.relaxed_slabs.term_0
energy_term_0 = wg.outputs.slab_energies.term_0
remote_term_0 = wg.outputs.slab_remote.term_0
```

### Debug Issues

```bash
# Check daemon
verdi daemon status

# Check VASP calculation
verdi process show <VASP_PK>
verdi calcjob outputcat <VASP_PK> --path OUTCAR

# Check remote folder
verdi node show <REMOTE_PK>
```

## Version History

- **v1.0** (2025-01-09): Initial implementation
  - Expose RemoteData nodes
  - Basic restart functionality
  
- **v2.0** (2025-01-09): Complete feature parity
  - Thermodynamics support
  - Cleavage calculations support
  - All outputs exposed
  - Collector task for proper output handling

## References

- **Main PS-TEROS Documentation**: `../../docs/`
- **AiiDA Documentation**: https://aiida.readthedocs.io
- **aiida-vasp Documentation**: https://aiida-vasp.readthedocs.io
- **VASP Restart Guide**: VASP manual section on ISTART

---

**Need more help?** Check `TUTORIAL.md` for a step-by-step guide with examples.
