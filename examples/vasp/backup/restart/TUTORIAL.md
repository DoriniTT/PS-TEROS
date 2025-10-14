# PS-TEROS Restart Feature Tutorial

## Step-by-Step Guide with Examples

This tutorial walks you through using the restart feature with a complete example.

## Prerequisites

1. AiiDA profile 'psteros' configured
2. VASP code installed and configured in AiiDA
3. PS-TEROS with restart feature installed
4. Daemon running: `verdi daemon start`

## Tutorial: Restarting a Failed Ag2O Slab Relaxation

### Step 1: Run Initial Calculation

Let's start with a calculation that won't converge (intentionally low NSW):

**File: `tutorial_step1.py`**

```python
#!/home/thiagotd/envs/psteros/bin/python
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

load_profile('psteros')

# Parameters that won't converge (NSW too low)
slab_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 20,          # TOO LOW - won't converge!
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

wg = build_core_workgraph(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    nonmetal_name="Ag.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={"Ag": "Ag", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"Ag": "Ag"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters={"ENCUT": 500, "NSW": 100, "EDIFFG": -0.1},
    bulk_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    metal_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    metal_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    nonmetal_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    nonmetal_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    oxygen_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    oxygen_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    clean_workdir=False,  # IMPORTANT: Keep files for restart!
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    slab_potential_mapping={"Ag": "Ag", "O": "O"},
    slab_kpoints_spacing=0.3,
    compute_thermodynamics=True,
    thermodynamics_sampling=100,
    name="Ag2O_Tutorial_Initial_LowNSW",
)

wg.submit(wait=False)
print(f"Submitted workflow: PK {wg.pk}")
print(f"Monitor with: verdi process show {wg.pk}")
```

**Run it:**

```bash
python tutorial_step1.py
```

**Output:**
```
Submitted workflow: PK 30001
Monitor with: verdi process show 30001
```

### Step 2: Wait and Check Status

Wait for the calculation to complete (will not converge):

```bash
# Check status
verdi process show 30001

# Check if slabs didn't converge
verdi process report 30001 | grep -i "converge"
```

**Expected output:**
```
Property     Value
-----------  ------------------------------------
type         WorkGraph<Ag2O_Tutorial_Initial_LowNSW>
state        Finished [0]
...

Outputs             PK     Type
------------------  -----  -------------
slab_remote
    term_0          30050  RemoteData  ← THESE ARE KEY!
    term_1          30051  RemoteData
relaxed_slabs
    term_0          30048  StructureData
    term_1          30049  StructureData
...
```

✓ Key observation: We have `slab_remote` outputs, so we can restart!

### Step 3: Examine Why It Didn't Converge

Check the VASP output:

```bash
# Find the VASP calculation PK
verdi process show 30001 | grep VaspWorkChain_slab

# Check the OUTCAR
verdi calcjob outputcat 30055 --path OUTCAR | tail -50
```

**You'll see:**
```
...
  20 F= -.12345E+03 E0= -.12346E+03  d E =-.123456E+00
 reached required accuracy - stopping structural energy minimisation
```

The calculation stopped at NSW=20 (didn't reach convergence).

### Step 4: Prepare Restart with More Steps

**File: `tutorial_step2_restart.py`**

```python
#!/home/thiagotd/envs/psteros/bin/python
from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph

load_profile('psteros')

# ===== RESTART CONFIGURATION =====
PREVIOUS_RUN_PK = 30001  # From Step 1

# Check previous run
prev_node = orm.load_node(PREVIOUS_RUN_PK)
print(f"Restarting from: {prev_node.label}")
print(f"  State: {prev_node.process_state}")
print(f"  Slabs: {list(prev_node.outputs.slab_structures.keys())}")
print(f"  RemoteData: {list(prev_node.outputs.slab_remote.keys())}")
print()

# NEW parameters - increase NSW
slab_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 100,         # INCREASED from 20!
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}

wg = build_core_workgraph(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    nonmetal_name="Ag.cif",
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={"Ag": "Ag", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"Ag": "Ag"},
    oxygen_potential_mapping={"O": "O"},
    kpoints_spacing=0.3,
    bulk_parameters={"ENCUT": 500, "NSW": 100, "EDIFFG": -0.1},
    bulk_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    metal_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    metal_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    nonmetal_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    nonmetal_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    oxygen_parameters={"ENCUT": 520, "NSW": 100, "EDIFFG": -0.1},
    oxygen_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    clean_workdir=False,  # Keep for potential future restart
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options={"resources": {"num_machines": 1, "num_cores_per_machine": 40}, "queue_name": "par40"},
    slab_potential_mapping={"Ag": "Ag", "O": "O"},
    slab_kpoints_spacing=0.3,
    compute_thermodynamics=True,
    thermodynamics_sampling=100,
    # ===== KEY PARAMETER =====
    restart_from_node=PREVIOUS_RUN_PK,
    name=f"Ag2O_Tutorial_Restart_from_{PREVIOUS_RUN_PK}",
)

wg.submit(wait=False)
print(f"\n✓ Restart submitted: PK {wg.pk}")
print(f"Monitor with: verdi process show {wg.pk}")
```

**Run it:**

```bash
python tutorial_step2_restart.py
```

**Output:**
```
Restarting from: Ag2O_Tutorial_Initial_LowNSW
  State: ProcessState.FINISHED
  Slabs: ['term_0', 'term_1']
  RemoteData: ['term_0', 'term_1']

======================================================================
RESTART MODE: Loading data from node 30001
======================================================================
  ✓ Extracted restart folders: ['term_0', 'term_1']
  ✓ Extracted slab structures: ['term_0', 'term_1']
  → Using slabs from previous run
======================================================================

  → Building workgraph with restart_folder for each slab
    → term_0: using restart_folder PK 30050
    → term_1: using restart_folder PK 30051
  ✓ Created 2 VASP tasks with restart capability
  ✓ All slab outputs connected via collector task
  ✓ Thermodynamics calculation enabled (surface energies)
  ✓ Cleavage energies calculation enabled

✓ Restart submitted: PK 30100
Monitor with: verdi process show 30100
```

### Step 5: Monitor Restart

```bash
# Watch progress
watch -n 10 'verdi process show 30100 | head -20'

# Check specific VASP task
verdi process show 30100 | grep VaspWorkChain_slab
# VaspWorkChain_slab_term_0  30105  VaspWorkChain

verdi process show 30105 | grep restart_folder
# restart_folder     30050  RemoteData  ← Using previous run!
```

### Step 6: Verify Results

After completion:

```bash
verdi process show 30100
```

**Output:**
```
Property     Value
-----------  ---------------------------------------------
type         WorkGraph<Ag2O_Tutorial_Restart_from_30001>
state        Finished [0]  ← SUCCESS!

Outputs             PK     Type
------------------  -----  -------------
cleavage_energies
    pair_0_1        30180  Dict
relaxed_slabs
    term_0          30170  StructureData  ← CONVERGED!
    term_1          30171  StructureData
slab_energies
    term_0          30172  Float
    term_1          30173  Float
slab_remote
    term_0          30174  RemoteData  ← NEW! Can restart again
    term_1          30175  RemoteData
surface_energies
    term_0          30181  Dict
    term_1          30182  Dict
...
```

### Step 7: Analyze Results

```python
from aiida import orm, load_profile
load_profile('psteros')

# Load both runs
initial = orm.load_node(30001)
restart = orm.load_node(30100)

print("=" * 70)
print("COMPARISON: Initial vs Restart")
print("=" * 70)

# Compare structures
print("\nInitial run slab_remote:")
for label in ['term_0', 'term_1']:
    remote = initial.outputs.slab_remote[label]
    print(f"  {label}: PK {remote.pk}")

print("\nRestart run slab_remote:")
for label in ['term_0', 'term_1']:
    remote = restart.outputs.slab_remote[label]
    print(f"  {label}: PK {remote.pk} ← NEW!")

# Compare energies
print("\nSlab energies:")
for label in ['term_0', 'term_1']:
    e_init = initial.outputs.slab_energies[label].value
    e_restart = restart.outputs.slab_energies[label].value
    print(f"  {label}: {e_init:.4f} eV → {e_restart:.4f} eV (Δ = {e_restart-e_init:.4f} eV)")

print("\n✓ Restart successful with improved convergence!")
```

## Advanced Tutorial: Chain Multiple Restarts

### Scenario: Progressive Convergence

Start loose → medium → tight convergence.

**Step 1: Loose (Quick)**
```python
# NSW=50, EDIFFG=-0.1
wg1 = build_core_workgraph(..., slab_parameters={"NSW": 50, "EDIFFG": -0.1})
# PK: 40001
```

**Step 2: Medium (Restart from Step 1)**
```python
# NSW=100, EDIFFG=-0.05
wg2 = build_core_workgraph(
    ..., 
    slab_parameters={"NSW": 100, "EDIFFG": -0.05},
    restart_from_node=40001,
)
# PK: 40002
```

**Step 3: Tight (Restart from Step 2)**
```python
# NSW=200, EDIFFG=-0.01
wg3 = build_core_workgraph(
    ..., 
    slab_parameters={"NSW": 200, "EDIFFG": -0.01},
    restart_from_node=40002,
)
# PK: 40003
```

Each restart uses WAVECAR from previous, converging faster!

## Common Patterns

### Pattern 1: Increase NSW Incrementally

```python
# Run 1
wg1 = build_core_workgraph(..., slab_parameters={"NSW": 50})

# Check convergence
# Not converged? Restart with more NSW

# Run 2
wg2 = build_core_workgraph(
    ..., 
    slab_parameters={"NSW": 100},
    restart_from_node=wg1.pk,
)
```

### Pattern 2: Switch Algorithms

```python
# Run 1: CG algorithm
wg1 = build_core_workgraph(..., slab_parameters={"IBRION": 2})

# Slow convergence? Try RMM-DIIS

# Run 2: RMM-DIIS algorithm
wg2 = build_core_workgraph(
    ..., 
    slab_parameters={"IBRION": 1, "POTIM": 0.5},
    restart_from_node=wg1.pk,
)
```

### Pattern 3: Emergency Recovery

Calculation crashed? Just restart:

```bash
# Check last good state
verdi process show <CRASHED_PK>

# Restart from last checkpoint
python restart_example.py  # Set PREVIOUS_RUN_PK to crashed run
```

## Tips and Tricks

### Tip 1: Always Set clean_workdir=False

```python
wg = build_core_workgraph(
    # ...
    clean_workdir=False,  # Keep for restart!
)
```

### Tip 2: Name Your Restarts Descriptively

```python
wg = build_core_workgraph(
    # ...
    restart_from_node=30001,
    name=f"Ag2O_Restart_{30001}_NSW200_tight",
)
```

### Tip 3: Check Convergence Before Restarting

```bash
# Check if forces are decreasing
verdi calcjob outputcat <VASP_PK> --path OUTCAR | grep "FORCES" -A 10
```

### Tip 4: Monitor Remote Storage

```bash
# Check remote folder size
du -sh ~/.aiida/<computer>/<uuid>
```

## Troubleshooting Tutorial Problems

### Problem: "No slab_remote outputs"

**Solution**: Previous run was with old PS-TEROS. Update and re-run initial calculation.

### Problem: Restart still doesn't converge

**Diagnosis**:
```bash
verdi calcjob outputcat <VASP_PK> --path OUTCAR | tail -100
```

**Solutions**:
- Increase NSW more
- Relax EDIFFG
- Try different algorithm (IBRION)
- Check if structure is fundamentally unstable

### Problem: Remote files deleted

**Cause**: `clean_workdir=True` was used.

**Solution**: Cannot restart. Run from scratch with `clean_workdir=False`.

## Next Steps

1. **Read README.md**: Full API documentation
2. **Read ADVANCED.md**: Advanced patterns and optimization
3. **Try your system**: Apply to your own materials
4. **Explore outputs**: Use AiiDA tools to analyze results

## Questions?

Check `README.md` for detailed API reference and troubleshooting guide.

---

**Congratulations!** You now know how to use the PS-TEROS restart feature. 🎉
