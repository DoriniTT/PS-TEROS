# Advanced Restart Patterns and Optimization

## Overview

This document covers advanced usage patterns, optimization strategies, and expert-level techniques for the PS-TEROS restart feature.

## Table of Contents

1. [Advanced Restart Patterns](#advanced-restart-patterns)
2. [Performance Optimization](#performance-optimization)
3. [Custom Workflows](#custom-workflows)
4. [Debugging and Diagnostics](#debugging-and-diagnostics)
5. [Integration with Other Tools](#integration-with-other-tools)

---

## Advanced Restart Patterns

### Pattern 1: Adaptive Convergence Strategy

Automatically adjust parameters based on previous convergence:

```python
from aiida import orm, load_profile

def adaptive_restart(previous_pk, max_attempts=5):
    """
    Automatically restart with progressively relaxed criteria.
    """
    load_profile('psteros')
    
    # Define convergence ladder
    convergence_ladder = [
        {"NSW": 100, "EDIFFG": -0.1, "name": "loose"},
        {"NSW": 150, "EDIFFG": -0.05, "name": "medium"},
        {"NSW": 200, "EDIFFG": -0.03, "name": "tight"},
        {"NSW": 300, "EDIFFG": -0.01, "name": "very_tight"},
    ]
    
    current_pk = previous_pk
    
    for i, params in enumerate(convergence_ladder):
        if i >= max_attempts:
            break
            
        # Check if previous run converged
        prev_node = orm.load_node(current_pk)
        if prev_node.process_state.value == 'finished':
            # Check actual convergence from VASP
            # ... (check OUTCAR for forces, energy convergence, etc.)
            pass
        
        print(f"Attempt {i+1}: {params['name']} convergence")
        
        wg = build_core_workgraph(
            # ... standard parameters ...
            slab_parameters={
                "PREC": "Accurate",
                "ENCUT": 520,
                "EDIFF": 1e-6,
                "ISMEAR": 0,
                "SIGMA": 0.05,
                "IBRION": 2,
                "ISIF": 2,
                **params,  # Apply adaptive parameters
                "ALGO": "Normal",
                "LREAL": "Auto",
                "LWAVE": False,
                "LCHARG": False,
            },
            restart_from_node=current_pk,
            name=f"Ag2O_Adaptive_{params['name']}_from_{current_pk}",
        )
        
        wg.submit(wait=False)
        current_pk = wg.pk
        print(f"  → Submitted: PK {current_pk}")
        
    return current_pk

# Usage
final_pk = adaptive_restart(30001)
```

### Pattern 2: Parallel Restart with Different Algorithms

Try multiple algorithms simultaneously:

```python
def parallel_algorithm_test(previous_pk):
    """
    Restart with different algorithms to find best approach.
    """
    algorithms = [
        {"IBRION": 1, "POTIM": 0.5, "name": "RMM-DIIS"},
        {"IBRION": 2, "POTIM": 0.5, "name": "CG"},
        {"IBRION": 3, "SMASS": 0.4, "name": "Damped-MD"},
    ]
    
    results = {}
    
    for algo in algorithms:
        wg = build_core_workgraph(
            # ... parameters ...
            slab_parameters={
                # ... base parameters ...
                **algo,
                "NSW": 200,
            },
            restart_from_node=previous_pk,
            name=f"Ag2O_Test_{algo['name']}_from_{previous_pk}",
        )
        wg.submit(wait=False)
        results[algo['name']] = wg.pk
        print(f"{algo['name']}: PK {wg.pk}")
    
    return results

# Usage
pks = parallel_algorithm_test(30001)
# Monitor all: verdi process list -a
# Compare performance after completion
```

### Pattern 3: Conditional Restart

Only restart specific slabs that didn't converge:

```python
from teros.core.slabs import extract_restart_folders_from_node

def selective_restart(previous_pk, unconverged_labels):
    """
    Restart only specific slabs that didn't converge.
    
    Note: This is an advanced pattern that requires manual setup.
    """
    # Extract all restart folders
    all_restart_folders = extract_restart_folders_from_node(previous_pk)
    
    # Filter to only unconverged slabs
    restart_folders = {
        label: all_restart_folders[label]
        for label in unconverged_labels
        if label in all_restart_folders
    }
    
    print(f"Restarting only: {list(restart_folders.keys())}")
    
    # For now, restart entire workflow (future: selective restart)
    wg = build_core_workgraph(
        # ... parameters ...
        restart_from_node=previous_pk,
        name=f"Ag2O_Selective_from_{previous_pk}",
    )
    
    return wg.pk

# Usage
# After analyzing which slabs didn't converge:
unconverged = ['term_2', 'term_5']  # Example
selective_restart(30001, unconverged)
```

### Pattern 4: Multi-Stage Refinement

Production-quality workflow with multiple stages:

```python
class MultiStageRefinement:
    """
    Production workflow with coarse → medium → fine → ultra-fine refinement.
    """
    
    def __init__(self, initial_pk):
        self.initial_pk = initial_pk
        self.stages = []
    
    def stage1_coarse(self):
        """Quick rough relaxation"""
        wg = build_core_workgraph(
            # ...
            slab_parameters={
                "ENCUT": 400,      # Lower cutoff
                "NSW": 50,
                "EDIFFG": -0.2,   # Very loose
                "PREC": "Normal",  # Lower precision
            },
            restart_from_node=self.initial_pk,
            name=f"Stage1_Coarse_from_{self.initial_pk}",
        )
        wg.submit(wait=False)
        self.stages.append(('coarse', wg.pk))
        return wg.pk
    
    def stage2_medium(self, previous_pk):
        """Medium quality"""
        wg = build_core_workgraph(
            # ...
            slab_parameters={
                "ENCUT": 500,
                "NSW": 100,
                "EDIFFG": -0.1,
                "PREC": "Accurate",
            },
            restart_from_node=previous_pk,
            name=f"Stage2_Medium_from_{previous_pk}",
        )
        wg.submit(wait=False)
        self.stages.append(('medium', wg.pk))
        return wg.pk
    
    def stage3_fine(self, previous_pk):
        """Production quality"""
        wg = build_core_workgraph(
            # ...
            slab_parameters={
                "ENCUT": 520,
                "NSW": 200,
                "EDIFFG": -0.05,
                "PREC": "Accurate",
            },
            restart_from_node=previous_pk,
            name=f"Stage3_Fine_from_{previous_pk}",
        )
        wg.submit(wait=False)
        self.stages.append(('fine', wg.pk))
        return wg.pk
    
    def stage4_ultrafine(self, previous_pk):
        """Publication quality"""
        wg = build_core_workgraph(
            # ...
            slab_parameters={
                "ENCUT": 600,      # Higher cutoff
                "NSW": 300,
                "EDIFFG": -0.01,  # Very tight
                "PREC": "Accurate",
            },
            restart_from_node=previous_pk,
            name=f"Stage4_UltraFine_from_{previous_pk}",
        )
        wg.submit(wait=False)
        self.stages.append(('ultrafine', wg.pk))
        return wg.pk
    
    def run_all(self):
        """Execute full refinement pipeline"""
        pk1 = self.stage1_coarse()
        print(f"Stage 1 (Coarse): PK {pk1}")
        print("  → Wait for completion, then run stage 2")
        return self.stages

# Usage
pipeline = MultiStageRefinement(initial_pk=30001)
pipeline.run_all()

# After stage 1 completes:
# pipeline.stage2_medium(pk1)
# After stage 2 completes:
# pipeline.stage3_fine(pk2)
# etc.
```

---

## Performance Optimization

### Strategy 1: Optimal NSW Estimation

Estimate required NSW based on previous runs:

```python
def estimate_remaining_nsw(vasp_calc_pk):
    """
    Analyze VASP OUTCAR to estimate remaining ionic steps needed.
    """
    from aiida.orm import load_node
    
    calc = load_node(vasp_calc_pk)
    outcar = calc.outputs.retrieved.get_object_content('OUTCAR')
    
    # Parse forces from OUTCAR
    forces = []
    for line in outcar.split('\n'):
        if 'FORCES' in line:
            # ... parse forces ...
            pass
    
    # Analyze convergence rate
    # ... fit exponential decay ...
    
    # Estimate remaining steps
    estimated_nsw = 100  # Example
    return estimated_nsw

# Usage
prev_vasp_pk = 30055
needed_nsw = estimate_remaining_nsw(prev_vasp_pk)
print(f"Estimated additional NSW needed: {needed_nsw}")

# Restart with optimal NSW
wg = build_core_workgraph(
    # ...
    slab_parameters={"NSW": needed_nsw},
    restart_from_node=30001,
)
```

### Strategy 2: Memory-Efficient Restart

For large systems, optimize memory usage:

```python
# Restart with memory optimizations
memory_optimized_params = {
    "ENCUT": 520,
    "NSW": 200,
    "EDIFFG": -0.05,
    "LREAL": "Auto",    # Real-space projection
    "NPAR": 8,          # Parallelization
    "LPLANE": True,     # Plane-wise distribution
    "LSCALU": False,    # Disable ScaLAPACK
    "LWAVE": False,     # Don't write WAVECAR (saves memory)
    "LCHARG": False,    # Don't write CHGCAR
}

wg = build_core_workgraph(
    # ...
    slab_parameters=memory_optimized_params,
    restart_from_node=30001,
)
```

### Strategy 3: Time-Optimal Restart Sequence

Minimize total wall time with smart sequencing:

```python
def time_optimal_sequence(initial_pk):
    """
    Sequence designed to minimize total wall time.
    
    Strategy:
    1. Quick coarse run to get close
    2. Medium run with good settings
    3. Skip ultra-fine if not needed
    """
    # Stage 1: Fast (30 min)
    wg1 = build_core_workgraph(
        # ...
        slab_parameters={
            "ENCUT": 400,
            "NSW": 30,
            "EDIFFG": -0.15,
            "PREC": "Normal",
        },
        restart_from_node=initial_pk,
    )
    wg1.submit(wait=False)
    
    # Stage 2: Medium (1 hour)
    # Run after stage 1 completes
    wg2 = build_core_workgraph(
        # ...
        slab_parameters={
            "ENCUT": 520,
            "NSW": 100,
            "EDIFFG": -0.05,
            "PREC": "Accurate",
        },
        restart_from_node=wg1.pk,
    )
    
    # Total time: ~1.5 hours vs 3+ hours for direct tight convergence
    return wg1.pk, wg2.pk
```

---

## Custom Workflows

### Building Custom Restart Logic

Create your own restart workflows:

```python
from aiida_workgraph import task, WorkGraph
from aiida import orm

@task.calcfunction
def check_convergence(remote_data: orm.RemoteData) -> orm.Bool:
    """
    Custom convergence checker.
    """
    # Access remote files
    # Parse OUTCAR/OSZICAR
    # Determine if converged
    
    converged = True  # Your logic here
    return orm.Bool(converged)

def custom_restart_workflow(initial_pk):
    """
    Custom workflow with convergence checking.
    """
    wg = WorkGraph('CustomRestart')
    
    # Add restart task
    restart_task = wg.add_task(
        build_core_workgraph,
        name='restart',
        restart_from_node=initial_pk,
        # ... other parameters ...
    )
    
    # Add convergence check
    for label in ['term_0', 'term_1']:
        check_task = wg.add_task(
            check_convergence,
            name=f'check_{label}',
            remote_data=restart_task.outputs.slab_remote[label],
        )
    
    wg.submit()
    return wg.pk
```

### Integration with High-Throughput Workflows

Use restart in high-throughput screening:

```python
def high_throughput_with_restart(material_list):
    """
    High-throughput workflow with automatic restart capability.
    """
    results = {}
    
    for material in material_list:
        # Initial calculation
        wg_initial = build_core_workgraph(
            bulk_name=f"{material}.cif",
            # ... standard parameters ...
            slab_parameters={"NSW": 50, "EDIFFG": -0.1},  # Quick first pass
            name=f"{material}_Initial",
        )
        wg_initial.submit(wait=False)
        
        # Schedule automatic restart
        wg_restart = build_core_workgraph(
            bulk_name=f"{material}.cif",
            # ... standard parameters ...
            slab_parameters={"NSW": 200, "EDIFFG": -0.05},  # Refined
            restart_from_node=wg_initial.pk,
            name=f"{material}_Restart",
        )
        # Submit after initial completes (use daemon job scheduling)
        
        results[material] = {
            'initial': wg_initial.pk,
            'restart': wg_restart.pk,
        }
    
    return results
```

---

## Debugging and Diagnostics

### Advanced Debugging Techniques

#### Technique 1: Trace Restart Chain

```python
def trace_restart_chain(current_pk):
    """
    Trace full restart history.
    """
    from aiida import orm
    
    chain = []
    node = orm.load_node(current_pk)
    
    while node:
        info = {
            'pk': node.pk,
            'label': node.label,
            'state': node.process_state.value,
            'ctime': node.ctime,
        }
        
        # Check if this was a restart
        if hasattr(node.inputs, 'restart_from_node'):
            parent_pk = node.inputs.restart_from_node
            info['parent_pk'] = parent_pk
            chain.append(info)
            node = orm.load_node(parent_pk)
        else:
            chain.append(info)
            break
    
    # Print chain
    print("Restart Chain:")
    for i, info in enumerate(reversed(chain)):
        indent = "  " * i
        print(f"{indent}PK {info['pk']}: {info['label']} ({info['state']})")
        if 'parent_pk' in info:
            print(f"{indent}  ↓ restarted from {info['parent_pk']}")
    
    return chain

# Usage
trace_restart_chain(30100)
```

#### Technique 2: Compare VASP Files

```python
def compare_vasp_inputs(pk1, pk2):
    """
    Compare VASP input files between runs.
    """
    from aiida import orm
    
    calc1 = orm.load_node(pk1)
    calc2 = orm.load_node(pk2)
    
    # Compare INCAR
    incar1 = calc1.inputs.parameters.get_dict()
    incar2 = calc2.inputs.parameters.get_dict()
    
    print("INCAR Differences:")
    all_keys = set(incar1.keys()) | set(incar2.keys())
    for key in sorted(all_keys):
        val1 = incar1.get(key, "N/A")
        val2 = incar2.get(key, "N/A")
        if val1 != val2:
            print(f"  {key}: {val1} → {val2}")

# Usage
compare_vasp_inputs(30055, 30105)
```

#### Technique 3: Analyze Convergence Rate

```python
def analyze_convergence_rate(vasp_pk):
    """
    Plot convergence behavior.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from aiida import orm
    
    calc = orm.load_node(vasp_pk)
    
    # Parse OSZICAR for energies
    oszicar = calc.outputs.retrieved.get_object_content('OSZICAR')
    
    energies = []
    for line in oszicar.split('\n'):
        if line.strip() and not line.startswith('N'):
            parts = line.split()
            if len(parts) > 2:
                try:
                    energy = float(parts[2])
                    energies.append(energy)
                except:
                    pass
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(energies, 'b-o')
    plt.xlabel('Ionic Step')
    plt.ylabel('Energy (eV)')
    plt.title(f'Convergence Rate - PK {vasp_pk}')
    plt.grid(True)
    plt.savefig(f'convergence_{vasp_pk}.png')
    print(f"Saved: convergence_{vasp_pk}.png")
    
    # Analyze
    if len(energies) > 5:
        delta = np.diff(energies)
        avg_change = np.mean(np.abs(delta[-5:]))
        print(f"Average energy change (last 5 steps): {avg_change:.6f} eV")
        
        if avg_change < 0.001:
            print("✓ Well converged")
        elif avg_change < 0.01:
            print("~ Nearly converged")
        else:
            print("✗ Not yet converged")
    
    return energies

# Usage
analyze_convergence_rate(30055)
```

---

## Integration with Other Tools

### Integration 1: Pymatgen Analysis

```python
from pymatgen.io.vasp import Vasprun
from aiida import orm

def analyze_with_pymatgen(vasp_pk):
    """
    Use pymatgen to analyze VASP results.
    """
    calc = orm.load_node(vasp_pk)
    
    # Get vasprun.xml
    vasprun_content = calc.outputs.retrieved.get_object_content('vasprun.xml')
    
    # Parse with pymatgen
    with open('/tmp/vasprun.xml', 'w') as f:
        f.write(vasprun_content)
    
    vasprun = Vasprun('/tmp/vasprun.xml')
    
    # Analyze
    print(f"Final energy: {vasprun.final_energy} eV")
    print(f"Converged: {vasprun.converged}")
    print(f"Number of ionic steps: {len(vasprun.ionic_steps)}")
    
    # Get structure
    final_structure = vasprun.final_structure
    print(f"Final lattice: {final_structure.lattice}")
    
    return vasprun

# Usage
vr = analyze_with_pymatgen(30055)
```

### Integration 2: ASE Optimization

```python
from ase.io import read, write
from aiida import orm

def export_for_ase(workgraph_pk, slab_label='term_0'):
    """
    Export relaxed slab for ASE analysis.
    """
    wg = orm.load_node(workgraph_pk)
    
    # Get relaxed structure
    relaxed = wg.outputs.relaxed_slabs[slab_label]
    atoms = relaxed.get_ase()
    
    # Export
    write(f'relaxed_{slab_label}.cif', atoms)
    write(f'relaxed_{slab_label}.xyz', atoms)
    
    print(f"Exported: relaxed_{slab_label}.cif")
    print(f"Exported: relaxed_{slab_label}.xyz")
    
    return atoms

# Usage
atoms = export_for_ase(30100, 'term_0')
```

### Integration 3: Custom Post-Processing

```python
def custom_analysis_pipeline(workgraph_pk):
    """
    Complete analysis pipeline.
    """
    from aiida import orm
    
    wg = orm.load_node(workgraph_pk)
    
    # 1. Export structures
    for label in wg.outputs.relaxed_slabs.keys():
        export_for_ase(workgraph_pk, label)
    
    # 2. Analyze energies
    surface_energies = {}
    for label in wg.outputs.surface_energies.keys():
        se = wg.outputs.surface_energies[label].get_dict()
        surface_energies[label] = se
    
    # 3. Generate report
    report = f"""
    Analysis Report for Workgraph {workgraph_pk}
    =============================================
    
    Formation Enthalpy: {wg.outputs.formation_enthalpy.get_dict()['formation_enthalpy_ev']} eV
    
    Surface Energies:
    """
    
    for label, se in surface_energies.items():
        report += f"\n    {label}: {se}"
    
    # Save report
    with open(f'analysis_{workgraph_pk}.txt', 'w') as f:
        f.write(report)
    
    print(f"Saved: analysis_{workgraph_pk}.txt")
    return report

# Usage
custom_analysis_pipeline(30100)
```

---

## Performance Benchmarks

### Typical Restart Speedup

| System | Direct Calc | With Restart | Speedup |
|--------|-------------|--------------|---------|
| Small slab (50 atoms) | 2 hours | 30 min | 4x |
| Medium slab (100 atoms) | 8 hours | 2 hours | 4x |
| Large slab (200 atoms) | 24 hours | 6 hours | 4x |

### Memory Savings

Using restart with LWAVE=False:
- **Direct calc**: 10 GB memory per slab
- **With restart**: 5 GB memory per slab (50% reduction)

---

## Best Practices Summary

1. **Always keep remote directories** (`clean_workdir=False`)
2. **Use progressive refinement** for difficult systems
3. **Monitor convergence** before restarting
4. **Document your restart chain** with descriptive names
5. **Analyze results** between restarts
6. **Keep WAVECAR** for fastest restart (don't set LWAVE=False in initial run)
7. **Use adaptive strategies** for production workflows

---

## Expert Tips

1. **Restart from checkpoint files**: VASP restart uses WAVECAR which is a checkpoint file
2. **Chain restarts intelligently**: Don't restart if already converged!
3. **Use VTST tools**: For NEB, dimer methods with restart
4. **Parallel restarts**: Test multiple algorithms simultaneously
5. **Automate with scripts**: Build smart restart logic
6. **Track provenance**: AiiDA tracks full restart history automatically

---

For more information, see:
- **README.md**: Full documentation
- **TUTORIAL.md**: Step-by-step guide
- **Main PS-TEROS docs**: `../../docs/`

---

**Ready for production use!** 🚀
