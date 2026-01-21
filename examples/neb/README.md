# NEB (Nudged Elastic Band) Calculations

This module provides tools for calculating minimum energy paths (MEP) and activation barriers using the Nudged Elastic Band method with VASP.

## Overview

The NEB method finds the minimum energy path between two stable configurations (initial and final states). The activation barrier is the energy difference between the highest-energy point along the path (transition state) and the initial state.

### Workflow Stages

```
┌─────────────────┐     ┌─────────────────┐
│ Relax Initial   │     │ Relax Final     │  (parallel, optional)
│ (vasp.v2.vasp)  │     │ (vasp.v2.vasp)  │
└────────┬────────┘     └────────┬────────┘
         │                       │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │ Interpolate Structures │
         │ (IDPP/linear)          │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │ Stage 1: NEB          │
         │ (vasp.neb, LCLIMB=F)  │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │ Stage 2: CI-NEB       │  (optional, climb=True)
         │ (vasp.neb, LCLIMB=T)  │
         └───────────┬───────────┘
                     ▼
         ┌───────────────────────┐
         │ Calculate Barrier     │
         └───────────────────────┘
```

## Quick Start

```python
from aiida import load_profile, orm
from teros.core.neb import build_neb_workgraph, print_neb_summary

load_profile()

# Load your initial and final structures
initial = orm.load_node(...)  # or create from file
final = orm.load_node(...)

# Configure calculation
builder_inputs = {
    'parameters': {
        'incar': {
            'encut': 520,
            'ediff': 1e-6,
            'ismear': 0,
            'sigma': 0.05,
        }
    },
    'kpoints_spacing': 0.03,
    'potential_family': 'PBE',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    'options': {
        'resources': {'num_machines': 3, 'num_cores_per_machine': 40},
        'queue_name': 'par120',
    },
}

# Build and submit
wg = build_neb_workgraph(
    initial_structure=initial,
    final_structure=final,
    n_images=5,
    code_label='VASP-6.4.3@bohr',
    builder_inputs=builder_inputs,
    relax_endpoints=True,
    climb=True,  # Two-stage NEB → CI-NEB
)

wg.submit(wait=False)
print(f'Submitted: {wg.pk}')

# After completion:
print_neb_summary(wg.pk)
```

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `initial_structure` | StructureData/int | Initial state (reactant) |
| `final_structure` | StructureData/int | Final state (product) |
| `n_images` | int | Number of intermediate images (3-7 typical) |
| `code_label` | str | VASP code label (VTST-compiled recommended) |
| `builder_inputs` | dict | VASP configuration |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `relax_endpoints` | True | Relax initial/final before NEB |
| `interpolation_method` | 'idpp' | Image generation method ('idpp' or 'linear') |
| `climb` | True | Enable two-stage NEB→CI-NEB workflow |
| `spring_constant` | -5.0 | NEB spring constant (negative = variable) |
| `neb_optimizer` | 1 | VTST optimizer (1=LBFGS, 3=FIRE) |
| `force_convergence` | -0.05 | EDIFFG value (negative = force-based) |
| `max_steps` | 500 | Max ionic steps per NEB stage |
| `restart_folder` | None | RemoteData for continuing interrupted NEB |

## Interpolation Methods

### IDPP (Recommended)
Image Dependent Pair Potential generates initial paths that avoid atomic clashes by optimizing a surrogate potential. Better for reactive systems.

### Linear
Simple linear interpolation of atomic positions. Faster but may produce poor initial paths for complex reactions.

## Climbing Image NEB (CI-NEB)

When `climb=True` (default), the workflow runs two stages:

1. **Stage 1 (NEB)**: Regular NEB to find approximate MEP
2. **Stage 2 (CI-NEB)**: Climbing image refinement for accurate transition state

The climbing image is allowed to move uphill along the reaction coordinate, providing a more accurate saddle point location.

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| `barrier` | Dict | Forward/reverse barriers, reaction energy |
| `summary` | Dict | Comprehensive calculation summary |
| `energies` | Dict | Per-image energies |
| `misc` | Dict | Raw VASP NEB output |
| `remote_folder` | RemoteData | Final NEB files (for restart/analysis) |
| `stage1_remote_folder` | RemoteData | Stage 1 files (for analysis) |
| `relaxed_initial` | StructureData | Relaxed initial (if `relax_endpoints=True`) |
| `relaxed_final` | StructureData | Relaxed final (if `relax_endpoints=True`) |

## Accessing Results

```python
from teros.core.neb import get_neb_results, print_neb_summary

# Get all results as a dictionary
results = get_neb_results(wg.pk)

# Access barrier information
barrier = results['barrier']
print(f"Forward barrier: {barrier['forward_barrier']:.3f} eV")
print(f"Reverse barrier: {barrier['reverse_barrier']:.3f} eV")
print(f"Reaction energy: {barrier['reaction_energy']:.3f} eV")
print(f"Saddle point at image: {barrier['saddle_point_index']}")

# Energy profile along MEP
for i, energy in enumerate(barrier['energies_list']):
    rel_energy = energy - barrier['energies_list'][0]
    print(f"Image {i}: {rel_energy:+.4f} eV")

# Or use the formatted summary
print_neb_summary(wg.pk)
```

## VASP Requirements

### VTST Extensions (Recommended)
For optimal NEB performance, VASP should be compiled with VTST extensions which provide:
- Improved optimizers (IOPT=1 LBFGS, IOPT=3 FIRE)
- Better convergence for transition state searches

Without VTST, the workflow will still run but may use less efficient optimizers.

### Key INCAR Parameters

The workflow automatically sets NEB-specific parameters:

```
IMAGES = n_images      # Number of intermediate images
IBRION = 3             # BFGS optimizer (required for NEB)
POTIM = 0.0            # Use IOPT optimizer
SPRING = -5.0          # Variable spring constant
LCLIMB = .TRUE./.FALSE.  # Climbing image (Stage 2 only)
IOPT = 1               # LBFGS optimizer (VTST)
NSW = 500              # Max ionic steps
EDIFFG = -0.05         # Force convergence
```

## Troubleshooting

### Structure Validation Error
```
ValueError: Initial and final structures have different atom counts
```
**Solution**: Ensure both structures have identical composition and atom count. NEB requires the same atoms in both endpoints.

### SCF Convergence Issues
- Increase ENCUT
- Reduce EDIFF for tighter electronic convergence
- Check k-point sampling

### Poor Initial Path
- Switch from 'linear' to 'idpp' interpolation
- Increase number of images
- Manually adjust intermediate images

### Slow Convergence
- Try IOPT=3 (FIRE) instead of LBFGS
- Reduce force convergence threshold initially
- Ensure endpoints are well-relaxed

## Example Applications

### Vacancy Migration
Calculate the activation barrier for vacancy diffusion in a crystal.

### Surface Diffusion
Find the energy barrier for adatom migration on surfaces.

### Chemical Reactions
Determine transition state energies for bond breaking/forming processes.

### Phase Transitions
Map the minimum energy path for structural transformations.

## References

- Henkelman, G., Uberuaga, B. P., & Jónsson, H. (2000). A climbing image nudged elastic band method for finding saddle points and minimum energy paths. J. Chem. Phys., 113(22), 9901-9904.
- Henkelman, G., & Jónsson, H. (2000). Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points. J. Chem. Phys., 113(22), 9978-9985.
- VASP VTST Tools: http://theory.cm.utexas.edu/vtsttools/
