# psteros

`psteros` is a reproducible AiiDA workflow library for ab-initio surface
thermodynamics.  Quantum ESPRESSO, through `aiida-quantumespresso`, is the
primary backend.  VASP remains available as a secondary, tested backend for
existing studies.

The public API is deliberately small:

```python
import psteros

recipe = psteros.SurfaceWorkflowConfig(
    backend="qe",
    calculation=psteros.QeCalculationConfig(
        code_label="QE-7.6-PW-GPU-A100@bohr",
        pseudo_family="SSSP/1.3/PBE/precision",
        parameters={
            "CONTROL": {"calculation": "relax"},
            "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
            "ELECTRONS": {"conv_thr": 1.0e-8},
        },
    ),
)
bulk = psteros.rutile_sno2_bulk()
graph = psteros.build_surface_workgraph({"bulk_sno2": bulk}, recipe)
```

`ExecutionPolicy` defaults to the maintained Bohr A100 contract: queue
`gpu_a100`, one GPU job at a time, one machine, and one MPI rank.  The library
does not start an external controller; a Tessera Odyssey Auto project owns
campaign execution, retries, provenance, and archival.

## SnO2(110) reference campaign

The included structure generator produces symmetric rutile SnO2(110) starting
models for the O, SnO, and Sn2O terminations.  For nine triple layers their
formulas are Sn18O36, Sn18O34, and Sn18O32, respectively.  The analysis uses
the explicit oxygen chemical-potential convention

`mu_O = E(O2)/2 + Delta mu_O`.

For slabs in bulk SnO2 equilibrium,

`gamma = [E_slab - N_Sn E_bulk - (N_O - 2 N_Sn) mu_O] / (2 A)`.

The corresponding pure analysis helper is
`psteros.surface_energy_oxide_equilibrium`; it returns both eV/A2 and J/m2.
The campaign record in Tessera documents the exact convergence, structure,
hardware, and result evidence for any published value.

## Installation and verification

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -c constraints/aiida-qe-2026-08.txt \
  aiida-core aiida-workgraph aiida-quantumespresso aiida-pseudo pytest
python -m pip install --no-deps -e .
pytest -q tests/unit/test_public_api.py
pytest -q
```

For real QE calculations, register an AiiDA `pw.x` code and an SSSP family in
the active profile.  The maintained Bohr campaign uses QE 7.6 GPU executables
and `SSSP/1.3/PBE/precision`; those are campaign configuration, not hidden
library defaults.

The constraint file is intentional: psteros uses the WorkGraph task-socket
API and therefore verifies a compatible AiiDA/WorkGraph stack rather than
silently replacing packages in an existing scientific profile.

For an established VASP project, install the separately optional adapter with
`python -m pip install -e '.[vasp]'`; it is retained as the secondary backend.

## Compatibility

The broad pre-1.0 VASP builder is intentionally not imported at package level.
It remains available only as `psteros.compat.build_core_workgraph` while an
established VASP project migrates.  New work should use the typed public API.
