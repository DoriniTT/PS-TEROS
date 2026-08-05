# psteros

`psteros` is a reproducible AiiDA workflow library for ab-initio surface
thermodynamics.  Quantum ESPRESSO, through `aiida-quantumespresso`, is the
primary backend.  VASP remains available as a secondary, tested backend for
existing studies.

The public API is deliberately small. Use the code and pseudopotential-family
identifiers registered in your own AiiDA profile, and set an execution policy
that matches your scheduler and resource allocation:

```python
import psteros

recipe = psteros.SurfaceWorkflowConfig(
    backend="qe",
    calculation=psteros.QeCalculationConfig(
        code_label="your-qe-code@your-computer",
        pseudo_family="your-pseudo-family",
        parameters={
            "CONTROL": {"calculation": "relax"},
            "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
            "ELECTRONS": {"conv_thr": 1.0e-8},
        },
    ),
    execution=psteros.ExecutionPolicy(
        computer="your-computer",
        queue="your-scheduler-queue",
        resources={
            "num_machines": 1,  # adjust for your scheduler and calculation
            "num_mpiprocs_per_machine": 1,
        },
        max_wallclock_seconds=86_400,  # adjust for your scheduler policy
        with_mpi=True,  # adjust for your executable and scheduler
        max_concurrent_jobs=1,  # public graph-local bound
    ),
)
bulk = psteros.rutile_sno2_bulk()
graph = psteros.build_surface_workgraph({"bulk_sno2": bulk}, recipe)
```

`ExecutionPolicy` records the graph-local execution bound and the scheduler
settings for your environment. The public builder limits a graph to one active
calculation. The library does not start an external controller; submission,
monitoring, retries, provenance, and archival remain the responsibility of the
user or operator in the selected AiiDA environment.

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
Record the exact convergence, structure, hardware, and result evidence for any
published value in the provenance records of the project that produced it.

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

For real QE calculations, register an AiiDA `pw.x` code and a pseudopotential
family in the active profile. Choose the executable, pseudopotential family,
and scientific parameters for your own project, then record those choices in
its AiiDA provenance.

The constraint file is intentional: psteros uses the WorkGraph task-socket
API and therefore verifies a compatible AiiDA/WorkGraph stack rather than
silently replacing packages in an existing scientific profile.

For an established VASP project, install the separately optional adapter with
`python -m pip install -e '.[vasp]'`; it is retained as the secondary backend.

## Compatibility

The broad pre-1.0 VASP builder is intentionally not imported at package level.
It remains available only as `psteros.compat.build_core_workgraph` while an
established VASP project migrates.  New work should use the typed public API.
