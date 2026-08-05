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

## Guided SnO2(110) example

The documentation includes a beginner-friendly walkthrough of the rutile
SnO2(110) surface example. It explains what a surface termination and a
symmetric slab are, then develops the thermodynamic calculation with rendered
equations rather than placing it in this overview.

Read the [SnO2(110) surface-thermodynamics walkthrough](docs/source/examples.rst)
when you are ready for the worked example.

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
