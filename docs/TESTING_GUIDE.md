# psteros testing

The supported validation path is intentionally short.

```bash
pytest -q tests/unit/test_public_api.py
pytest -q
```

Tier 1 exercises configuration validation, SnO2(110) structure generation,
and thermodynamic analysis without an AiiDA database.  Tier 2 builds QE or
VASP WorkGraphs against an AiiDA test profile.  Tier 3 is the real Bohr A100
campaign and is recorded by the Tessera project rather than represented by
static node identifiers in this repository.

Before a QE campaign is widened, verify all of the following from retrieved
artifacts rather than from scheduler state alone:

- the `pw.x` executable and GPU runtime used by the CalcJob;
- the `aiida-quantumespresso` parser completed successfully;
- the requested SSSP pseudopotentials and input parameters were used;
- SCF, forces, stress, and relaxation criteria meet the plan's acceptance
  criteria.

After changing installed workflow Python, restart the AiiDA daemon before
submitting new work.

