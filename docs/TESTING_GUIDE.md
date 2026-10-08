# psteros testing

The supported validation path is intentionally short.

```bash
pytest -q tests/unit/test_public_api.py
pytest -q
```

Tier 1 exercises configuration validation, structure generation and
thermodynamic analysis without an AiiDA database.  Tier 2 builds VASP or QE
WorkGraphs against an AiiDA test profile.  Tier 3 runs real calculations on
your own computer and is not part of the automated suite.

Before a campaign is widened, verify all of the following from retrieved
artifacts rather than from scheduler state alone:

- the executable and runtime used by the CalcJob;
- the aiida-vasp (or aiida-quantumespresso) parser completed successfully;
- the requested POTCARs (or pseudopotentials) and input parameters were used;
- SCF, forces, stress, and relaxation criteria meet the plan's acceptance
  criteria.

After changing installed workflow Python, restart the AiiDA daemon before
submitting new work.

