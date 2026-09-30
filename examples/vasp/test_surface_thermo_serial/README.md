# Serial Surface Thermodynamics Test

This directory contains test scripts for the experimental serial surface thermodynamics preset.

## Purpose

Test that the flat-graph implementation correctly limits concurrent VASP jobs using `max_concurrent_jobs`.

## Prerequisites

1. Pre-generated slab structures (or modify script to use miller_indices)
2. VASP code configured: `VASP-6.4.1@cluster02`
3. Potential family: `PBE.54`
4. Structure files in a directory

## Running the Test

1. Update `structures_dir` path in `test_serial_preset.py`
2. Ensure AiiDA daemon is running: `verdi daemon start`
3. Run the test:
   ```bash
   source ~/envs/aiida/bin/activate
   python test_serial_preset.py
   ```

## Verification

Monitor execution:
```bash
verdi process list
```

Expected behavior:
- Maximum 2 VASP jobs running simultaneously (due to `max_concurrent_jobs=2`)
- All jobs complete successfully
- Main workgraph returns exit code [0]

Check results:
```bash
verdi process show <PK>
verdi process report <PK>
```

## Success Criteria

1. `verdi process list` never shows more than 2 running VASP jobs
2. All calculations complete successfully
3. Surface energies calculated for all slabs
4. Main node exits with [0]
