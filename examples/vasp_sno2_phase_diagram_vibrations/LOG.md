# Log: SnO2(110) with vibrational references (VASP, Lovelace par128)

Plan: `PLAN.md`. Branch: `features/references-vibrational-contributions`.

## Environment

| Item | Value |
|---|---|
| Python / aiida-core / aiida-vasp / aiida-workgraph | 3.13.13 / 2.8.0 / 5.0.0 / 0.8.1 (`~/envs/aiida`, shared) |
| psteros | this worktree, through `PYTHONPATH` for my shell and my daemon; the editable install of the shared env points to another checkout and was not modified |
| AiiDA profile | `psteros_sno2_vibrations` (new, sqlite_dos + RabbitMQ). The default profile was not changed. `psteros_vibrations_lovelace` belongs to another session working on `features/vibrational-contributions` and was not touched |
| Computer / code | `lovelace` and `VASP-6.5.1@lovelace`, copied from profile `presto` (`verdi computer export`, `verdi code export`) |
| POTCARs | family `PBE` (329 POTCARs) exported from `presto` and imported with `verdi archive`; `Sn -> Sn_d`, `O -> O` resolve |
| Queue | `par128`: walltime max 168 h, max 6 queued / 6 running per user. The user already had other jobs there |
| Daemon | `verdi -p psteros_sno2_vibrations daemon start 1` with `PYTHONPATH=<worktree>` |
| Tests before any change | `python -m pytest`: 328 passed, 7 skipped |

Every `verdi` call uses `-p psteros_sno2_vibrations`: the default profile of this machine is somebody else's.

## Timeline

- Scripts: `campaign.py` (phases `o2`, `refs`, `slabs`), built and printed with `submit=False`; INCAR,
  k-points, options, settings, supercells (2x2x3 for SnO2 = 72 atoms, 2x2x2 for alpha-Sn = 64 atoms),
  Gamma-only O2 (`kpoints_spacing=5.0`), `IBRION` 5/6, `NSW=1`, `import_sys_environment=False` checked
  against the tables of `PLAN.md`.
