# Results: SrTiO3(001) ternary smoke test of the vibrational contributions on Lovelace

Status: **in progress** (see PLAN.md for the design and the checks T1 to T8).

Same profile, computer, code, queue and settings as `validation/vibrations_lovelace/RESULTS.md` (section 1), including aiida-workgraph 0.9.0
through the overlay and the daemon started from this worktree. `kpoints_spacing=0.04`. POTCARs `Sr_sv`, `Ti_pv`, `O`.

## Process keys

| Graph | PK | State |
|---|---:|---|
| refs (5 solids, 10 jobs) | 6819 | **Finished [0]**, 10 of 10 calculations exit 0 (47 min, mostly queue) |
| slabs (2 terminations, 4 jobs; lattice of 6819, a = 3.908 A) | 6965 | **Finished [0]**, 4 of 4 exit 0 (about 5 h of queue and run time on 12 Oct) |
| vibrations, attempt 1 (282 jobs, one at a time) | 9561 | submitted 12 Oct 13:56; 11 calculations finished by 13 Oct 07:23, with queue waits of 22, 70, 114, 173, 184, 23, 4, 4, 3, 61 and 17 min (mean 61 min; `par128` had 118 waiting and 32 running); the host restarted on 13 Oct about 07:25 and the graph then stalled. Cancelled by me on 13 Oct at 10:00 (the kill of a graph that is still loading its checkpoint has to be forced, and it left one work chain, 10253, and its calculation, 10258, which I cancelled too) |
| vibrations, attempt 2 (282 jobs, `--max-concurrent 8`) | 12807 | submitted 13 Oct 10:14, running; I chose 8 at a time after attempt 1 showed about 12 days for one at a time |

O2 is taken from the SnO2 test: energy of graph 977, vibrations of graph 3400.

## References and slabs (relaxed, static energies)

| Label | E (eV) | cell |
|---|---:|---|
| srtio3_bulk (Sr Ti O3) | -40.0271 | a = 3.908 A (PBE about 3.94: 0.8 % small at ENCUT 400) |
| sr_metal (Sr4, fcc) | -6.5489 | a = 6.017 A |
| ti_metal (Ti2, hcp) | -15.5992 | a = 2.942 A, c = 4.635 A |
| sro_bulk (Sr4 O4) | -48.2419 | a = 5.151 A |
| tio2_bulk (Ti2 O4, rutile) | -53.5035 | a = 4.596 A, c = 2.955 A |
| slab_tio2 (Sr2 Ti3 O8) | -105.0667 | 3.908 x 3.908 x 17.816 A |
| slab_sro (Sr3 Ti2 O7) | -91.2045 | 3.908 x 3.908 x 17.816 A |

k-meshes written: SrTiO3 7x7x7, Sr 5x5x5, SrO 5x5x5, Ti 10x10x6, TiO2 6x6x9, slabs 7x7x2. Slab relaxations use selective dynamics; the
statics do not. Delta H_f(SrTiO3) from these energies and E(O2) = -9.8568 eV: -15.8 eV (raw PBE; experiment about -17.3 eV).
