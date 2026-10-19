# RuO₂(110) with a fine-tuned MACE model

The MLIP route of psteros (see `docs/source/mlip.rst`) with your own slabs: CP2K static calculations label a
replay fine-tune of MACE-MH-1, which then relaxes the terminations on a GPU. RuO₂ is a metal, so hand-cut slabs are
fine and the graph checks the bulk; for a semiconductor see `examples/mlip_routes/`.

| File | What it does |
|---|---|
| `ruo2_mlip.py labels` | Builds (and with `--submit` submits) the CP2K static graph of bulk RuO₂, hcp Ru, O₂ and the six RuO₂(110) terminations. |
| `ruo2_mlip.py jobs` | Reads the finished graph, checks that the bulk is a metal, and writes the fine-tune and relaxation job folders. |
| `ruo2_mlip_graph.py build` | The whole route as one WorkGraph (`build_mlip_surface_workgraph`): CP2K statics, fine-tune and relaxations on a GPU code, one active-learning round, a CP2K static at each relaxed geometry. |
| `ruo2_mlip_graph.py read` | Reads the energies of that graph from one source (`read_mlip_energies`, `"mlip_dft"`) and the relaxations. |

On the GPU node, in each folder: `python finetune.py`, then `python relax.py` (with the fine-tuned `models/`
next to it). The scripts need only ase, torch and mace-torch.

The CP2K settings are a smoke test. Converge cutoff, basis (DZVP against TZV2P), k-points and vacuum before
interpreting numbers, and take every energy of one phase diagram from one source (`check_one_energy_source`).
