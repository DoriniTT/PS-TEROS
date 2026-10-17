# RuO₂(110) with a fine-tuned MACE model

The opt-in MLIP route of psteros (metallic systems only; see `docs/source/mlip.rst`): CP2K static calculations
label a replay fine-tune of MACE-MH-1, which then relaxes the terminations on a GPU.

| File | What it does |
|---|---|
| `ruo2_mlip.py labels` | Builds (and with `--submit` submits) the CP2K static graph of bulk RuO₂, hcp Ru, O₂ and the six RuO₂(110) terminations. |
| `ruo2_mlip.py jobs` | Reads the finished graph, checks that the bulk is a metal, and writes the fine-tune and relaxation job folders. |

On the GPU node, in each folder: `python finetune.py`, then `python relax.py` (with the fine-tuned `models/`
next to it). The scripts need only ase, torch and mace-torch.

The CP2K settings are a smoke test. Converge cutoff, basis (DZVP against TZV2P), k-points and vacuum before
interpreting numbers, and take every energy of one phase diagram from one source (`check_one_energy_source`).
