# The MLIP route with slabs cut by psteros

`psteros.MlipRoute` runs a study of psteros on the MLIP route (see `docs/MLIP_SURFACES.md`): DFT statics label a
replay fine-tune of MACE-MH-1, the model relaxes every slab on a GPU, one active-learning round adds DFT labels near
the relaxed geometries, and a DFT static at each relaxed geometry gives the energies of the diagram. The studies cut
slabs with the right electron count, so these examples are semiconductors.

| File | What it does |
|---|---|
| `sno2_charge_neutral_mlip.py show` | Cuts the charge-neutral SnO₂(110) slabs (`ChargeNeutralSurfaceStudy`) and prints the route; no AiiDA needed. |
| `sno2_charge_neutral_mlip.py build` | The route as one WorkGraph with VASP statics; `--vibrations` adds the model's harmonic modes (large bulk supercells, no extra DFT). |
| `sno2_charge_neutral_mlip.py read` | The phase diagram from the final VASP energies; `--temperature T` uses the free energies. |
| `gaas_polar_mlip.py show / build / read` | The polar GaAs (111) and (-1-1-1) faces on a pseudo-hydrogen passivated bottom (`PolarSurfaceStudy`); `read` gives the absolute surface energies and the bottom check of the MACE-relaxed slabs. |

Each slab of the SnO₂ example keeps a centre of whole SnO₂ units in place (`fixed_half_width_angstrom` per slab):
the vibrations count that frozen centre as bulk. Polar slabs relax every atom, as in the polar method.

The bulk and the references must be relaxed with the same DFT settings first, and the VASP settings here are not
converged: converge them for your system and check a new system against a few DFT relaxations.
