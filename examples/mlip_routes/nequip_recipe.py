"""The MLIP route with a NequIP or Allegro model (nequip.net) instead of MACE: the recipe and its job, no AiiDA needed.

``psteros.NequipFineTuneConfig`` goes where a ``MaceFineTuneConfig`` goes (``MlipRoute.build``,
``MlipSurfaceStudy.build``, ``build_mlip_surface_workgraph``), with an ``MlipJobConfig`` whose code is the Python of
an environment with nequip (and nequip-allegro) on the GPU computer. Download a package once on a computer with
internet (``load_saved_model("nequip.net:mir-group/NequIP-OAM-L:0.1")`` caches it in ``~/.nequip/model_cache``).

Usage
-----
    python nequip_recipe.py show                  # the nequip.net models and the level-of-theory checks
    python nequip_recipe.py job <folder>          # the fine-tune job folder for a few hand-made frames
"""

from __future__ import annotations

import argparse
from pathlib import Path

import psteros


def show() -> None:
    for model in psteros.NEQUIP_FOUNDATION_MODELS.values():
        print(f"{model.name:14s} {model.architecture:8s} {model.download:42s} {model.licence}  {model.use}")
    recipe = psteros.VaspCalculationConfig("vasp@cluster", incar={"ENCUT": 520, "ISMEAR": 0, "NSW": 0})
    for elements in (["Mg", "O"], ["Co", "O"]):
        messages = psteros.nequip_level_warnings("nequip-oam-l", recipe, elements)
        print(f"\nPBE recipe, {' '.join(elements)}: {messages or 'no warning'}")


def job(folder: Path) -> None:
    from pymatgen.core import Lattice, Structure

    pair = Structure(Lattice.cubic(6.0), ["Mg", "O"], [[0, 0, 0], [0.3, 0, 0]])
    triple = Structure(Lattice.cubic(6.0), ["Mg", "O", "O"], [[0, 0, 0], [0.3, 0, 0], [0.7, 0, 0]])
    frames = [psteros.LabelledStructure("pair", pair, -10.0, [[0.1, 0, 0], [-0.1, 0, 0]]),
              psteros.LabelledStructure("triple", triple, -15.0, [[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]])]
    fine_tune = psteros.NequipFineTuneConfig.for_foundation(
        "nequip-oam-l", "/work/me/models/NequIP-OAM-L-0.1.nequip.zip", replay=False, epochs=30)
    written = psteros.prepare_nequip_finetune_job(folder, frames, fine_tune, name="mgo_round0", device="cuda")
    print(f"{written}: {sorted(path.name for path in written.iterdir())}")
    print("run it on the GPU computer with the nequip Python:  python -u finetune.py")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("show")
    making = sub.add_parser("job")
    making.add_argument("folder", type=Path)
    args = parser.parse_args()
    show() if args.command == "show" else job(args.folder)


if __name__ == "__main__":
    main()
