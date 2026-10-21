"""The MLIP route with a UMA model (fairchem) instead of MACE: the recipe and its fine-tune job, no AiiDA needed.

``psteros.UmaFineTuneConfig`` goes where a ``MaceFineTuneConfig`` goes (``MlipRoute.build``, ``MlipSurfaceStudy.build``,
``build_mlip_surface_workgraph``), with an ``MlipJobConfig`` whose code is the Python of a fairchem environment on the
GPU computer. UMA's weights are gated: accept the licence at https://huggingface.co/facebook/UMA, then download the
checkpoint to the GPU computer.

Usage
-----
    python uma_recipe.py show                  # the UMA checkpoints, their tasks and the level-of-theory checks
    python uma_recipe.py job <folder>          # the fine-tune job folder for a few hand-made frames
"""

from __future__ import annotations

import argparse
from pathlib import Path

import psteros


def show() -> None:
    for model in psteros.UMA_FOUNDATION_MODELS.values():
        print(f"{model.name}  ({model.licence})  {model.use}")
        print(f"    {model.url}")
        for task in model.tasks:
            mark = "*" if task.name == model.default_task else " "
            print(f"  {mark} {task.name:6s} {task.functional:8s} {task.training_data}")
    recipe = psteros.VaspCalculationConfig("vasp@cluster", incar={"ENCUT": 520, "ISMEAR": 0, "NSW": 0})
    for elements in (["Ru", "O"], ["Fe", "O"]):
        messages = psteros.uma_level_warnings("uma-s-1p2p1", "omat", recipe, elements)
        print(f"\nPBE recipe, {' '.join(elements)}: {messages or 'no warning'}")


def job(folder: Path) -> None:
    from pymatgen.core import Lattice, Structure

    pair = Structure(Lattice.cubic(6.0), ["Ru", "O"], [[0, 0, 0], [0.3, 0, 0]])
    triple = Structure(Lattice.cubic(6.0), ["Ru", "O", "O"], [[0, 0, 0], [0.3, 0, 0], [0.7, 0, 0]])
    frames = [psteros.LabelledStructure("pair", pair, -10.0, [[0.1, 0, 0], [-0.1, 0, 0]]),
              psteros.LabelledStructure("triple", triple, -15.0, [[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]])]
    fine_tune = psteros.UmaFineTuneConfig.for_foundation("uma-s-1p2p1", "/work/me/models/uma-s-1p2p1.pt",
                                                         replay=False, epochs=20)
    written = psteros.prepare_uma_finetune_job(folder, frames, fine_tune, name="ruo2_round0", device="cuda")
    print(f"{written}: {sorted(path.name for path in written.iterdir())}")
    print("run it on the GPU computer with the fairchem Python:  python -u finetune.py")


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
