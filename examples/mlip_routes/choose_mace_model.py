"""Which MACE foundation model for the MLIP route? The catalogue, and the heads that match a DFT recipe.

``psteros.MACE_FOUNDATION_MODELS`` lists the released MACE models the route has been run with: their files,
heads (as the files name them), the DFT each head was trained on, licence and replay data. A recipe built with
``MaceFineTuneConfig.for_foundation`` is checked against the DFT of the study: the builders warn when the
functional or the +U differs (``psteros.foundation_level_warnings``). No AiiDA and no GPU needed here.

Usage
-----
    python choose_mace_model.py list                                  # every model and head
    python choose_mace_model.py check --elements Fe O                 # PBE without +U (the default VASP recipe)
    python choose_mace_model.py check --elements Fe O --ldau          # PBE+U
    python choose_mace_model.py check --elements Sn O --metagga R2SCAN
"""

from __future__ import annotations

import argparse

import psteros


def show_catalogue() -> None:
    for model in psteros.MACE_FOUNDATION_MODELS.values():
        print(f"{model.name}  ({model.licence}, mace-torch >= {model.min_mace_torch})  {model.use}")
        print(f"    file   {model.url}")
        print(f"    replay {model.replay_url}")
        for head in model.heads:
            plus_u = {"mp": "Materials Project +U", "none": "no +U", None: "+U not stated"}[head.hubbard]
            mark = "*" if head.name == model.default_head else " "
            slabs = "" if head.for_slabs else "  (molecules only)"
            print(f"  {mark} {head.name:18s} {head.functional:13s} {plus_u:22s} {head.training_data}{slabs}")


def check(elements: list[str], ldau: bool, metagga: str | None) -> None:
    incar = {"ENCUT": 520, "ISMEAR": 0, "SIGMA": 0.05, "IBRION": -1, "NSW": 0}
    if ldau:
        incar.update(LDAU=True, LDAUTYPE=2)
    if metagga:
        incar["METAGGA"] = metagga
    recipe = psteros.VaspCalculationConfig("vasp@cluster", incar=incar)
    functional, uses_u = psteros.recipe_level_of_theory(recipe)
    print(f"recipe: {functional}, {'with' if uses_u else 'without'} +U; elements {' '.join(elements)}\n")
    for model in psteros.MACE_FOUNDATION_MODELS.values():
        for head in model.heads:
            messages = psteros.foundation_level_warnings(model.name, head.name, recipe, elements)
            print(f"{'ok  ' if not messages else 'warn'} {model.name} ({head.name})")
            for message in messages:
                print(f"       {message}")
    config = psteros.MaceFineTuneConfig.for_foundation("mace-matpes-pbe-0", "/work/me/models/MACE-matpes-pbe-omat-ft.model",
                                                       replay_file="matpes-pbe-replay-data.xyz")
    print(f"\nfor example: MaceFineTuneConfig.for_foundation(...) -> head {config.foundation_head!r}, "
          f"checksum {config.foundation_model_sha256[:12]}...")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    checking = sub.add_parser("check")
    checking.add_argument("--elements", nargs="+", required=True)
    checking.add_argument("--ldau", action="store_true")
    checking.add_argument("--metagga")
    args = parser.parse_args()
    if args.command == "list":
        show_catalogue()
    else:
        check(args.elements, args.ldau, args.metagga)


if __name__ == "__main__":
    main()
