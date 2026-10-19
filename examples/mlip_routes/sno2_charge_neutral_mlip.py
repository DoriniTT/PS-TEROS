"""SnO2(110), a semiconductor, on the MLIP route: charge-neutral slabs cut by psteros, VASP statics, MACE on a GPU.

``psteros.ChargeNeutralSurfaceStudy`` cuts the charge-neutral (110) slabs and names every calculation;
``psteros.MlipRoute`` runs them on the MLIP route (docs/MLIP_SURFACES.md): VASP statics label a replay fine-tune of
MACE-MH-1, the model relaxes the slabs on a GPU, one active-learning round, a VASP static at each relaxed geometry.
With ``--vibrations`` the model also gives the harmonic modes of every structure, and ``read --temperature T`` draws
the diagram with free energies at T.

The bulk, alpha-Sn and O2 must be relaxed with the same VASP settings first (``build_relax_static_workgraph``); the
built-in starting structures used here are not. Converge the VASP settings for your system.

Usage
-----
    python sno2_charge_neutral_mlip.py show                      # slabs, fixed sites, DFT statics, GPU jobs; no AiiDA
    python sno2_charge_neutral_mlip.py build --profile P --code vasp@cpu-cluster --mace-code mace-python@gpu-cluster \\
        --model /path/on/gpu/mace-mh-1.model --replay replay-data-mh-1-omat-pbe.xyz [--vibrations] [--submit]
    python sno2_charge_neutral_mlip.py read --profile P --pk <PK> [--temperature 800]
"""

from __future__ import annotations

import argparse

import psteros

#: VASP static of every calculation (forces are written by default).
INCAR = {"ENCUT": 520, "PREC": "Accurate", "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05, "LREAL": False,
         "IBRION": -1, "NSW": 0}

#: Half-width (A) of the fixed centre of each slab: the smallest that freezes whole SnO2 units, which the
#: vibrations need (the frozen centre is counted as bulk).
HALF_WIDTHS = (1.5, 2.0, 2.5, 3.0, 3.5)


def study() -> psteros.ChargeNeutralSurfaceStudy:
    return psteros.ChargeNeutralSurfaceStudy(
        psteros.rutile_sno2_bulk(), [(1, 1, 0)],
        references={"Sn": psteros.alpha_sn_bulk(), "O": psteros.triplet_o2_cell()},
    )


def stoichiometric_half_widths(calculations: psteros.ChargeNeutralSurfaceStudy) -> dict[str, float]:
    from collections import Counter

    from pymatgen.core import Composition

    bulk = calculations.bulk.composition.reduced_composition
    widths = {}
    for label, role in calculations.roles.items():
        if role != "slab":
            continue
        slab = calculations.structures[label]
        for width in HALF_WIDTHS:
            frozen = Counter(slab[index].specie.symbol for index in psteros.central_sites(slab, width))
            if frozen and Composition(frozen).reduced_composition == bulk:
                widths[label] = width
                break
        else:
            raise ValueError(f"{label}: no half-width in {HALF_WIDTHS} freezes whole SnO2 units")
    return widths


def route() -> psteros.MlipRoute:
    calculations = study()
    return psteros.MlipRoute(calculations, fixed_half_width_angstrom=stoichiometric_half_widths(calculations))


def recipes(args) -> tuple:
    calculations = route().study
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="sno2_mlip",
        calculation=psteros.VaspCalculationConfig(code_label=args.code, incar=INCAR, potential_family="PBE",
                                                  potential_mapping={"Sn": "Sn_d"}, kpoints_spacing=0.25),
        execution=psteros.ExecutionPolicy(resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
                                          max_wallclock_seconds=4 * 3600, max_concurrent_jobs=2),
        role_overrides=calculations.vasp_overrides("static"),  # triplet O2, Gamma-only molecule
    )
    fine_tune = psteros.MaceFineTuneConfig(args.model, replay_file=args.replay)
    job = psteros.MlipJobConfig(args.mace_code, psteros.ExecutionPolicy(
        queue=args.queue, with_mpi=False, max_wallclock_seconds=12 * 3600,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1},
        extra_options={"environment_variables": {"OMP_NUM_THREADS": "8", "PYTHONUNBUFFERED": "1"}}))
    return static, fine_tune, job


def vibrations() -> psteros.VibrationsConfig:
    # MACE is cheap: a bulk supercell much larger than a DFT one costs little
    return psteros.VibrationsConfig(supercells={"bulk": (3, 3, 4), "ref_Sn": (2, 2, 2)})


def show(args) -> None:
    print(route())


def build(args) -> None:
    static, fine_tune, job = recipes(args)
    sno2 = route()
    print(sno2)
    graph = sno2.build(static, fine_tune, job=job, vibrations=vibrations() if args.vibrations else None,
                       submit=args.submit)
    print(graph.name, len(graph.tasks), "tasks", f"submitted as {graph.pk}" if args.submit else "not submitted")


def read(args) -> None:
    from aiida import orm

    result = route().analyse(orm.load_node(args.pk), temperature_k=args.temperature)
    print(result)
    suffix = "" if args.temperature is None else f"_{args.temperature:g}K"
    result.diagram.plot(f"sno2_110_mlip{suffix}.png", title="SnO$_2$(110), MLIP route")
    result.diagram.to_csv(f"sno2_110_mlip{suffix}.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="phase", required=True)
    sub.add_parser("show")
    first = sub.add_parser("build")
    for name in ("--profile", "--code", "--mace-code", "--model", "--replay"):
        first.add_argument(name, required=True)
    first.add_argument("--queue", default="gpu")
    first.add_argument("--vibrations", action="store_true")
    first.add_argument("--submit", action="store_true")
    second = sub.add_parser("read")
    second.add_argument("--profile", required=True)
    second.add_argument("--pk", type=int, required=True)
    second.add_argument("--temperature", type=float, default=None)
    args = parser.parse_args()
    if args.phase != "show":
        from aiida import load_profile

        load_profile(args.profile)
    {"show": show, "build": build, "read": read}[args.phase](args)


if __name__ == "__main__":
    main()
