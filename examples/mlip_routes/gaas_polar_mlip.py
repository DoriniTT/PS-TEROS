"""GaAs (111) and (-1-1-1), polar faces, on the MLIP route: passivated slabs cut by psteros, VASP statics, MACE on a GPU.

``psteros.PolarSurfaceStudy`` cuts the polar slabs on a pseudo-hydrogen passivated bottom and adds the
pseudo-molecules; ``psteros.MlipRoute`` runs every calculation on the MLIP route (docs/MLIP_SURFACES.md): VASP statics
label a replay fine-tune of MACE-MH-1, the model relaxes every structure on a GPU (every atom free, as in the polar
method), one active-learning round, a VASP static at each relaxed geometry. ``read`` gives the absolute surface
energies with the bottom check of the MACE-relaxed slabs.

The POTCAR family must contain H.75 and H1.25. The bulk and the elemental references must be relaxed with the same
VASP settings first; the starting structures used here are not. Converge the VASP settings for your system.

Usage
-----
    python gaas_polar_mlip.py show                               # slabs and calculations; no AiiDA
    python gaas_polar_mlip.py build --profile P --code vasp@cpu-cluster --mace-code mace-python@gpu-cluster \\
        --model /path/on/gpu/mace-mh-1.model --replay replay-data-mh-1-omat-pbe.xyz [--submit]
    python gaas_polar_mlip.py read --profile P --pk <PK>
"""

from __future__ import annotations

import argparse

from pymatgen.core import Lattice, Structure

import psteros

#: VASP static of every calculation; the dipole correction of the asymmetric slabs comes from the study.
INCAR = {"ENCUT": 400, "PREC": "Accurate", "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05, "LREAL": False,
         "IBRION": -1, "NSW": 0}


def study() -> psteros.PolarSurfaceStudy:
    """Replace these structures by your relaxed ones."""

    bulk = Structure.from_spacegroup("F-43m", Lattice.cubic(5.75), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])
    gallium = Structure.from_spacegroup("Cmce", Lattice.orthorhombic(4.52, 7.66, 4.53), ["Ga"], [[0, 0.1549, 0.081]])
    arsenic = Structure.from_spacegroup("R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]])
    return psteros.PolarSurfaceStudy(bulk, faces=[(1, 1, 1), (-1, -1, -1)], references={"Ga": gallium, "As": arsenic},
                                     bilayers=9)


def route() -> psteros.MlipRoute:
    return psteros.MlipRoute(study())


def recipes(args, calculations: psteros.PolarSurfaceStudy) -> tuple:
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="gaas_mlip",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code, incar=INCAR, potential_family="PBE", kpoints_spacing=0.2,
            potential_mapping=calculations.potential_mapping({"Ga": "Ga_d", "As": "As"})),
        execution=psteros.ExecutionPolicy(resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
                                          max_wallclock_seconds=4 * 3600, max_concurrent_jobs=2),
        role_overrides=calculations.vasp_overrides("static"),
    )
    fine_tune = psteros.MaceFineTuneConfig(args.model, replay_file=args.replay)
    job = psteros.MlipJobConfig(args.mace_code, psteros.ExecutionPolicy(
        queue=args.queue, with_mpi=False, max_wallclock_seconds=12 * 3600,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1},
        extra_options={"environment_variables": {"OMP_NUM_THREADS": "8", "PYTHONUNBUFFERED": "1"}}))
    return static, fine_tune, job


def show(args) -> None:
    gaas = route()
    print(gaas)
    print(f"{len(gaas.structures)} calculations: {', '.join(gaas.structures)}")


def build(args) -> None:
    gaas = route()
    static, fine_tune, job = recipes(args, gaas.study)
    print(gaas)
    graph = gaas.build(static, fine_tune, job=job, submit=args.submit)
    print(graph.name, len(graph.tasks), "tasks", f"submitted as {graph.pk}" if args.submit else "not submitted")


def read(args) -> None:
    from aiida import orm

    result = route().analyse(orm.load_node(args.pk))  # PolarStudyResult, bottom checks of the MACE-relaxed slabs
    print(result)
    for face, report in result.bottom_checks.items():
        print(f"\n{face}\n{report}")
    result.diagram.plot("gaas_111_mlip.png", title="GaAs(111), MLIP route")
    result.diagram.to_csv("gaas_111_mlip.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="phase", required=True)
    sub.add_parser("show")
    first = sub.add_parser("build")
    for name in ("--profile", "--code", "--mace-code", "--model", "--replay"):
        first.add_argument(name, required=True)
    first.add_argument("--queue", default="gpu")
    first.add_argument("--submit", action="store_true")
    second = sub.add_parser("read")
    second.add_argument("--profile", required=True)
    second.add_argument("--pk", type=int, required=True)
    args = parser.parse_args()
    if args.phase != "show":
        from aiida import load_profile

        load_profile(args.profile)
    {"show": show, "build": build, "read": read}[args.phase](args)


if __name__ == "__main__":
    main()
