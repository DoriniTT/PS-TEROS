"""End-to-end test of the vibrational contributions on Lovelace (queue par128).

Builds, and with ``--submit`` submits, one graph at a time with the settings of
``PLAN.md``. Only the public ``psteros`` API is used.

    python run_test.py refs                         # build and count, nothing is submitted
    python run_test.py refs --submit
    python run_test.py slabs --refs-pk <REFS> --submit
    python run_test.py vibrations --refs-pk <REFS> --slabs-pk <SLABS> --submit
    python run_test.py ibrion5 --refs-pk <REFS> --submit     # optional O2 cross-check

Every graph runs with ``max_concurrent_jobs=1`` on one Lovelace ``par128`` node
(128 MPI ranks). ``import_sys_environment=False`` is passed through
``CalculationOverride.metadata`` of every label, because ``ExecutionPolicy`` has
no field for it and Lmod functions in the job environment break the ``qstat``
parser.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import psteros

HERE = Path(__file__).resolve().parent

PROFILE = "psteros_vibrations_lovelace"
CODE = "VASP-6.5.1@lovelace"
POTENTIAL_FAMILY = "PBE"
COMPUTER = "lovelace"
QUEUE = "par128"
RANKS = 128
NCORE = 16  # 128 ranks = 8 band groups x 16 cores per band

TERMINATIONS = ("o", "sn2o")
SLABS = tuple(f"slab_{termination}" for termination in TERMINATIONS)
REFERENCES = ("sno2_bulk", "alpha_sn", "o2")
TRIPLE_LAYERS = 3
VACUUM = 10.0
POTCARS = {"Sn": "Sn_d", "O": "O"}
BULK_SUPERCELL = (1, 1, 2)

RELAX_WALLTIME = 7200  # refs and slabs: relaxation and static share one ExecutionPolicy
VIBRATIONS_WALLTIME = 3600

#: The electronic settings of examples/vasp_surface_phase_diagram/campaign.py plus NCORE.
ELECTRONIC = {
    "ENCUT": 400,
    "PREC": "Accurate",
    "EDIFF": 1.0e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
    "NCORE": NCORE,
}
RELAX = {"IBRION": 2, "NSW": 100, "ISIF": 2, "EDIFFG": -0.01}
STATIC = {"IBRION": -1, "NSW": 0}
VIBRATIONS_EDIFF = 1.0e-7

#: Not an ExecutionPolicy field: merged into the job options of the label's override.
NO_SYS_ENVIRONMENT = {"import_sys_environment": False}

DISPLACEMENT = 0.01
EXPECTED_JOBS = {
    "refs": 6,
    "slabs": 4,
    "vibrations": 252,  # 72 + 48 + 72 + 48 + 12
    "ibrion5": 1,
}
EXPECTED_VIBRATION_JOBS = {"slab_o": 72, "slab_sn2o": 48, "sno2_bulk": 72, "alpha_sn": 48, "o2": 12}


def override(**kwargs) -> psteros.CalculationOverride:
    """A CalculationOverride that also switches off ``import_sys_environment``."""

    return psteros.CalculationOverride(metadata=NO_SYS_ENVIRONMENT, **kwargs)


def triplet_o2(**kwargs) -> psteros.CalculationOverride:
    """Triplet O2, Gamma only in the 12 A box."""

    parameters = {"INCAR": {"ISPIN": 2, "MAGMOM": [1.0, 1.0], **kwargs.pop("incar", {})}}
    return override(parameters=parameters, kpoints_distance=10.0, **kwargs)


def execution(walltime: int) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=COMPUTER,
        queue=QUEUE,
        resources={"num_machines": 1, "num_cores_per_machine": RANKS, "num_mpiprocs_per_machine": RANKS},
        max_wallclock_seconds=walltime,
        with_mpi=True,
        max_concurrent_jobs=1,
        custom_scheduler_commands="#PBS -j oe",
    )


def vasp(incar: dict, args) -> psteros.VaspCalculationConfig:
    return psteros.VaspCalculationConfig(
        code_label=args.code,
        incar={**ELECTRONIC, **incar},
        potential_family=args.potential_family,
        potential_mapping=POTCARS,
        kpoints_spacing=0.3,
        max_iterations=3,  # restarts continue a relaxation stopped by the walltime
    )


def recipe(name, incar, walltime, args, overrides) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp", calculation=vasp(incar, args), execution=execution(walltime),
        name=name, role_overrides=overrides,
    )


# ----------------------------------------------------------------------------- graphs


def build_refs(args):
    structures = {
        "sno2_bulk": psteros.rutile_sno2_bulk(),
        "alpha_sn": psteros.alpha_sn_bulk(),
        "o2": psteros.triplet_o2_cell(cell_length=12.0),
    }
    cell_relaxation = override(parameters={"INCAR": {"ISIF": 3}})
    relax = recipe(
        "sno2_refs", RELAX, RELAX_WALLTIME, args,
        {"sno2_bulk": cell_relaxation, "alpha_sn": cell_relaxation, "o2": triplet_o2()},
    )
    static = recipe("sno2_refs_static", STATIC, RELAX_WALLTIME, args, {
        "sno2_bulk": override(), "alpha_sn": override(), "o2": triplet_o2(),
    })
    return psteros.build_relax_static_workgraph(structures, relax, static, submit=args.submit)


def bulk_lattice(refs_graph) -> tuple[float, float]:
    """``(a, c)`` of the relaxed rutile cell, from a finished refs graph."""

    relaxed = refs_graph.outputs.sno2_bulk_relaxed_structure.get_pymatgen_structure()
    lengths = sorted(relaxed.lattice.abc)  # rutile: a = b > c
    return (lengths[1] + lengths[2]) / 2.0, lengths[0]


def unrelaxed_slabs(a: float, c: float) -> dict:
    return {
        f"slab_{termination}": psteros.sno2_110_slab(
            termination=termination, triple_layers=TRIPLE_LAYERS, vacuum_angstrom=VACUUM, a=a, c=c
        )[0]
        for termination in TERMINATIONS
    }


def fixed_sites(slabs: dict) -> dict[str, tuple[int, ...]]:
    """The central triple layer (one Sn2O4 bulk cell) that is frozen in every slab."""

    return {label: tuple(psteros.central_sites(slab, half_width=1.5)) for label, slab in slabs.items()}


def build_slabs(args):
    from aiida import orm

    a, c = bulk_lattice(orm.load_node(args.refs_pk)) if args.refs_pk is not None else (args.a, args.c)
    print(f"slabs on the bulk lattice a={a:.4f} A, c={c:.4f} A")
    slabs = unrelaxed_slabs(a, c)
    fixed = {label: override(fixed_sites=sites) for label, sites in fixed_sites(slabs).items()}
    relax = recipe("sno2_110_slabs", RELAX, RELAX_WALLTIME, args, fixed)
    static = recipe("sno2_110_slabs_static", STATIC, RELAX_WALLTIME, args, {label: override() for label in slabs})
    return psteros.build_relax_static_workgraph(slabs, relax, static, submit=args.submit)


def build_vibrations(args):
    from aiida import orm

    refs = orm.load_node(args.refs_pk) if args.refs_pk is not None else None
    a, c = bulk_lattice(refs) if refs is not None else (args.a, args.c)
    fixed = fixed_sites(unrelaxed_slabs(a, c))  # rebuilt as for the relaxation, on the same lattice
    if args.unrelaxed:  # job count only: the unrelaxed structures of the campaign
        if args.submit:
            raise SystemExit("--unrelaxed only counts the jobs; it cannot be submitted")
        structures = {
            "sno2_bulk": psteros.rutile_sno2_bulk(), "alpha_sn": psteros.alpha_sn_bulk(),
            "o2": psteros.triplet_o2_cell(cell_length=12.0), **unrelaxed_slabs(a, c),
        }
    else:
        _, references = psteros.read_vasp_results(refs, REFERENCES)
        _, slabs = psteros.read_vasp_results(orm.load_node(args.slabs_pk), SLABS)
        structures = {**references, **slabs}
    overrides = {label: override(fixed_sites=sites) for label, sites in fixed.items()}
    overrides.update({"sno2_bulk": override(), "alpha_sn": override(), "o2": triplet_o2()})
    static = recipe(
        "sno2_110_vibrations", {**STATIC, "EDIFF": VIBRATIONS_EDIFF}, VIBRATIONS_WALLTIME, args, overrides
    )
    return psteros.build_vibrations_workgraph(
        structures,
        static,
        psteros.VibrationsConfig(
            displacement_angstrom=DISPLACEMENT,
            supercells={"sno2_bulk": BULK_SUPERCELL},
            molecules=("o2",),
        ),
        submit=args.submit,
    )


def build_ibrion5(args):
    """One VASP IBRION=5 job (NFREE=2, POTIM=0.01) on the relaxed O2: an independent frequency reference."""

    from aiida import orm

    _, references = psteros.read_vasp_results(orm.load_node(args.refs_pk), ("o2",))
    incar = {"IBRION": 5, "NFREE": 2, "POTIM": DISPLACEMENT, "NSW": 1, "EDIFF": VIBRATIONS_EDIFF}
    config = recipe("sno2_o2_ibrion5", incar, VIBRATIONS_WALLTIME, args, {"o2": triplet_o2()})
    return psteros.build_surface_workgraph({"o2": references["o2"]}, config, submit=args.submit)


BUILDERS = {"refs": build_refs, "slabs": build_slabs, "vibrations": build_vibrations, "ibrion5": build_ibrion5}


def count_jobs(graph) -> dict[str, int]:
    """VASP jobs per structure label, from the task names ``<label>_..._vasp``."""

    counts: dict[str, int] = {}
    for task in graph.tasks:
        if not task.name.endswith("_vasp"):
            continue
        label = next((l for l in sorted((*REFERENCES, *SLABS), key=len, reverse=True) if task.name.startswith(l)), "?")
        counts[label] = counts.get(label, 0) + 1
    return counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=tuple(BUILDERS))
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--code", default=CODE)
    parser.add_argument("--potential-family", default=POTENTIAL_FAMILY)
    parser.add_argument("--refs-pk", type=int, help="PK of the finished refs graph (slabs, vibrations, ibrion5)")
    parser.add_argument("--slabs-pk", type=int, help="PK of the finished slabs graph (vibrations)")
    parser.add_argument("--a", type=float, default=4.737, help="bulk a when --refs-pk is not given (build check only)")
    parser.add_argument("--c", type=float, default=3.186, help="bulk c when --refs-pk is not given (build check only)")
    parser.add_argument("--unrelaxed", action="store_true", help="vibrations: count jobs on unrelaxed structures")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and count")
    args = parser.parse_args(argv)

    needs = {"slabs": (), "vibrations": ("slabs_pk",) if not args.unrelaxed else (), "ibrion5": ("refs_pk",)}
    for name in needs.get(args.phase, ()):
        if getattr(args, name) is None:
            parser.error(f"{args.phase} needs --{name.replace('_', '-')}")
    if args.phase == "vibrations" and not args.unrelaxed and args.refs_pk is None:
        parser.error("vibrations needs --refs-pk (or --unrelaxed for a job count)")

    from aiida import load_profile

    load_profile(args.profile)
    graph = BUILDERS[args.phase](args)
    names = [task.name for task in graph.tasks if task.name not in ("graph_inputs", "graph_outputs", "graph_ctx")]
    counts = count_jobs(graph)
    jobs = sum(counts.values())
    print(f"graph {graph.name!r}: {len(names)} tasks, {jobs} VASP jobs, max_concurrent_jobs={graph.max_number_jobs}")
    print("  jobs per label:", dict(sorted(counts.items())))
    if args.phase != "vibrations":
        print("  tasks:", ", ".join(names))
    expected = EXPECTED_JOBS[args.phase]
    status = "matches" if jobs == expected else "DIFFERS FROM"
    print(f"  {status} the plan ({expected} jobs)")
    if args.phase == "vibrations" and counts != EXPECTED_VIBRATION_JOBS:
        print(f"  DIFFERS FROM the plan per label: {EXPECTED_VIBRATION_JOBS}")
        return 1
    if jobs != expected:
        return 1
    if args.submit:
        print(f"submitted: PK={graph.pk}")
        results = HERE / "results"
        results.mkdir(exist_ok=True)
        with (results / "submissions.log").open("a") as log:
            log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {args.phase} PK={graph.pk} profile={args.profile}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
