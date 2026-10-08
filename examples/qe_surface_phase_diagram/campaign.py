"""SnO2(110) surface phase diagram with Quantum ESPRESSO: the calculations.

A small, end-to-end campaign built only on the public psteros API.  It is
sized as a smoke test for a short debug queue (<= 5 cores, <= 5 GB,
<= 20 min per job, one queued job per user): 3 triple-layer slabs,
10 A vacuum, 40/320 Ry and coarse k-points.  Converge these before drawing
physical conclusions.

Phase ``refs``   bulk rutile SnO2 (vc-relax), alpha-Sn (vc-relax) and triplet
                 O2 (relax), each followed by a static SCF.
Phase ``slabs``  the o / sno / sn2o terminations of SnO2(110) built on the
                 *relaxed* bulk lattice, central triple layer frozen,
                 relax -> static SCF.

Usage
-----
    python campaign.py refs  --profile P --code pw@my-cluster --pseudo-family SSSP/1.3/PBE/efficiency --submit
    python campaign.py slabs --profile P --code pw@my-cluster --pseudo-family SSSP/1.3/PBE/efficiency \\
        --refs-pk <PK of the finished refs graph> --submit

Then run ``phase_diagram.py`` on the two graph PKs.
"""

from __future__ import annotations

import argparse

import psteros

TERMINATIONS = ("o", "sno", "sn2o")
TRIPLE_LAYERS = 3
VACUUM = 10.0

SYSTEM = {
    "ecutwfc": 40.0,
    "ecutrho": 320.0,
    "occupations": "smearing",
    "smearing": "mv",
    "degauss": 0.01,
}
ELECTRONS = {"conv_thr": 1.0e-7, "mixing_beta": 0.4, "electron_maxstep": 150}
RELAX_CONTROL = {
    "calculation": "relax",
    "forc_conv_thr": 1.0e-3,
    "etot_conv_thr": 1.0e-4,
    "nstep": 80,
    "tprnfor": True,
    "tstress": True,
}
STATIC_CONTROL = {"calculation": "scf", "tprnfor": True, "tstress": True}

VC_RELAX = psteros.CalculationOverride(
    parameters={"CONTROL": {"calculation": "vc-relax"}, "CELL": {"press_conv_thr": 0.5}}
)
TRIPLET_O2 = psteros.CalculationOverride(
    parameters={"SYSTEM": {"nspin": 2, "tot_magnetization": 2, "starting_magnetization": {"O": 0.5}}},
    kpoints_distance=2.0,  # Gamma only in the 12 A box
)


def execution(args) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=args.computer,
        queue=args.queue,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": 4, "num_cores_per_machine": 4},
        max_wallclock_seconds=args.walltime,
        with_mpi=True,
        max_concurrent_jobs=1,
    )


def recipes(args, name, relax_overrides=None, static_overrides=None, max_iterations=3):
    """Relaxation and static recipes sharing code, pseudopotentials, cutoffs and k-points."""

    def qe(control):
        return psteros.QeCalculationConfig(
            code_label=args.code,
            pseudo_family=args.pseudo_family,
            parameters={"CONTROL": control, "SYSTEM": SYSTEM, "ELECTRONS": ELECTRONS},
            kpoints_distance=0.5,
            # Restarts let PwBaseWorkChain continue a relaxation stopped by the walltime.
            max_iterations=max_iterations,
        )

    relax = psteros.SurfaceWorkflowConfig(
        backend="qe", calculation=qe(RELAX_CONTROL), execution=execution(args),
        name=name, role_overrides=relax_overrides or {},
    )
    static = psteros.SurfaceWorkflowConfig(
        backend="qe", calculation=qe(STATIC_CONTROL), execution=execution(args),
        name=f"{name}_static", role_overrides=static_overrides or {},
    )
    return relax, static


def central_layer_fixed(slab, half_width: float = 1.5) -> psteros.CalculationOverride:
    """Freeze the atoms within ``half_width`` A of the slab mid-plane (the central triple layer)."""

    heights = [site.coords.dot(slab.normal) for site in slab]
    middle = (max(heights) + min(heights)) / 2.0
    fixed = [index for index, height in enumerate(heights) if abs(height - middle) < half_width]
    return psteros.CalculationOverride(
        settings={"FIXED_COORDS": psteros.qe_fixed_coordinate_flags(len(slab), fixed)}
    )


def build_refs(args):
    structures = {
        "sno2_bulk": psteros.rutile_sno2_bulk(),
        "alpha_sn": psteros.alpha_sn_bulk(),
        "o2": psteros.triplet_o2_cell(cell_length=12.0),
    }
    relax, static = recipes(
        args, "sno2_refs",
        relax_overrides={"sno2_bulk": VC_RELAX, "alpha_sn": VC_RELAX, "o2": TRIPLET_O2},
        static_overrides={"o2": TRIPLET_O2},
    )
    return psteros.build_qe_relax_static_workgraph(structures, relax, static, submit=args.submit)


def relaxed_bulk_lattice(refs_pk: int) -> tuple[float, float]:
    """Return (a, c) of the vc-relaxed rutile cell from a finished refs graph."""
    from aiida import orm

    relaxed = orm.load_node(refs_pk).outputs.sno2_bulk_relaxed_structure.get_pymatgen_structure()
    lengths = sorted(relaxed.lattice.abc)  # rutile: a = b > c
    return (lengths[1] + lengths[2]) / 2.0, lengths[0]


def build_slabs(args):
    a, c = relaxed_bulk_lattice(args.refs_pk) if args.refs_pk is not None else (args.a, args.c)
    print(f"building slabs on bulk lattice a={a:.4f} A, c={c:.4f} A")
    structures = {
        f"slab_{termination}": psteros.sno2_110_slab(
            termination=termination, triple_layers=TRIPLE_LAYERS, vacuum_angstrom=VACUUM, a=a, c=c
        )[0]
        for termination in TERMINATIONS
    }
    fixed = {label: central_layer_fixed(slab) for label, slab in structures.items()}
    relax, static = recipes(args, "sno2_110_slabs", relax_overrides=fixed, max_iterations=4)
    return psteros.build_qe_relax_static_workgraph(structures, relax, static, submit=args.submit)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("refs", "slabs"))
    parser.add_argument("--profile", required=True)
    parser.add_argument("--code", required=True, help="AiiDA label of a quantumespresso.pw code")
    parser.add_argument("--pseudo-family", required=True)
    parser.add_argument("--computer", help="name of your AiiDA computer (for the record)")
    parser.add_argument("--queue", help="queue or partition; the scheduler's default when left out")
    parser.add_argument("--walltime", type=int, default=1140, help="seconds per job (a 20-minute queue by default)")
    parser.add_argument("--refs-pk", type=int, help="slabs phase: PK of the finished refs graph")
    parser.add_argument("--a", type=float, default=4.737, help="bulk a when --refs-pk is not given")
    parser.add_argument("--c", type=float, default=3.186, help="bulk c when --refs-pk is not given")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and list the graph")
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    graph = build_refs(args) if args.phase == "refs" else build_slabs(args)
    tasks = [task.name for task in graph.tasks if task.name not in ("graph_inputs", "graph_outputs", "graph_ctx")]
    print(f"graph {graph.name!r}: {len(tasks)} tasks, one active job at a time")
    print("  " + ", ".join(tasks))
    if args.submit:
        print(f"submitted: PK={graph.pk}")
    return graph


if __name__ == "__main__":
    main()
