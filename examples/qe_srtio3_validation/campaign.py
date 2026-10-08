"""SrTiO3(001) validation campaign with Quantum ESPRESSO (PBE): the calculations.

Validates the psteros ternary phase diagram against experiment and literature
(see README.md).  Three graphs, run one after the other on one node:

``refs``       SrTiO3, alpha-Sr, alpha-Ti, SrO, rutile and anatase TiO2
               (vc-relax -> static SCF) and triplet O2 (relax -> static SCF).
``slabs``      7-layer symmetric SrO- and TiO2-terminated SrTiO3(001) 1x1 slabs
               on the relaxed bulk lattice, central layer fixed, relax -> static.
``unrelaxed``  static SCF of the same slabs before relaxation (cleavage energy).

Usage
-----
    python campaign.py refs      --profile P --code pw@my-cluster --queue my-queue --submit
    python campaign.py slabs     --profile P --code pw@my-cluster --refs-pk <REFS_PK> --submit
    python campaign.py unrelaxed --profile P --code pw@my-cluster --refs-pk <REFS_PK> --submit

Set the MPI ranks, wall time and queue for your computer. For a QE wrapper that
launches MPI itself, pass ``--no-mpi`` and the lines it needs with ``--prepend``.
"""

from __future__ import annotations

import argparse

import psteros

LAYERS = 7
VACUUM = 22.0  # A; with kpoints_distance 0.2 the vacuum direction gets one k-point
TERMINATIONS = ("SrO", "TiO2")

SYSTEM = {
    "ecutwfc": 60.0,  # SSSP 1.3 efficiency recommends 50/400 Ry for Sr, Ti, O
    "ecutrho": 480.0,
    "occupations": "smearing",
    "smearing": "mv",
    "degauss": 0.01,
}
ELECTRONS = {"conv_thr": 1.0e-9, "mixing_beta": 0.4, "electron_maxstep": 200}
RELAX_CONTROL = {
    "calculation": "relax",
    "forc_conv_thr": 1.0e-3,
    "etot_conv_thr": 1.0e-5,
    "nstep": 150,
    "tprnfor": True,
    "tstress": True,
}
STATIC_CONTROL = {"calculation": "scf", "tprnfor": True, "tstress": True}


def reference_structures() -> dict:
    """Experimental-neighbourhood starting cells; every one is vc-relaxed."""
    from pymatgen.core import Lattice, Structure

    return {
        "srtio3": Structure.from_spacegroup(
            "Pm-3m", Lattice.cubic(3.905), ["Sr", "Ti", "O"], [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0]]
        ),
        "sr_metal": Structure.from_spacegroup("Fm-3m", Lattice.cubic(6.08), ["Sr"], [[0, 0, 0]]).get_primitive_structure(),
        "ti_metal": Structure.from_spacegroup(
            "P6_3/mmc", Lattice.hexagonal(2.951, 4.686), ["Ti"], [[1 / 3, 2 / 3, 0.25]]
        ),
        "sro": Structure.from_spacegroup(
            "Fm-3m", Lattice.cubic(5.16), ["Sr", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        ).get_primitive_structure(),
        "tio2_rutile": Structure.from_spacegroup(
            "P4_2/mnm", Lattice.tetragonal(4.594, 2.959), ["Ti", "O"], [[0, 0, 0], [0.3048, 0.3048, 0]]
        ),
        # I4_1/amd, origin choice 1: Ti 4a (0, 0, 0), O 8e (0, 0, z)
        "tio2_anatase": Structure.from_spacegroup(
            "I4_1/amd", Lattice.tetragonal(3.785, 9.514), ["Ti", "O"], [[0, 0, 0], [0, 0, 0.2081]]
        ).get_primitive_structure(),
        "o2": psteros.triplet_o2_cell(cell_length=14.0),
    }


def srtio3_001_slab(a: float, termination: str, *, layers: int = LAYERS, vacuum: float = VACUUM):
    """Symmetric 1x1 SrTiO3(001) slab of alternating SrO and TiO2 planes.

    ``SrO`` gives Sr(n+1)Ti(n)O(3n+1) and ``TiO2`` gives Sr(n)Ti(n+1)O(3n+2) for
    ``layers = 2n + 1``, the models of Eglitis and Vanderbilt, PRB 77, 195408.
    """
    from pymatgen.core import Lattice, Structure

    if termination not in TERMINATIONS:
        raise ValueError(f"termination must be one of {TERMINATIONS}")
    if layers < 3 or layers % 2 == 0:
        raise ValueError("layers must be odd and at least 3")
    spacing = a / 2.0
    sro = [("Sr", 0.0, 0.0), ("O", 0.5, 0.5)]
    tio2 = [("Ti", 0.5, 0.5), ("O", 0.5, 0.0), ("O", 0.0, 0.5)]
    outer, inner = (sro, tio2) if termination == "SrO" else (tio2, sro)
    species, coords = [], []
    for index in range(layers):
        for element, fx, fy in outer if index % 2 == 0 else inner:
            species.append(element)
            coords.append([fx * a, fy * a, vacuum / 2.0 + index * spacing])
    return Structure(
        Lattice.tetragonal(a, (layers - 1) * spacing + vacuum), species, coords, coords_are_cartesian=True
    )


def execution(args) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=args.computer,
        queue=args.queue,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": args.ranks},
        max_wallclock_seconds=args.walltime,
        with_mpi=not args.no_mpi,
        prepend_text=args.prepend,
        max_concurrent_jobs=1,
    )


def recipe(args, name, control, overrides) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="qe",
        calculation=psteros.QeCalculationConfig(
            code_label=args.code,
            pseudo_family=args.pseudo_family,
            parameters={"CONTROL": control, "SYSTEM": SYSTEM, "ELECTRONS": ELECTRONS},
            kpoints_distance=0.2,
            max_iterations=3,
        ),
        execution=execution(args),
        name=name,
        role_overrides=overrides,
    )


def pools(count: int, **extra) -> dict:
    return {"CMDLINE": ["-nk", str(count)], **extra}


def build_refs(args):
    structures = reference_structures()
    vc_relax = {"CONTROL": {"calculation": "vc-relax"}, "CELL": {"press_conv_thr": 0.2, "cell_dofree": "all"}}
    triplet = {"SYSTEM": {"nspin": 2, "tot_magnetization": 2, "starting_magnetization": {"O": 0.5}}}
    relax, static = {}, {}
    for label in structures:
        if label == "o2":
            relax[label] = static[label] = psteros.CalculationOverride(
                parameters=triplet, kpoints_distance=2.0, settings=pools(2)  # Gamma, two spin channels
            )
            continue
        dense = 0.15 if label.endswith("_metal") else None  # metals need denser sampling
        relax[label] = psteros.CalculationOverride(parameters=vc_relax, kpoints_distance=dense, settings=pools(4))
        static[label] = psteros.CalculationOverride(kpoints_distance=dense, settings=pools(4))
    return psteros.build_qe_relax_static_workgraph(
        structures,
        recipe(args, "srtio3_refs", RELAX_CONTROL, relax),
        recipe(args, "srtio3_refs_static", STATIC_CONTROL, static),
        submit=args.submit,
    )


def relaxed_lattice_constant(refs_pk: int) -> float:
    from aiida import orm

    structure = orm.load_node(refs_pk).outputs.srtio3_relaxed_structure.get_pymatgen_structure()
    return sum(structure.lattice.abc) / 3.0


def slabs(args) -> dict:
    a = relaxed_lattice_constant(args.refs_pk) if args.refs_pk is not None else args.a
    print(f"building {LAYERS}-layer slabs on a = {a:.4f} A")
    return {f"slab_{t.lower()}": srtio3_001_slab(a, t) for t in TERMINATIONS}


def central_layer_fixed(slab):
    heights = [site.coords[2] for site in slab]
    middle = (max(heights) + min(heights)) / 2.0
    fixed = [index for index, height in enumerate(heights) if abs(height - middle) < 0.3]
    return psteros.qe_fixed_coordinate_flags(len(slab), fixed)


def build_slabs(args):
    structures = slabs(args)
    relax = {
        label: psteros.CalculationOverride(settings=pools(4, FIXED_COORDS=central_layer_fixed(slab)))
        for label, slab in structures.items()
    }
    static = {label: psteros.CalculationOverride(settings=pools(4)) for label in structures}
    return psteros.build_qe_relax_static_workgraph(
        structures,
        recipe(args, "srtio3_001_slabs", RELAX_CONTROL, relax),
        recipe(args, "srtio3_001_slabs_static", STATIC_CONTROL, static),
        submit=args.submit,
    )


def build_unrelaxed(args):
    structures = slabs(args)
    overrides = {label: psteros.CalculationOverride(settings=pools(4)) for label in structures}
    return psteros.build_surface_workgraph(
        structures, recipe(args, "srtio3_001_unrelaxed", STATIC_CONTROL, overrides), submit=args.submit
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("refs", "slabs", "unrelaxed"))
    parser.add_argument("--profile", required=True)
    parser.add_argument("--code", required=True, help="AiiDA label of your quantumespresso.pw code")
    parser.add_argument("--pseudo-family", default="SSSP/1.3/PBE/efficiency")
    parser.add_argument("--computer", help="name of your AiiDA computer (for the record)")
    parser.add_argument("--queue", help="queue or partition; the scheduler's default when left out")
    parser.add_argument("--ranks", type=int, default=4, help="MPI ranks; must be divisible by the pools (4, 2)")
    parser.add_argument("--walltime", type=int, default=4 * 3600, help="seconds per job")
    parser.add_argument("--no-mpi", action="store_true", help="the code (a wrapper) launches MPI itself")
    parser.add_argument("--prepend", default="", help="shell lines before the executable, e.g. module loads")
    parser.add_argument("--refs-pk", type=int, help="slabs/unrelaxed: PK of the finished refs graph")
    parser.add_argument("--a", type=float, default=3.94, help="lattice constant when --refs-pk is not given")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    graph = {"refs": build_refs, "slabs": build_slabs, "unrelaxed": build_unrelaxed}[args.phase](args)
    tasks = [t.name for t in graph.tasks if t.name not in ("graph_inputs", "graph_outputs", "graph_ctx")]
    print(f"graph {graph.name!r}: {len(tasks)} tasks, one active job at a time")
    print("  " + ", ".join(tasks))
    if args.submit:
        print(f"submitted: PK={graph.pk}")
    return graph


if __name__ == "__main__":
    main()
