"""Ternary-oxide smoke test of the vibrational contributions on Lovelace: SrTiO3(001). See PLAN.md.

    python run_ternary.py refs [--submit]
    python run_ternary.py slabs --refs-pk <REFS> [--submit]
    python run_ternary.py vibrations --refs-pk <REFS> --slabs-pk <SLABS> [--submit]

Same cluster, code, queue, electronic settings and ``import_sys_environment=False`` as
``validation/vibrations_lovelace/run_test.py``, whose helpers are imported. Structures are built with pymatgen (psteros has no
SrTiO3 builders); everything that is calculated goes through the public ``psteros`` API.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure

import psteros

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "vibrations_lovelace"))
import run_test as sno2  # noqa: E402  (settings and helpers of the SnO2 test)

POTCARS = {"Sr": "Sr_sv", "Ti": "Ti_pv", "O": "O"}
BULK_LABELS = ("srtio3_bulk", "sr_metal", "ti_metal", "sro_bulk", "tio2_bulk")
SLAB_LABELS = ("slab_tio2", "slab_sro")
VACUUM = 10.0
BULK_SUPERCELLS = {"srtio3_bulk": (1, 1, 2), "ti_metal": (1, 1, 2)}  # the others use their own cell
EXPECTED_JOBS = {"refs": 10, "slabs": 4, "vibrations": 282}
EXPECTED_VIBRATION_JOBS = {
    "srtio3_bulk": 60, "tio2_bulk": 36, "sro_bulk": 48, "sr_metal": 24, "ti_metal": 24, "slab_tio2": 48, "slab_sro": 42,
}


# ----------------------------------------------------------------------------- structures


def reference_structures() -> dict[str, Structure]:
    """Starting structures near the experimental lattices; the cells are relaxed (ISIF=3)."""

    return {
        "srtio3_bulk": Structure.from_spacegroup(
            "Pm-3m", Lattice.cubic(3.905), ["Sr", "Ti", "O"], [[0.5, 0.5, 0.5], [0, 0, 0], [0.5, 0, 0]]
        ),
        "sr_metal": Structure.from_spacegroup("Fm-3m", Lattice.cubic(6.085), ["Sr"], [[0, 0, 0]]),
        "ti_metal": Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.951, 4.686), ["Ti"], [[1 / 3, 2 / 3, 0.25]]),
        "sro_bulk": Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.16), ["Sr", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
        "tio2_bulk": Structure.from_spacegroup(
            "P4_2/mnm", Lattice.tetragonal(4.594, 2.959), ["Ti", "O"], [[0, 0, 0], [0.305, 0.305, 0]]
        ),
    }


#: In-plane positions (fractions of a) of each layer type of the perovskite with Ti at the origin.
_LAYER = {
    "TiO2": [("Ti", (0.0, 0.0)), ("O", (0.5, 0.0)), ("O", (0.0, 0.5))],
    "SrO": [("Sr", (0.5, 0.5)), ("O", (0.0, 0.0))],
}
_STACKS = {
    "slab_tio2": ("TiO2", "SrO", "TiO2", "SrO", "TiO2"),
    "slab_sro": ("SrO", "TiO2", "SrO", "TiO2", "SrO"),
}
FROZEN_LAYERS = (1, 2)  # one SrO and one TiO2 layer in the middle: one SrTiO3 formula unit, same in both stacks


def perovskite_slab(label: str, a: float, vacuum: float = VACUUM) -> tuple[Structure, tuple[int, ...]]:
    """SrTiO3(001) slab (1x1, five layers of spacing a/2) centred in the cell, and the indices of the frozen sites."""

    stack = _STACKS[label]
    height = (len(stack) - 1) * a / 2.0
    c = height + vacuum
    species, coords, layer_of = [], [], []
    for index, kind in enumerate(stack):
        z = (vacuum / 2.0 + index * a / 2.0) / c
        for element, (x, y) in _LAYER[kind]:
            species.append(element)
            coords.append([x, y, z])
            layer_of.append(index)
    slab = Structure(Lattice.from_parameters(a, a, c, 90, 90, 90), species, coords)
    frozen = tuple(i for i, layer in enumerate(layer_of) if layer in FROZEN_LAYERS)
    return slab, frozen


def relaxed_lattice(refs_graph) -> float:
    """Cubic lattice parameter a of the relaxed SrTiO3 bulk."""

    structure = refs_graph.outputs.srtio3_bulk_relaxed_structure.get_pymatgen_structure()
    return float(np.mean(structure.lattice.abc))


# ----------------------------------------------------------------------------- recipes


def recipe(name: str, incar: dict, walltime: int, args, overrides) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code, incar={**sno2.ELECTRONIC, **incar}, potential_family=args.potential_family,
            potential_mapping=POTCARS, kpoints_spacing=args.kpoints_spacing, max_iterations=3,
        ),
        execution=dataclasses.replace(sno2.execution(walltime), max_concurrent_jobs=args.max_concurrent),
        name=name, role_overrides=overrides,
    )


def build_refs(args):
    structures = reference_structures()
    cell = sno2.override(parameters={"INCAR": {"ISIF": 3}})
    relax = recipe("srtio3_refs", sno2.RELAX, sno2.RELAX_WALLTIME, args, {label: cell for label in structures})
    static = recipe("srtio3_refs_static", sno2.STATIC, sno2.RELAX_WALLTIME, args, {l: sno2.override() for l in structures})
    return psteros.build_relax_static_workgraph(structures, relax, static, submit=args.submit)


def build_slabs(args):
    from aiida import orm

    a = relaxed_lattice(orm.load_node(args.refs_pk)) if args.refs_pk is not None else args.a
    print(f"slabs on the SrTiO3 lattice a = {a:.4f} A")
    slabs, fixed = {}, {}
    for label in SLAB_LABELS:
        slabs[label], frozen = perovskite_slab(label, a)
        fixed[label] = sno2.override(fixed_sites=frozen)
    relax = recipe("srtio3_slabs", sno2.RELAX, sno2.RELAX_WALLTIME, args, fixed)
    static = recipe("srtio3_slabs_static", sno2.STATIC, sno2.RELAX_WALLTIME, args, {l: sno2.override() for l in slabs})
    return psteros.build_relax_static_workgraph(slabs, relax, static, submit=args.submit)


def build_vibrations(args):
    from aiida import orm

    refs = orm.load_node(args.refs_pk) if args.refs_pk is not None else None
    a = relaxed_lattice(refs) if refs is not None else args.a
    slabs0 = {label: perovskite_slab(label, a) for label in SLAB_LABELS}
    if args.unrelaxed:  # job count only
        if args.submit:
            raise SystemExit("--unrelaxed only counts the jobs; it cannot be submitted")
        structures = {**reference_structures(), **{label: s for label, (s, _) in slabs0.items()}}
    else:
        _, references = psteros.read_vasp_results(refs, BULK_LABELS)
        _, slabs = psteros.read_vasp_results(orm.load_node(args.slabs_pk), SLAB_LABELS)
        structures = {**references, **slabs}
    overrides = {label: sno2.override(fixed_sites=frozen) for label, (_, frozen) in slabs0.items()}
    overrides.update({label: sno2.override() for label in BULK_LABELS})
    static = recipe(
        "srtio3_vibrations", {**sno2.STATIC, "EDIFF": sno2.VIBRATIONS_EDIFF}, sno2.VIBRATIONS_WALLTIME, args, overrides
    )
    return psteros.build_vibrations_workgraph(
        structures, static,
        psteros.VibrationsConfig(displacement_angstrom=sno2.DISPLACEMENT, supercells=BULK_SUPERCELLS),
        submit=args.submit,
    )


BUILDERS = {"refs": build_refs, "slabs": build_slabs, "vibrations": build_vibrations}


def count_jobs(graph) -> dict[str, int]:
    counts: dict[str, int] = {}
    labels = sorted((*BULK_LABELS, *SLAB_LABELS), key=len, reverse=True)
    for task in graph.tasks:
        if task.name.endswith("_vasp"):
            label = next((l for l in labels if task.name.startswith(l)), "?")
            counts[label] = counts.get(label, 0) + 1
    return counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=tuple(BUILDERS))
    parser.add_argument("--profile", default=sno2.PROFILE)
    parser.add_argument("--code", default=sno2.CODE)
    parser.add_argument("--potential-family", default=sno2.POTENTIAL_FAMILY)
    parser.add_argument("--kpoints-spacing", type=float, default=sno2.KPOINTS_SPACING, help="aiida-vasp units (x 2 pi)")
    parser.add_argument("--max-concurrent", type=int, default=1,
                        help="calculations of the graph running at once (1: the one-at-a-time limit of the SnO2 test)")
    parser.add_argument("--refs-pk", type=int)
    parser.add_argument("--slabs-pk", type=int)
    parser.add_argument("--a", type=float, default=3.94, help="SrTiO3 a when --refs-pk is not given (build check only)")
    parser.add_argument("--unrelaxed", action="store_true", help="vibrations: count jobs on unrelaxed structures")
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    if args.phase == "vibrations" and not args.unrelaxed and (args.refs_pk is None or args.slabs_pk is None):
        parser.error("vibrations needs --refs-pk and --slabs-pk (or --unrelaxed for a job count)")
    if args.phase == "slabs" and args.refs_pk is None:
        print("note: no --refs-pk, using --a")

    from aiida import load_profile

    load_profile(args.profile)
    graph = BUILDERS[args.phase](args)
    counts = count_jobs(graph)
    jobs = sum(counts.values())
    print(f"graph {graph.name!r}: {jobs} VASP jobs, max_concurrent_jobs={graph.max_number_jobs}")
    print("  jobs per label:", dict(sorted(counts.items())))
    expected = EXPECTED_JOBS[args.phase]
    print(f"  {'matches' if jobs == expected else 'DIFFERS FROM'} the plan ({expected} jobs)")
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
