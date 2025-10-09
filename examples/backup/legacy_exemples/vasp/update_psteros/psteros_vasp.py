#!/usr/bin/env python3
"""Minimal runner for the PS-TEROS workflow with VASP.

This script keeps the configuration compact while remaining explicit about the
inputs that matter for the Ag3PO4 example distributed with the repository.
Adjust the constants in :func:`main` as needed for different materials or
computing environments.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

from aiida import load_profile
from aiida.orm import (
    Bool,
    Dict,
    Float,
    Int,
    KpointsData,
    List,
    StructureData,
    Str,
    load_code,
)
from aiida.plugins import WorkflowFactory
from ase.io import read

from teros.core.workgraph import create_teros_workgraph

# Ensure the active profile is available before touching AiiDA objects.
load_profile()

VASP_WORKCHAIN = WorkflowFactory("vasp.vasp")


@dataclass(frozen=True)
class RunnerConfig:
    """Container with the few knobs we typically change."""

    code_label: str = "VASP-VTST-6.4.3@bohr"
    queue: str | None = "par40"
    num_machines: int = 1
    num_cores_per_machine: int = 40
    wallclock_seconds: int = 72 * 3600
    max_restarts: int = 3
    potential_family: str = "PBE"
    potential_mapping: Mapping[str, str] = field(
        default_factory=lambda: {"Ag": "Ag", "P": "P", "O": "O"}
    )

    def scheduler_options(self) -> dict[str, object]:
        resources: dict[str, object] = {"num_machines": self.num_machines}
        if self.num_cores_per_machine:
            resources["num_cores_per_machine"] = self.num_cores_per_machine
        return {
            "resources": resources,
            "max_wallclock_seconds": self.wallclock_seconds,
            "withmpi": True,
            **({"queue_name": self.queue} if self.queue else {}),
        }


@dataclass(frozen=True)
class StructurePaths:
    """Locations of the example structures shipped with the project."""

    root: Path = Path(__file__).resolve().parent
    base: Path = root / "structures"
    bulk: Path = base / "bulk/Ag6O8P2_optimized.cif"
    ag: Path = base / "pure_elements/Ag.cif"
    p: Path = base / "pure_elements/P.cif"
    o2: Path = base / "pure_elements/O2.cif"


def read_structure(path: Path) -> StructureData:
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")
    return StructureData(ase=read(path))


def mesh_kpoints(mesh: tuple[int, int, int]) -> KpointsData:
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(mesh)
    return kpoints


def make_incar(base: Mapping[str, object], **updates: object) -> dict[str, object]:
    incar = dict(base)
    incar.update(updates)
    return incar


def make_builder(
    *,
    code,
    config: RunnerConfig,
    structure: StructureData | None,
    incar: Mapping[str, object],
    k_mesh: tuple[int, int, int],
    label: str,
    description: str,
) -> object:
    builder = VASP_WORKCHAIN.get_builder()
    builder.code = code
    if structure is not None:
        builder.structure = structure
    builder.potential_family = Str(config.potential_family)
    builder.potential_mapping = Dict(dict(config.potential_mapping))
    builder.parameters = Dict(dict({"incar": dict(incar)}))
    builder.kpoints = mesh_kpoints(k_mesh)
    builder.options = Dict(dict(config.scheduler_options()))
    builder.max_iterations = Int(config.max_restarts)
    builder.clean_workdir = Bool(False)
    builder.metadata.label = label
    builder.metadata.description = description
    return builder


def main() -> None:
    config = RunnerConfig()
    paths = StructurePaths()
    code = load_code(config.code_label)

    common_incar = {
        "ENCUT": 500,
        "EDIFF": 1.0e-5,
        "ISMEAR": 0,
        "SIGMA": 0.01,
        "ISPIN": 2,
        "PREC": "Accurate",
        "ALGO": "Fast",
        "NCORE": 2,
        "NELM": 80,
        "NELMIN": 5,
        "LREAL": "Auto",
    }

    bulk_structure = read_structure(paths.bulk)
    bulk_builder = make_builder(
        code=code,
        config=config,
        structure=bulk_structure,
        incar=make_incar(
            common_incar,
            ISIF=3,
            IBRION=2,
            NSW=200,
            EDIFFG=-0.01,
        ),
        k_mesh=(3, 3, 3),
        label="Ag3PO4 bulk relaxation",
        description="VASP bulk relaxation for PS-TEROS",
    )

    slab_builder = make_builder(
        code=code,
        config=config,
        structure=None,
        incar=make_incar(
            common_incar,
            ISIF=2,
            IBRION=2,
            NSW=300,
            EDIFFG=-0.05,
            LDIPOL=True,
            IDIPOL=3,
        ),
        k_mesh=(2, 2, 1),
        label="Ag3PO4 slab relaxation",
        description="VASP slab relaxations driven by PS-TEROS",
    )

    reference_builders = {
        "Ag": make_builder(
            code=code,
            config=config,
            structure=read_structure(paths.ag),
            incar=make_incar(
                common_incar,
                ISIF=3,
                IBRION=2,
                NSW=120,
                EDIFFG=-0.01,
                ISMEAR=1,
                SIGMA=0.2,
            ),
            k_mesh=(7, 7, 7),
            label="Ag reference",
            description="Silver bulk reference for PS-TEROS",
        ),
        "P": make_builder(
            code=code,
            config=config,
            structure=read_structure(paths.p),
            incar=make_incar(
                common_incar,
                ISIF=3,
                IBRION=2,
                NSW=150,
                EDIFFG=-0.01,
            ),
            k_mesh=(4, 4, 4),
            label="P reference",
            description="Phosphorus reference for PS-TEROS",
        ),
        "O2": make_builder(
            code=code,
            config=config,
            structure=read_structure(paths.o2),
            incar=make_incar(
                common_incar,
                ISIF=2,
                IBRION=2,
                NSW=80,
                EDIFFG=-0.02,
            ),
            k_mesh=(1, 1, 1),
            label="O2 reference",
            description="Molecular oxygen reference for PS-TEROS",
        ),
    }

    wg = create_teros_workgraph(
        dft_workchain=VASP_WORKCHAIN,
        builder_bulk=bulk_builder,
        builder_slab=slab_builder,
        reference_builders=reference_builders,
        # workgraph_name="psteros_vasp_ag3po4",
        sampling=Int(100),
        miller_indices=List(list=[[1, 0, 0]]),
        min_slab_thickness=Float(15.0),
        min_vacuum_thickness=Float(15.0),
        lll_reduce=Bool(True),
        center_slab=Bool(True),
        symmetrize=Bool(True),
        primitive=Bool(True),
        in_unit_planes=Bool(False),
    )

    print("Submitting PS-TEROS VASP workflow...")
    wg.submit(wait=False)
    wg.to_html()
    print(f"Workflow submitted. WorkGraph PK: {wg.pk}")


if __name__ == "__main__":
    main()
