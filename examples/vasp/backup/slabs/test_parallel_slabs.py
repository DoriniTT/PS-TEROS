#!/usr/bin/env python
"""Test script for parallel slab relaxation using the Map zone pattern."""

from __future__ import annotations

from pathlib import Path

from aiida import load_profile, orm
from ase.io import read

from teros.workgraph_parallel_slabs import build_parallel_slab_workgraph


def load_bulk_structure(structure_path: Path) -> orm.StructureData:
    """Create a ``StructureData`` directly from a structure file."""
    atoms = read(structure_path)
    return orm.StructureData(ase=atoms)


def main() -> None:
    """Configure parameters and launch the parallel slab relaxation workflow."""
    load_profile()

    structures_dir = Path("/home/thiagotd/git/PS-TEROS/teros/structures")
    bulk_path = structures_dir / "ag3po4.cif"

    if not bulk_path.exists():
        raise FileNotFoundError(f"Bulk structure file not found: {bulk_path}")

    bulk_structure = load_bulk_structure(bulk_path)

    code_label = "VASP-VTST-6.4.3@bohr"
    potential_family = "PBE"
    potential_mapping = {"Ag": "Ag", "P": "P", "O": "O"}

    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 100,
        "EDIFFG": -0.02,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }

    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    workflow = build_parallel_slab_workgraph(
        bulk_structure=bulk_structure,
        miller_indices=(1, 0, 0),
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,
        code_label=code_label,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
        parameters=slab_parameters,
        options=slab_options,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    # Submit workflow
    workflow.submit(wait=False)

    print("\nWorkflow submitted:")
    print(f"  PK: {workflow.pk}")
    print(f"\nMonitor with: verdi process show {workflow.pk}")
    print(f"Check report: verdi process report {workflow.pk}")


if __name__ == "__main__":
    main()
