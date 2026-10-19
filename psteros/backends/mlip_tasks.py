"""Calcfunctions of :func:`psteros.build_mlip_surface_workgraph`.

This module imports AiiDA at import time, so psteros loads it only while
building a graph; the daemon imports it when it runs one, which is why
psteros (with ase) must be installed in the environment of the AiiDA daemon.
"""

from __future__ import annotations

import io
import json
import tempfile
from pathlib import Path

import numpy as np
from aiida import orm
from aiida_workgraph import task

from psteros import mlip

# Files of a retrieved DFT folder that the MLIP route reads (forces, occupations).
NEEDED_FILES = {
    "cp2k": ("psteros_forces.dat", "psteros.dos"),
    "qe": ("aiida.out", "data-file-schema.xml"),
    "vasp": ("vasprun.xml",),
}


def retrieved_files(folder: orm.FolderData, backend: str) -> dict[str, str]:
    """Text of the files of ``NEEDED_FILES[backend]`` that ``folder`` holds."""

    names = set(folder.base.repository.list_object_names())
    return {name: folder.base.repository.get_object_content(name) for name in NEEDED_FILES[backend] if name in names}


def _read_extxyz(content: str) -> list:
    from ase.io import read

    return read(io.StringIO(content), index=":", format="extxyz")


@task.calcfunction
def mlip_training_folder(settings: orm.Dict, **inputs) -> orm.FolderData:
    """The folder of one fine-tune job (:func:`psteros.mlip.prepare_finetune_job`) from finished DFT statics.

    ``settings`` holds ``backend``, ``name``, ``device``, the
    :class:`~psteros.mlip.MaceFineTuneConfig` fields (``fine_tune``), the
    labels of the training frames (``frames``) and the label of the bulk that
    must be metallic (``metallic_bulk``). ``inputs`` holds, for every frame,
    ``structure_<label>``, ``output_<label>`` (output parameters, VASP
    ``misc``) and ``retrieved_<label>``, and for a replay fine-tune
    ``replay_train`` and ``replay_valid`` (extended XYZ of the replay
    geometries chosen when the graph was built).
    """

    options = settings.get_dict()
    backend = options["backend"]
    bulk = options["metallic_bulk"]
    files = retrieved_files(inputs[f"retrieved_{bulk}"], backend)
    mlip.require_metallic(bulk, partially_occupied=mlip.dft_partially_occupied(backend, label=bulk, files=files))
    frames = [
        mlip.dft_label(
            backend,
            label=label,
            structure=inputs[f"structure_{label}"].get_pymatgen_structure(),
            parameters=inputs[f"output_{label}"].get_dict(),
            files=retrieved_files(inputs[f"retrieved_{label}"], backend),
        )
        for label in options["frames"]
    ]
    config = mlip.MaceFineTuneConfig(**options["fine_tune"])
    replay = {}
    if config.replay:
        replay = {
            "replay_train": _read_extxyz(inputs["replay_train"].get_content()),
            "replay_valid": _read_extxyz(inputs["replay_valid"].get_content()),
        }
    with tempfile.TemporaryDirectory() as scratch:
        folder = mlip.prepare_finetune_job(
            Path(scratch) / "job", frames, config, name=options["name"], device=options["device"], **replay
        )
        return orm.FolderData(tree=str(folder))


@task.calcfunction(outputs=["structure", "energy"])
def mlip_relaxed(
    structure: orm.StructureData, results: orm.SinglefileData, relaxed: orm.FolderData, settings: orm.Dict
) -> dict:
    """Structure and energy (eV) of ``settings["label"]`` relaxed with the model ``settings["model"]``.

    ``structure`` is the starting structure, whose kinds and cell the relaxed
    structure keeps; ``results`` and ``relaxed`` are the ``results.json`` and
    ``relaxed/`` folder of the relaxation job. Raises unless that relaxation
    converged.
    """

    options = settings.get_dict()
    label, model = options["label"], options["model"]
    data = json.loads(results.get_content())
    row = data.get("models", {}).get(model, {}).get("structures", {}).get(label)
    if row is None:
        raise ValueError(f"{label}: results.json has no relaxation with the model {model} (status {data.get('status')})")
    if not row["converged"]:
        raise ValueError(
            f"{label}: the relaxation with {model} did not converge (largest force "
            f"{row['max_force_ev_per_angstrom']:.4f} eV/A after {row['steps']} steps); raise MlipRelaxConfig.max_steps"
        )
    (atoms,) = _read_extxyz(relaxed.base.repository.get_object_content(f"{model}/{label}.extxyz"))
    symbols = [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    if atoms.get_chemical_symbols() != symbols:
        raise ValueError(f"{label}: the relaxed atoms are not those of the starting structure")
    if not np.allclose(atoms.cell.array, np.array(structure.cell), atol=1e-6):
        raise ValueError(f"{label}: the relaxation changed the cell")
    copy = orm.StructureData(cell=structure.cell, pbc=structure.pbc)
    for kind in structure.kinds:
        copy.append_kind(kind)
    for site, position in zip(structure.sites, atoms.positions):
        copy.append_site(orm.Site(kind_name=site.kind_name, position=[float(value) for value in position]))
    return {"structure": copy, "energy": orm.Float(float(row["energy_ev"]))}


@task.calcfunction
def mlip_dft_energy(output: orm.Dict, settings: orm.Dict) -> orm.Float:
    """Energy (eV) of a DFT static on a relaxed geometry (:func:`psteros.mlip.dft_diagram_energy`)."""

    options = settings.get_dict()
    return orm.Float(mlip.dft_diagram_energy(options["backend"], label=options["label"], parameters=output.get_dict()))
