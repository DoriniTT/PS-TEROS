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
    :class:`~psteros.mlip.MaceFineTuneConfig` fields (``fine_tune``; with
    ``family="uma"`` those of :class:`~psteros.uma.UmaFineTuneConfig`), the
    labels of the training frames (``frames``) and the label of the bulk that
    must be metallic (``metallic_bulk``; not checked when ``metal_check`` is
    false). ``inputs`` holds, for every frame,
    ``structure_<label>``, ``output_<label>`` (output parameters, VASP
    ``misc``) and ``retrieved_<label>``, and for a replay fine-tune
    ``replay_train`` and ``replay_valid`` (extended XYZ of the replay
    geometries chosen when the graph was built).
    """

    options = settings.get_dict()
    backend = options["backend"]
    bulk = options["metallic_bulk"]
    if options.get("metal_check", True):
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
    if options.get("family", "mace") == "uma":
        from psteros.uma import UmaFineTuneConfig, prepare_uma_finetune_job

        config, prepare = UmaFineTuneConfig(**options["fine_tune"]), prepare_uma_finetune_job
    elif options.get("family") == "nequip":
        from psteros.nequip import NequipFineTuneConfig, prepare_nequip_finetune_job

        config, prepare = NequipFineTuneConfig(**options["fine_tune"]), prepare_nequip_finetune_job
    else:
        config, prepare = mlip.MaceFineTuneConfig(**options["fine_tune"]), mlip.prepare_finetune_job
    replay = {}
    if config.replay:
        replay = {
            "replay_train": _read_extxyz(inputs["replay_train"].get_content()),
            "replay_valid": _read_extxyz(inputs["replay_valid"].get_content()),
        }
    with tempfile.TemporaryDirectory() as scratch:
        folder = prepare(Path(scratch) / "job", frames, config, name=options["name"], device=options["device"], **replay)
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


@task.calcfunction
def mlip_vibrations_folder(settings: orm.Dict, **structures) -> orm.FolderData:
    """The folder of the vibrations job (:func:`psteros.mlip.prepare_vibrations_job`) on the relaxed structures.

    ``settings`` holds ``model`` (``[path, head]``), ``displaced_sites``
    (label -> sites), ``displacement_angstrom``, ``default_dtype`` and
    ``device``; ``structures`` holds ``structure_<label>``, the relaxed
    structure of each label (a bulk as its supercell).
    """

    options = settings.get_dict()
    plain = {label: structures[f"structure_{label}"].get_pymatgen_structure() for label in options["displaced_sites"]}
    with tempfile.TemporaryDirectory() as scratch:
        folder = mlip.prepare_vibrations_job(
            Path(scratch) / "job", plain, tuple(options["model"]), options["displaced_sites"],
            displacement_angstrom=options["displacement_angstrom"], default_dtype=options["default_dtype"],
            device=options["device"], family=options.get("family", "mace"),
        )
        return orm.FolderData(tree=str(folder))


@task.calcfunction
def mlip_harmonic_modes(structure: orm.StructureData, settings: orm.Dict, forces: orm.FolderData) -> orm.Dict:
    """Harmonic modes of ``settings["label"]`` from the forces of the MACE vibrations job.

    ``forces`` is the job's ``vibrations/`` folder; ``structure`` the
    displaced structure (a bulk as its supercell); ``settings`` the options of
    :func:`psteros.backends.vibration_tasks.harmonic_modes`. The output reads
    like a ``<label>_vibrations`` of :func:`psteros.build_vibrations_workgraph`.
    """

    from psteros.backends.vibration_tasks import modes_from_forces

    options = settings.get_dict()
    label = options["label"]
    names = set(forces.base.repository.list_object_names())
    if f"{label}.npz" not in names:
        raise ValueError(f"{label}: the vibrations job wrote no {label}.npz")
    with np.load(io.BytesIO(forces.base.repository.get_object_content(f"{label}.npz", mode="rb"))) as data:
        plus, minus = data["forces_plus"], data["forces_minus"]
    expected = (3 * len(options["displaced_sites"]), len(structure.sites), 3)
    if plus.shape != expected or minus.shape != expected:
        raise ValueError(f"{label}: the vibrations job wrote forces of shape {plus.shape}, expected {expected}")
    return orm.Dict({**modes_from_forces(structure, options, plus, minus), "model": options["model"]})


@task.calcfunction
def mlip_displacements_folder(settings: orm.Dict, **structures) -> orm.FolderData:
    """The folder of the forces job (:func:`psteros.mlip.prepare_forces_job`) of phonopy vibrations.

    ``settings`` holds ``model`` (``[path, head]``), ``device``,
    ``default_dtype``, ``displacement_angstrom``, ``plus_minus``, ``symprec``,
    ``phonopy`` (label -> phonopy supercell) and ``molecules`` (label -> sites
    displaced by ``+d`` and ``-d`` along x, y and z, the Gamma-point modes);
    ``structures`` holds ``structure_<label>``, the relaxed structure of each
    label. The phonopy displacements of every supercell are generated here,
    from the relaxed structure, and written to ``phonopy.json`` next to the job.
    """

    from psteros.phonons import phonopy_displacements

    options = settings.get_dict()
    step = float(options["displacement_angstrom"])
    plain, displacements, plans = {}, {}, {}
    for label, supercell in options["phonopy"].items():
        structure = structures[f"structure_{label}"].get_pymatgen_structure()
        try:
            plan = phonopy_displacements(structure, supercell, displacement_angstrom=step,
                                         plus_minus=bool(options["plus_minus"]), symprec=float(options["symprec"]))
        except ValueError as error:
            raise ValueError(f"{label}: {error}") from error
        plans[label] = plan.to_dict()
        plain[label] = plan.structures(structure)[0]
        displacements[label] = list(plan.displacements)
    for label, sites in options["molecules"].items():
        plain[label] = structures[f"structure_{label}"].get_pymatgen_structure()
        displacements[label] = [(site, [sign * step if a == axis else 0.0 for a in range(3)])
                                for site in sites for axis in range(3) for sign in (1.0, -1.0)]
    with tempfile.TemporaryDirectory() as scratch:
        folder = mlip.prepare_forces_job(
            Path(scratch) / "job", plain, tuple(options["model"]), displacements,
            default_dtype=options["default_dtype"], device=options["device"], family=options.get("family", "mace"),
        )
        (folder / "phonopy.json").write_text(json.dumps(plans) + "\n")
        return orm.FolderData(tree=str(folder))


@task.calcfunction
def mlip_phonon_modes(structure: orm.StructureData, settings: orm.Dict, job: orm.FolderData,
                      forces: orm.FolderData) -> orm.Dict:
    """Modes of ``settings["label"]`` from the forces job of :func:`mlip_displacements_folder`.

    ``job`` is the job's folder (its ``phonopy.json``), ``forces`` the job's
    ``forces/`` folder and ``structure`` the relaxed structure. With
    ``settings["method"] == "phonopy"`` the settings are those of
    :func:`psteros.vibrations_workflow.phonopy_settings` (without
    ``displacements``, read from ``phonopy.json``); with ``"gamma"`` (a
    molecule) those of :func:`psteros.backends.vibration_tasks.harmonic_modes`.
    The output reads like a ``<label>_vibrations`` of :func:`psteros.build_vibrations_workgraph`.
    """

    from psteros.backends.vibration_tasks import modes_from_forces, phonopy_modes_from_forces

    options = settings.get_dict()
    label = options["label"]
    if f"{label}.npz" not in set(forces.base.repository.list_object_names()):
        raise ValueError(f"{label}: the forces job wrote no {label}.npz")
    with np.load(io.BytesIO(forces.base.repository.get_object_content(f"{label}.npz", mode="rb"))) as data:
        values, reference = data["forces"], data["reference_forces"]
    if options["method"] == "gamma":
        sites = list(options["displaced_sites"])
        if values.shape != (6 * len(sites), len(structure.sites), 3):
            raise ValueError(f"{label}: the forces job wrote forces of shape {values.shape} for {len(sites)} sites")
        return orm.Dict({**modes_from_forces(structure, options, values[0::2], values[1::2]), "model": options["model"]})
    plans = json.loads(job.base.repository.get_object_content("phonopy.json"))
    try:
        result = phonopy_modes_from_forces(structure, {**options, "displacements": plans[label]}, values, reference)
    except ValueError as error:
        raise ValueError(f"{label}: {error}") from error
    return orm.Dict({**result, "model": options["model"]})
