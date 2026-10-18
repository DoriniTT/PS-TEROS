"""The AiiDA-free MLIP layer and the rutile (110) slabs it is validated on."""

from __future__ import annotations

import io

import numpy as np
import pytest

import psteros
from psteros import mlip
from psteros.structures import RUTILE_110_TERMINATIONS, rutile_110_slab

pymatgen = pytest.importorskip("pymatgen")


def ruo2_bulk():
    return psteros.structures.rutile_bulk("Ru", a=4.492, c=3.107)


def frame(label="bulk", structure=None, energy=-100.0, role="train"):
    structure = structure or ruo2_bulk()
    forces = np.arange(3 * len(structure), dtype=float).reshape(-1, 3) / 100.0
    return mlip.LabelledStructure(label, structure, energy, forces, role=role)


def test_labelled_structure_checks_shapes_and_roles() -> None:
    example = frame()
    assert example.composition == {"Ru": 2, "O": 4}
    assert len(example.forces_ev_per_angstrom) == 6
    with pytest.raises(ValueError, match="one row of 3 per site"):
        mlip.LabelledStructure("bulk", ruo2_bulk(), -1.0, [[0.0, 0.0, 0.0]])
    with pytest.raises(ValueError, match="role"):
        frame(role="holdout")
    with pytest.raises(ValueError, match="finite"):
        frame(energy=float("nan"))
    with pytest.raises(ValueError, match="non-empty word"):
        frame(label="two words")


def test_extxyz_round_trips_through_ase() -> None:
    ase_io = pytest.importorskip("ase.io")
    stream = io.StringIO()
    frames = [frame(), frame("ru", psteros.structures.Structure.from_spacegroup(
        "P6_3/mmc", psteros.structures.Lattice.hexagonal(2.71, 4.28), ["Ru"], [[1 / 3, 2 / 3, 0.25]]), -9.5, "valid")]
    mlip.write_extxyz(frames, stream)
    stream.seek(0)
    atoms = ase_io.read(stream, ":", format="extxyz")
    assert [a.info["label"] for a in atoms] == ["bulk", "ru"] and atoms[1].info["role"] == "valid"
    assert atoms[0].info["REF_energy"] == pytest.approx(-100.0)
    assert atoms[0].arrays["REF_forces"] == pytest.approx(np.asarray(frames[0].forces_ev_per_angstrom))
    assert atoms[0].get_positions() == pytest.approx(frames[0].structure.cart_coords, abs=1e-9)
    assert atoms[0].cell.array == pytest.approx(frames[0].structure.lattice.matrix)
    assert all(atoms[0].pbc)
    with pytest.raises(ValueError, match="duplicate frame labels"):
        mlip.write_extxyz([frame(), frame()], io.StringIO())


def test_element_shifts_need_separating_compositions() -> None:
    shifts = {"Ru": -2.0, "O": 0.5}
    compositions = [{"Ru": 2, "O": 4}, {"Ru": 2}, {"O": 2}, {"Ru": 6, "O": 14}]
    residuals = [sum(n * shifts[e] for e, n in c.items()) for c in compositions]
    assert mlip.fit_element_shifts(compositions, residuals) == {
        "O": pytest.approx(0.5), "Ru": pytest.approx(-2.0)}
    with pytest.raises(ValueError, match="do not separate"):
        mlip.fit_element_shifts([{"Ru": 2, "O": 4}, {"Ru": 4, "O": 8}], [1.0, 2.0])


def test_fine_tune_recipe_and_training_options() -> None:
    config = mlip.MaceFineTuneConfig("/models/mace-mh-1.model", replay_file="/data/replay.xyz")
    options = mlip.mace_training_options(config, name="ruo2_m2", seed=20260920, atomic_numbers=[44, 8, 8])
    assert options["atomic_numbers"] == "[8,44]"
    assert options["multiheads_finetuning"] is True and options["force_mh_ft_lr"] is True
    assert options["pt_train_file"] == "replay_labels/train.extxyz"
    assert options["max_num_epochs"] == 100 and options["patience"] == 101
    assert (options["lr"], options["ema_decay"], options["huber_delta"]) == (1.0e-3, 0.999, 0.01)
    assert options["scheduler"] == "ExponentialLR" and options["lr_scheduler_gamma"] == 1.0
    ordinary = mlip.mace_training_options(
        mlip.MaceFineTuneConfig("/models/m.model", replay=False), name="m1", seed=20260921, atomic_numbers=[8])
    assert ordinary["multiheads_finetuning"] is False and "pt_train_file" not in ordinary
    with pytest.raises(ValueError, match="replay_file"):
        mlip.MaceFineTuneConfig("/models/m.model")
    with pytest.raises(ValueError, match="seeds"):
        mlip.MaceFineTuneConfig("/models/m.model", replay=False, seeds=(1, 1))
    with pytest.raises(ValueError, match="not one of the recipe's seeds"):
        mlip.mace_training_options(config, name="x", seed=7, atomic_numbers=[8])
    with pytest.raises(ValueError, match="optimizer"):
        mlip.MlipRelaxConfig(optimizer="CG")


def test_learning_rate_decay_is_opt_in() -> None:
    config = mlip.MaceFineTuneConfig("/models/m.model", replay_file="/data/replay.xyz")
    default = mlip.mace_training_options(config, name="m2", seed=20260920, atomic_numbers=[8, 44])
    assert default["lr_scheduler_gamma"] == 1.0  # constant learning rate, as before
    decaying = mlip.MaceFineTuneConfig("/models/m.model", replay_file="/data/replay.xyz", epochs=500,
                                       lr_scheduler_gamma=0.990832)
    options = mlip.mace_training_options(decaying, name="m2l", seed=20260920, atomic_numbers=[8, 44])
    assert options["lr_scheduler_gamma"] == 0.990832 and options["max_num_epochs"] == 500
    assert {k: v for k, v in options.items() if k not in ("lr_scheduler_gamma", "max_num_epochs", "patience", "name")} \
        == {k: v for k, v in default.items() if k not in ("lr_scheduler_gamma", "max_num_epochs", "patience", "name")}
    with pytest.raises(ValueError, match="lr_scheduler_gamma"):
        mlip.MaceFineTuneConfig("/models/m.model", replay=False, lr_scheduler_gamma=1.5)


def test_relaxation_keeps_fixed_sites_and_converges() -> None:
    pytest.importorskip("ase")
    from ase.calculators.emt import EMT
    from pymatgen.core import Lattice, Structure

    lattice = Lattice.cubic(4.05)
    cell = Structure.from_spacegroup("Fm-3m", lattice, ["Al"], [[0, 0, 0]])
    cell.make_supercell([2, 2, 2])
    start = cell.copy()
    start.translate_sites([0], [0.02, 0.0, 0.0])
    start.translate_sites([1], [0.0, 0.02, 0.0])
    result = mlip.relax_with_calculator(start, EMT(), config=mlip.MlipRelaxConfig(fmax_ev_per_angstrom=0.005),
                                        fixed_sites=(1,))
    assert result.converged and result.max_force_ev_per_angstrom < 0.005 and result.steps > 0
    assert result.structure[1].coords == pytest.approx(start[1].coords)
    assert not np.allclose(result.structure[0].coords, start[0].coords)
    with pytest.raises(ValueError, match="fixed_sites"):
        mlip.relax_with_calculator(start, EMT(), fixed_sites=(99,))


def test_metallic_check_counts_partly_filled_states() -> None:
    assert mlip.partially_occupied_count([[2.0, 2.0, 1.2, 0.0]]) == 1
    assert mlip.partially_occupied_count([[1.0, 1.0, 1e-6, 0.0], [1.0, 0.999999, 0.0, 0.0]]) == 0
    mlip.require_metallic("bulk", partially_occupied=3)
    with pytest.raises(ValueError, match="bulk: the bulk has no partly filled states .* metallic systems only"):
        mlip.require_metallic("bulk", partially_occupied=0)


def test_one_energy_source_per_diagram() -> None:
    assert mlip.check_one_energy_source({"bulk": "mlip", "slab": "mlip", "o2": "mlip"}) == "mlip"
    with pytest.raises(ValueError, match=r"dft: \['bulk', 'o2'\]; mlip: \['slab'\]"):
        mlip.check_one_energy_source({"bulk": "dft", "slab": "mlip", "o2": "dft"})
    with pytest.raises(ValueError, match="unknown energy source"):
        mlip.check_one_energy_source({"slab": "vasp"})


@pytest.mark.parametrize("termination", RUTILE_110_TERMINATIONS)
def test_rutile_110_terminations(termination) -> None:
    slab, identity = rutile_110_slab(metal="Ru", termination=termination, a=4.492, c=3.107)
    cells = 2 if termination.endswith("_half") else 1
    excess_per_cell_and_surface = {"o_bridge": 0, "o_cus": 1, "metal": -1, "o_cus_half": 0.5,
                                   "o_bridge_half": -0.5, "o_cus_no_bridge": 0}[termination]
    composition = slab.composition
    assert composition["Ru"] == 6 * cells
    assert (composition["O"] - 2 * composition["Ru"]) / (2 * cells) == excess_per_cell_and_surface
    assert slab.is_symmetric() and identity.termination == termination
    assert slab.lattice.abc[0] == pytest.approx(3.107 * cells) or slab.lattice.abc[1] == pytest.approx(3.107 * cells)
    distances = slab.distance_matrix + 10 * np.eye(len(slab))
    assert distances.min() > 1.65


def test_rutile_110_slab_rejects_bad_requests_and_keeps_sno2() -> None:
    with pytest.raises(ValueError, match="termination must be one of"):
        rutile_110_slab(metal="Ru", termination="ru", a=4.492, c=3.107)
    with pytest.raises(ValueError, match="odd integer"):
        rutile_110_slab(metal="Ru", triple_layers=4, a=4.492, c=3.107)
    sno2, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    same, _ = rutile_110_slab(metal="Sn", a=4.737, c=3.186, triple_layers=3, vacuum_angstrom=20.0)
    assert sno2.composition == same.composition


def _replay_file(tmp_path):
    lines = []
    for count, symbols in ((2, "O O"), (3, "Ru O O"), (2, "O O"), (2, "Ag O"), (1, "Ru"), (3, "Ru O O")):
        lines.append(str(count))
        lines.append('Lattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3 pbc="T T T"')
        for i, symbol in enumerate(symbols.split()):
            lines.append(f"{symbol} {i * 1.2:.3f} 0.0 0.0")
    path = tmp_path / "replay.xyz"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_replay_selection_keeps_distinct_allowed_structures(tmp_path) -> None:
    pytest.importorskip("ase")
    train, valid, record = mlip.select_replay_structures(
        _replay_file(tmp_path), elements=["Ru", "O"], min_atoms=2, validation_every=2)
    # Ag-O is not allowed, the lone Ru is too small, the repeated O2 and RuO2 are duplicates
    assert (record["source_frames"], record["kept"]) == (6, 2)
    assert [a.get_chemical_formula() for a in train + valid] == ["O2", "O2Ru"]
    assert valid[0].info["replay_source_frame"] == 1 and len(record["source_sha256"]) == 64
    with pytest.raises(ValueError, match="no replay structure"):
        mlip.select_replay_structures(_replay_file(tmp_path), elements=["Ti"])


def test_job_folders_hold_inputs_and_scripts(tmp_path) -> None:
    import json

    pytest.importorskip("ase")
    frames = [frame(), frame("o2", psteros.triplet_o2_cell(), -9.0), frame("test", role="test")]
    config = mlip.MaceFineTuneConfig("/models/mh1.model", replay_file="/data/replay.xyz", epochs=3)
    train, valid, _ = mlip.select_replay_structures(_replay_file(tmp_path), elements=["Ru", "O"], validation_every=2)
    folder = mlip.prepare_finetune_job(tmp_path / "ft", frames, config, name="ruo2", replay_train=train,
                                       replay_valid=valid, seeds=(20260921,))
    job = json.loads((folder / "job.json").read_text())
    assert job["replay"] is True and [run["name"] for run in job["runs"]] == ["ruo2_seed20260921"]
    assert job["runs"][0]["options"]["atomic_numbers"] == "[8,44]"
    assert job["runs"][0]["options"]["max_num_epochs"] == 3
    assert (folder / "train.extxyz").read_text().count("label=") == 2  # the test frame stays out
    assert (folder / "finetune.py").read_text().startswith('"""psteros MACE fine-tune job')
    assert (folder / "replay_train_geometries.extxyz").is_file()
    with pytest.raises(FileExistsError):
        mlip.prepare_finetune_job(tmp_path / "ft", frames, config, name="ruo2", replay_train=train, replay_valid=valid)
    with pytest.raises(ValueError, match="needs replay_train"):
        mlip.prepare_finetune_job(tmp_path / "other", frames, config, name="ruo2")

    slab, _ = rutile_110_slab(metal="Ru", termination="o_cus", a=4.492, c=3.107)
    relax = mlip.prepare_relax_job(
        tmp_path / "relax", {"slab": slab}, {"m2": ("/models/m2.model", "Default"), "m0": ("/models/mh1.model", "omat_pbe")},
        fixed_sites={"slab": psteros.central_sites(slab)})
    job = json.loads((relax / "job.json").read_text())
    assert job["models"]["m0"] == {"path": "/models/mh1.model", "head": "omat_pbe"}
    assert job["structures"][0]["fixed"] == sorted(psteros.central_sites(slab))
    assert job["fmax_ev_per_angstrom"] == 0.01 and (relax / "relax.py").is_file()
    with pytest.raises(ValueError, match="not in structures"):
        mlip.prepare_relax_job(tmp_path / "relax2", {"slab": slab}, {"m": ("/m", "Default")}, fixed_sites={"bulk": (0,)})
