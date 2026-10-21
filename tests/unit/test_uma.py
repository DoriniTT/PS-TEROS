"""UMA (fairchem) on the MLIP route: the recipe, its job folder and the graph built with it (no fairchem needed)."""

from __future__ import annotations

import json
import warnings

import pytest

import psteros
from psteros import mlip


def frames():
    from pymatgen.core import Lattice, Structure

    first = Structure(Lattice.cubic(6.0), ["Ru", "O"], [[0, 0, 0], [0.3, 0, 0]])
    second = Structure(Lattice.cubic(6.0), ["Ru", "O", "O"], [[0, 0, 0], [0.3, 0, 0], [0.7, 0, 0]])
    return [psteros.LabelledStructure("pair", first, -10.0, [[0.1, 0, 0], [-0.1, 0, 0]]),
            psteros.LabelledStructure("triple", second, -15.0, [[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]])]


def test_uma_catalogue_and_recipe() -> None:
    models = psteros.UMA_FOUNDATION_MODELS
    assert set(models) == {"uma-s-1p2p1", "uma-m-1p1"}
    small = psteros.uma_foundation_model("uma-s-1p2p1")
    assert small.url == "https://huggingface.co/facebook/UMA/resolve/main/checkpoints/uma-s-1p2p1.pt"
    assert small.default_task == "omat" and small.task().functional == "PBE" and small.task().hubbard == "mp"
    with pytest.raises(ValueError, match=r"uma-s-1p2p1 has no task 'omc'; its tasks are \['omat', 'oc20', 'oc22', 'omol'\]"):
        small.task("omc")
    with pytest.raises(ValueError, match="unknown UMA foundation model 'uma-s-9'"):
        psteros.uma_foundation_model("uma-s-9")

    config = psteros.UmaFineTuneConfig.for_foundation("uma-s-1p2p1", "/work/me/uma-s-1p2p1.pt", replay_file="r.xyz")
    assert config.task == "omat" and config.foundation_name == "uma-s-1p2p1" and config.seeds == (20260920,)
    with pytest.raises(ValueError, match="replay=True needs replay_file"):
        psteros.UmaFineTuneConfig("/m.pt")
    with pytest.raises(ValueError, match="epochs must be a positive integer"):
        psteros.UmaFineTuneConfig("/m.pt", replay=False, epochs=0)
    with pytest.raises(ValueError, match="ema_decay"):
        psteros.UmaFineTuneConfig("/m.pt", replay=False, ema_decay=1.0)
    with pytest.raises(ValueError, match="has no task 'omc'"):
        psteros.UmaFineTuneConfig("/m.pt", replay=False, task="omc", foundation_name="uma-s-1p2p1")


def test_uma_level_warnings() -> None:
    vasp = psteros.VaspCalculationConfig("vasp@cluster", incar={"ENCUT": 520})
    assert psteros.uma_level_warnings("uma-s-1p2p1", None, vasp, ["Ru", "O"]) == []
    (rpbe,) = psteros.uma_level_warnings("uma-s-1p2p1", "oc20", vasp, ["Pt"])
    assert rpbe.startswith("the DFT recipe uses PBE but uma-s-1p2p1 (task 'oc20') was trained on RPBE data")
    assert "(trained on PBE: uma-s-1p2p1 (task 'omat'), uma-s-1p2p1 (task 'oc22'), uma-m-1p1 (task 'omat'))" in rpbe
    assert psteros.uma_foundation_model("uma-m-1p1").task_names == ("omat", "oc20", "omol")
    (hubbard,) = psteros.uma_level_warnings("uma-s-1p2p1", "omat", vasp, ["Fe", "O"])
    assert hubbard == ("uma-s-1p2p1 (task 'omat') was trained with the Materials Project +U on Fe in oxides and "
                       "fluorides, but the DFT recipe has no +U")


def test_uma_job_folder(tmp_path) -> None:
    from ase import Atoms

    config = psteros.UmaFineTuneConfig("/work/me/uma-s-1p2p1.pt", foundation_model_sha256="a" * 64,
                                       replay_file="r.xyz", epochs=7, seeds=(1, 2))
    replay = [Atoms("RuO", positions=[[0, 0, 0], [1.9, 0, 0]], cell=[8, 8, 8], pbc=True)]
    folder = psteros.prepare_uma_finetune_job(tmp_path / "job", frames(), config, name="ruo2_round0",
                                              replay_train=replay, replay_valid=replay, seeds=(2,))
    job = json.loads((folder / "job.json").read_text())
    assert job["family"] == "uma" and job["task"] == "omat" and job["replay"] is True
    assert job["foundation_model"] == "/work/me/uma-s-1p2p1.pt" and job["foundation_model_sha256"] == "a" * 64
    assert job["training"] == {"epochs": 7, "learning_rate": 1e-4, "weight_decay": 1e-3, "batch_size": 2,
                               "energy_weight": 20.0, "forces_weight": 2.0, "ema_decay": 0.999}
    assert job["runs"] == [{"name": "ruo2_round0_seed2", "seed": 2}]
    script = (folder / "finetune.py").read_text()
    from psteros.mlip_jobs import script_path

    assert script == script_path("uma_finetune.py").read_text()
    assert "fairchem" in script and "mace" not in script.lower().replace("macecalculator", "")
    assert {path.name for path in folder.iterdir()} == {"job.json", "finetune.py", "train.extxyz",
                                                       "replay_train_geometries.extxyz",
                                                       "replay_valid_geometries.extxyz"}
    with pytest.raises(ValueError, match=r"seeds \[3\] are not seeds of the recipe"):
        psteros.prepare_uma_finetune_job(tmp_path / "other", frames(), config, name="x", replay_train=replay,
                                         replay_valid=replay, seeds=(3,))


def test_relax_and_vibrations_jobs_name_the_family_only_for_uma(tmp_path) -> None:
    structure = frames()[0].structure
    mace = mlip.prepare_relax_job(tmp_path / "mace", {"pair": structure}, {"seed1": ("m.model", "Default")})
    assert json.loads((mace / "job.json").read_text())["models"] == {"seed1": {"path": "m.model", "head": "Default"}}
    uma = mlip.prepare_relax_job(tmp_path / "uma", {"pair": structure}, {"seed1": ("m.model", "omat")}, family="uma")
    assert json.loads((uma / "job.json").read_text())["models"] == {
        "seed1": {"path": "m.model", "head": "omat", "family": "uma"}}
    vib = mlip.prepare_vibrations_job(tmp_path / "vib", {"pair": structure}, ("m.model", "omat"), {"pair": [0]},
                                      family="uma")
    assert json.loads((vib / "job.json").read_text())["model"] == {"path": "m.model", "head": "omat", "family": "uma"}
    with pytest.raises(ValueError, match=r"family must be one of \('mace', 'uma', 'nequip'\)"):
        mlip.prepare_relax_job(tmp_path / "bad", {"pair": structure}, {"seed1": ("m", "h")}, family="orb")


def test_relax_job_loads_uma_through_fairchem(tmp_path, monkeypatch) -> None:
    """relax.py with a UMA model asks fairchem for the calculator (fairchem stubbed by EMT)."""

    import subprocess
    import sys
    import textwrap

    pytest.importorskip("ase")
    stub = tmp_path / "stub" / "fairchem" / "core"
    (stub / "units" / "mlip_unit").mkdir(parents=True)
    (stub.parent / "__init__.py").write_text("")
    (stub / "__init__.py").write_text(textwrap.dedent('''
        from ase.calculators.emt import EMT
        CALLS = []
        def FAIRChemCalculator(unit, task_name):
            assert unit == ("predict-unit", task_name == "omat")
            return EMT()
    '''))
    (stub / "units" / "__init__.py").write_text("")
    (stub / "units" / "mlip_unit" / "__init__.py").write_text(
        "def load_predict_unit(path, inference_settings=None, device=None):\n"
        "    assert path.endswith('m.pt') and device == 'cpu'\n"
        "    assert inference_settings.compile is False and inference_settings.merge_mole is True\n"
        "    assert inference_settings.base_precision_dtype == 'float64'\n"
        "    return ('predict-unit', True)\n")
    (stub / "units" / "mlip_unit" / "api").mkdir()
    (stub / "units" / "mlip_unit" / "api" / "__init__.py").write_text("")
    (stub / "units" / "mlip_unit" / "api" / "inference.py").write_text(
        "from dataclasses import dataclass\n"
        "@dataclass\nclass InferenceSettings:\n    tf32: bool = False\n    activation_checkpointing: bool = True\n"
        "    merge_mole: bool = False\n    compile: bool = False\n    base_precision_dtype: object = None\n")
    from pymatgen.core import Lattice, Structure

    copper = Structure(Lattice.cubic(3.61), ["Cu"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.52]])
    folder = mlip.prepare_relax_job(tmp_path / "job", {"cu": copper}, {"seed": (str(tmp_path / "m.pt"), "omat")},
                                    psteros.MlipRelaxConfig(fmax_ev_per_angstrom=0.05), device="cpu", family="uma")
    environment = {"PYTHONPATH": str(tmp_path / "stub"), "PATH": "/usr/bin:/bin"}
    completed = subprocess.run([sys.executable, "-u", "relax.py"], cwd=folder, env=environment, capture_output=True,
                               text=True, timeout=300)
    assert completed.returncode == 0, completed.stderr[-2000:]
    results = json.loads((folder / "results.json").read_text())
    assert results["status"] == "completed" and results["models"]["seed"]["structures"]["cu"]["converged"]


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_graph_with_uma(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("ase")
    from test_mlip_workflow import _computer, _recipes, _structures

    computer = _computer(tmp_path)
    static, mace, job = _recipes(tmp_path, computer)
    structures = _structures()
    uma = psteros.UmaFineTuneConfig.for_foundation("uma-s-1p2p1", "/models/uma-s-1p2p1.pt",
                                                   replay_file=mace.replay_file, epochs=3, seeds=(11, 12))
    vibrations = psteros.VibrationsConfig(molecules=("o2",))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        graph = psteros.build_mlip_surface_workgraph(structures, static, uma, job=job, metallic_bulk="bulk",
                                                     vibrations=vibrations)
    assert not [item for item in caught if "UmaFineTuneConfig" in str(item.message)]  # PBE recipe, no +U oxide
    settings = graph.tasks["mlip_round1_training"].inputs.settings.value.get_dict()
    assert settings["family"] == "uma" and settings["fine_tune"]["task"] == "omat"
    assert psteros.UmaFineTuneConfig(**settings["fine_tune"]) == uma  # what the training calcfunction rebuilds
    plan = json.loads(graph.tasks["mlip_round0_relax"].inputs.nodes.job.value.get_object_content("job.json"))
    assert plan["models"] == {
        "seed11": {"path": "finetune/models/ruo2_round0_seed11.model", "head": "omat", "family": "uma"},
        "seed12": {"path": "finetune/models/ruo2_round0_seed12.model", "head": "omat", "family": "uma"}}
    vib = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert vib["family"] == "uma" and vib["model"] == ["finetune/models/ruo2_round1_seed11.model", "omat"]

    plain = psteros.build_mlip_surface_workgraph(structures, static, mace, job=job, metallic_bulk="bulk",
                                                 vibrations=vibrations)
    assert "family" not in plain.tasks["mlip_round1_training"].inputs.settings.value.get_dict()
    assert "family" not in plain.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    relax = json.loads(plain.tasks["mlip_round0_relax"].inputs.nodes.job.value.get_object_content("job.json"))
    assert relax["models"]["seed11"] == {"path": "finetune/models/ruo2_round0_seed11.model", "head": "Default"}
    with pytest.raises(TypeError, match="MaceFineTuneConfig, a UmaFineTuneConfig or a NequipFineTuneConfig"):
        psteros.build_mlip_surface_workgraph(structures, static, {"model": "uma"}, job=job, metallic_bulk="bulk")


def test_uma_recipe_example(tmp_path, capsys) -> None:
    import importlib.util
    import sys
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "examples" / "mlip_routes" / "uma_recipe.py"
    spec = importlib.util.spec_from_file_location("example_uma_recipe", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.show()
    shown = capsys.readouterr().out
    assert "uma-s-1p2p1" in shown and "* omat" in shown and "PBE recipe, Ru O: no warning" in shown
    assert "Materials Project +U on Fe" in shown
    module.job(tmp_path / "job")
    assert json.loads((tmp_path / "job" / "job.json").read_text())["replay"] is False
