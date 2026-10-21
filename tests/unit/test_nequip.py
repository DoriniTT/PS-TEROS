"""NequIP and Allegro (the NequIP framework) on the MLIP route: recipe, job folder and graph (no nequip needed)."""

from __future__ import annotations

import json
import warnings

import pytest

import psteros
from psteros import mlip


def frames():
    from pymatgen.core import Lattice, Structure

    first = Structure(Lattice.cubic(6.0), ["Mg", "O"], [[0, 0, 0], [0.3, 0, 0]])
    second = Structure(Lattice.cubic(6.0), ["Mg", "O", "O"], [[0, 0, 0], [0.3, 0, 0], [0.7, 0, 0]])
    return [psteros.LabelledStructure("pair", first, -10.0, [[0.1, 0, 0], [-0.1, 0, 0]]),
            psteros.LabelledStructure("triple", second, -15.0, [[0.1, 0, 0], [-0.1, 0, 0], [0, 0, 0]])]


def test_nequip_catalogue_and_recipe() -> None:
    models = psteros.NEQUIP_FOUNDATION_MODELS
    assert set(models) == {"nequip-oam-s", "nequip-oam-m", "nequip-oam-l", "nequip-oam-xl", "allegro-oam-l"}
    large = psteros.nequip_foundation_model("nequip-oam-l")
    assert large.model_id == "mir-group/NequIP-OAM-L:0.1" and large.download == "nequip.net:mir-group/NequIP-OAM-L:0.1"
    assert large.url == "https://www.nequip.net/models/mir-group/NequIP-OAM-L:0.1" and large.licence == "CC BY 4.0"
    assert models["allegro-oam-l"].architecture == "Allegro" and large.level.functional == "PBE"
    with pytest.raises(ValueError, match="unknown NequIP foundation model 'nequip-oam-xxl'"):
        psteros.nequip_foundation_model("nequip-oam-xxl")

    config = psteros.NequipFineTuneConfig.for_foundation("nequip-oam-l", "/work/me/NequIP-OAM-L-0.1.nequip.zip",
                                                         replay_file="r.xyz")
    assert config.foundation_name == "nequip-oam-l" and config.epochs == 30 and config.seeds == (20260920,)
    with pytest.raises(ValueError, match="must be a packaged model"):
        psteros.NequipFineTuneConfig("/work/me/model.ckpt", replay=False)
    with pytest.raises(ValueError, match="replay=True needs replay_file"):
        psteros.NequipFineTuneConfig("/m.nequip.zip")
    with pytest.raises(ValueError, match="ema_decay"):
        psteros.NequipFineTuneConfig("/m.nequip.zip", replay=False, ema_decay=0.0)


def test_nequip_level_warnings() -> None:
    vasp = psteros.VaspCalculationConfig("vasp@cluster", incar={"ENCUT": 520})
    assert psteros.nequip_level_warnings("nequip-oam-l", vasp, ["Mg", "O"]) == []
    (hubbard,) = psteros.nequip_level_warnings("allegro-oam-l", vasp, ["Co", "O"])
    assert hubbard == ("allegro-oam-l was trained with the Materials Project +U on Co in oxides and fluorides, "
                       "but the DFT recipe has no +U")
    (functional,) = psteros.nequip_level_warnings("nequip-oam-s", psteros.VaspCalculationConfig(
        "vasp@cluster", incar={"METAGGA": "R2SCAN"}), ["Mg", "O"])
    assert functional.startswith("the DFT recipe uses r2SCAN but nequip-oam-s was trained on PBE data")


def test_nequip_job_folder(tmp_path) -> None:
    from ase import Atoms

    from psteros.mlip_jobs import script_path

    config = psteros.NequipFineTuneConfig("/work/me/NequIP-OAM-L-0.1.nequip.zip", replay_file="r.xyz", epochs=4,
                                          seeds=(1, 2))
    replay = [Atoms("MgO", positions=[[0, 0, 0], [2.1, 0, 0]], cell=[8, 8, 8], pbc=True)]
    folder = psteros.prepare_nequip_finetune_job(tmp_path / "job", frames(), config, name="mgo_round0",
                                                 replay_train=replay, replay_valid=replay)
    job = json.loads((folder / "job.json").read_text())
    assert job["family"] == "nequip" and job["foundation_model"] == "/work/me/NequIP-OAM-L-0.1.nequip.zip"
    assert job["training"] == {"epochs": 4, "learning_rate": 1e-4, "batch_size": 2, "energy_weight": 1.0,
                               "forces_weight": 1.0, "ema_decay": 0.999}
    assert job["runs"] == [{"name": "mgo_round0_seed1", "seed": 1}, {"name": "mgo_round0_seed2", "seed": 2}]
    assert (folder / "finetune.py").read_text() == script_path("nequip_finetune.py").read_text()
    with pytest.raises(ValueError, match="replay structures were given to a fine-tune without replay"):
        psteros.prepare_nequip_finetune_job(tmp_path / "other", frames(),
                                            psteros.NequipFineTuneConfig("/m.nequip.zip", replay=False), name="x",
                                            replay_train=replay, replay_valid=replay)


def test_relax_job_loads_nequip_packages(tmp_path) -> None:
    """relax.py with a NequIP-framework model asks nequip for the calculator (nequip stubbed by EMT)."""

    import subprocess
    import sys
    import textwrap

    pytest.importorskip("ase")
    stub = tmp_path / "stub" / "nequip" / "integrations"
    stub.mkdir(parents=True)
    (stub.parent / "__init__.py").write_text("")
    (stub / "__init__.py").write_text("")
    (stub / "ase.py").write_text(textwrap.dedent('''
        from ase.calculators.emt import EMT
        class NequIPCalculator:
            @classmethod
            def _from_saved_model(cls, path, device=None, chemical_species_to_atom_type_map=None):
                assert path.endswith(".nequip.zip") and device == "cpu" and chemical_species_to_atom_type_map is True
                return EMT()
    '''))
    from pymatgen.core import Lattice, Structure

    copper = Structure(Lattice.cubic(3.61), ["Cu"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.52]])
    folder = mlip.prepare_relax_job(tmp_path / "job", {"cu": copper}, {"seed": ("m.nequip.zip", "Default")},
                                    psteros.MlipRelaxConfig(fmax_ev_per_angstrom=0.05), device="cpu", family="nequip")
    assert json.loads((folder / "job.json").read_text())["models"]["seed"]["family"] == "nequip"
    environment = {"PYTHONPATH": str(tmp_path / "stub"), "PATH": "/usr/bin:/bin"}
    completed = subprocess.run([sys.executable, "-u", "relax.py"], cwd=folder, env=environment, capture_output=True,
                               text=True, timeout=300)
    assert completed.returncode == 0, completed.stderr[-2000:]
    assert json.loads((folder / "results.json").read_text())["models"]["seed"]["structures"]["cu"]["converged"]


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_graph_with_nequip(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("ase")
    from test_mlip_workflow import _computer, _recipes, _structures

    computer = _computer(tmp_path)
    static, mace, job = _recipes(tmp_path, computer)
    structures = _structures()
    nequip = psteros.NequipFineTuneConfig.for_foundation("allegro-oam-l", "/models/Allegro-OAM-L-0.1.nequip.zip",
                                                         replay_file=mace.replay_file, epochs=3, seeds=(11,))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        graph = psteros.build_mlip_surface_workgraph(structures, static, nequip, job=job, metallic_bulk="bulk",
                                                     vibrations=psteros.VibrationsConfig(molecules=("o2",)))
    assert not [item for item in caught if "NequipFineTuneConfig" in str(item.message)]
    settings = graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()
    assert settings["family"] == "nequip" and psteros.NequipFineTuneConfig(**settings["fine_tune"]) == nequip
    plan = json.loads(graph.tasks["mlip_round1_relax"].inputs.nodes.job.value.get_object_content("job.json"))
    assert plan["models"] == {"seed11": {"path": "finetune/models/ruo2_round1_seed11.nequip.zip", "head": "Default",
                                         "family": "nequip"}}
    vib = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert vib["family"] == "nequip" and vib["model"][0] == "finetune/models/ruo2_round1_seed11.nequip.zip"


def test_nequip_recipe_example(tmp_path, capsys) -> None:
    import importlib.util
    import sys
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "examples" / "mlip_routes" / "nequip_recipe.py"
    spec = importlib.util.spec_from_file_location("example_nequip_recipe", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.show()
    shown = capsys.readouterr().out
    assert "allegro-oam-l" in shown and "nequip.net:mir-group/NequIP-OAM-L:0.1" in shown
    assert "PBE recipe, Mg O: no warning" in shown and "Materials Project +U on Co" in shown
    module.job(tmp_path / "job")
    assert json.loads((tmp_path / "job" / "job.json").read_text())["family"] == "nequip"
