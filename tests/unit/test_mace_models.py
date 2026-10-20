"""The MACE foundation models of the MLIP route: catalogue, recipe and level-of-theory checks."""

from __future__ import annotations

import json
import re
import warnings
from dataclasses import asdict

import pytest

import psteros
from psteros import mlip

MH1_SHA256 = "a522eb7f59c7879963d41586528f4980baf33e086c94aa92e3eafdeccad3be47"


def vasp(incar=None, **changes):
    return psteros.VaspCalculationConfig("vasp@cluster", incar=incar or {"ENCUT": 520, "ISMEAR": 0}, **changes)


def qe(system=None, pseudo_family="SSSP/1.3/PBE/efficiency"):
    parameters = {"CONTROL": {"calculation": "scf"}, "SYSTEM": {"ecutwfc": 60.0, **(system or {})},
                  "ELECTRONS": {"conv_thr": 1e-8}}
    return psteros.QeCalculationConfig("pw@cluster", pseudo_family, parameters)


def cp2k(xc, **dft):
    return psteros.Cp2kCalculationConfig(
        "cp2k@cluster", {"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": {"XC": {"XC_FUNCTIONAL": xc}, **dft}}},
        basis_sets={"Fe": "DZVP-MOLOPT-SR-GTH"}, potentials={"Fe": "GTH-PBE"})


def test_catalogue_entries_are_consistent() -> None:
    models = psteros.MACE_FOUNDATION_MODELS
    assert set(models) >= {"mace-mh-1", "mace-omat-0-medium", "mace-omat-0-small", "mace-matpes-pbe-0",
                           "mace-matpes-r2scan-0", "mace-mpa-0-medium", "mace-mp-0b3-medium", "mace-mp-0-medium"}
    for name, model in models.items():
        assert model.name == name
        assert model.url.startswith("https://github.com/ACEsuit/mace-foundations/releases/download/")
        assert model.filename.endswith(".model") and model.url.endswith("/" + model.filename)
        assert re.fullmatch("[0-9a-f]{64}", model.sha256)
        assert model.licence in ("MIT", "ASL")
        assert model.default_head in model.head_names and model.head().for_slabs
        assert re.fullmatch(r"\d+\.\d+\.\d+", model.min_mace_torch)
        assert model.replay_url.startswith("https://") and model.use
    assert models["mace-mh-1"].sha256 == MH1_SHA256  # the model the route was validated with
    # the heads as the files name them (read from the released files with mace-torch 0.3.16)
    assert set(models["mace-mh-1"].head_names) == {"matpes_r2scan", "mp_pbe_refit_add", "spice_wB97M", "oc20_usemppbe",
                                                   "omol", "omat_pbe"}
    single = {name: model.head_names for name, model in models.items() if name != "mace-mh-1"}
    assert single == {"mace-omat-0-medium": ("default",), "mace-omat-0-small": ("omat_pbe",),
                      "mace-matpes-pbe-0": ("default",), "mace-matpes-r2scan-0": ("default",),
                      "mace-mpa-0-medium": ("default",), "mace-mp-0b3-medium": ("default",),
                      "mace-mp-0-medium": ("Default",)}
    assert [model.name for model in models.values() if model.licence == "MIT"] == [
        "mace-mpa-0-medium", "mace-mp-0b3-medium", "mace-mp-0-medium"]
    with pytest.raises(TypeError):
        models["mine"] = models["mace-mh-1"]  # read-only


def test_unknown_model_and_head_are_errors() -> None:
    with pytest.raises(ValueError, match=r"unknown MACE foundation model 'mace-mp-9'; known: \['mace-matpes"):
        psteros.mace_foundation_model("mace-mp-9")
    with pytest.raises(ValueError, match=r"mace-omat-0-medium has no head 'omat_pbe'; its heads are \['default'\]"):
        psteros.mace_foundation_model("mace-omat-0-medium").head("omat_pbe")
    with pytest.raises(ValueError, match="hubbard"):
        psteros.MaceFoundationHead("x", "PBE", "dudarev", "data")
    with pytest.raises(ValueError, match="default_head"):
        psteros.MaceFoundationModel("m", "https://x/m.model", "0" * 64, "MIT", "0.3.10",
                                    psteros.mace_foundation_model("mace-omat-0-medium").heads, "omat_pbe", "u", "use")


def test_for_foundation_fills_in_head_and_checksum() -> None:
    config = psteros.MaceFineTuneConfig.for_foundation("mace-omat-0-medium", "/work/m/mace-omat-0-medium.model",
                                                       replay_file="mp_traj_combined_omat.xyz", epochs=40)
    assert config.foundation_name == "mace-omat-0-medium" and config.foundation_head == "default"
    assert config.foundation_model_sha256 == psteros.MACE_FOUNDATION_MODELS["mace-omat-0-medium"].sha256
    assert config.epochs == 40 and config.replay_file == "mp_traj_combined_omat.xyz"
    head = psteros.MaceFineTuneConfig.for_foundation("mace-mh-1", "/m.model", head="matpes_r2scan", replay=False)
    assert head.foundation_head == "matpes_r2scan" and head.foundation_model_sha256 == MH1_SHA256
    unchecked = psteros.MaceFineTuneConfig.for_foundation("mace-mh-1", "/m.model", replay=False,
                                                          foundation_model_sha256=None)
    assert unchecked.foundation_model_sha256 is None
    with pytest.raises(ValueError, match="has no head 'omol25'"):
        psteros.MaceFineTuneConfig.for_foundation("mace-mh-1", "/m.model", head="omol25", replay=False)


def test_a_named_recipe_is_checked_against_the_catalogue() -> None:
    with pytest.raises(ValueError, match=r"mace-mpa-0-medium has no head 'omat_pbe'; its heads are \['default'\]"):
        psteros.MaceFineTuneConfig("/m.model", replay=False, foundation_name="mace-mpa-0-medium")
    with pytest.raises(ValueError, match="foundation_model_sha256 is not the checksum of mace-mh-1"):
        psteros.MaceFineTuneConfig("/m.model", replay=False, foundation_model_sha256="0" * 64, foundation_name="mace-mh-1")
    with pytest.raises(ValueError, match="unknown MACE foundation model"):
        psteros.MaceFineTuneConfig("/m.model", replay=False, foundation_name="mh-1")


def test_recipes_without_a_name_are_unchanged(tmp_path) -> None:
    """Leaving foundation_name out: no check, the same training options and the same job folder as before."""

    pymatgen = pytest.importorskip("pymatgen.core")
    config = psteros.MaceFineTuneConfig("/models/any.model", foundation_head="whatever", replay=False, seeds=(7,))
    assert config.foundation_name is None and asdict(config)["foundation_name"] is None
    options = mlip.mace_training_options(config, name="m_seed7", seed=7, atomic_numbers=[8, 44])
    assert options["foundation_model"] == "/models/any.model" and options["foundation_head"] == "whatever"
    assert "foundation_name" not in options
    structure = pymatgen.Structure(pymatgen.Lattice.cubic(6.0), ["Ru", "O"], [[0, 0, 0], [0.3, 0, 0]])
    frame = psteros.LabelledStructure("pair", structure, -10.0, [[0.1, 0, 0], [-0.1, 0, 0]])
    folder = psteros.prepare_finetune_job(tmp_path / "job", [frame], config, name="m")
    job = json.loads((folder / "job.json").read_text())
    assert set(job) == {"foundation_model", "foundation_model_sha256", "foundation_head", "replay", "device", "runs"}
    named = psteros.MaceFineTuneConfig.for_foundation("mace-mh-1", "/models/any.model", replay=False, seeds=(7,))
    assert mlip.mace_training_options(named, name="m_seed7", seed=7, atomic_numbers=[8, 44]) == \
        mlip.mace_training_options(psteros.MaceFineTuneConfig("/models/any.model", foundation_model_sha256=MH1_SHA256,
                                                              replay=False, seeds=(7,)),
                                   name="m_seed7", seed=7, atomic_numbers=[8, 44])


def test_level_of_theory_of_vasp_recipes() -> None:
    level = psteros.recipe_level_of_theory
    assert level(vasp()) == ("PBE", False)
    assert level(psteros.VaspCalculationConfig("vasp@cluster")) == ("PBE", False)  # precision level, PBE POTCARs
    assert level(vasp({"GGA": "PS"})) == ("PBEsol", False)
    assert level(vasp({"gga": "rp", "LDAU": ".TRUE."})) == ("RPBE", True)
    assert level(vasp({"METAGGA": "R2SCAN", "LDAU": True})) == ("r2SCAN", True)
    assert level(vasp({"METAGGA": "SCAN"})) == ("SCAN", False)
    assert level(vasp({"LHFCALC": True, "HFSCREEN": 0.2})) == ("hybrid", False)
    assert level(vasp({"incar": {"METAGGA": "R2SCAN"}, "dynamics": {}})) == ("r2SCAN", False)
    assert level(vasp({"GGA": "AM"})) == (None, False)  # AM05: unknown to the catalogue
    assert level(vasp(potential_family="LDA")) == (None, False)


def test_level_of_theory_of_qe_and_cp2k_recipes() -> None:
    level = psteros.recipe_level_of_theory
    assert level(qe()) == ("PBE", False)
    assert level(qe(pseudo_family="SSSP/1.3/PBEsol/efficiency")) == ("PBEsol", False)
    assert level(qe({"input_dft": "r2scan"})) == ("r2SCAN", False)
    assert level(qe({"lda_plus_u": True})) == ("PBE", True)
    assert level(qe(pseudo_family="PseudoDojo/0.4/LDA/SR/standard/upf")) == (None, False)
    assert level(cp2k({"_": "PBE"})) == ("PBE", False)
    assert level(cp2k({"PBE": {"PARAMETRIZATION": "PBESOL"}})) == ("PBEsol", False)
    assert level(cp2k({"PBE": {}}, PLUS_U_METHOD="LOWDIN")) == ("PBE", True)
    assert level(cp2k({"LIBXC": {"FUNCTIONAL": "MGGA_X_R2SCAN MGGA_C_R2SCAN"}})) == ("r2SCAN", False)
    assert level(cp2k({"BLYP": {}})) == (None, False)
    with pytest.raises(TypeError, match="QE, VASP or CP2K"):
        level({"GGA": "PE"})


def test_level_of_theory_warnings() -> None:
    warn = psteros.foundation_level_warnings
    assert warn("mace-mh-1", None, vasp(), ["Ru", "O"]) == []  # PBE, no +U element: nothing to say
    assert warn("mace-mh-1", "omat_pbe", vasp({"LDAU": True}), ["Fe", "O"]) == []
    assert warn("mace-matpes-r2scan-0", None, vasp({"METAGGA": "R2SCAN"}), ["Fe", "O"]) == []
    assert warn("mace-omat-0-medium", None, vasp(), ["Fe", "Pt"]) == []  # no oxide or fluoride: no +U in MP either
    assert warn("mace-mh-1", "oc20_usemppbe", vasp(), ["Fe", "O"]) == []  # +U scheme of the head not known

    (functional,) = warn("mace-mh-1", None, vasp({"METAGGA": "R2SCAN"}), ["Ru", "O"])
    assert "the DFT recipe uses r2SCAN but mace-mh-1 (head 'omat_pbe') was trained on PBE" in functional
    assert "mace-mh-1 (head 'matpes_r2scan'), mace-matpes-r2scan-0 (head 'default')" in functional
    (hubbard,) = warn("mace-omat-0-medium", None, vasp(), ["Fe", "Ni", "O"])
    assert "Materials Project +U on Fe, Ni in oxides and fluorides, but the DFT recipe has no +U" in hubbard
    assert hubbard.endswith("(trained on PBE without +U: mace-matpes-pbe-0 (head 'default'))")
    (without,) = warn("mace-matpes-pbe-0", None, vasp({"LDAU": True}), ["Mn", "O"])
    assert without == "mace-matpes-pbe-0 (head 'default') was trained without +U, but the DFT recipe uses +U"
    molecules = warn("mace-mh-1", "omol", vasp(), ["Ru", "O"])
    assert molecules[0].startswith("mace-mh-1 (head 'omol') was trained on molecules")
    assert "was trained on wB97M-VV10 data" in molecules[1]
    assert warn("mace-mh-1", None, vasp(potential_family="LDA"), ["Fe", "O"]) == []  # functional unknown: silent


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_builder_warns_only_for_a_named_model(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("ase")
    from test_mlip_workflow import _computer, _recipes, _structures

    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        plain = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                                     active_learning_rounds=0, dft_check=False)
    assert not [item for item in caught if "MaceFineTuneConfig" in str(item.message)]
    settings = plain.tasks["mlip_round0_training"].inputs.settings.value.get_dict()
    assert "foundation_name" not in settings["fine_tune"]  # the stored recipe is the one of earlier graphs

    named = psteros.MaceFineTuneConfig.for_foundation(
        "mace-mh-1", "/models/mace-mh-1.model", head="matpes_r2scan", replay_file=fine_tune.replay_file, epochs=2,
        seeds=(11,))
    with pytest.warns(UserWarning, match="MaceFineTuneConfig: the DFT recipe uses PBE but mace-mh-1 "
                                         r"\(head 'matpes_r2scan'\) was trained on r2SCAN"):
        graph = psteros.build_mlip_surface_workgraph(structures, static, named, job=job, metallic_bulk="bulk",
                                                     active_learning_rounds=0, dft_check=False)
    recipe = graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()["fine_tune"]
    assert recipe["foundation_name"] == "mace-mh-1" and recipe["foundation_head"] == "matpes_r2scan"
    assert recipe["foundation_model_sha256"] == MH1_SHA256
    assert mlip.MaceFineTuneConfig(**recipe) == named  # what the training calcfunction rebuilds


def test_choose_mace_model_example(capsys) -> None:
    import importlib.util
    import sys
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "examples" / "mlip_routes" / "choose_mace_model.py"
    spec = importlib.util.spec_from_file_location("example_choose_mace_model", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.show_catalogue()
    listed = capsys.readouterr().out
    assert "* omat_pbe" in listed and "mp_pbe_refit_add" in listed and "(molecules only)" in listed
    module.check(["Sn", "O"], ldau=False, metagga="R2SCAN")
    checked = capsys.readouterr().out
    assert checked.startswith("recipe: r2SCAN, without +U; elements Sn O")
    assert "ok   mace-mh-1 (matpes_r2scan)" in checked and "ok   mace-matpes-r2scan-0 (default)" in checked
    assert "warn mace-mh-1 (omat_pbe)" in checked and "head 'default', checksum" in checked
