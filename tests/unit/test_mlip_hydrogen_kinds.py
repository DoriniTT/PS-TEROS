"""Hydrogen kinds on the MLIP route (issue #5).

The models of the route see elements, not kinds: a pseudo-hydrogen (``H0p75``, ``H1p25``) is hydrogen for them.
"""

from __future__ import annotations

import pytest

from psteros.mlip_workflow import _hydrogen_kinds, _replay_elements

pymatgen = pytest.importorskip("pymatgen.core")


def structure(species, kinds=None):
    lattice = pymatgen.Lattice.tetragonal(4.0, 20.0)
    coords = [[0.1 * i, 0.1 * i, 0.05 * i] for i in range(len(species))]
    properties = {"kind_name": kinds} if kinds else None
    return pymatgen.Structure(lattice, species, coords, site_properties=properties)


def gaas_slab(*hydrogen_kinds):
    species = ["Ga", "As", "Ga", "As", *["H"] * len(hydrogen_kinds)]
    return structure(species, ["Ga", "As", "Ga", "As", *hydrogen_kinds])


def plain(structures):
    return [pymatgen.Structure(item.lattice, item.species, item.frac_coords) for item in structures.values()]


def test_hydrogen_kinds() -> None:
    assert _hydrogen_kinds(gaas_slab("H0p75", "H0p75")) == {"H0p75"}
    assert _hydrogen_kinds(gaas_slab("H0p75", "H1p25")) == {"H0p75", "H1p25"}
    assert _hydrogen_kinds(structure(["H", "H"])) == {"H"}  # no kind names: hydrogen
    assert _hydrogen_kinds(structure(["H", "H"], ["H1", "H2"])) == {"H"}  # other kind names: hydrogen
    assert _hydrogen_kinds(structure(["Ga", "As"], ["Ga", "As"])) == set()


def test_replay_of_a_passivated_study_has_no_hydrogen() -> None:
    passivated = {"bulk": structure(["Ga", "As"]), "slab": gaas_slab("H0p75", "H1p25")}
    assert _replay_elements(passivated, plain(passivated)) == ["As", "Ga"]
    hydrogenated = {"bulk": structure(["Ga", "As"]), "ref_H": structure(["H", "H"])}
    assert _replay_elements(hydrogenated, plain(hydrogenated)) == ["As", "Ga", "H"]  # real hydrogen: as before
    oxide = {"bulk": structure(["Ru", "O"])}
    assert _replay_elements(oxide, plain(oxide)) == ["O", "Ru"]


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_polar_route_replays_no_real_hydrogen(tmp_path) -> None:
    """Regression: the replay set of a polar study held real-hydrogen structures (GaH, H2, AsH) next to its
    pseudo-hydrogens, all one species for the model."""

    pytest.importorskip("aiida_vasp")
    pytest.importorskip("aiida_shell")
    import argparse

    from aiida import orm
    from ase import Atoms
    from ase.io import write
    from test_mlip_route import _example

    computer = orm.Computer(label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
                            scheduler_type="core.direct", workdir=str(tmp_path)).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    orm.InstalledCode(computer=computer, filepath_executable="/opt/mace/bin/python", label="mace",
                      default_calc_job_plugin="core.shell").store()
    frames = [Atoms(symbols, positions=[[0, 0, 0], [1.6 + 0.01 * i, 0, 0]], cell=[8, 8, 8], pbc=True)
              for i in range(6) for symbols in ("GaAs", "GaH", "H2", "AsH")]
    write(tmp_path / "replay.xyz", frames, format="extxyz")
    args = argparse.Namespace(code=f"vasp@{computer.label}", mace_code=f"mace@{computer.label}",
                              model="/models/mace-mh-1.model", replay=str(tmp_path / "replay.xyz"), queue="gpu",
                              vibrations=False, submit=False)
    example = _example("gaas_polar_mlip")
    route = example.route()
    static, fine_tune, job = example.recipes(args, route.study)
    graph = route.build(static, fine_tune, job=job)
    selection = graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()["replay_selection"]
    assert selection["elements"] == ["As", "Ga"] and selection["kept"] == 6  # the GaAs pairs only
