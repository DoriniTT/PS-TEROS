"""Tier-1 tests for calculation blocks and reference descriptions (no AiiDA needed)."""

from __future__ import annotations

import pytest

import psteros
from psteros.blocks import check_blocks


def test_blocks_lower_case_their_incar_and_check_names() -> None:
    assert psteros.Relax(incar={"NSW": 50}).incar == {"nsw": 50}
    with pytest.raises(ValueError, match="block name"):
        psteros.Static(name="static scf")
    with pytest.raises(ValueError, match="ibrion must be 5 or 6"):
        psteros.Vibrations(ibrion=8)
    with pytest.raises(ValueError, match="through the Vibrations fields"):
        psteros.Vibrations(incar={"POTIM": 0.02})
    with pytest.raises(ValueError, match="nfree"):
        psteros.Vibrations(nfree=3)


def test_vibrations_choose_finite_differences_by_phase() -> None:
    assert psteros.Vibrations().required("gas")["ibrion"] == 5
    assert psteros.Vibrations().required("solid")["ibrion"] == 6
    assert psteros.Vibrations(ibrion=5).required("solid") == {"ibrion": 5, "potim": 0.015, "nfree": 2, "nsw": 1}
    assert psteros.Static().required("solid") == {"nsw": 0}


def test_vibrations_avoid_a_k_point_change_that_band_parallelisation_forbids() -> None:
    # VASP stops with "requested a change of the k-point set ... remove NPAR" when displaced
    # (lower-symmetry) cells run with NCORE > 1; found on Lovelace (par128, 128 ranks).
    # A gas drops the symmetry (NCORE = 1 on 128 ranks gives 128 bands for 6 occupied ones and
    # "EDDDAV: Call to ZHEGV failed"); a solid keeps its symmetry-reduced displacements and uses NCORE = 1.
    assert psteros.Vibrations().defaults("gas") == {"isif": 2, "isym": 0}
    assert psteros.Vibrations().defaults("solid") == {"isif": 2, "ncore": 1}
    assert psteros.Vibrations(incar={"NCORE": 4}).incar == {"ncore": 4}
    # Other blocks keep their defaults: the fix adds nothing to them.
    assert psteros.Relax().defaults("solid") == {}
    assert psteros.Static().defaults("solid") == {"ibrion": -1}


def test_block_sequence_links_only_to_earlier_blocks() -> None:
    assert len(check_blocks([psteros.Relax(), psteros.Vibrations(structure_from="relax")])) == 2
    with pytest.raises(ValueError, match="earlier block"):
        check_blocks([psteros.Static(structure_from="relax"), psteros.Relax()])
    with pytest.raises(ValueError, match="duplicate"):
        check_blocks([psteros.Relax(), psteros.Relax()])
    with pytest.raises(ValueError, match="at least one"):
        check_blocks([])
    with pytest.raises(TypeError, match="Relax, Static or Vibrations"):
        check_blocks(["relax"])


def test_reference_system_requires_gas_thermochemistry_metadata() -> None:
    o2 = psteros.ReferenceSystem(psteros.triplet_o2_cell(), "gas", symmetry_number=2, spin=1.0)
    assert o2.supercell == (1, 1, 1)
    with pytest.raises(ValueError, match="symmetry_number and spin"):
        psteros.ReferenceSystem(psteros.triplet_o2_cell(), "gas")
    with pytest.raises(ValueError, match="own box"):
        psteros.ReferenceSystem(psteros.triplet_o2_cell(), "gas", supercell=(2, 2, 2), symmetry_number=2, spin=1.0)
    with pytest.raises(ValueError, match="leave them unset"):
        psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid", spin=0.0)
    with pytest.raises(ValueError, match="phase"):
        psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "bulk")
    with pytest.raises(ValueError, match="supercell"):
        psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid", supercell=(2, 2))
    with pytest.raises(ValueError, match="under 'INCAR'"):
        psteros.ReferenceSystem(
            psteros.alpha_sn_bulk(),
            "solid",
            override=psteros.CalculationOverride(parameters={"SYSTEM": {"ecutwfc": 40}}),
        )


def test_reference_builder_rejects_a_qe_recipe_before_touching_aiida() -> None:
    qe = psteros.SurfaceWorkflowConfig(
        backend="qe",
        calculation=psteros.QeCalculationConfig(
            "qe@x", "sssp", {"CONTROL": {}, "SYSTEM": {}, "ELECTRONS": {}}
        ),
    )
    reference = psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")
    with pytest.raises(TypeError, match="requires a VASP recipe"):
        psteros.build_vasp_reference_workgraph({"sn": reference}, qe)
