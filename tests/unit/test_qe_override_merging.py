"""Regression tests for per-role QE namelist overrides."""

from __future__ import annotations

from psteros.backends.qe import _merged_parameters
from psteros.config import CalculationOverride


def test_qe_role_override_deep_merges_only_the_requested_namelist_values() -> None:
    base = {
        "SYSTEM": {"ecutwfc": 90.0, "occupations": "smearing", "nbnd": 280},
        "ELECTRONS": {"conv_thr": 1.0e-10},
    }
    merged = _merged_parameters(
        base,
        CalculationOverride(parameters={"SYSTEM": {"degauss": 0.01}}),
    )
    assert merged == {
        "SYSTEM": {
            "ecutwfc": 90.0,
            "occupations": "smearing",
            "nbnd": 280,
            "degauss": 0.01,
        },
        "ELECTRONS": {"conv_thr": 1.0e-10},
    }
    assert base["SYSTEM"] == {"ecutwfc": 90.0, "occupations": "smearing", "nbnd": 280}
