# Test for calculate_surface_energy_ternary in ternary.py
# This test dynamically disables the @task.calcfunction() decorator during import.

import types
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms


def _load_module_without_decorator(module_path: str):
    """
    Load the ternary module after commenting the @task.calcfunction() decorator.
    This avoids requiring an AiiDA engine context during unit testing.
    """
    src_path = Path(module_path)
    source = src_path.read_text()
    source = source.replace("@task.calcfunction()", "# @task.calcfunction()  # disabled for tests")
    # Stub aiida.orm imports if present
    if "from aiida.orm import Dict, load_node, Int" in source:
        source = source.replace(
            "from aiida.orm import Dict, load_node, Int",
            (
                "# Stubbed aiida.orm for tests\n"
                "class _DummyAiiDADict:\n"
                "    def __init__(self, dict=None):\n"
                "        self._d = dict or {}\n"
                "    def get_dict(self):\n"
                "        return self._d\n"
                "class _DummyAiiDAInt:\n"
                "    def __init__(self, value):\n"
                "        self.value = int(value)\n"
                "def load_node(*args, **kwargs):\n"
                "    raise RuntimeError('load_node is stubbed in tests')\n"
                "Dict = _DummyAiiDADict\n"
                "Int = _DummyAiiDAInt\n"
            ),
        )
    # If aiida_workgraph isn't available, stub the import inline to avoid ImportError
    if "from aiida_workgraph import task" in source:
        source = source.replace(
            "from aiida_workgraph import task",
            (
                "# Stubbed aiida_workgraph.task for tests\n"
                "class _DummyTask:\n"
                "    def calcfunction(self):\n"
                "        def _decorator(fn):\n"
                "            return fn\n"
                "        return _decorator\n"
                "task = _DummyTask()\n"
            ),
        )
    module = types.ModuleType("ternary_testable")
    compiled = compile(source, str(src_path), "exec")
    exec(compiled, module.__dict__)
    return module


class DummyDict:
    def __init__(self, dict):
        self._d = dict

    def get_dict(self):
        return self._d


class DummyInt:
    def __init__(self, value: int):
        self.value = int(value)


class DummyStructure:
    def __init__(self, atoms: Atoms):
        self._atoms = atoms

    def get_ase(self) -> Atoms:
        return self._atoms


def _build_ag3po4_bulk_structure():
    """Create a minimal Ag3PO4 bulk-like StructureData (one formula unit) for testing.
    Positions are arbitrary; only composition and counts are relevant for this unit test.
    """
    # Simple cubic-ish cell
    cell = np.diag([6.0, 6.0, 6.0])

    # One formula unit: Ag3PO4 (3 Ag, 1 P, 4 O)
    symbols = ["Ag", "Ag", "Ag", "P", "O", "O", "O", "O"]
    positions = [
        [1.0, 1.0, 1.0],  # Ag
        [4.5, 1.5, 1.5],  # Ag
        [1.5, 4.5, 1.5],  # Ag
        [3.0, 3.0, 3.0],  # P
        [3.0, 3.0, 1.8],  # O
        [3.0, 1.8, 3.0],  # O
        [1.8, 3.0, 3.0],  # O
        [4.2, 3.0, 3.0],  # O
    ]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return DummyStructure(atoms)


def _build_ag3po4_slab_structure():
    """Create a minimal Ag3PO4 slab-like StructureData for testing surface area and counts.
    Uses a large c lattice vector to emulate vacuum; lateral area is 10x10.
    """
    cell = np.array(
        [
            [10.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 30.0],  # vacuum along c
        ]
    )

    # One formula unit placed roughly in the middle along c
    z0 = 15.0
    symbols = ["Ag", "Ag", "Ag", "P", "O", "O", "O", "O"]
    positions = [
        [2.0, 2.0, z0],
        [7.5, 2.5, z0],
        [2.5, 7.5, z0],
        [5.0, 5.0, z0],
        [5.0, 5.0, z0 - 1.2],
        [5.0, 3.8, z0],
        [3.8, 5.0, z0],
        [6.2, 5.0, z0],
    ]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=[True, True, False])
    return DummyStructure(atoms)


@pytest.mark.parametrize("code", ["VASP"])  # extend if needed: QUANTUM_ESPRESSO, CP2K
def test_calculate_surface_energy_ternary_ag3po4(code):
    # Load module with decorator disabled
    ternary_path = "/home/thiagotd/projects/aiida_teros/teros/functions/thermodynamics/ternary.py"
    mod = _load_module_without_decorator(ternary_path)
    # Monkeypatch module Dict to our dummy so the function returns DummyDict
    mod.Dict = DummyDict
    calc = mod.calculate_surface_energy_ternary

    # Bulk structure and parameters (tentative values for testing only)
    bulk_structure = _build_ag3po4_bulk_structure()
    if code in ["QUANTUM_ESPRESSO", "CP2K"]:
        bulk_parameters = DummyDict({"energy": -50.0})
    else:  # VASP
        bulk_parameters = DummyDict({"total_energies": {"energy_extrapolated": -50.0}})

    # Formation enthalpy and reference energies per atom (tentative)
    formation_enthalpy = DummyDict(
        {
            "ag_energy_per_atom": -2.6,
            "p_energy_per_atom": -5.0,
            "o_energy_per_atom": -4.5,
            "formation_enthalpy_ev": -6.0,
        }
    )

    # One simple slab
    slab_structures = {"slab1": _build_ag3po4_slab_structure()}
    if code in ["QUANTUM_ESPRESSO", "CP2K"]:
        slab_parameters = {"slab1": DummyDict({"energy": -49.0})}
    else:  # VASP
        slab_parameters = {"slab1": DummyDict({"total_energies": {"energy_extrapolated": -49.0}})}

    # Small sampling for quick test
    sampling = DummyInt(3)

    # Execute
    results = calc(
        bulk_structure,
        bulk_parameters,
        sampling=sampling,
        formation_enthalpy=formation_enthalpy,
        code=code,
        slab_structures=slab_structures,
        slab_parameters=slab_parameters,
    )

    # Basic validations
    assert isinstance(results, dict)
    assert "s_0" in results

    s0 = results["s_0"]
    assert hasattr(s0, "get_dict"), "Expected a dict-like node in results"

    out = s0.get_dict()
    print("\n--- Output for s_0 ---")
    for k, v in out.items():
        if isinstance(v, dict) and len(v) > 10:
            print(f"{k}: <dict with {len(v)} entries>")
        else:
            print(f"{k}: {v}")
    for key in [
        "phi",
        "Gamma_M_vs_Nref",
        "Gamma_O_vs_Nref",
        "gamma_values_grid",
        "gamma_values_fixed_muM_zero",
        "area_A2",
        "element_M_independent",
        "element_N_reference",
        "bulk_stoichiometry_MxNyOz",
        "slab_atom_counts",
        "reference_energies_per_atom",
        "E_slab_eV",
        "E_bulk_fu_eV",
    ]:
        assert key in out, f"Missing key in output: {key}"

    assert isinstance(out["phi"], float)
    assert isinstance(out["Gamma_M_vs_Nref"], float)
    assert isinstance(out["Gamma_O_vs_Nref"], float)
    assert isinstance(out["gamma_values_grid"], dict)
    assert isinstance(out["gamma_values_fixed_muM_zero"], dict)
    assert out["element_M_independent"] in {"Ag", "P"}  # metals
    assert out["element_N_reference"] in {"Ag", "P"}

    # Ensure some gamma values were computed
    assert len(out["gamma_values_grid"]) == int(sampling.value) ** 2
    assert len(out["gamma_values_fixed_muM_zero"]) == int(sampling.value)
