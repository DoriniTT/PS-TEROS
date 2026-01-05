"""
Import Validation Tests

These tests verify that all modules can be imported without errors.
They catch syntax errors and missing dependencies early.
"""

import pytest
import sys
import importlib


class TestCoreModuleImports:
    """Test that core modules can be imported."""

    def test_import_workflow_presets(self):
        """Test importing workflow_presets module."""
        from teros.core import workflow_presets
        assert workflow_presets is not None

    def test_import_workflow_presets_functions(self):
        """Test importing specific functions from workflow_presets."""
        from teros.core.workflow_presets import (
            WORKFLOW_PRESETS,
            DEFAULT_PRESET,
            list_workflow_presets,
            get_preset_config,
            get_preset_summary,
            resolve_preset,
            validate_preset_inputs,
            validate_flag_dependencies,
        )

        assert callable(list_workflow_presets)
        assert callable(get_preset_config)
        assert callable(resolve_preset)

    def test_import_helper_functions(self):
        """Test importing helper_functions module."""
        from teros.core import helper_functions
        assert helper_functions is not None

    def test_import_slabs(self):
        """Test importing slabs module."""
        from teros.core import slabs
        assert slabs is not None

    def test_import_thermodynamics(self):
        """Test importing thermodynamics module."""
        from teros.core import thermodynamics
        assert thermodynamics is not None

    def test_import_hf(self):
        """Test importing hf module."""
        from teros.core import hf
        assert hf is not None

    def test_import_adsorption_energy(self):
        """Test importing adsorption_energy module."""
        from teros.core import adsorption_energy
        assert adsorption_energy is not None

    def test_import_workgraph(self):
        """Test importing workgraph module."""
        from teros.core import workgraph
        assert workgraph is not None


class TestDependencyImports:
    """Test that external dependencies can be imported."""

    def test_import_numpy(self):
        """Test numpy import."""
        import numpy as np
        assert np is not None

    def test_import_pymatgen_core(self):
        """Test pymatgen core imports."""
        from pymatgen.core import Composition, Structure
        from pymatgen.core.surface import SlabGenerator
        assert Composition is not None
        assert SlabGenerator is not None

    def test_import_ase(self):
        """Test ASE imports."""
        from ase import Atoms
        from ase.io import read, write
        assert Atoms is not None

    def test_import_collections(self):
        """Test standard library imports."""
        from collections import Counter
        from functools import reduce
        from math import gcd
        assert Counter is not None


class TestPythonSyntax:
    """Test Python files for syntax errors."""

    def get_python_files(self):
        """Get list of Python files in teros/core."""
        from pathlib import Path

        core_dir = Path(__file__).parent.parent / 'teros' / 'core'
        if not core_dir.exists():
            pytest.skip("teros/core directory not found")

        python_files = list(core_dir.glob('*.py'))
        return python_files

    def test_all_files_have_valid_syntax(self):
        """Test that all Python files have valid syntax."""
        import py_compile

        python_files = self.get_python_files()

        for filepath in python_files:
            try:
                py_compile.compile(str(filepath), doraise=True)
            except py_compile.PyCompileError as e:
                pytest.fail(f"Syntax error in {filepath}: {e}")


class TestModuleDocstrings:
    """Test that modules have docstrings."""

    def test_workflow_presets_has_docstring(self):
        """Test workflow_presets module has docstring."""
        from teros.core import workflow_presets
        assert workflow_presets.__doc__ is not None
        assert len(workflow_presets.__doc__) > 10
