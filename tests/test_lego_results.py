"""Unit tests for lego result extraction and print_stage_results() functions.

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import pytest


# ---------------------------------------------------------------------------
# TestExtractEnergyFromMisc
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestExtractEnergyFromMisc:
    """Tests for teros.core.lego.results._extract_energy_from_misc()."""

    def _extract(self, misc):
        from teros.core.lego.results import _extract_energy_from_misc
        return _extract_energy_from_misc(misc)

    def test_energy_extrapolated(self):
        misc = {'total_energies': {'energy_extrapolated': -123.4}}
        assert self._extract(misc) == pytest.approx(-123.4)

    def test_energy_no_entropy_fallback(self):
        misc = {'total_energies': {'energy_no_entropy': -99.0}}
        assert self._extract(misc) == pytest.approx(-99.0)

    def test_energy_key_fallback(self):
        misc = {'total_energies': {'energy': -50.0}}
        assert self._extract(misc) == pytest.approx(-50.0)

    def test_returns_none_no_keys(self):
        misc = {'total_energies': {'foo': 1}}
        assert self._extract(misc) is None

    def test_flat_dict_top_level(self):
        misc = {'energy_extrapolated': -99.9}
        assert self._extract(misc) == pytest.approx(-99.9)

    def test_returns_float_type(self):
        misc = {'total_energies': {'energy_extrapolated': -123.4}}
        result = self._extract(misc)
        assert isinstance(result, float)

    def test_empty_dict_returns_none(self):
        assert self._extract({}) is None

    def test_empty_total_energies(self):
        misc = {'total_energies': {}}
        assert self._extract(misc) is None


# ---------------------------------------------------------------------------
# TestPrintVaspStageResults
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintVaspStageResults:
    """Tests for teros.core.lego.bricks.vasp.print_stage_results()."""

    def _print(self, index, stage_name, stage_result):
        from teros.core.lego.bricks.vasp import print_stage_results
        print_stage_results(index, stage_name, stage_result)

    def _base_result(self, **overrides):
        result = {
            'energy': None, 'structure': None, 'misc': None,
            'remote': None, 'files': None, 'pk': 1, 'stage': 'relax', 'type': 'vasp',
        }
        result.update(overrides)
        return result

    def test_prints_energy(self, capsys):
        self._print(1, 'relax', self._base_result(energy=-123.456789))
        out = capsys.readouterr().out
        assert '-123.456789' in out

    def test_prints_none_energy(self, capsys):
        self._print(1, 'relax', self._base_result(energy=None))
        out = capsys.readouterr().out
        # Should not contain "Energy:" line at all when None
        assert 'Energy:' not in out

    def test_prints_index_and_name(self, capsys):
        self._print(1, 'relax', self._base_result())
        out = capsys.readouterr().out
        assert '[1]' in out
        assert 'relax' in out

    def test_prints_force_from_misc(self, capsys):
        misc = {'run_status': 'finished', 'maximum_force': 0.0123}
        self._print(1, 'relax', self._base_result(misc=misc))
        out = capsys.readouterr().out
        assert '0.0123' in out


# ---------------------------------------------------------------------------
# TestPrintDosStageResults
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintDosStageResults:
    """Tests for teros.core.lego.bricks.dos.print_stage_results()."""

    def _print(self, index, stage_name, stage_result):
        from teros.core.lego.bricks.dos import print_stage_results
        print_stage_results(index, stage_name, stage_result)

    def _base_result(self, **overrides):
        result = {
            'energy': None, 'scf_misc': None, 'scf_remote': None,
            'scf_retrieved': None, 'dos_misc': None, 'dos_remote': None,
            'files': None, 'pk': 1, 'stage': 'dos_calc', 'type': 'dos',
        }
        result.update(overrides)
        return result

    def test_prints_dos_label(self, capsys):
        self._print(2, 'dos_calc', self._base_result())
        out = capsys.readouterr().out
        assert '(DOS)' in out

    def test_prints_scf_energy(self, capsys):
        self._print(2, 'dos_calc', self._base_result(energy=-50.123456))
        out = capsys.readouterr().out
        assert '-50.123456' in out

    def test_prints_band_gap(self, capsys):
        dos_misc = {
            'band_properties': {'band_gap': 1.2345, 'is_direct_gap': False},
            'fermi_level': 5.678,
        }
        self._print(2, 'dos_calc', self._base_result(dos_misc=dos_misc))
        out = capsys.readouterr().out
        assert '1.2345' in out
        assert 'indirect' in out


# ---------------------------------------------------------------------------
# TestPrintBatchStageResults
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintBatchStageResults:
    """Tests for teros.core.lego.bricks.batch.print_stage_results()."""

    def _print(self, index, stage_name, stage_result):
        from teros.core.lego.bricks.batch import print_stage_results
        print_stage_results(index, stage_name, stage_result)

    def _base_result(self, **overrides):
        result = {
            'calculations': {}, 'pk': 1, 'stage': 'fukui', 'type': 'batch',
        }
        result.update(overrides)
        return result

    def test_prints_batch_label(self, capsys):
        self._print(3, 'fukui', self._base_result())
        out = capsys.readouterr().out
        assert '(BATCH)' in out

    def test_prints_per_calc_energies(self, capsys):
        calcs = {
            'neutral': {'energy': -100.0, 'misc': None, 'remote': None, 'files': None},
            'cation': {'energy': -95.5, 'misc': None, 'remote': None, 'files': None},
        }
        self._print(3, 'fukui', self._base_result(calculations=calcs))
        out = capsys.readouterr().out
        assert 'neutral' in out
        assert '-100.000000' in out
        assert 'cation' in out
        assert '-95.500000' in out

    def test_prints_empty_calculations(self, capsys):
        self._print(3, 'fukui', self._base_result(calculations={}))
        out = capsys.readouterr().out
        assert 'No calculation results' in out


# ---------------------------------------------------------------------------
# TestPrintBaderStageResults
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintBaderStageResults:
    """Tests for teros.core.lego.bricks.bader.print_stage_results()."""

    def _print(self, index, stage_name, stage_result):
        from teros.core.lego.bricks.bader import print_stage_results
        print_stage_results(index, stage_name, stage_result)

    def _base_result(self, **overrides):
        result = {
            'charges': None, 'dat_files': {}, 'pk': 1,
            'stage': 'bader', 'type': 'bader',
        }
        result.update(overrides)
        return result

    def test_prints_bader_label(self, capsys):
        self._print(4, 'bader', self._base_result())
        out = capsys.readouterr().out
        assert '(BADER)' in out

    def test_prints_atom_count_and_charge(self, capsys):
        charges = {
            'atoms': [
                {'index': 1, 'x': 0, 'y': 0, 'z': 0, 'charge': 6.5, 'min_dist': 1.0, 'volume': 10.0,
                 'element': 'Sn', 'valence': 14.0, 'bader_charge': 7.5},
                {'index': 2, 'x': 1, 'y': 1, 'z': 1, 'charge': 8.0, 'min_dist': 0.9, 'volume': 12.0,
                 'element': 'O', 'valence': 6.0, 'bader_charge': -2.0},
            ],
            'total_charge': 14.5,
            'vacuum_charge': 0.0,
            'vacuum_volume': 0.0,
        }
        self._print(4, 'bader', self._base_result(charges=charges))
        out = capsys.readouterr().out
        assert 'Atoms analyzed: 2' in out
        assert '14.50000' in out

    def test_prints_per_atom_charges(self, capsys):
        charges = {
            'atoms': [
                {'index': 1, 'x': 0, 'y': 0, 'z': 0, 'charge': 6.5, 'min_dist': 1.0, 'volume': 10.0,
                 'element': 'Sn', 'valence': 14.0, 'bader_charge': 7.5},
            ],
            'total_charge': 6.5,
            'vacuum_charge': 0.0,
            'vacuum_volume': 0.0,
        }
        self._print(4, 'bader', self._base_result(charges=charges))
        out = capsys.readouterr().out
        assert 'Sn' in out
        assert 'valence=14.0' in out


# ---------------------------------------------------------------------------
# TestPrintThicknessStageResults
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintThicknessStageResults:
    """Tests for teros.core.lego.bricks.thickness.print_stage_results()."""

    def _print(self, index, stage_name, stage_result):
        from teros.core.lego.bricks.thickness import print_stage_results
        print_stage_results(index, stage_name, stage_result)

    def _base_result(self, **overrides):
        result = {
            'convergence_results': None,
            'recommended_structure': None,
            'recommended_layers': None,
            'converged': False,
            'surface_energies': {},
            'pk': 1,
            'stage': 'thickness_test',
            'type': 'thickness',
        }
        result.update(overrides)
        return result

    def test_prints_thickness_label(self, capsys):
        self._print(2, 'thickness_test', self._base_result())
        out = capsys.readouterr().out
        assert '(THICKNESS CONVERGENCE)' in out

    def test_prints_index_and_name(self, capsys):
        self._print(2, 'thickness_test', self._base_result())
        out = capsys.readouterr().out
        assert '[2]' in out
        assert 'thickness_test' in out

    def test_prints_converged_with_layers(self, capsys):
        conv = {
            'miller_indices': [1, 1, 0],
            'summary': {
                'thicknesses': [3, 5, 7, 9],
                'surface_energies_J_m2': [1.50, 1.35, 1.32, 1.31],
                'converged': True,
                'recommended_layers': 7,
                'max_tested_layers': 9,
                'convergence_threshold': 0.01,
            },
        }
        self._print(2, 'thickness_test', self._base_result(convergence_results=conv))
        out = capsys.readouterr().out
        assert 'Converged: YES at 7 layers' in out
        assert '1, 1, 0' in out

    def test_prints_not_converged(self, capsys):
        conv = {
            'miller_indices': [1, 1, 1],
            'summary': {
                'thicknesses': [3, 5],
                'surface_energies_J_m2': [2.0, 1.5],
                'converged': False,
                'recommended_layers': 5,
                'max_tested_layers': 5,
                'convergence_threshold': 0.01,
            },
        }
        self._print(2, 'thickness_test', self._base_result(convergence_results=conv))
        out = capsys.readouterr().out
        assert 'Converged: NO' in out
        assert '5 layers' in out

    def test_prints_surface_energy_table(self, capsys):
        conv = {
            'miller_indices': [1, 1, 0],
            'summary': {
                'thicknesses': [3, 5, 7],
                'surface_energies_J_m2': [1.50, 1.35, 1.32],
                'converged': True,
                'recommended_layers': 5,
                'max_tested_layers': 7,
                'convergence_threshold': 0.05,
            },
        }
        self._print(2, 'thickness_test', self._base_result(convergence_results=conv))
        out = capsys.readouterr().out
        assert '1.5000' in out
        assert '1.3500' in out
        assert '1.3200' in out
        # Delta column: |1.35 - 1.50| * 1000 = 150.0
        assert '150.0' in out

    def test_prints_miller_indices(self, capsys):
        conv = {
            'miller_indices': [1, 1, 0],
            'summary': {
                'thicknesses': [3, 5],
                'surface_energies_J_m2': [1.5, 1.3],
                'converged': False,
                'recommended_layers': 5,
                'max_tested_layers': 5,
                'convergence_threshold': 0.01,
            },
        }
        self._print(2, 'thickness_test', self._base_result(convergence_results=conv))
        out = capsys.readouterr().out
        assert 'Miller indices: (1, 1, 0)' in out

    def test_prints_threshold(self, capsys):
        conv = {
            'miller_indices': [1, 1, 0],
            'summary': {
                'thicknesses': [3, 5],
                'surface_energies_J_m2': [1.5, 1.3],
                'converged': False,
                'recommended_layers': 5,
                'max_tested_layers': 5,
                'convergence_threshold': 0.02,
            },
        }
        self._print(2, 'thickness_test', self._base_result(convergence_results=conv))
        out = capsys.readouterr().out
        assert '0.02' in out

    def test_no_convergence_results(self, capsys):
        """When convergence_results is None, should still print the label."""
        self._print(2, 'thickness_test', self._base_result())
        out = capsys.readouterr().out
        assert 'THICKNESS CONVERGENCE' in out
        # No crash expected
