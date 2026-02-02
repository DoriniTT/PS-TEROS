"""Unit tests for the lego brick registry and brick-level validate_stage() functions.

Also covers Bader file parsers: _parse_acf_dat() and _enrich_acf_dat().

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import os
import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def valid_vasp_stage():
    return {'name': 'relax', 'type': 'vasp', 'incar': {'NSW': 100}, 'restart': None}


@pytest.fixture
def valid_dos_stage():
    return {
        'name': 'dos_calc', 'type': 'dos',
        'scf_incar': {'encut': 400},
        'dos_incar': {'nedos': 2000},
        'structure_from': 'relax',
    }


@pytest.fixture
def valid_batch_stage():
    return {
        'name': 'fukui', 'type': 'batch',
        'structure_from': 'relax',
        'base_incar': {'NSW': 0},
        'calculations': {'neutral': {'incar': {'NELECT': 100}}},
    }


@pytest.fixture
def valid_bader_stage():
    return {'name': 'bader', 'type': 'bader', 'charge_from': 'relax'}


@pytest.fixture
def valid_convergence_stage():
    return {'name': 'conv', 'type': 'convergence'}


@pytest.fixture
def sample_acf_dat_content():
    """Standard ACF.dat content with 2 atoms."""
    return (
        "    #         X           Y           Z        CHARGE     MIN DIST   ATOMIC VOL\n"
        " -----------------------------------------------------------------------\n"
        "    1      0.0000      0.0000      1.2345    6.5432      1.234     12.345\n"
        "    2      2.5000      2.5000      3.6789    8.7654      0.987     15.678\n"
        " -----------------------------------------------------------------------\n"
        "    NUMBER OF ELECTRONS:      15.30860\n"
        "    VACUUM CHARGE:             0.00000\n"
        "    VACUUM VOLUME:             0.00000\n"
    )


# ---------------------------------------------------------------------------
# TestBrickRegistry
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestBrickRegistry:
    """Tests for BRICK_REGISTRY, VALID_BRICK_TYPES, get_brick_module()."""

    def test_valid_brick_types_tuple(self):
        from teros.core.lego.bricks import VALID_BRICK_TYPES
        assert VALID_BRICK_TYPES == ('vasp', 'dos', 'batch', 'bader', 'convergence')

    def test_registry_has_five_entries(self):
        from teros.core.lego.bricks import BRICK_REGISTRY
        assert len(BRICK_REGISTRY) == 5

    def test_get_brick_module_valid_types(self):
        from teros.core.lego.bricks import get_brick_module, VALID_BRICK_TYPES
        for brick_type in VALID_BRICK_TYPES:
            mod = get_brick_module(brick_type)
            assert mod is not None

    def test_each_brick_has_five_functions(self):
        from teros.core.lego.bricks import get_brick_module, VALID_BRICK_TYPES
        required = (
            'validate_stage', 'create_stage_tasks', 'expose_stage_outputs',
            'get_stage_results', 'print_stage_results',
        )
        for brick_type in VALID_BRICK_TYPES:
            mod = get_brick_module(brick_type)
            for fn_name in required:
                assert callable(getattr(mod, fn_name, None)), \
                    f"{brick_type} module missing callable '{fn_name}'"

    def test_get_brick_module_invalid_raises(self):
        from teros.core.lego.bricks import get_brick_module
        with pytest.raises(ValueError, match="Unknown brick type"):
            get_brick_module('invalid')

    def test_get_brick_module_empty_raises(self):
        from teros.core.lego.bricks import get_brick_module
        with pytest.raises(ValueError):
            get_brick_module('')


# ---------------------------------------------------------------------------
# TestVaspValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestVaspValidateStage:
    """Tests for teros.core.lego.bricks.vasp.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.vasp import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_minimal_passes(self, valid_vasp_stage):
        self._validate(valid_vasp_stage)

    def test_missing_incar_raises(self):
        stage = {'name': 'relax', 'restart': None}
        with pytest.raises(ValueError, match="incar"):
            self._validate(stage)

    def test_missing_restart_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}}
        with pytest.raises(ValueError, match="restart"):
            self._validate(stage)

    def test_restart_none_passes(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None}
        self._validate(stage)

    def test_restart_valid_stage_passes(self):
        stage = {'name': 'scf', 'incar': {'NSW': 0}, 'restart': 'relax'}
        self._validate(stage, stage_names={'relax'})

    def test_restart_unknown_stage_raises(self):
        stage = {'name': 'scf', 'incar': {'NSW': 0}, 'restart': 'nope'}
        with pytest.raises(ValueError, match="unknown"):
            self._validate(stage, stage_names={'relax'})

    def test_structure_from_previous_passes(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'structure_from': 'previous'}
        self._validate(stage)

    def test_structure_from_input_passes(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'structure_from': 'input'}
        self._validate(stage)

    def test_structure_from_valid_name_passes(self):
        stage = {'name': 'scf', 'incar': {'NSW': 0}, 'restart': None, 'structure_from': 'relax'}
        self._validate(stage, stage_names={'relax'})

    def test_structure_from_invalid_raises(self):
        stage = {'name': 'scf', 'incar': {'NSW': 0}, 'restart': None, 'structure_from': 'unknown'}
        with pytest.raises(ValueError, match="structure_from"):
            self._validate(stage)

    def test_supercell_valid_passes(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': [2, 2, 1]}
        self._validate(stage)

    def test_supercell_wrong_length_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': [2, 2]}
        with pytest.raises(ValueError, match="supercell"):
            self._validate(stage)

    def test_supercell_not_list_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': '2x2x1'}
        with pytest.raises(ValueError, match="supercell"):
            self._validate(stage)

    def test_supercell_zero_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': [2, 0, 1]}
        with pytest.raises(ValueError, match="positive integers"):
            self._validate(stage)

    def test_supercell_negative_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': [2, -1, 1]}
        with pytest.raises(ValueError, match="positive integers"):
            self._validate(stage)

    def test_supercell_float_raises(self):
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None, 'supercell': [2.0, 2, 1]}
        with pytest.raises(ValueError, match="positive integers"):
            self._validate(stage)

    def test_fix_type_valid_passes(self):
        for fix_type in ('bottom', 'center', 'top'):
            stage = {
                'name': 'relax', 'incar': {'NSW': 100}, 'restart': None,
                'fix_type': fix_type, 'fix_thickness': 3.0,
            }
            self._validate(stage)

    def test_fix_type_invalid_raises(self):
        stage = {
            'name': 'relax', 'incar': {'NSW': 100}, 'restart': None,
            'fix_type': 'left', 'fix_thickness': 3.0,
        }
        with pytest.raises(ValueError, match="fix_type"):
            self._validate(stage)

    def test_fix_type_zero_thickness_raises(self):
        stage = {
            'name': 'relax', 'incar': {'NSW': 100}, 'restart': None,
            'fix_type': 'bottom', 'fix_thickness': 0,
        }
        with pytest.raises(ValueError, match="fix_thickness"):
            self._validate(stage)

    def test_fix_type_negative_thickness_raises(self):
        stage = {
            'name': 'relax', 'incar': {'NSW': 100}, 'restart': None,
            'fix_type': 'top', 'fix_thickness': -1,
        }
        with pytest.raises(ValueError, match="fix_thickness"):
            self._validate(stage)

    def test_fix_type_none_no_check(self):
        """When fix_type is not set, no thickness validation happens."""
        stage = {'name': 'relax', 'incar': {'NSW': 100}, 'restart': None}
        self._validate(stage)  # no fix_thickness required


# ---------------------------------------------------------------------------
# TestDosValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestDosValidateStage:
    """Tests for teros.core.lego.bricks.dos.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.dos import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_passes(self, valid_dos_stage):
        self._validate(valid_dos_stage, stage_names={'relax'})

    def test_missing_scf_incar_raises(self):
        stage = {'name': 'dos', 'dos_incar': {'nedos': 2000}, 'structure_from': 'relax'}
        with pytest.raises(ValueError, match="scf_incar"):
            self._validate(stage, stage_names={'relax'})

    def test_missing_dos_incar_raises(self):
        stage = {'name': 'dos', 'scf_incar': {'encut': 400}, 'structure_from': 'relax'}
        with pytest.raises(ValueError, match="dos_incar"):
            self._validate(stage, stage_names={'relax'})

    def test_missing_structure_from_raises(self):
        stage = {'name': 'dos', 'scf_incar': {'encut': 400}, 'dos_incar': {'nedos': 2000}}
        with pytest.raises(ValueError, match="structure_from"):
            self._validate(stage)

    def test_structure_from_unknown_raises(self):
        stage = {
            'name': 'dos', 'scf_incar': {'encut': 400},
            'dos_incar': {'nedos': 2000}, 'structure_from': 'nope',
        }
        with pytest.raises(ValueError, match="previous stage"):
            self._validate(stage, stage_names={'relax'})

    def test_structure_from_valid_passes(self):
        stage = {
            'name': 'dos', 'scf_incar': {'encut': 400},
            'dos_incar': {'nedos': 2000}, 'structure_from': 'relax',
        }
        self._validate(stage, stage_names={'relax'})


# ---------------------------------------------------------------------------
# TestBatchValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestBatchValidateStage:
    """Tests for teros.core.lego.bricks.batch.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.batch import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_passes(self, valid_batch_stage):
        self._validate(valid_batch_stage, stage_names={'relax'})

    def test_missing_structure_from_raises(self):
        stage = {
            'name': 'fukui',
            'base_incar': {'NSW': 0},
            'calculations': {'neutral': {'incar': {}}},
        }
        with pytest.raises(ValueError, match="structure_from"):
            self._validate(stage, stage_names={'relax'})

    def test_missing_base_incar_raises(self):
        stage = {
            'name': 'fukui', 'structure_from': 'relax',
            'calculations': {'neutral': {'incar': {}}},
        }
        with pytest.raises(ValueError, match="base_incar"):
            self._validate(stage, stage_names={'relax'})

    def test_missing_calculations_raises(self):
        stage = {
            'name': 'fukui', 'structure_from': 'relax',
            'base_incar': {'NSW': 0},
        }
        with pytest.raises(ValueError, match="calculations"):
            self._validate(stage, stage_names={'relax'})

    def test_empty_calculations_raises(self):
        stage = {
            'name': 'fukui', 'structure_from': 'relax',
            'base_incar': {'NSW': 0}, 'calculations': {},
        }
        with pytest.raises(ValueError, match="calculations"):
            self._validate(stage, stage_names={'relax'})

    def test_structure_from_unknown_raises(self):
        stage = {
            'name': 'fukui', 'structure_from': 'nope',
            'base_incar': {'NSW': 0},
            'calculations': {'neutral': {'incar': {}}},
        }
        with pytest.raises(ValueError, match="previous stage"):
            self._validate(stage, stage_names={'relax'})

    def test_structure_from_valid_passes(self):
        stage = {
            'name': 'fukui', 'structure_from': 'relax',
            'base_incar': {'NSW': 0},
            'calculations': {'neutral': {'incar': {}}},
        }
        self._validate(stage, stage_names={'relax'})


# ---------------------------------------------------------------------------
# TestBaderValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestBaderValidateStage:
    """Tests for teros.core.lego.bricks.bader.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.bader import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_passes(self, valid_bader_stage):
        self._validate(valid_bader_stage, stage_names={'relax'})

    def test_missing_charge_from_raises(self):
        stage = {'name': 'bader'}
        with pytest.raises(ValueError, match="charge_from"):
            self._validate(stage, stage_names={'relax'})

    def test_charge_from_unknown_raises(self):
        stage = {'name': 'bader', 'charge_from': 'nope'}
        with pytest.raises(ValueError, match="previous stage"):
            self._validate(stage, stage_names={'relax'})

    def test_charge_from_valid_passes(self):
        stage = {'name': 'bader', 'charge_from': 'relax'}
        self._validate(stage, stage_names={'relax'})


# ---------------------------------------------------------------------------
# TestConvergenceValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestConvergenceValidateStage:
    """Tests for teros.core.lego.bricks.convergence.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.convergence import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_minimal_passes(self):
        stage = {'name': 'conv', 'type': 'convergence'}
        self._validate(stage)

    def test_valid_full_passes(self, valid_convergence_stage):
        self._validate(valid_convergence_stage)

    def test_structure_from_valid_passes(self):
        stage = {'name': 'conv', 'type': 'convergence', 'structure_from': 'relax'}
        self._validate(stage, stage_names={'relax'})

    def test_structure_from_unknown_raises(self):
        stage = {'name': 'conv', 'type': 'convergence', 'structure_from': 'nope'}
        with pytest.raises(ValueError, match="previous stage"):
            self._validate(stage, stage_names={'relax'})

    def test_conv_settings_not_dict_raises(self):
        stage = {'name': 'conv', 'type': 'convergence', 'conv_settings': 'bad'}
        with pytest.raises(ValueError, match="conv_settings"):
            self._validate(stage)

    def test_threshold_negative_raises(self):
        stage = {'name': 'conv', 'type': 'convergence', 'convergence_threshold': -0.1}
        with pytest.raises(ValueError, match="positive"):
            self._validate(stage)

    def test_threshold_zero_raises(self):
        stage = {'name': 'conv', 'type': 'convergence', 'convergence_threshold': 0}
        with pytest.raises(ValueError, match="positive"):
            self._validate(stage)

    def test_threshold_string_raises(self):
        stage = {'name': 'conv', 'type': 'convergence', 'convergence_threshold': 'bad'}
        with pytest.raises(ValueError, match="number"):
            self._validate(stage)

    def test_conv_settings_dict_passes(self):
        stage = {'name': 'conv', 'type': 'convergence',
                 'conv_settings': {'cutoff_start': 300}}
        self._validate(stage)


# ---------------------------------------------------------------------------
# TestParseAcfDat
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestParseAcfDat:
    """Tests for teros.core.lego.bricks.bader._parse_acf_dat()."""

    def _parse(self, filepath):
        from teros.core.lego.bricks.bader import _parse_acf_dat
        return _parse_acf_dat(filepath)

    def test_parses_two_atoms(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        result = self._parse(str(acf))
        assert len(result['atoms']) == 2

    def test_atom_fields_present(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        result = self._parse(str(acf))
        expected_keys = {'index', 'x', 'y', 'z', 'charge', 'min_dist', 'volume'}
        for atom in result['atoms']:
            assert expected_keys.issubset(atom.keys())

    def test_atom_values_correct(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        result = self._parse(str(acf))
        atom0 = result['atoms'][0]
        assert atom0['charge'] == pytest.approx(6.5432)
        assert atom0['x'] == pytest.approx(0.0)
        assert atom0['y'] == pytest.approx(0.0)
        assert atom0['z'] == pytest.approx(1.2345)

    def test_total_electrons(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        result = self._parse(str(acf))
        assert result['total_charge'] == pytest.approx(15.30860)

    def test_vacuum_values(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        result = self._parse(str(acf))
        assert result['vacuum_charge'] == pytest.approx(0.0)
        assert result['vacuum_volume'] == pytest.approx(0.0)

    def test_empty_file(self, tmp_path):
        acf = tmp_path / 'ACF.dat'
        acf.write_text('')
        result = self._parse(str(acf))
        assert result['atoms'] == []
        assert result['total_charge'] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# TestEnrichAcfDat
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestEnrichAcfDat:
    """Tests for teros.core.lego.bricks.bader._enrich_acf_dat()."""

    def _enrich(self, filepath, elements, zval_dict):
        from teros.core.lego.bricks.bader import _enrich_acf_dat
        _enrich_acf_dat(filepath, elements, zval_dict)

    def test_header_has_new_columns(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        self._enrich(str(acf), ['Sn', 'O'], {'Sn': 14.0, 'O': 6.0})
        lines = acf.read_text().splitlines()
        header = lines[0]
        assert 'ELEMENT' in header
        assert 'VALENCE' in header
        assert 'BADER_CHARGE' in header

    def test_correct_element_assignment(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        self._enrich(str(acf), ['Sn', 'O'], {'Sn': 14.0, 'O': 6.0})
        text = acf.read_text()
        assert 'Sn' in text
        assert 'O' in text

    def test_bader_charge_calculation(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        self._enrich(str(acf), ['Sn', 'O'], {'Sn': 14.0, 'O': 6.0})
        lines = acf.read_text().splitlines()
        # Atom 1: VALENCE=14.0, CHARGE=6.5432 → BADER_CHARGE = 14.0 - 6.5432 = 7.4568
        atom1_line = lines[2]  # skip header + separator
        assert '7.4568' in atom1_line
        # Atom 2: VALENCE=6.0, CHARGE=8.7654 → BADER_CHARGE = 6.0 - 8.7654 = -2.7654
        atom2_line = lines[3]
        assert '-2.7654' in atom2_line

    def test_summary_lines_preserved(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        self._enrich(str(acf), ['Sn', 'O'], {'Sn': 14.0, 'O': 6.0})
        text = acf.read_text()
        assert 'NUMBER OF ELECTRONS' in text

    def test_empty_zval_dict(self, tmp_path, sample_acf_dat_content):
        acf = tmp_path / 'ACF.dat'
        acf.write_text(sample_acf_dat_content)
        self._enrich(str(acf), ['Sn', 'O'], {})
        lines = acf.read_text().splitlines()
        # valence=0.0, bader_charge = 0.0 - charge = -charge
        # Atom 1: bader_charge = 0 - 6.5432 = -6.5432
        atom1_line = lines[2]
        assert '-6.5432' in atom1_line


# ---------------------------------------------------------------------------
# TestConvergenceValidateStage
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestConvergenceValidateStage:
    """Tests for teros.core.lego.bricks.convergence.validate_stage()."""

    def _validate(self, stage, stage_names=None):
        from teros.core.lego.bricks.convergence import validate_stage
        if stage_names is None:
            stage_names = set()
        validate_stage(stage, stage_names)

    def test_valid_minimal_passes(self, valid_convergence_stage):
        self._validate(valid_convergence_stage)

    def test_structure_from_valid_passes(self):
        stage = {'name': 'conv', 'structure_from': 'relax'}
        self._validate(stage, stage_names={'relax'})

    def test_structure_from_unknown_raises(self):
        stage = {'name': 'conv', 'structure_from': 'nope'}
        with pytest.raises(ValueError, match="previous stage"):
            self._validate(stage, stage_names={'relax'})

    def test_conv_settings_not_dict_raises(self):
        stage = {'name': 'conv', 'conv_settings': 'bad'}
        with pytest.raises(ValueError, match="conv_settings"):
            self._validate(stage)

    def test_conv_settings_dict_passes(self):
        stage = {'name': 'conv', 'conv_settings': {'cutoff_start': 300}}
        self._validate(stage)

    def test_threshold_positive_passes(self):
        stage = {'name': 'conv', 'convergence_threshold': 0.005}
        self._validate(stage)

    def test_threshold_zero_raises(self):
        stage = {'name': 'conv', 'convergence_threshold': 0}
        with pytest.raises(ValueError, match="positive"):
            self._validate(stage)

    def test_threshold_negative_raises(self):
        stage = {'name': 'conv', 'convergence_threshold': -0.001}
        with pytest.raises(ValueError, match="positive"):
            self._validate(stage)

    def test_threshold_non_numeric_raises(self):
        stage = {'name': 'conv', 'convergence_threshold': 'bad'}
        with pytest.raises(ValueError, match="number"):
            self._validate(stage)
