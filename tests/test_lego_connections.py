"""Unit tests for the lego brick connection / port validation system.

Tests the PORTS declarations, validate_connections(), port type registry,
prerequisites checking, conditional output warnings, and the full
validation pipeline.

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import warnings

import pytest

# Import from connections.py directly (no AiiDA dependency).
# We use importlib to avoid triggering teros.core.__init__ which pulls in AiiDA.
import importlib.util
import os
_connections_path = os.path.join(
    os.path.dirname(__file__), os.pardir,
    'teros', 'core', 'lego', 'bricks', 'connections.py',
)
_connections_path = os.path.normpath(_connections_path)
_spec = importlib.util.spec_from_file_location(
    'teros.core.lego.bricks.connections', _connections_path,
)
_connections = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_connections)

PORT_TYPES = _connections.PORT_TYPES
ALL_PORTS = _connections.ALL_PORTS
VASP_PORTS = _connections.VASP_PORTS
DOS_PORTS = _connections.DOS_PORTS
BATCH_PORTS = _connections.BATCH_PORTS
BADER_PORTS = _connections.BADER_PORTS
CONVERGENCE_PORTS = _connections.CONVERGENCE_PORTS
validate_connections = _connections.validate_connections
_validate_port_types = _connections._validate_port_types
_evaluate_conditional = _connections._evaluate_conditional


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def relax_stage():
    """Valid VASP relaxation stage (nsw=100)."""
    return {
        'name': 'relax',
        'type': 'vasp',
        'incar': {'encut': 520, 'nsw': 100, 'ibrion': 2, 'isif': 3},
        'restart': None,
        'retrieve': ['CONTCAR', 'OUTCAR'],
    }


@pytest.fixture
def scf_stage():
    """Valid VASP SCF stage (nsw=0) with Bader prerequisites."""
    return {
        'name': 'scf',
        'type': 'vasp',
        'incar': {
            'encut': 520, 'nsw': 0, 'ibrion': -1,
            'lcharg': True, 'laechg': True,
        },
        'restart': None,
        'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
    }


@pytest.fixture
def scf_stage_no_bader():
    """VASP SCF stage without Bader prerequisites."""
    return {
        'name': 'scf',
        'type': 'vasp',
        'incar': {'encut': 520, 'nsw': 0, 'ibrion': -1},
        'restart': None,
        'retrieve': ['OUTCAR'],
    }


@pytest.fixture
def dos_stage():
    """Valid DOS stage pointing at relax."""
    return {
        'name': 'dos',
        'type': 'dos',
        'structure_from': 'relax',
        'scf_incar': {'encut': 520},
        'dos_incar': {'nedos': 3000},
    }


@pytest.fixture
def batch_stage():
    """Valid batch stage pointing at relax."""
    return {
        'name': 'charge_scan',
        'type': 'batch',
        'structure_from': 'relax',
        'base_incar': {'encut': 520, 'nsw': 0},
        'calculations': {
            'neutral': {},
            'plus1': {'incar': {'nelect': 47}},
        },
    }


@pytest.fixture
def bader_stage():
    """Valid bader stage pointing at scf."""
    return {
        'name': 'bader',
        'type': 'bader',
        'charge_from': 'scf',
    }


@pytest.fixture
def convergence_stage():
    """Valid convergence stage (no structure_from)."""
    return {
        'name': 'conv',
        'type': 'convergence',
        'conv_settings': {'cutoff_start': 300},
    }


# ---------------------------------------------------------------------------
# TestPortTypeRegistry (issue #16)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPortTypeRegistry:
    """Every port type used in PORTS declarations must be in PORT_TYPES."""

    def test_all_vasp_output_types_recognized(self):
        for port_name, port in VASP_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"VASP output '{port_name}' has unrecognized type '{port['type']}'"

    def test_all_dos_output_types_recognized(self):
        for port_name, port in DOS_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"DOS output '{port_name}' has unrecognized type '{port['type']}'"

    def test_all_batch_output_types_recognized(self):
        for port_name, port in BATCH_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"Batch output '{port_name}' has unrecognized type '{port['type']}'"

    def test_all_bader_output_types_recognized(self):
        for port_name, port in BADER_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"Bader output '{port_name}' has unrecognized type '{port['type']}'"

    def test_all_convergence_output_types_recognized(self):
        for port_name, port in CONVERGENCE_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"Convergence output '{port_name}' has unrecognized type '{port['type']}'"

    def test_typo_in_port_type_caught(self):
        """A misspelled type should be caught by _validate_port_types."""
        bad_ports = {
            'inputs': {},
            'outputs': {
                'data': {'type': 'retrived', 'description': 'typo'},
            },
        }
        with pytest.raises(ValueError, match="Unknown port type 'retrived'"):
            _validate_port_types(bad_ports, 'test_brick')

    def test_all_input_types_recognized(self):
        """All input port types across all bricks must be in PORT_TYPES."""
        for brick_name, ports in ALL_PORTS.items():
            for port_name, port in ports['inputs'].items():
                assert port['type'] in PORT_TYPES, \
                    f"{brick_name} input '{port_name}' has unrecognized type '{port['type']}'"


# ---------------------------------------------------------------------------
# TestPortDeclarations (issues #5, #6, #7)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPortDeclarations:
    """Verify PORTS dicts are complete and consistent with actual code."""

    def test_vasp_has_five_outputs(self):
        assert len(VASP_PORTS['outputs']) == 5

    def test_vasp_outputs_include_structure(self):
        assert 'structure' in VASP_PORTS['outputs']

    def test_vasp_structure_is_conditional(self):
        assert 'conditional' in VASP_PORTS['outputs']['structure']

    def test_dos_has_no_structure_output(self):
        """DOS brick must NOT declare a structure output (decision #4)."""
        assert 'structure' not in DOS_PORTS['outputs']

    def test_dos_has_nine_outputs(self):
        """DOS must declare all 9 outputs including remote/retrieved (issue #5)."""
        assert len(DOS_PORTS['outputs']) == 9

    def test_dos_has_scf_remote(self):
        assert 'scf_remote' in DOS_PORTS['outputs']

    def test_dos_has_scf_retrieved(self):
        assert 'scf_retrieved' in DOS_PORTS['outputs']

    def test_dos_has_dos_remote(self):
        assert 'dos_remote' in DOS_PORTS['outputs']

    def test_dos_has_dos_retrieved(self):
        assert 'dos_retrieved' in DOS_PORTS['outputs']

    def test_batch_has_no_structure_output(self):
        assert 'structure' not in BATCH_PORTS['outputs']

    def test_batch_outputs_are_per_calculation(self):
        for port in BATCH_PORTS['outputs'].values():
            assert port.get('per_calculation') is True

    def test_bader_has_four_outputs(self):
        """Bader must declare charges + acf + bcf + avf (issue #6)."""
        assert len(BADER_PORTS['outputs']) == 4
        for key in ('charges', 'acf', 'bcf', 'avf'):
            assert key in BADER_PORTS['outputs']

    def test_bader_both_inputs_have_compatible_bricks(self):
        """Both bader inputs should have compatible_bricks (issue #15)."""
        for input_name, port in BADER_PORTS['inputs'].items():
            assert 'compatible_bricks' in port, \
                f"Bader input '{input_name}' missing compatible_bricks"

    def test_convergence_has_no_structure_output(self):
        """Convergence brick must not declare a structure output (issue #7)."""
        assert 'structure' not in CONVERGENCE_PORTS['outputs']

    def test_convergence_has_three_outputs(self):
        assert len(CONVERGENCE_PORTS['outputs']) == 3

    def test_convergence_structure_input_is_optional(self):
        """Convergence structure input is optional (falls back to initial)."""
        assert CONVERGENCE_PORTS['inputs']['structure']['required'] is False


# ---------------------------------------------------------------------------
# TestConditionalEvaluation (issue #9)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestConditionalEvaluation:
    """Test the conditional output evaluation mechanism."""

    def test_none_conditional_is_true(self):
        assert _evaluate_conditional(None, {}) is True

    def test_string_conditional_raises(self):
        """String conditionals must be rejected (no eval)."""
        with pytest.raises(ValueError, match="must be a dict"):
            _evaluate_conditional('nsw > 0', {'incar': {'nsw': 100}})

    def test_nsw_greater_than_zero_true(self):
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {'incar': {'nsw': 100}}) is True

    def test_nsw_greater_than_zero_false(self):
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {'incar': {'nsw': 0}}) is False

    def test_missing_incar_key_defaults_to_zero(self):
        """VASP defaults most keys to 0; missing key → 0."""
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {'incar': {}}) is False

    def test_missing_incar_dict_defaults_to_zero(self):
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {}) is False

    def test_equals_operator(self):
        cond = {'incar_key': 'ismear', 'operator': '==', 'value': -5}
        assert _evaluate_conditional(cond, {'incar': {'ismear': -5}}) is True
        assert _evaluate_conditional(cond, {'incar': {'ismear': 0}}) is False

    def test_unknown_operator_raises(self):
        cond = {'incar_key': 'nsw', 'operator': '~', 'value': 0}
        with pytest.raises(ValueError, match="Unknown operator"):
            _evaluate_conditional(cond, {'incar': {'nsw': 0}})


# ---------------------------------------------------------------------------
# TestValidateConnectionsBasic
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBasic:
    """Basic connection validation: happy paths and simple errors."""

    def test_single_vasp_stage_passes(self, relax_stage):
        warnings = validate_connections([relax_stage])
        assert warnings == []

    def test_two_vasp_stages_chain_passes(self, relax_stage, scf_stage):
        warnings = validate_connections([relax_stage, scf_stage])
        # scf has nsw=0, but relax→scf is auto(previous) with relax having structure
        assert warnings == []

    def test_vasp_then_dos_passes(self, relax_stage, dos_stage):
        warnings = validate_connections([relax_stage, dos_stage])
        assert warnings == []

    def test_vasp_then_batch_passes(self, relax_stage, batch_stage):
        warnings = validate_connections([relax_stage, batch_stage])
        assert warnings == []

    def test_vasp_then_convergence_passes(self, relax_stage, convergence_stage):
        conv = {**convergence_stage, 'structure_from': 'relax'}
        warnings = validate_connections([relax_stage, conv])
        assert warnings == []

    def test_convergence_without_structure_from_passes(self, convergence_stage):
        """Convergence with no structure_from is valid (optional input)."""
        warnings = validate_connections([convergence_stage])
        assert warnings == []

    def test_full_pipeline_passes(self, relax_stage, scf_stage, dos_stage,
                                  batch_stage, bader_stage):
        stages = [relax_stage, scf_stage, dos_stage, batch_stage, bader_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)


# ---------------------------------------------------------------------------
# TestValidateConnectionsOutputExists (decision #4, #5)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsOutputExists:
    """Validate that source stages produce the required output type."""

    def test_dos_structure_from_vasp_passes(self, relax_stage, dos_stage):
        validate_connections([relax_stage, dos_stage])

    def test_dos_structure_from_dos_rejected(self, relax_stage, dos_stage):
        """DOS has no structure output → structure_from='dos' must fail."""
        dos1 = {**dos_stage, 'name': 'dos1'}
        dos2 = {
            'name': 'dos2', 'type': 'dos', 'structure_from': 'dos1',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos1, dos2])

    def test_batch_structure_from_dos_rejected(self, relax_stage, dos_stage):
        """Batch needs structure, DOS doesn't produce one."""
        batch = {
            'name': 'batch1', 'type': 'batch', 'structure_from': 'dos',
            'base_incar': {'encut': 520}, 'calculations': {'a': {}},
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos_stage, batch])

    def test_error_message_suggests_valid_stages(self, relax_stage, dos_stage):
        """Error message should list stages that DO produce the needed type."""
        batch = {
            'name': 'batch1', 'type': 'batch', 'structure_from': 'dos',
            'base_incar': {'encut': 520}, 'calculations': {'a': {}},
        }
        with pytest.raises(ValueError, match="relax") as exc_info:
            validate_connections([relax_stage, dos_stage, batch])
        assert 'relax' in str(exc_info.value)

    def test_structure_from_unknown_stage_rejected(self, relax_stage):
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'nonexistent',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, dos])


# ---------------------------------------------------------------------------
# TestValidateConnectionsBrickCompat
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBrickCompat:
    """Test compatible_bricks constraints."""

    def test_bader_from_vasp_passes(self, relax_stage, scf_stage, bader_stage):
        validate_connections([relax_stage, scf_stage, bader_stage])

    def test_bader_from_dos_rejected(self, relax_stage, dos_stage):
        """Bader requires charge_from to point at a vasp brick."""
        bader = {'name': 'bader', 'type': 'bader', 'charge_from': 'dos'}
        with pytest.raises(ValueError, match="compatible with bricks.*vasp"):
            validate_connections([relax_stage, dos_stage, bader])

    def test_bader_from_batch_rejected(self, relax_stage, batch_stage):
        bader = {'name': 'bader', 'type': 'bader', 'charge_from': 'charge_scan'}
        with pytest.raises(ValueError, match="compatible with bricks.*vasp"):
            validate_connections([relax_stage, batch_stage, bader])

    def test_bader_from_convergence_rejected(self, relax_stage, convergence_stage):
        """Convergence doesn't produce 'retrieved' → type check catches it
        before compatible_bricks check even fires."""
        conv = {**convergence_stage, 'structure_from': 'relax'}
        bader = {'name': 'bader', 'type': 'bader', 'charge_from': 'conv'}
        with pytest.raises(ValueError, match="doesn't produce it"):
            validate_connections([relax_stage, conv, bader])


# ---------------------------------------------------------------------------
# TestValidateConnectionsPrerequisites (issue #7, #11)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsPrerequisites:
    """Test that prerequisite INCAR/retrieve requirements are enforced."""

    def test_bader_with_all_prereqs_passes(self, relax_stage, scf_stage,
                                           bader_stage):
        """SCF has laechg=True, lcharg=True, retrieves AECCAR files."""
        validate_connections([relax_stage, scf_stage, bader_stage])

    def test_bader_missing_laechg_rejected(self, relax_stage, bader_stage):
        """SCF without laechg=True should fail."""
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'lcharg': True},
            'restart': None,
            'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
        }
        with pytest.raises(ValueError, match="laechg"):
            validate_connections([relax_stage, scf_bad, bader_stage])

    def test_bader_missing_lcharg_rejected(self, relax_stage, bader_stage):
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'laechg': True},
            'restart': None,
            'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
        }
        with pytest.raises(ValueError, match="lcharg"):
            validate_connections([relax_stage, scf_bad, bader_stage])

    def test_bader_missing_retrieve_files_rejected(self, relax_stage, bader_stage):
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'laechg': True, 'lcharg': True},
            'restart': None,
            'retrieve': ['OUTCAR'],  # missing AECCAR0, AECCAR2, CHGCAR
        }
        with pytest.raises(ValueError, match="Missing retrieve"):
            validate_connections([relax_stage, scf_bad, bader_stage])

    def test_bader_missing_both_incar_and_retrieve(self, relax_stage, bader_stage):
        """Both INCAR and retrieve are wrong → error mentions both."""
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'retrieve': ['OUTCAR'],
        }
        with pytest.raises(ValueError) as exc_info:
            validate_connections([relax_stage, scf_bad, bader_stage])
        msg = str(exc_info.value)
        assert 'laechg' in msg
        assert 'AECCAR0' in msg

    def test_bader_partial_retrieve_rejected(self, relax_stage, bader_stage):
        """Having some but not all required files should fail."""
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'laechg': True, 'lcharg': True},
            'restart': None,
            'retrieve': ['AECCAR0', 'OUTCAR'],  # missing AECCAR2, CHGCAR
        }
        with pytest.raises(ValueError, match="Missing retrieve"):
            validate_connections([relax_stage, scf_bad, bader_stage])


# ---------------------------------------------------------------------------
# TestValidateConnectionsConditional (decision #6, issue #13, #14)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsConditional:
    """Test conditional output warnings."""

    def test_structure_from_relaxation_no_warning(self, relax_stage, dos_stage):
        """Relaxation (nsw=100) → no warning about structure."""
        warnings = validate_connections([relax_stage, dos_stage])
        assert len(warnings) == 0

    def test_structure_from_static_warns(self, relax_stage, scf_stage):
        """A VASP stage after static SCF (nsw=0) should warn."""
        # Third stage gets structure auto from scf (nsw=0)
        third = {
            'name': 'rerelax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100, 'ibrion': 2},
            'restart': None,
        }
        warnings = validate_connections([relax_stage, scf_stage, third])
        assert len(warnings) == 1
        assert 'nsw=0' in warnings[0]

    def test_explicit_structure_from_static_warns(self, relax_stage, scf_stage):
        """Explicit structure_from pointing at nsw=0 stage should warn."""
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'scf',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        warnings = validate_connections([relax_stage, scf_stage, dos])
        assert len(warnings) == 1
        assert 'scf' in warnings[0]

    def test_first_stage_static_no_warning(self):
        """First stage with nsw=0 doesn't warn (no downstream auto yet)."""
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        warnings = validate_connections([scf])
        assert len(warnings) == 0

    def test_structure_from_input_no_warning(self, relax_stage, scf_stage):
        """structure_from='input' bypasses conditional check."""
        third = {
            'name': 'fresh', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
            'structure_from': 'input',
        }
        warnings = validate_connections([relax_stage, scf_stage, third])
        assert len(warnings) == 0


# ---------------------------------------------------------------------------
# TestValidateConnectionsAutoResolution (issue #2, #12)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsAutoResolution:
    """Test VASP 'auto' structure resolution and its edge cases."""

    def test_first_stage_auto_always_valid(self):
        """First stage gets initial structure → always valid."""
        stage = {
            'name': 'relax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        warnings = validate_connections([stage])
        assert warnings == []

    def test_auto_previous_with_vasp_passes(self, relax_stage):
        """Second VASP stage gets structure from previous VASP → valid."""
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        validate_connections([relax_stage, scf])

    def test_auto_previous_after_dos_fails(self, relax_stage, dos_stage):
        """VASP after DOS with structure_from='previous' (default) should fail
        because DOS has no structure output."""
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            # structure_from defaults to 'previous' → dos
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_auto_previous_after_batch_fails(self, relax_stage, batch_stage):
        """Batch has no structure output → auto(previous) should fail."""
        vasp_after_batch = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, batch_stage, vasp_after_batch])

    def test_auto_previous_after_bader_fails(self, relax_stage, scf_stage,
                                             bader_stage):
        """Bader has no structure output → auto(previous) should fail."""
        vasp_after_bader = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, scf_stage, bader_stage,
                                  vasp_after_bader])

    def test_explicit_structure_from_bypasses_previous(self, relax_stage,
                                                       dos_stage):
        """VASP after DOS with explicit structure_from='relax' should pass."""
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'relax',
        }
        validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_structure_from_input_always_valid(self, relax_stage, dos_stage):
        """structure_from='input' uses initial structure → always valid."""
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'input',
        }
        validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_structure_from_nonexistent_rejected(self, relax_stage):
        """Explicit structure_from pointing to unknown stage should fail."""
        vasp = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'nonexistent',
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, vasp])

    def test_error_suggests_stages_with_structure(self, relax_stage, dos_stage):
        """When auto fails, error should list stages that have structure."""
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError) as exc_info:
            validate_connections([relax_stage, dos_stage, vasp_after_dos])
        assert 'relax' in str(exc_info.value)


# ---------------------------------------------------------------------------
# TestValidateConnectionsRestart (issue #3)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsRestart:
    """Test restart connection validation."""

    def test_restart_none_passes(self, relax_stage):
        validate_connections([relax_stage])

    def test_restart_from_vasp_passes(self, relax_stage):
        """VASP produces remote_folder → restart is valid."""
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': 'relax',
        }
        validate_connections([relax_stage, scf])

    # NOTE: restart validation through port system requires 'restart' to be
    # a non-'auto' source. Currently VASP PORTS declares restart_folder
    # with source='restart' and required=False. When restart=None, it's
    # skipped. When restart='relax', the port system should check that
    # 'relax' produces 'remote_folder'. This test verifies that.


# ---------------------------------------------------------------------------
# TestValidateConnectionsBaderMultipleInputs (issue #4)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBaderMultipleInputs:
    """Test that bader's two inputs from the same source are both validated."""

    def test_both_inputs_satisfied_by_vasp(self, relax_stage, scf_stage,
                                           bader_stage):
        """VASP produces both 'retrieved' and 'structure' → both pass."""
        validate_connections([relax_stage, scf_stage, bader_stage])

    def test_missing_charge_from_rejected(self, relax_stage):
        """Bader without charge_from field should fail."""
        bader = {'name': 'bader', 'type': 'bader'}
        with pytest.raises(ValueError, match="charge_from.*missing"):
            validate_connections([relax_stage, bader])


# ---------------------------------------------------------------------------
# TestValidateConnectionsBatchOutputs (issue #8)
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBatchOutputs:
    """Test that batch template outputs are handled correctly."""

    def test_batch_registers_outputs(self, relax_stage, batch_stage):
        """Batch stage should register its outputs without error."""
        validate_connections([relax_stage, batch_stage])

    def test_batch_has_no_structure_output(self, relax_stage, batch_stage):
        """A stage after batch with auto(previous) should fail."""
        post = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, batch_stage, post])


# ---------------------------------------------------------------------------
# TestValidateConnectionsMissingSourceField
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsMissingSourceField:
    """Test that missing required connection fields produce clear errors."""

    def test_dos_missing_structure_from_rejected(self, relax_stage):
        dos = {
            'name': 'dos', 'type': 'dos',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="structure_from.*missing"):
            validate_connections([relax_stage, dos])

    def test_batch_missing_structure_from_rejected(self, relax_stage):
        batch = {
            'name': 'batch', 'type': 'batch',
            'base_incar': {'encut': 520}, 'calculations': {'a': {}},
        }
        with pytest.raises(ValueError, match="structure_from.*missing"):
            validate_connections([relax_stage, batch])

    def test_bader_missing_charge_from_rejected(self, relax_stage, scf_stage):
        bader = {'name': 'bader', 'type': 'bader'}
        with pytest.raises(ValueError, match="charge_from.*missing"):
            validate_connections([relax_stage, scf_stage, bader])


# ---------------------------------------------------------------------------
# TestValidateConnectionsFullPipeline
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsFullPipeline:
    """Integration tests: full pipeline from the plan's example."""

    def test_example_pipeline_passes(self, relax_stage, scf_stage, dos_stage,
                                     batch_stage, bader_stage):
        """The exact pipeline from PLAN_CONNECTIONS.md should validate."""
        stages = [relax_stage, scf_stage, dos_stage, batch_stage, bader_stage]
        warnings = validate_connections(stages)
        # Only the scf→dos (structure_from='relax') conditional doesn't fire;
        # but scf (nsw=0) auto-chains → the vasp after relax doesn't warn
        # because auto(previous) checks relax (nsw=100).
        # However, scf is nsw=0 and nothing references its structure output
        # via auto(previous), so no warnings expected for this pipeline.
        assert isinstance(warnings, list)

    def test_reversed_order_fails(self, relax_stage, bader_stage):
        """Bader before its source stage should fail."""
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([bader_stage, relax_stage])

    def test_circular_reference_impossible(self, relax_stage):
        """A stage can't reference itself (it's not in available_outputs yet)."""
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'dos',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, dos])

    def test_forward_reference_fails(self, relax_stage, dos_stage):
        """DOS can't reference a stage that comes after it."""
        dos_forward = {**dos_stage, 'structure_from': 'late_relax'}
        late_relax = {
            'name': 'late_relax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, dos_forward, late_relax])


# ---------------------------------------------------------------------------
# TestValidateConnectionsEdgeCases
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsEdgeCases:
    """Edge cases and corner scenarios."""

    def test_empty_stages_no_error(self):
        """Empty list should return empty warnings (upstream validates emptiness)."""
        warnings = validate_connections([])
        assert warnings == []

    def test_single_dos_stage_fails(self):
        """DOS needs structure_from, and with no previous stages it fails."""
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'relax',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([dos])

    def test_convergence_with_unknown_structure_from_fails(self, relax_stage):
        conv = {
            'name': 'conv', 'type': 'convergence',
            'structure_from': 'nonexistent',
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, conv])

    def test_multiple_baders_from_different_sources(self, relax_stage):
        """Two bader stages pointing at different VASP stages."""
        scf1 = {
            'name': 'scf1', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'laechg': True, 'lcharg': True},
            'restart': None,
            'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
        }
        scf2 = {
            'name': 'scf2', 'type': 'vasp',
            'incar': {'encut': 400, 'nsw': 0, 'laechg': True, 'lcharg': True},
            'restart': None,
            'structure_from': 'relax',
            'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
        }
        bader1 = {'name': 'bader1', 'type': 'bader', 'charge_from': 'scf1'}
        bader2 = {'name': 'bader2', 'type': 'bader', 'charge_from': 'scf2'}
        validate_connections([relax_stage, scf1, scf2, bader1, bader2])

    def test_dos_after_convergence_with_explicit_structure_from(
        self, relax_stage, convergence_stage
    ):
        """DOS after convergence must use structure_from pointing at relax."""
        conv = {**convergence_stage, 'structure_from': 'relax'}
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'relax',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        validate_connections([relax_stage, conv, dos])

    def test_many_dos_stages_from_same_source(self, relax_stage):
        """Multiple DOS stages can all point at the same relaxation."""
        stages = [relax_stage]
        for i in range(5):
            stages.append({
                'name': f'dos_{i}', 'type': 'dos', 'structure_from': 'relax',
                'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
            })
        validate_connections(stages)

    def test_default_type_is_vasp(self):
        """Stage without 'type' field should be treated as vasp."""
        stage = {
            'name': 'relax',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        validate_connections([stage])
