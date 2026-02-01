"""Unit tests for the lego brick connection / port validation system.

Tests the PORTS declarations, validate_connections(), port type registry,
prerequisites checking, conditional output warnings, and the full
validation pipeline.

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

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
HUBBARD_RESPONSE_PORTS = _connections.HUBBARD_RESPONSE_PORTS
HUBBARD_ANALYSIS_PORTS = _connections.HUBBARD_ANALYSIS_PORTS
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
def ground_state_stage():
    """Valid VASP ground state stage for Hubbard U (lorbit=11, lwave, lcharg)."""
    return {
        'name': 'gs',
        'type': 'vasp',
        'incar': {
            'encut': 520, 'nsw': 0, 'ibrion': -1,
            'lorbit': 11, 'lwave': True, 'lcharg': True,
            'ldau': False, 'lmaxmix': 4,
        },
        'restart': None,
        'retrieve': ['OUTCAR'],
    }


@pytest.fixture
def hubbard_response_stage():
    """Valid hubbard_response stage pointing at gs."""
    return {
        'name': 'response',
        'type': 'hubbard_response',
        'ground_state_from': 'gs',
        'structure_from': 'input',
        'target_species': 'Ni',
        'potential_values': [-0.2, -0.1, 0.1, 0.2],
        'ldaul': 2,
    }


@pytest.fixture
def hubbard_analysis_stage():
    """Valid hubbard_analysis stage pointing at response."""
    return {
        'name': 'analysis',
        'type': 'hubbard_analysis',
        'response_from': 'response',
        'structure_from': 'input',
        'target_species': 'Ni',
        'ldaul': 2,
    }


# ---------------------------------------------------------------------------
# TestPortTypeRegistry
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

    def test_all_hubbard_response_output_types_recognized(self):
        for port_name, port in HUBBARD_RESPONSE_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"Hubbard response output '{port_name}' has unrecognized type '{port['type']}'"

    def test_all_hubbard_analysis_output_types_recognized(self):
        for port_name, port in HUBBARD_ANALYSIS_PORTS['outputs'].items():
            assert port['type'] in PORT_TYPES, \
                f"Hubbard analysis output '{port_name}' has unrecognized type '{port['type']}'"

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
# TestPortDeclarations
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
        assert 'structure' not in DOS_PORTS['outputs']

    def test_dos_has_nine_outputs(self):
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
        assert len(BADER_PORTS['outputs']) == 4
        for key in ('charges', 'acf', 'bcf', 'avf'):
            assert key in BADER_PORTS['outputs']

    def test_bader_both_inputs_have_compatible_bricks(self):
        for input_name, port in BADER_PORTS['inputs'].items():
            assert 'compatible_bricks' in port, \
                f"Bader input '{input_name}' missing compatible_bricks"

    def test_hubbard_response_has_no_structure_output(self):
        assert 'structure' not in HUBBARD_RESPONSE_PORTS['outputs']

    def test_hubbard_response_has_two_outputs(self):
        assert len(HUBBARD_RESPONSE_PORTS['outputs']) == 2
        assert 'responses' in HUBBARD_RESPONSE_PORTS['outputs']
        assert 'ground_state_occupation' in HUBBARD_RESPONSE_PORTS['outputs']

    def test_hubbard_response_requires_ground_state(self):
        inputs = HUBBARD_RESPONSE_PORTS['inputs']
        assert 'ground_state_remote' in inputs
        assert inputs['ground_state_remote']['required'] is True

    def test_hubbard_response_ground_state_has_prerequisites(self):
        gs_input = HUBBARD_RESPONSE_PORTS['inputs']['ground_state_remote']
        assert 'prerequisites' in gs_input
        prereqs = gs_input['prerequisites']
        assert 'incar' in prereqs
        assert 'retrieve' in prereqs

    def test_hubbard_response_ground_state_prereqs_require_lorbit(self):
        prereqs = HUBBARD_RESPONSE_PORTS['inputs']['ground_state_remote']['prerequisites']
        incar_prereqs = prereqs['incar']
        assert 'lorbit' in incar_prereqs
        assert incar_prereqs['lorbit'] == 11

    def test_hubbard_response_ground_state_prereqs_require_lwave(self):
        prereqs = HUBBARD_RESPONSE_PORTS['inputs']['ground_state_remote']['prerequisites']
        incar_prereqs = prereqs['incar']
        assert 'lwave' in incar_prereqs
        assert incar_prereqs['lwave'] is True

    def test_hubbard_response_ground_state_prereqs_require_lcharg(self):
        prereqs = HUBBARD_RESPONSE_PORTS['inputs']['ground_state_remote']['prerequisites']
        incar_prereqs = prereqs['incar']
        assert 'lcharg' in incar_prereqs
        assert incar_prereqs['lcharg'] is True

    def test_hubbard_response_ground_state_prereqs_require_outcar(self):
        prereqs = HUBBARD_RESPONSE_PORTS['inputs']['ground_state_remote']['prerequisites']
        assert 'OUTCAR' in prereqs['retrieve']

    def test_hubbard_analysis_has_no_structure_output(self):
        assert 'structure' not in HUBBARD_ANALYSIS_PORTS['outputs']

    def test_hubbard_analysis_has_two_outputs(self):
        assert len(HUBBARD_ANALYSIS_PORTS['outputs']) == 2
        assert 'summary' in HUBBARD_ANALYSIS_PORTS['outputs']
        assert 'hubbard_u_result' in HUBBARD_ANALYSIS_PORTS['outputs']

    def test_hubbard_analysis_requires_responses(self):
        inputs = HUBBARD_ANALYSIS_PORTS['inputs']
        assert 'responses' in inputs
        assert inputs['responses']['required'] is True

    def test_hubbard_analysis_responses_compatible_with_hubbard_response(self):
        responses_input = HUBBARD_ANALYSIS_PORTS['inputs']['responses']
        assert 'hubbard_response' in responses_input['compatible_bricks']

    def test_hubbard_analysis_requires_ground_state_occupation(self):
        inputs = HUBBARD_ANALYSIS_PORTS['inputs']
        assert 'ground_state_occupation' in inputs
        assert inputs['ground_state_occupation']['required'] is True


# ---------------------------------------------------------------------------
# TestConditionalEvaluation
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestConditionalEvaluation:
    """Test the conditional output evaluation mechanism."""

    def test_none_conditional_is_true(self):
        assert _evaluate_conditional(None, {}) is True

    def test_string_conditional_raises(self):
        with pytest.raises(ValueError, match="must be a dict"):
            _evaluate_conditional('nsw > 0', {'incar': {'nsw': 100}})

    def test_nsw_greater_than_zero_true(self):
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {'incar': {'nsw': 100}}) is True

    def test_nsw_greater_than_zero_false(self):
        cond = {'incar_key': 'nsw', 'operator': '>', 'value': 0}
        assert _evaluate_conditional(cond, {'incar': {'nsw': 0}}) is False

    def test_missing_incar_key_defaults_to_zero(self):
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
        assert warnings == []

    def test_vasp_then_dos_passes(self, relax_stage, dos_stage):
        warnings = validate_connections([relax_stage, dos_stage])
        assert warnings == []

    def test_vasp_then_batch_passes(self, relax_stage, batch_stage):
        warnings = validate_connections([relax_stage, batch_stage])
        assert warnings == []

    def test_full_pipeline_passes(self, relax_stage, scf_stage, dos_stage,
                                  batch_stage, bader_stage):
        stages = [relax_stage, scf_stage, dos_stage, batch_stage, bader_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)

    def test_hubbard_pipeline_passes(self, relax_stage, ground_state_stage,
                                     hubbard_response_stage,
                                     hubbard_analysis_stage):
        """Full Hubbard U pipeline: relax → gs → response → analysis."""
        gs = {**ground_state_stage, 'structure_from': 'relax'}
        stages = [relax_stage, gs, hubbard_response_stage,
                  hubbard_analysis_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)


# ---------------------------------------------------------------------------
# TestValidateConnectionsOutputExists
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsOutputExists:
    """Validate that source stages produce the required output type."""

    def test_dos_structure_from_vasp_passes(self, relax_stage, dos_stage):
        validate_connections([relax_stage, dos_stage])

    def test_dos_structure_from_dos_rejected(self, relax_stage, dos_stage):
        dos1 = {**dos_stage, 'name': 'dos1'}
        dos2 = {
            'name': 'dos2', 'type': 'dos', 'structure_from': 'dos1',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos1, dos2])

    def test_batch_structure_from_dos_rejected(self, relax_stage, dos_stage):
        batch = {
            'name': 'batch1', 'type': 'batch', 'structure_from': 'dos',
            'base_incar': {'encut': 520}, 'calculations': {'a': {}},
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos_stage, batch])

    def test_error_message_suggests_valid_stages(self, relax_stage, dos_stage):
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
        bader = {'name': 'bader', 'type': 'bader', 'charge_from': 'dos'}
        with pytest.raises(ValueError, match="compatible with bricks.*vasp"):
            validate_connections([relax_stage, dos_stage, bader])

    def test_bader_from_batch_rejected(self, relax_stage, batch_stage):
        bader = {'name': 'bader', 'type': 'bader', 'charge_from': 'charge_scan'}
        with pytest.raises(ValueError, match="compatible with bricks.*vasp"):
            validate_connections([relax_stage, batch_stage, bader])

    def test_hubbard_analysis_from_vasp_rejected(
        self, relax_stage, ground_state_stage
    ):
        """Analysis responses input must come from hubbard_response, not vasp."""
        gs = {**ground_state_stage, 'structure_from': 'relax'}
        analysis = {
            'name': 'analysis', 'type': 'hubbard_analysis',
            'response_from': 'gs',
            'structure_from': 'input',
            'target_species': 'Ni', 'ldaul': 2,
        }
        with pytest.raises(ValueError, match="doesn't produce it"):
            validate_connections([relax_stage, gs, analysis])


# ---------------------------------------------------------------------------
# TestValidateConnectionsPrerequisites
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsPrerequisites:
    """Test that prerequisite INCAR/retrieve requirements are enforced."""

    def test_bader_with_all_prereqs_passes(self, relax_stage, scf_stage,
                                           bader_stage):
        validate_connections([relax_stage, scf_stage, bader_stage])

    def test_bader_missing_laechg_rejected(self, relax_stage, bader_stage):
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
            'retrieve': ['OUTCAR'],
        }
        with pytest.raises(ValueError, match="Missing retrieve"):
            validate_connections([relax_stage, scf_bad, bader_stage])

    def test_bader_missing_both_incar_and_retrieve(self, relax_stage, bader_stage):
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
        scf_bad = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'laechg': True, 'lcharg': True},
            'restart': None,
            'retrieve': ['AECCAR0', 'OUTCAR'],
        }
        with pytest.raises(ValueError, match="Missing retrieve"):
            validate_connections([relax_stage, scf_bad, bader_stage])

    def test_hubbard_response_missing_lorbit_rejected(
        self, relax_stage, hubbard_response_stage
    ):
        """Ground state without lorbit=11 should fail prerequisites."""
        gs_bad = {
            'name': 'gs', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'lwave': True, 'lcharg': True},
            'restart': None,
            'retrieve': ['OUTCAR'],
        }
        with pytest.raises(ValueError, match="lorbit"):
            validate_connections([relax_stage, gs_bad, hubbard_response_stage])

    def test_hubbard_response_missing_lwave_rejected(
        self, relax_stage, hubbard_response_stage
    ):
        gs_bad = {
            'name': 'gs', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'lorbit': 11, 'lcharg': True},
            'restart': None,
            'retrieve': ['OUTCAR'],
        }
        with pytest.raises(ValueError, match="lwave"):
            validate_connections([relax_stage, gs_bad, hubbard_response_stage])

    def test_hubbard_response_missing_lcharg_rejected(
        self, relax_stage, hubbard_response_stage
    ):
        gs_bad = {
            'name': 'gs', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0, 'lorbit': 11, 'lwave': True},
            'restart': None,
            'retrieve': ['OUTCAR'],
        }
        with pytest.raises(ValueError, match="lcharg"):
            validate_connections([relax_stage, gs_bad, hubbard_response_stage])

    def test_hubbard_response_missing_outcar_retrieve_rejected(
        self, relax_stage, hubbard_response_stage
    ):
        gs_bad = {
            'name': 'gs', 'type': 'vasp',
            'incar': {
                'encut': 520, 'nsw': 0, 'lorbit': 11,
                'lwave': True, 'lcharg': True,
            },
            'restart': None,
            'retrieve': [],
        }
        with pytest.raises(ValueError, match="Missing retrieve"):
            validate_connections([relax_stage, gs_bad, hubbard_response_stage])

    def test_hubbard_response_with_all_prereqs_passes(
        self, relax_stage, ground_state_stage, hubbard_response_stage
    ):
        validate_connections([relax_stage, ground_state_stage,
                              hubbard_response_stage])


# ---------------------------------------------------------------------------
# TestValidateConnectionsConditional
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsConditional:
    """Test conditional output warnings."""

    def test_structure_from_relaxation_no_warning(self, relax_stage, dos_stage):
        warnings = validate_connections([relax_stage, dos_stage])
        assert len(warnings) == 0

    def test_structure_from_static_warns(self, relax_stage, scf_stage):
        third = {
            'name': 'rerelax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100, 'ibrion': 2},
            'restart': None,
        }
        warnings = validate_connections([relax_stage, scf_stage, third])
        assert len(warnings) == 1
        assert 'nsw=0' in warnings[0]

    def test_explicit_structure_from_static_warns(self, relax_stage, scf_stage):
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'scf',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        warnings = validate_connections([relax_stage, scf_stage, dos])
        assert len(warnings) == 1
        assert 'scf' in warnings[0]

    def test_first_stage_static_no_warning(self):
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        warnings = validate_connections([scf])
        assert len(warnings) == 0

    def test_structure_from_input_no_warning(self, relax_stage, scf_stage):
        third = {
            'name': 'fresh', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
            'structure_from': 'input',
        }
        warnings = validate_connections([relax_stage, scf_stage, third])
        assert len(warnings) == 0


# ---------------------------------------------------------------------------
# TestValidateConnectionsAutoResolution
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsAutoResolution:
    """Test VASP 'auto' structure resolution and its edge cases."""

    def test_first_stage_auto_always_valid(self):
        stage = {
            'name': 'relax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        warnings = validate_connections([stage])
        assert warnings == []

    def test_auto_previous_with_vasp_passes(self, relax_stage):
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        validate_connections([relax_stage, scf])

    def test_auto_previous_after_dos_fails(self, relax_stage, dos_stage):
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_auto_previous_after_batch_fails(self, relax_stage, batch_stage):
        vasp_after_batch = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, batch_stage, vasp_after_batch])

    def test_auto_previous_after_bader_fails(self, relax_stage, scf_stage,
                                             bader_stage):
        vasp_after_bader = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, scf_stage, bader_stage,
                                  vasp_after_bader])

    def test_auto_previous_after_hubbard_response_fails(
        self, relax_stage, ground_state_stage, hubbard_response_stage
    ):
        """Hubbard response has no structure output → auto(previous) should fail."""
        vasp_after_response = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, ground_state_stage,
                                  hubbard_response_stage, vasp_after_response])

    def test_auto_previous_after_hubbard_analysis_fails(
        self, relax_stage, ground_state_stage, hubbard_response_stage,
        hubbard_analysis_stage
    ):
        """Hubbard analysis has no structure output → auto(previous) should fail."""
        vasp_after_analysis = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, ground_state_stage,
                                  hubbard_response_stage,
                                  hubbard_analysis_stage,
                                  vasp_after_analysis])

    def test_explicit_structure_from_bypasses_previous(self, relax_stage,
                                                       dos_stage):
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'relax',
        }
        validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_structure_from_input_always_valid(self, relax_stage, dos_stage):
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'input',
        }
        validate_connections([relax_stage, dos_stage, vasp_after_dos])

    def test_structure_from_nonexistent_rejected(self, relax_stage):
        vasp = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
            'structure_from': 'nonexistent',
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, vasp])

    def test_error_suggests_stages_with_structure(self, relax_stage, dos_stage):
        vasp_after_dos = {
            'name': 'post', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': None,
        }
        with pytest.raises(ValueError) as exc_info:
            validate_connections([relax_stage, dos_stage, vasp_after_dos])
        assert 'relax' in str(exc_info.value)


# ---------------------------------------------------------------------------
# TestValidateConnectionsRestart
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsRestart:
    """Test restart connection validation."""

    def test_restart_none_passes(self, relax_stage):
        validate_connections([relax_stage])

    def test_restart_from_vasp_passes(self, relax_stage):
        scf = {
            'name': 'scf', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 0},
            'restart': 'relax',
        }
        validate_connections([relax_stage, scf])


# ---------------------------------------------------------------------------
# TestValidateConnectionsHubbard
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsHubbard:
    """Hubbard U specific connection validation tests."""

    def test_full_hubbard_pipeline(self, relax_stage, ground_state_stage,
                                   hubbard_response_stage,
                                   hubbard_analysis_stage):
        """Complete relax → gs → response → analysis pipeline."""
        stages = [relax_stage, ground_state_stage,
                  hubbard_response_stage, hubbard_analysis_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)

    def test_hubbard_response_structure_from_vasp_passes(
        self, relax_stage, ground_state_stage
    ):
        response = {
            'name': 'response', 'type': 'hubbard_response',
            'ground_state_from': 'gs', 'structure_from': 'relax',
            'target_species': 'Ni',
        }
        validate_connections([relax_stage, ground_state_stage, response])

    def test_hubbard_response_structure_from_input_passes(
        self, ground_state_stage
    ):
        response = {
            'name': 'response', 'type': 'hubbard_response',
            'ground_state_from': 'gs', 'structure_from': 'input',
            'target_species': 'Ni',
        }
        validate_connections([ground_state_stage, response])

    def test_hubbard_response_structure_from_hubbard_rejected(
        self, relax_stage, ground_state_stage, hubbard_response_stage
    ):
        """Can't get structure from a hubbard_response stage."""
        response2 = {
            'name': 'response2', 'type': 'hubbard_response',
            'ground_state_from': 'gs', 'structure_from': 'response',
            'target_species': 'Fe',
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, ground_state_stage,
                                  hubbard_response_stage, response2])

    def test_hubbard_analysis_from_wrong_brick_rejected(
        self, relax_stage, ground_state_stage
    ):
        """Analysis response_from must point to hubbard_response, not vasp."""
        analysis = {
            'name': 'analysis', 'type': 'hubbard_analysis',
            'response_from': 'gs', 'structure_from': 'input',
            'target_species': 'Ni', 'ldaul': 2,
        }
        with pytest.raises(ValueError, match="doesn't produce it"):
            validate_connections([relax_stage, ground_state_stage, analysis])

    def test_hubbard_analysis_structure_from_input_passes(
        self, ground_state_stage, hubbard_response_stage,
        hubbard_analysis_stage
    ):
        stages = [ground_state_stage, hubbard_response_stage,
                  hubbard_analysis_stage]
        validate_connections(stages)

    def test_hubbard_analysis_structure_from_analysis_rejected(
        self, relax_stage, ground_state_stage, hubbard_response_stage,
        hubbard_analysis_stage
    ):
        """Can't get structure from a hubbard_analysis stage."""
        analysis2 = {
            'name': 'analysis2', 'type': 'hubbard_analysis',
            'response_from': 'response', 'structure_from': 'analysis',
            'target_species': 'Ni', 'ldaul': 2,
        }
        with pytest.raises(ValueError, match="doesn't produce.*structure"):
            validate_connections([relax_stage, ground_state_stage,
                                  hubbard_response_stage,
                                  hubbard_analysis_stage, analysis2])


# ---------------------------------------------------------------------------
# TestValidateConnectionsBaderMultipleInputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBaderMultipleInputs:
    """Test that bader's two inputs from the same source are both validated."""

    def test_both_inputs_satisfied_by_vasp(self, relax_stage, scf_stage,
                                           bader_stage):
        validate_connections([relax_stage, scf_stage, bader_stage])

    def test_missing_charge_from_rejected(self, relax_stage):
        bader = {'name': 'bader', 'type': 'bader'}
        with pytest.raises(ValueError, match="charge_from.*missing"):
            validate_connections([relax_stage, bader])


# ---------------------------------------------------------------------------
# TestValidateConnectionsBatchOutputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsBatchOutputs:
    """Test that batch template outputs are handled correctly."""

    def test_batch_registers_outputs(self, relax_stage, batch_stage):
        validate_connections([relax_stage, batch_stage])

    def test_batch_has_no_structure_output(self, relax_stage, batch_stage):
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
        stages = [relax_stage, scf_stage, dos_stage, batch_stage, bader_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)

    def test_reversed_order_fails(self, relax_stage, bader_stage):
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([bader_stage, relax_stage])

    def test_circular_reference_impossible(self, relax_stage):
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'dos',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, dos])

    def test_forward_reference_fails(self, relax_stage, dos_stage):
        dos_forward = {**dos_stage, 'structure_from': 'late_relax'}
        late_relax = {
            'name': 'late_relax', 'type': 'vasp',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([relax_stage, dos_forward, late_relax])

    def test_full_hubbard_pipeline_passes(self, relax_stage, ground_state_stage,
                                          hubbard_response_stage,
                                          hubbard_analysis_stage):
        """Full pipeline: relax → gs → response → analysis."""
        stages = [relax_stage, ground_state_stage,
                  hubbard_response_stage, hubbard_analysis_stage]
        warnings = validate_connections(stages)
        assert isinstance(warnings, list)


# ---------------------------------------------------------------------------
# TestValidateConnectionsEdgeCases
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestValidateConnectionsEdgeCases:
    """Edge cases and corner scenarios."""

    def test_empty_stages_no_error(self):
        warnings = validate_connections([])
        assert warnings == []

    def test_single_dos_stage_fails(self):
        dos = {
            'name': 'dos', 'type': 'dos', 'structure_from': 'relax',
            'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
        }
        with pytest.raises(ValueError, match="unknown stage"):
            validate_connections([dos])

    def test_multiple_baders_from_different_sources(self, relax_stage):
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

    def test_many_dos_stages_from_same_source(self, relax_stage):
        stages = [relax_stage]
        for i in range(5):
            stages.append({
                'name': f'dos_{i}', 'type': 'dos', 'structure_from': 'relax',
                'scf_incar': {'encut': 520}, 'dos_incar': {'nedos': 3000},
            })
        validate_connections(stages)

    def test_default_type_is_vasp(self):
        stage = {
            'name': 'relax',
            'incar': {'encut': 520, 'nsw': 100},
            'restart': None,
        }
        validate_connections([stage])
