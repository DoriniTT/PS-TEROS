"""Unit tests for the lego brick connection / port validation system.

Tests the PORTS declarations, validate_connections(), port type registry,
prerequisites checking, conditional output warnings, and the full
validation pipeline. These tests define the expected behavior from the
PLAN_CONNECTIONS.md design document and serve as a specification that the
implementation must satisfy.

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import warnings

import pytest


# ---------------------------------------------------------------------------
# Helpers – minimal PORTS used across many tests
# ---------------------------------------------------------------------------

VASP_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',
            'required': True,
            'source': 'auto',
            'description': 'Atomic structure',
        },
        'restart_folder': {
            'type': 'remote_folder',
            'required': False,
            'source': 'restart',
            'description': 'Remote folder for WAVECAR/CHGCAR restart',
        },
    },
    'outputs': {
        'structure': {
            'type': 'structure',
            'conditional': {'incar_key': 'nsw', 'operator': '>', 'value': 0},
            'description': 'Relaxed structure (only meaningful if nsw > 0)',
        },
        'energy': {
            'type': 'energy',
            'description': 'Total energy (eV)',
        },
        'misc': {
            'type': 'misc',
            'description': 'Parsed VASP results dict',
        },
        'remote_folder': {
            'type': 'remote_folder',
            'description': 'Remote calculation directory',
        },
        'retrieved': {
            'type': 'retrieved',
            'description': 'Retrieved files from cluster',
        },
    },
}

DOS_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',
            'required': True,
            'source': 'structure_from',
            'description': 'Structure to compute DOS for',
        },
    },
    'outputs': {
        'energy': {
            'type': 'energy',
            'description': 'SCF energy',
        },
        'scf_misc': {
            'type': 'misc',
            'description': 'SCF parsed results',
        },
        'dos_misc': {
            'type': 'misc',
            'description': 'DOS parsed results',
        },
        'dos': {
            'type': 'dos_data',
            'description': 'DOS ArrayData',
        },
        'projectors': {
            'type': 'projectors',
            'description': 'Projected DOS data',
        },
        'scf_remote': {
            'type': 'remote_folder',
            'description': 'SCF remote folder',
        },
        'scf_retrieved': {
            'type': 'retrieved',
            'description': 'SCF retrieved files',
        },
        'dos_remote': {
            'type': 'remote_folder',
            'description': 'DOS remote folder',
        },
        'dos_retrieved': {
            'type': 'retrieved',
            'description': 'DOS retrieved files',
        },
        # NOTE: no 'structure' output
    },
}

BATCH_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',
            'required': True,
            'source': 'structure_from',
            'description': 'Structure for all sub-calculations',
        },
    },
    'outputs': {
        '{label}_energy': {
            'type': 'energy',
            'description': 'Energy per sub-calculation',
            'per_calculation': True,
        },
        '{label}_misc': {
            'type': 'misc',
            'description': 'Parsed results per sub-calculation',
            'per_calculation': True,
        },
        '{label}_remote_folder': {
            'type': 'remote_folder',
            'description': 'Remote folder per sub-calculation',
            'per_calculation': True,
        },
        '{label}_retrieved': {
            'type': 'retrieved',
            'description': 'Retrieved files per sub-calculation',
            'per_calculation': True,
        },
        # NOTE: no 'structure' output
    },
}

BADER_PORTS = {
    'inputs': {
        'charge_files': {
            'type': 'retrieved',
            'required': True,
            'source': 'charge_from',
            'compatible_bricks': ['vasp'],
            'prerequisites': {
                'incar': {'laechg': True, 'lcharg': True},
                'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR', 'OUTCAR'],
            },
            'description': 'Retrieved folder with AECCAR0, AECCAR2, CHGCAR',
        },
        'structure': {
            'type': 'structure',
            'required': True,
            'source': 'charge_from',
            'compatible_bricks': ['vasp'],
            'description': 'Structure for element-to-atom mapping',
        },
    },
    'outputs': {
        'charges': {
            'type': 'bader_charges',
            'description': 'Per-atom Bader charges (Dict)',
        },
        'acf': {
            'type': 'file',
            'description': 'ACF.dat (SinglefileData)',
        },
        'bcf': {
            'type': 'file',
            'description': 'BCF.dat (SinglefileData)',
        },
        'avf': {
            'type': 'file',
            'description': 'AVF.dat (SinglefileData)',
        },
    },
}

CONVERGENCE_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',
            'required': False,
            'source': 'structure_from',
            'description': 'Structure (optional, falls back to initial)',
        },
    },
    'outputs': {
        'cutoff_analysis': {
            'type': 'convergence',
            'description': 'ENCUT convergence analysis',
        },
        'kpoints_analysis': {
            'type': 'convergence',
            'description': 'K-points convergence analysis',
        },
        'recommendations': {
            'type': 'convergence',
            'description': 'Recommended parameters',
        },
        # NOTE: no 'structure' output
    },
}

# Registry mapping brick type name → PORTS dict
ALL_PORTS = {
    'vasp': VASP_PORTS,
    'dos': DOS_PORTS,
    'batch': BATCH_PORTS,
    'bader': BADER_PORTS,
    'convergence': CONVERGENCE_PORTS,
}

# Recognized port types (PORT_TYPES registry)
PORT_TYPES = {
    'structure', 'energy', 'misc', 'remote_folder', 'retrieved',
    'dos_data', 'projectors', 'bader_charges', 'trajectory',
    'convergence', 'file',
}


# ---------------------------------------------------------------------------
# Validation functions under test
#
# These implement the validate_connections() logic from the plan.
# When the real implementation lands, these imports should be replaced
# with the actual module. For now, we implement the spec here so that
# the tests are self-contained and runnable immediately.
# ---------------------------------------------------------------------------

def _get_ports(brick_type):
    """Look up PORTS dict for a brick type."""
    if brick_type not in ALL_PORTS:
        raise ValueError(f"Unknown brick type '{brick_type}'")
    return ALL_PORTS[brick_type]


def _validate_port_types(ports, brick_name):
    """Validate that all port type strings are recognized (issue #16)."""
    for section in ('inputs', 'outputs'):
        for port_name, port in ports.get(section, {}).items():
            ptype = port['type']
            if ptype not in PORT_TYPES:
                raise ValueError(
                    f"Unknown port type '{ptype}' in {brick_name}.{section}.{port_name}"
                )


def _evaluate_conditional(conditional, stage_config):
    """Evaluate a conditional dict against a stage's config (issue #9).

    Args:
        conditional: dict with 'incar_key', 'operator', 'value'
        stage_config: the raw stage dict

    Returns:
        True if the condition is met (output is available), False otherwise.
    """
    if conditional is None:
        return True
    if isinstance(conditional, str):
        raise ValueError(
            f"Conditional must be a dict, not a string: '{conditional}'. "
            f"Use {{'incar_key': 'nsw', 'operator': '>', 'value': 0}} instead."
        )
    incar = stage_config.get('incar', {})
    key = conditional['incar_key']
    op = conditional['operator']
    threshold = conditional['value']
    actual = incar.get(key, 0)  # VASP defaults most to 0
    if op == '>':
        return actual > threshold
    elif op == '>=':
        return actual >= threshold
    elif op == '==':
        return actual == threshold
    elif op == '!=':
        return actual != threshold
    elif op == '<':
        return actual < threshold
    elif op == '<=':
        return actual <= threshold
    else:
        raise ValueError(f"Unknown operator '{op}' in conditional")


def validate_connections(stages):
    """Validate all inter-stage connections before submission.

    This is the core function from the plan. It checks:
    - Port type compatibility between source and destination
    - Brick compatibility constraints (compatible_bricks)
    - Prerequisites (INCAR settings, retrieve lists)
    - Conditional output warnings
    - Auto-resolution for VASP structure input (issue #2, #12)

    Args:
        stages: list of stage config dicts

    Returns:
        list of warning strings (empty if no warnings)

    Raises:
        ValueError: if any connection is invalid
    """
    available_outputs = {}   # stage_name -> {port_name: port_type}
    stage_configs = {}       # stage_name -> raw stage dict (issue #11)
    stage_types = {}         # stage_name -> brick type string (issue #17)
    output_conditionals = {} # stage_name -> {port_name: conditional_dict}
    warn_list = []

    for i, stage in enumerate(stages):
        name = stage['name']
        brick_type = stage.get('type', 'vasp')
        ports = _get_ports(brick_type)

        stage_configs[name] = stage
        stage_types[name] = brick_type

        # Validate port type strings are recognized (issue #16)
        _validate_port_types(ports, brick_type)

        # Check every required input can be satisfied
        for input_name, input_port in ports['inputs'].items():
            if not input_port.get('required', True) is False:
                # Port is required (default True)
                pass
            else:
                # Optional port — skip if source field not present
                source_key = input_port['source']
                if stage.get(source_key) is None:
                    continue

            source_key = input_port['source']

            # ── Handle 'auto' source (VASP structure) ── issue #2, #12
            if source_key == 'auto':
                structure_from = stage.get('structure_from', 'previous')

                if i == 0:
                    # First stage: always uses initial structure → valid
                    pass
                elif structure_from == 'input':
                    # Explicit initial structure → valid (issue #12d)
                    pass
                elif structure_from == 'previous':
                    # Check previous stage has a structure output (#12b)
                    prev_name = stages[i - 1]['name']
                    if prev_name in available_outputs:
                        prev_types = {
                            pname: ptype
                            for pname, ptype in available_outputs[prev_name].items()
                        }
                        has_structure = any(
                            t == 'structure' for t in prev_types.values()
                        )
                        if not has_structure:
                            # Collect stages that DO have structure
                            producers = [
                                sn for sn, outs in available_outputs.items()
                                if any(t == 'structure' for t in outs.values())
                            ]
                            raise ValueError(
                                f"Stage '{name}': previous stage '{prev_name}' "
                                f"(type: {stage_types[prev_name]}) doesn't produce "
                                f"a 'structure' output. Use 'structure_from' to "
                                f"reference a stage that does. "
                                f"Stages with structure: {producers}"
                            )
                        # Check conditional warning
                        prev_conds = output_conditionals.get(prev_name, {})
                        if 'structure' in prev_conds:
                            cond = prev_conds['structure']
                            if not _evaluate_conditional(cond, stage_configs[prev_name]):
                                warn_list.append(
                                    f"Warning: Stage '{prev_name}' has "
                                    f"{cond['incar_key']}="
                                    f"{stage_configs[prev_name].get('incar', {}).get(cond['incar_key'], 0)} "
                                    f"(static calculation). Its 'structure' output "
                                    f"may not be meaningful."
                                )
                else:
                    # structure_from is an explicit stage name (#12c)
                    if structure_from not in available_outputs:
                        raise ValueError(
                            f"Stage '{name}': structure_from='{structure_from}' "
                            f"references unknown stage"
                        )
                    ref_outputs = available_outputs[structure_from]
                    if not any(t == 'structure' for t in ref_outputs.values()):
                        producers = [
                            sn for sn, outs in available_outputs.items()
                            if any(t == 'structure' for t in outs.values())
                        ]
                        raise ValueError(
                            f"Stage '{name}': structure_from='{structure_from}' "
                            f"references stage '{structure_from}' "
                            f"(type: {stage_types[structure_from]}), which doesn't "
                            f"produce a 'structure' output. "
                            f"Stages with structure: {producers}"
                        )
                    # Check conditional warning
                    ref_conds = output_conditionals.get(structure_from, {})
                    if 'structure' in ref_conds:
                        cond = ref_conds['structure']
                        if not _evaluate_conditional(cond, stage_configs[structure_from]):
                            warn_list.append(
                                f"Warning: Stage '{structure_from}' has "
                                f"{cond['incar_key']}="
                                f"{stage_configs[structure_from].get('incar', {}).get(cond['incar_key'], 0)} "
                                f"(static calculation). Its 'structure' output "
                                f"may not be meaningful."
                            )
                continue  # auto handling done

            # ── Handle explicit source fields ──
            ref_stage_name = stage.get(source_key)
            if ref_stage_name is None:
                if input_port.get('required', True):
                    raise ValueError(
                        f"Stage '{name}': input '{input_name}' requires "
                        f"'{source_key}' field but it's missing"
                    )
                continue

            # Check referenced stage exists
            if ref_stage_name not in available_outputs:
                raise ValueError(
                    f"Stage '{name}': '{source_key}={ref_stage_name}' "
                    f"references unknown stage"
                )

            # Check type compatibility
            ref_outputs = available_outputs[ref_stage_name]
            matching = [
                pname for pname, ptype in ref_outputs.items()
                if ptype == input_port['type']
            ]
            if not matching:
                # Build helpful suggestion
                producers = [
                    sn for sn, outs in available_outputs.items()
                    if any(t == input_port['type'] for t in outs.values())
                ]
                raise ValueError(
                    f"Stage '{name}': input '{input_name}' needs type "
                    f"'{input_port['type']}' but stage '{ref_stage_name}' "
                    f"(type: {stage_types[ref_stage_name]}) doesn't produce it. "
                    f"Stages with '{input_port['type']}': {producers}"
                )

            # Check brick compatibility constraint
            if 'compatible_bricks' in input_port:
                ref_brick_type = stage_types[ref_stage_name]  # issue #17
                if ref_brick_type not in input_port['compatible_bricks']:
                    raise ValueError(
                        f"Stage '{name}': input '{input_name}' is only "
                        f"compatible with bricks: {input_port['compatible_bricks']}, "
                        f"but '{ref_stage_name}' is type '{ref_brick_type}'"
                    )

            # Check prerequisites (issue #7, #11)
            prereqs = input_port.get('prerequisites')
            if prereqs:
                ref_config = stage_configs[ref_stage_name]
                ref_incar = ref_config.get('incar', {})
                ref_retrieve = ref_config.get('retrieve', [])

                missing_incar = {}
                for key, required_val in prereqs.get('incar', {}).items():
                    if ref_incar.get(key) != required_val:
                        missing_incar[key] = required_val

                missing_retrieve = [
                    f for f in prereqs.get('retrieve', [])
                    if f not in ref_retrieve
                ]

                if missing_incar or missing_retrieve:
                    msg = (
                        f"Stage '{name}' connects to stage '{ref_stage_name}' "
                        f"via '{source_key}', but '{ref_stage_name}' is missing "
                        f"required settings:"
                    )
                    if missing_incar:
                        items = ', '.join(
                            f"{k}={v}" for k, v in missing_incar.items()
                        )
                        msg += f"\n  Missing INCAR: {items}"
                    if missing_retrieve:
                        msg += f"\n  Missing retrieve: {', '.join(missing_retrieve)}"
                    raise ValueError(msg)

            # Check conditional warning on the referenced output
            ref_conds = output_conditionals.get(ref_stage_name, {})
            for matched_port in matching:
                if matched_port in ref_conds:
                    cond = ref_conds[matched_port]
                    if not _evaluate_conditional(cond, stage_configs[ref_stage_name]):
                        # Issue #13: suppress for bricks that handle it
                        if not input_port.get('handles_conditional', False):
                            warn_list.append(
                                f"Warning: Stage '{ref_stage_name}' has "
                                f"{cond['incar_key']}="
                                f"{stage_configs[ref_stage_name].get('incar', {}).get(cond['incar_key'], 0)} "
                                f"(static calculation). Its '{matched_port}' output "
                                f"may not be meaningful."
                            )

        # ── Register this stage's outputs ──
        stage_outputs = {}
        stage_conds = {}
        for port_name, port in ports['outputs'].items():
            # Issue #8: batch template outputs
            if port.get('per_calculation'):
                # Register the base type as available (not the template key)
                stage_outputs[port_name] = port['type']
            else:
                stage_outputs[port_name] = port['type']
            # Track conditionals
            if 'conditional' in port:
                stage_conds[port_name] = port['conditional']

        available_outputs[name] = stage_outputs
        output_conditionals[name] = stage_conds

    return warn_list


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
