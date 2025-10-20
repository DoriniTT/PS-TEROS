"""Tests for VASP builder input construction in adsorption_energy module."""

import pytest
from aiida import orm, load_profile
from teros.core.adsorption_energy import _build_vasp_inputs


@pytest.fixture(scope="module", autouse=True)
def load_aiida_profile():
    """Load AiiDA profile for all tests."""
    load_profile('psteros')


@pytest.fixture
def create_test_structure_and_code():
    """Create reusable H2 structure and VASP code for tests."""
    from ase import Atoms
    ase_struct = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    structure = orm.StructureData(ase=ase_struct)
    code = orm.load_code('VASP-VTST-6.4.3@bohr')
    return structure, code


def test_build_vasp_inputs_from_builder_dict(create_test_structure_and_code):
    """Test building inputs from new-style builder_inputs dict."""
    structure, code = create_test_structure_and_code
    # Builder inputs
    builder_inputs = {
        'parameters': {'incar': {'PREC': 'Accurate', 'ENCUT': 520}},
        'options': {'resources': {'num_machines': 1}},
        'potential_family': 'PBE',
        'potential_mapping': {'H': 'H'},
        'kpoints_spacing': 0.3,
    }

    # Build inputs
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=builder_inputs,
    )

    # Verify structure and code are set
    assert result['structure'] == structure
    assert result['code'] == code

    # Verify builder_inputs are preserved
    assert result['parameters']['incar']['PREC'] == 'Accurate'
    assert result['parameters']['incar']['ENCUT'] == 520
    assert result['options']['resources']['num_machines'] == 1
    assert result['potential_family'] == 'PBE'
    assert result['potential_mapping']['H'] == 'H'
    assert result['kpoints_spacing'] == 0.3


def test_build_vasp_inputs_from_parameters_dict(create_test_structure_and_code):
    """Test building inputs from old-style parameters dict (backward compat)."""
    structure, code = create_test_structure_and_code

    # Old-style parameters
    parameters = {'PREC': 'Accurate', 'ENCUT': 520, 'NSW': 100}
    options = {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40}}
    potential_mapping = {'H': 'H'}

    # Build inputs
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=None,
        parameters=parameters,
        options=options,
        potential_family='PBE',
        potential_mapping=potential_mapping,
        kpoints_spacing=0.3,
        clean_workdir=True,
    )

    # Verify structure and code are set
    assert result['structure'] == structure
    assert result['code'] == code

    # Verify parameters are wrapped in 'incar' dict
    assert result['parameters']['incar']['PREC'] == 'Accurate'
    assert result['parameters']['incar']['ENCUT'] == 520
    assert result['parameters']['incar']['NSW'] == 100

    # Verify options
    assert result['options']['resources']['num_machines'] == 1
    assert result['options']['resources']['num_mpiprocs_per_machine'] == 40

    # Verify other fields
    assert result['potential_family'] == 'PBE'
    assert result['potential_mapping']['H'] == 'H'
    assert result['kpoints_spacing'] == 0.3
    assert result['clean_workdir'] == True

    # Verify settings dict is created
    assert 'settings' in result
    assert isinstance(result['settings'], orm.Dict)


def test_build_vasp_inputs_force_scf(create_test_structure_and_code):
    """Test that force_scf=True enforces NSW=0 and IBRION=-1."""
    structure, code = create_test_structure_and_code

    # Builder inputs with relaxation settings
    builder_inputs = {
        'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ENCUT': 520}},
        'options': {'resources': {'num_machines': 1}},
        'potential_family': 'PBE',
        'potential_mapping': {'H': 'H'},
    }

    # Build inputs with force_scf=True
    result = _build_vasp_inputs(
        structure=structure,
        code=code,
        builder_inputs=builder_inputs,
        force_scf=True,
    )

    # Verify NSW=0 and IBRION=-1 are enforced
    assert result['parameters']['incar']['NSW'] == 0
    assert result['parameters']['incar']['IBRION'] == -1

    # Verify other parameters are preserved
    assert result['parameters']['incar']['ENCUT'] == 520


def test_build_vasp_inputs_raises_on_missing_inputs(create_test_structure_and_code):
    """Test that ValueError is raised when neither builder_inputs nor parameters are provided."""
    structure, code = create_test_structure_and_code

    # Call without builder_inputs or parameters
    with pytest.raises(ValueError) as exc_info:
        _build_vasp_inputs(
            structure=structure,
            code=code,
            builder_inputs=None,
            parameters=None,
        )

    # Verify error message is helpful
    error_message = str(exc_info.value)
    assert "Must provide either 'builder_inputs' or 'parameters'" in error_message
    assert "builder_inputs for full control" in error_message
    assert "backward compatibility" in error_message
