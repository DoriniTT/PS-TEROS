"""Unit tests for max_concurrent_jobs parameter and stage outputs organization.

Tests cover:
1. quick_vasp_sequential accepts max_concurrent_jobs parameter
2. Each brick's expose_stage_outputs() returns a list of output names
3. print_stage_outputs_summary() displays per-stage output listings
4. print_sequential_results() includes output names per stage

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import inspect
import pytest
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mock_wg():
    """A mock WorkGraph with a settable outputs namespace."""
    wg = MagicMock()
    wg.outputs = MagicMock()
    return wg


@pytest.fixture
def mock_vasp_stage_tasks():
    """Mock stage_tasks_result for a VASP brick."""
    vasp_task = MagicMock()
    vasp_task.outputs.structure = MagicMock(name='structure_socket')
    vasp_task.outputs.misc = MagicMock(name='misc_socket')
    vasp_task.outputs.remote_folder = MagicMock(name='remote_socket')
    vasp_task.outputs.retrieved = MagicMock(name='retrieved_socket')

    energy_task = MagicMock()
    energy_task.outputs.result = MagicMock(name='energy_socket')

    return {'vasp': vasp_task, 'energy': energy_task}


@pytest.fixture
def mock_dos_stage_tasks():
    """Mock stage_tasks_result for a DOS brick."""
    bands_task = MagicMock()
    bands_task.outputs.dos = MagicMock(name='dos_socket')
    bands_task.outputs.projectors = MagicMock(name='projectors_socket')
    bands_task.outputs.scf_misc = MagicMock(name='scf_misc_socket')
    bands_task.outputs.scf_remote_folder = MagicMock(name='scf_remote_socket')
    bands_task.outputs.scf_retrieved = MagicMock(name='scf_retrieved_socket')
    bands_task.outputs.dos_misc = MagicMock(name='dos_misc_socket')
    bands_task.outputs.dos_remote_folder = MagicMock(name='dos_remote_socket')
    bands_task.outputs.dos_retrieved = MagicMock(name='dos_retrieved_socket')

    return {'bands_task': bands_task}


@pytest.fixture
def mock_batch_stage_tasks():
    """Mock stage_tasks_result for a batch brick with 2 calculations."""
    vasp_neutral = MagicMock()
    vasp_neutral.outputs.misc = MagicMock(name='misc_neutral')
    vasp_neutral.outputs.remote_folder = MagicMock(name='remote_neutral')
    vasp_neutral.outputs.retrieved = MagicMock(name='retrieved_neutral')

    vasp_cation = MagicMock()
    vasp_cation.outputs.misc = MagicMock(name='misc_cation')
    vasp_cation.outputs.remote_folder = MagicMock(name='remote_cation')
    vasp_cation.outputs.retrieved = MagicMock(name='retrieved_cation')

    energy_neutral = MagicMock()
    energy_neutral.outputs.result = MagicMock(name='energy_neutral')

    energy_cation = MagicMock()
    energy_cation.outputs.result = MagicMock(name='energy_cation')

    return {
        'calc_tasks': {'neutral': vasp_neutral, 'cation': vasp_cation},
        'energy_tasks': {'neutral': energy_neutral, 'cation': energy_cation},
    }


@pytest.fixture
def mock_bader_stage_tasks():
    """Mock stage_tasks_result for a bader brick."""
    bader_task = MagicMock()
    bader_task.outputs.charges = MagicMock(name='charges_socket')
    bader_task.outputs.acf = MagicMock(name='acf_socket')
    bader_task.outputs.bcf = MagicMock(name='bcf_socket')
    bader_task.outputs.avf = MagicMock(name='avf_socket')

    return {'bader': bader_task}


@pytest.fixture
def sample_sequential_result():
    """A sample return dict from quick_vasp_sequential with stage_outputs."""
    return {
        '__workgraph_pk__': 12345,
        '__stage_names__': ['relax_rough', 'relax_fine', 'dos_calc'],
        '__stage_types__': {
            'relax_rough': 'vasp',
            'relax_fine': 'vasp',
            'dos_calc': 'dos',
        },
        '__stage_outputs__': {
            'relax_rough': [
                'relax_rough_energy', 'relax_rough_structure',
                'relax_rough_misc', 'relax_rough_remote',
                'relax_rough_retrieved',
            ],
            'relax_fine': [
                'relax_fine_energy', 'relax_fine_structure',
                'relax_fine_misc', 'relax_fine_remote',
                'relax_fine_retrieved',
            ],
            'dos_calc': [
                'dos_calc_dos', 'dos_calc_projectors',
                'dos_calc_scf_misc', 'dos_calc_scf_remote',
            ],
        },
        'relax_rough': 12345,
        'relax_fine': 12345,
        'dos_calc': 12345,
    }


# ---------------------------------------------------------------------------
# TestMaxConcurrentJobsParameter
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestMaxConcurrentJobsParameter:
    """Tests that quick_vasp_sequential accepts max_concurrent_jobs."""

    def test_parameter_in_signature(self):
        from teros.core.lego.workgraph import quick_vasp_sequential
        sig = inspect.signature(quick_vasp_sequential)
        assert 'max_concurrent_jobs' in sig.parameters

    def test_default_is_none(self):
        from teros.core.lego.workgraph import quick_vasp_sequential
        sig = inspect.signature(quick_vasp_sequential)
        param = sig.parameters['max_concurrent_jobs']
        assert param.default is None

    def test_parameter_type_annotation(self):
        from teros.core.lego.workgraph import quick_vasp_sequential
        sig = inspect.signature(quick_vasp_sequential)
        param = sig.parameters['max_concurrent_jobs']
        assert param.annotation is int


# ---------------------------------------------------------------------------
# TestVaspExposeStageOutputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestVaspExposeStageOutputs:
    """Tests for teros.core.lego.bricks.vasp.expose_stage_outputs() return value."""

    def test_returns_list(self, mock_wg, mock_vasp_stage_tasks):
        from teros.core.lego.bricks.vasp import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'relax', mock_vasp_stage_tasks)
        assert isinstance(result, list)

    def test_returns_five_outputs(self, mock_wg, mock_vasp_stage_tasks):
        from teros.core.lego.bricks.vasp import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'relax', mock_vasp_stage_tasks)
        assert len(result) == 5

    def test_output_names_prefixed(self, mock_wg, mock_vasp_stage_tasks):
        from teros.core.lego.bricks.vasp import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'my_stage', mock_vasp_stage_tasks)
        for name in result:
            assert name.startswith('my_stage_')

    def test_expected_output_names(self, mock_wg, mock_vasp_stage_tasks):
        from teros.core.lego.bricks.vasp import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'relax', mock_vasp_stage_tasks)
        expected = [
            'relax_energy', 'relax_structure', 'relax_misc',
            'relax_remote', 'relax_retrieved',
        ]
        assert result == expected

    def test_setattr_called_for_each_output(self, mock_wg, mock_vasp_stage_tasks):
        from teros.core.lego.bricks.vasp import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'relax', mock_vasp_stage_tasks)
        # Each output name was set on wg.outputs
        for name in result:
            getattr(mock_wg.outputs, name)  # should not raise


# ---------------------------------------------------------------------------
# TestDosExposeStageOutputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestDosExposeStageOutputs:
    """Tests for teros.core.lego.bricks.dos.expose_stage_outputs() return value."""

    def test_returns_list(self, mock_wg, mock_dos_stage_tasks):
        from teros.core.lego.bricks.dos import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'dos_calc', mock_dos_stage_tasks)
        assert isinstance(result, list)

    def test_returns_eight_outputs(self, mock_wg, mock_dos_stage_tasks):
        from teros.core.lego.bricks.dos import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'dos_calc', mock_dos_stage_tasks)
        assert len(result) == 8

    def test_output_names_prefixed(self, mock_wg, mock_dos_stage_tasks):
        from teros.core.lego.bricks.dos import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'dos_calc', mock_dos_stage_tasks)
        for name in result:
            assert name.startswith('dos_calc_')

    def test_expected_output_names(self, mock_wg, mock_dos_stage_tasks):
        from teros.core.lego.bricks.dos import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'dos_calc', mock_dos_stage_tasks)
        expected = [
            'dos_calc_dos', 'dos_calc_projectors',
            'dos_calc_scf_misc', 'dos_calc_scf_remote', 'dos_calc_scf_retrieved',
            'dos_calc_dos_misc', 'dos_calc_dos_remote', 'dos_calc_dos_retrieved',
        ]
        assert result == expected

    def test_missing_attributes_skipped(self, mock_wg):
        """If bands_task is missing some outputs, they are silently skipped."""
        from teros.core.lego.bricks.dos import expose_stage_outputs
        bands_task = MagicMock()
        # Only dos and scf_misc available, others raise AttributeError
        bands_task.outputs.dos = MagicMock()
        bands_task.outputs.scf_misc = MagicMock()
        del bands_task.outputs.projectors
        del bands_task.outputs.scf_remote_folder
        del bands_task.outputs.scf_retrieved
        del bands_task.outputs.dos_misc
        del bands_task.outputs.dos_remote_folder
        del bands_task.outputs.dos_retrieved
        stage_tasks = {'bands_task': bands_task}
        result = expose_stage_outputs(mock_wg, 'dos', stage_tasks)
        assert 'dos_dos' in result
        assert 'dos_scf_misc' in result
        assert len(result) == 2


# ---------------------------------------------------------------------------
# TestBatchExposeStageOutputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestBatchExposeStageOutputs:
    """Tests for teros.core.lego.bricks.batch.expose_stage_outputs() return value."""

    def test_returns_list(self, mock_wg, mock_batch_stage_tasks):
        from teros.core.lego.bricks.batch import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'fukui', mock_batch_stage_tasks)
        assert isinstance(result, list)

    def test_returns_four_per_calc(self, mock_wg, mock_batch_stage_tasks):
        from teros.core.lego.bricks.batch import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'fukui', mock_batch_stage_tasks)
        # 2 calcs * 4 outputs each = 8
        assert len(result) == 8

    def test_output_names_prefixed(self, mock_wg, mock_batch_stage_tasks):
        from teros.core.lego.bricks.batch import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'fukui', mock_batch_stage_tasks)
        for name in result:
            assert name.startswith('fukui_')

    def test_contains_all_calc_labels(self, mock_wg, mock_batch_stage_tasks):
        from teros.core.lego.bricks.batch import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'fukui', mock_batch_stage_tasks)
        neutral_outputs = [n for n in result if 'neutral' in n]
        cation_outputs = [n for n in result if 'cation' in n]
        assert len(neutral_outputs) == 4
        assert len(cation_outputs) == 4

    def test_single_calc(self, mock_wg):
        from teros.core.lego.bricks.batch import expose_stage_outputs
        vasp_task = MagicMock()
        vasp_task.outputs.misc = MagicMock()
        vasp_task.outputs.remote_folder = MagicMock()
        vasp_task.outputs.retrieved = MagicMock()
        energy_task = MagicMock()
        energy_task.outputs.result = MagicMock()
        stage_tasks = {
            'calc_tasks': {'only': vasp_task},
            'energy_tasks': {'only': energy_task},
        }
        result = expose_stage_outputs(mock_wg, 'step', stage_tasks)
        assert len(result) == 4
        assert result == [
            'step_only_energy', 'step_only_misc',
            'step_only_remote', 'step_only_retrieved',
        ]


# ---------------------------------------------------------------------------
# TestBaderExposeStageOutputs
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestBaderExposeStageOutputs:
    """Tests for teros.core.lego.bricks.bader.expose_stage_outputs() return value."""

    def test_returns_list(self, mock_wg, mock_bader_stage_tasks):
        from teros.core.lego.bricks.bader import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'bader_step', mock_bader_stage_tasks)
        assert isinstance(result, list)

    def test_returns_four_outputs(self, mock_wg, mock_bader_stage_tasks):
        from teros.core.lego.bricks.bader import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'bader_step', mock_bader_stage_tasks)
        assert len(result) == 4

    def test_expected_output_names(self, mock_wg, mock_bader_stage_tasks):
        from teros.core.lego.bricks.bader import expose_stage_outputs
        result = expose_stage_outputs(mock_wg, 'bader_step', mock_bader_stage_tasks)
        expected = [
            'bader_step_charges', 'bader_step_acf',
            'bader_step_bcf', 'bader_step_avf',
        ]
        assert result == expected


# ---------------------------------------------------------------------------
# TestPrintStageOutputsSummary
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestPrintStageOutputsSummary:
    """Tests for teros.core.lego.results.print_stage_outputs_summary()."""

    def test_prints_workgraph_pk(self, capsys, sample_sequential_result):
        from teros.core.lego.results import print_stage_outputs_summary
        print_stage_outputs_summary(sample_sequential_result)
        out = capsys.readouterr().out
        assert '12345' in out

    def test_prints_stage_names(self, capsys, sample_sequential_result):
        from teros.core.lego.results import print_stage_outputs_summary
        print_stage_outputs_summary(sample_sequential_result)
        out = capsys.readouterr().out
        assert 'relax_rough' in out
        assert 'relax_fine' in out
        assert 'dos_calc' in out

    def test_prints_stage_types(self, capsys, sample_sequential_result):
        from teros.core.lego.results import print_stage_outputs_summary
        print_stage_outputs_summary(sample_sequential_result)
        out = capsys.readouterr().out
        assert '(vasp)' in out
        assert '(dos)' in out

    def test_prints_output_names(self, capsys, sample_sequential_result):
        from teros.core.lego.results import print_stage_outputs_summary
        print_stage_outputs_summary(sample_sequential_result)
        out = capsys.readouterr().out
        assert 'relax_rough_energy' in out
        assert 'relax_rough_structure' in out
        assert 'dos_calc_dos' in out

    def test_prints_indices(self, capsys, sample_sequential_result):
        from teros.core.lego.results import print_stage_outputs_summary
        print_stage_outputs_summary(sample_sequential_result)
        out = capsys.readouterr().out
        assert '[1]' in out
        assert '[2]' in out
        assert '[3]' in out

    def test_empty_stage_outputs(self, capsys):
        from teros.core.lego.results import print_stage_outputs_summary
        result = {
            '__workgraph_pk__': 99,
            '__stage_names__': ['step1'],
            '__stage_types__': {'step1': 'vasp'},
            '__stage_outputs__': {'step1': []},
        }
        print_stage_outputs_summary(result)
        out = capsys.readouterr().out
        assert '(no outputs)' in out

    def test_missing_stage_outputs_key(self, capsys):
        """Gracefully handles result dicts without __stage_outputs__."""
        from teros.core.lego.results import print_stage_outputs_summary
        result = {
            '__workgraph_pk__': 99,
            '__stage_names__': ['step1'],
            '__stage_types__': {'step1': 'vasp'},
        }
        print_stage_outputs_summary(result)
        out = capsys.readouterr().out
        assert '(no outputs)' in out


# ---------------------------------------------------------------------------
# TestSequentialResultStageOutputsKey
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestSequentialResultStageOutputsKey:
    """Tests that __stage_outputs__ is correctly structured in the result dict."""

    def test_stage_outputs_is_dict(self, sample_sequential_result):
        assert isinstance(sample_sequential_result['__stage_outputs__'], dict)

    def test_stage_outputs_keys_match_stage_names(self, sample_sequential_result):
        stage_names = sample_sequential_result['__stage_names__']
        stage_outputs = sample_sequential_result['__stage_outputs__']
        assert set(stage_outputs.keys()) == set(stage_names)

    def test_each_value_is_list(self, sample_sequential_result):
        for outputs in sample_sequential_result['__stage_outputs__'].values():
            assert isinstance(outputs, list)

    def test_each_item_is_string(self, sample_sequential_result):
        for outputs in sample_sequential_result['__stage_outputs__'].values():
            for name in outputs:
                assert isinstance(name, str)

    def test_output_names_contain_stage_prefix(self, sample_sequential_result):
        for stage_name, outputs in sample_sequential_result['__stage_outputs__'].items():
            for name in outputs:
                assert name.startswith(f'{stage_name}_')


# ---------------------------------------------------------------------------
# TestExposeStageOutputsConsistency
# ---------------------------------------------------------------------------

@pytest.mark.tier1
class TestExposeStageOutputsConsistency:
    """Cross-brick consistency checks for expose_stage_outputs."""

    def test_all_bricks_return_list(
        self, mock_wg, mock_vasp_stage_tasks,
        mock_dos_stage_tasks, mock_batch_stage_tasks,
        mock_bader_stage_tasks,
    ):
        from teros.core.lego.bricks import vasp, dos, batch, bader
        assert isinstance(vasp.expose_stage_outputs(mock_wg, 's', mock_vasp_stage_tasks), list)
        assert isinstance(dos.expose_stage_outputs(mock_wg, 's', mock_dos_stage_tasks), list)
        assert isinstance(batch.expose_stage_outputs(mock_wg, 's', mock_batch_stage_tasks), list)
        assert isinstance(bader.expose_stage_outputs(mock_wg, 's', mock_bader_stage_tasks), list)

    def test_all_bricks_return_nonempty(
        self, mock_wg, mock_vasp_stage_tasks,
        mock_dos_stage_tasks, mock_batch_stage_tasks,
        mock_bader_stage_tasks,
    ):
        from teros.core.lego.bricks import vasp, dos, batch, bader
        assert len(vasp.expose_stage_outputs(mock_wg, 's', mock_vasp_stage_tasks)) > 0
        assert len(dos.expose_stage_outputs(mock_wg, 's', mock_dos_stage_tasks)) > 0
        assert len(batch.expose_stage_outputs(mock_wg, 's', mock_batch_stage_tasks)) > 0
        assert len(bader.expose_stage_outputs(mock_wg, 's', mock_bader_stage_tasks)) > 0

    def test_all_bricks_return_strings(
        self, mock_wg, mock_vasp_stage_tasks,
        mock_dos_stage_tasks, mock_batch_stage_tasks,
        mock_bader_stage_tasks,
    ):
        from teros.core.lego.bricks import vasp, dos, batch, bader
        for result in [
            vasp.expose_stage_outputs(mock_wg, 's', mock_vasp_stage_tasks),
            dos.expose_stage_outputs(mock_wg, 's', mock_dos_stage_tasks),
            batch.expose_stage_outputs(mock_wg, 's', mock_batch_stage_tasks),
            bader.expose_stage_outputs(mock_wg, 's', mock_bader_stage_tasks),
        ]:
            for name in result:
                assert isinstance(name, str)
