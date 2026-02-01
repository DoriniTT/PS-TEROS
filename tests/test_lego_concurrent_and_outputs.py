"""Unit tests for max_concurrent_jobs parameter in quick_vasp_sequential.

All tests are tier1 (pure Python, no AiiDA profile needed).
"""

import inspect
import pytest


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
