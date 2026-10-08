"""The aiida-vasp work chain that runs the VASP tasks of psteros.

psteros recipes write the INCAR in the usual upper case (``{"ENCUT": 520}``)
and ``vasp_parameters`` keeps it that way.  aiida-vasp 5 lower-cases the tags
itself, but its input validator refuses upper-case keys on a *stored* ``Dict``,
which is what a WorkGraph hands to a task: the graph then fails at once with
"Case inconsistency found in the parameters dictionary".

:class:`PsterosVaspWorkChain` is ``VaspWorkChain`` with one change: the
``parameters`` validator lower-cases a copy of the keys before applying the
validation of aiida-vasp, so the tags, the namespaces and the values are
checked exactly as before.  The inputs, the outputs, the exit codes and the
handlers are those of ``VaspWorkChain``.

This module imports aiida-vasp at import time, so psteros loads it only while
building a VASP graph; the daemon imports it when it runs one, which is why
psteros must be installed in the environment of the AiiDA daemon.
"""

from __future__ import annotations

from typing import Any

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_vasp.common import parameters_validator
from aiida_vasp.utils.aiida_utils import convert_dict_case

VaspWorkChain = WorkflowFactory("vasp.v2.vasp")


def case_tolerant_parameters_validator(node: orm.Dict | None, port: Any = None) -> str | None:
    """The validator of aiida-vasp applied to a copy of ``node`` with lower-case keys."""

    if not node:
        return None
    return parameters_validator(orm.Dict(convert_dict_case(node.get_dict(), lower=True)), port)


class PsterosVaspWorkChain(VaspWorkChain):
    """``VaspWorkChain`` that accepts the upper-case INCAR tags of the psteros recipes."""

    @classmethod
    def define(cls, spec: Any) -> None:
        super().define(spec)
        spec.inputs["parameters"].validator = case_tolerant_parameters_validator
