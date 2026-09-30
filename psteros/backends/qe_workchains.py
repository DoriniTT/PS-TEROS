"""Quantum ESPRESSO work chains that psteros graphs run on the AiiDA daemon.

This module imports aiida-quantumespresso at import time, so psteros loads it
only while building a graph.  The daemon imports it by module path when it
runs such a graph, which is why psteros must be installed in the daemon's
Python environment.
"""

from __future__ import annotations

from aiida.engine import ProcessHandlerReport, process_handler
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain


class PwRelaxStageWorkChain(PwBaseWorkChain):
    """``PwBaseWorkChain`` for the relaxation stage of a relaxation-to-static graph.

    A ``vc-relax`` can end with exit status 501: the ionic cycle converged, but
    the final SCF, recomputed with the plane-wave basis of the new cell,
    exceeds the force or stress thresholds.  aiida-quantumespresso keeps that
    structure as final yet still reports 501, so aiida-workgraph marks the task
    failed and skips the static SCF that depends on it.  Here the relaxed
    structure is a successful outcome, because the following static stage
    recomputes the energy on exactly that geometry.
    """

    @process_handler(
        priority=570,
        exit_codes=[PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF],
    )
    def handle_vcrelax_converged_except_final_scf(self, calculation):
        """Accept the relaxed structure; the static stage evaluates its energy."""

        self.ctx.is_finished = True
        self.report_error_handled(
            calculation, "accept the relaxed structure; the static stage recomputes its energy."
        )
        return ProcessHandlerReport(True)
