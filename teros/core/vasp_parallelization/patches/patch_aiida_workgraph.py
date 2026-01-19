"""Patch for aiida-workgraph error_handler_manager.py.

This patch adds the AcceptError exception and _accept_task_as_finished method
to allow error handlers to signal that a "failed" task should be accepted as
successful, enabling child tasks to run.

Usage:
    from teros.core.vasp_parallelization.patches.patch_aiida_workgraph import (
        apply_patch,
        verify_patch,
    )
    apply_patch()
    verify_patch()  # Raises if patch not applied
"""

import importlib.util
from pathlib import Path


def get_error_handler_manager_path() -> Path:
    """Get the path to the error_handler_manager.py file."""
    spec = importlib.util.find_spec("aiida_workgraph.engine.error_handler_manager")
    if spec is None or spec.origin is None:
        raise ImportError("Could not find aiida_workgraph.engine.error_handler_manager")
    return Path(spec.origin)


# The complete patched content for error_handler_manager.py
PATCHED_CONTENT = '''from __future__ import annotations
import traceback
from node_graph.error_handler import ErrorHandlerSpec


class AcceptError(Exception):
    """Raise this exception in an error handler to accept a failed task as finished.

    When an error handler raises AcceptError, the task will be marked as FINISHED
    instead of being reset for retry. This is useful when a task's "failure" is
    actually an expected outcome (e.g., benchmarking calculations that intentionally
    don't converge).

    Example usage in an error handler:
        def my_error_handler(task, engine):
            from aiida_workgraph.engine.error_handler_manager import AcceptError
            # Accept the task as-is without retry
            raise AcceptError("Task accepted despite non-zero exit code")
    """

    pass


class ErrorHandlerManager:
    def __init__(self, process, ctx_manager, logger):
        self.process = process
        self.ctx_manager = ctx_manager
        self.ctx = ctx_manager.ctx
        self.logger = logger

    def run_error_handlers(self, task_name: str) -> None:
        """Run error handlers for a task."""

        self.process.report(f"Run error handlers for {task_name}")

        node = self.process.task_manager.state_manager.get_task_runtime_info(
            task_name, "process"
        )
        if not node or not node.exit_status:
            return
        # error_handlers from the task
        error_handlers = self.process.wg.tasks[task_name].get_error_handlers()
        for data in error_handlers.values():
            if node.exit_status in data.exit_codes:
                self.run_error_handler(data, task_name)
                return
        # error_handlers from the workgraph
        for data in self.process.wg._error_handlers.values():
            if node.exit_code.status in data["tasks"].get(task_name, {}).get(
                "exit_codes", []
            ):
                self.run_error_handler(data, task_name)
                return

    def run_error_handler(self, handler: ErrorHandlerSpec, task_name: str) -> None:
        """Run the error handler for a task."""
        from inspect import signature
        from node_graph.executor import RuntimeExecutor

        executor = RuntimeExecutor(**handler.executor.to_dict()).callable
        executor_sig = signature(executor)
        self.process.report(f"Run error handler: {executor.__name__}")
        if handler.retry < handler.max_retries:
            task = self.process.task_manager.get_task(task_name)
            try:
                # Run the error handler to update the inputs of the task
                if "engine" in executor_sig.parameters:
                    msg = executor(task, engine=self, **(handler.kwargs or {}))
                else:
                    msg = executor(task, **(handler.kwargs or {}))
                # Reset the task to rerun it
                self.process.task_manager.state_manager.reset_task(task.name)
                # Save the updated task into self.ctx._wgdata
                tdata = task.to_dict()
                self.ctx._wgdata["tasks"][task.name] = tdata
                if msg:
                    self.process.report(msg)
                handler.retry += 1
            except AcceptError as e:
                # AcceptError signals that the task should be accepted as FINISHED
                # despite the non-zero exit code. This is useful for tasks where
                # "failure" is an expected outcome (e.g., benchmark calculations).
                self.process.report(f"Accepting task {task_name} as finished: {e}")
                self._accept_task_as_finished(task_name)
            except Exception as e:
                error_traceback = traceback.format_exc()  # Capture the full traceback
                self.logger.error(
                    f"Error in running error handler for {task_name}: {e}\\n{error_traceback}"
                )
                self.process.report(
                    f"Error in running error handler for {task_name}: {e}\\n{error_traceback}"
                )

    def _accept_task_as_finished(self, task_name: str) -> None:
        """Mark a failed task as FINISHED and un-skip its child tasks.

        This is called when an error handler raises AcceptError to signal that
        a task's "failure" should be accepted as a successful completion.
        """
        state_manager = self.process.task_manager.state_manager

        # Mark the task as FINISHED instead of FAILED
        state_manager.set_task_runtime_info(task_name, "state", "FINISHED")
        self.process.report(f"Task {task_name} marked as FINISHED (accepted).")

        # Internal WorkGraph tasks that should not be reset manually
        # These are special infrastructure tasks with no executor
        internal_tasks = {"graph_inputs", "graph_outputs", "graph_ctx"}

        # Un-skip child tasks so they can run
        child_names = self.process.wg.connectivity["child_node"].get(task_name, [])
        for child_name in child_names:
            # Skip internal WorkGraph tasks - they are handled automatically
            if child_name in internal_tasks:
                continue
            child_state = state_manager.get_task_runtime_info(child_name, "state")
            if child_state == "SKIPPED":
                state_manager.set_task_runtime_info(child_name, "state", "PLANNED")
                self.process.report(
                    f"Child task {child_name} reset from SKIPPED to PLANNED."
                )
'''


def verify_patch() -> bool:
    """Verify the aiida-workgraph patch is applied.

    Returns:
        True if patch is applied.

    Raises:
        RuntimeError: If patch is not applied.
    """
    try:
        from aiida_workgraph.engine.error_handler_manager import AcceptError

        # Check that _accept_task_as_finished exists
        from aiida_workgraph.engine.error_handler_manager import ErrorHandlerManager

        if not hasattr(ErrorHandlerManager, "_accept_task_as_finished"):
            raise RuntimeError(
                "aiida-workgraph patch incomplete: _accept_task_as_finished method missing"
            )

        return True
    except ImportError:
        raise RuntimeError(
            "aiida-workgraph patch not applied: AcceptError class not found. "
            "Run: python -m teros.core.vasp_parallelization.patches.apply_patches"
        )


def apply_patch(force: bool = False) -> bool:
    """Apply the patch to aiida-workgraph error_handler_manager.py.

    Args:
        force: If True, apply patch even if verification succeeds.

    Returns:
        True if patch was applied, False if already patched.
    """
    # Check if already patched
    if not force:
        try:
            verify_patch()
            print("aiida-workgraph patch already applied.")
            return False
        except RuntimeError:
            pass  # Needs patching

    # Get file path
    file_path = get_error_handler_manager_path()

    # Backup original file
    backup_path = file_path.with_suffix(".py.backup")
    if not backup_path.exists():
        import shutil

        shutil.copy(file_path, backup_path)
        print(f"Backed up original to: {backup_path}")

    # Write patched content
    with open(file_path, "w") as f:
        f.write(PATCHED_CONTENT)

    print(f"Patched: {file_path}")

    # Verify
    # Need to reload the module
    import importlib
    import aiida_workgraph.engine.error_handler_manager

    importlib.reload(aiida_workgraph.engine.error_handler_manager)

    verify_patch()
    print("aiida-workgraph patch verification: OK")

    return True


if __name__ == "__main__":
    apply_patch()
