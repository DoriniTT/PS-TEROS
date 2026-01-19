"""Custom exceptions for PS-TEROS.

This module defines the exception hierarchy for handling errors in the
PS-TEROS workflow system. All custom exceptions inherit from TerosError
to allow catching all PS-TEROS-specific errors with a single except clause.

Example:
    >>> from teros.core.exceptions import ValidationError
    >>> if param < 0:
    ...     raise ValidationError(
    ...         "Parameter must be positive",
    ...         parameter="encut",
    ...         value=param,
    ...         expected="value > 0"
    ...     )
"""

import typing as t


class TerosError(Exception):
    """Base exception for all PS-TEROS errors.

    All custom exceptions in PS-TEROS inherit from this class, allowing
    users to catch all PS-TEROS-specific errors with a single except clause.

    Attributes:
        message: The error message.
        context: Optional dictionary of additional context information.
    """

    def __init__(self, message: str, **context: t.Any) -> None:
        """Initialize the exception.

        Args:
            message: The error message.
            **context: Additional context information (e.g., parameter names,
                values, expected ranges).
        """
        super().__init__(message)
        self.message = message
        self.context = context

    def __str__(self) -> str:
        """Return a formatted error message with context.

        Returns:
            Formatted error message including context if available.
        """
        if not self.context:
            return self.message

        context_str = ", ".join(f"{k}={v!r}" for k, v in self.context.items())
        return f"{self.message} ({context_str})"


class ValidationError(TerosError):
    """Invalid input parameters or structures.

    Raised when input validation fails, such as:
    - Invalid parameter values (e.g., negative ENCUT)
    - Missing required parameters
    - Parameter values outside expected ranges
    - Invalid parameter combinations

    Attributes:
        message: The error message.
        parameter: Name of the invalid parameter (optional).
        value: The invalid value (optional).
        expected: Description of expected value/range (optional).

    Example:
        >>> raise ValidationError(
        ...     "ENCUT must be positive",
        ...     parameter="encut",
        ...     value=-100,
        ...     expected="value > 0"
        ... )
    """

    def __init__(
        self,
        message: str,
        parameter: t.Optional[str] = None,
        value: t.Optional[t.Any] = None,
        expected: t.Optional[str] = None,
        **context: t.Any
    ) -> None:
        """Initialize the validation error.

        Args:
            message: The error message.
            parameter: Name of the invalid parameter.
            value: The invalid value.
            expected: Description of expected value/range.
            **context: Additional context information.
        """
        if parameter is not None:
            context["parameter"] = parameter
        if value is not None:
            context["value"] = value
        if expected is not None:
            context["expected"] = expected
        super().__init__(message, **context)


class StructureError(TerosError):
    """Problems with atomic structures.

    Raised when issues are encountered with atomic structures, such as:
    - Invalid structure formats
    - Missing atoms or incomplete structures
    - Incompatible structure types
    - Structure generation failures
    - Invalid Miller indices

    Attributes:
        message: The error message.
        structure_type: Type of structure (e.g., 'slab', 'bulk') (optional).
        miller_index: Miller indices for surfaces (optional).

    Example:
        >>> raise StructureError(
        ...     "Cannot generate slab from non-periodic structure",
        ...     structure_type="slab",
        ...     miller_index=(1, 1, 1)
        ... )
    """

    def __init__(
        self,
        message: str,
        structure_type: t.Optional[str] = None,
        miller_index: t.Optional[t.Tuple[int, int, int]] = None,
        **context: t.Any
    ) -> None:
        """Initialize the structure error.

        Args:
            message: The error message.
            structure_type: Type of structure.
            miller_index: Miller indices for surfaces.
            **context: Additional context information.
        """
        if structure_type is not None:
            context["structure_type"] = structure_type
        if miller_index is not None:
            context["miller_index"] = miller_index
        super().__init__(message, **context)


class ConvergenceError(TerosError):
    """DFT calculation convergence failures.

    Raised when DFT calculations fail to converge, such as:
    - Electronic SCF convergence failures
    - Ionic relaxation convergence failures
    - Exceeded maximum iteration limits
    - Energy/force convergence not achieved

    Attributes:
        message: The error message.
        calculation_type: Type of calculation (e.g., 'scf', 'relax') (optional).
        max_iterations: Maximum iterations allowed (optional).
        final_value: Final unconverged value (optional).
        threshold: Convergence threshold (optional).

    Example:
        >>> raise ConvergenceError(
        ...     "SCF failed to converge",
        ...     calculation_type="scf",
        ...     max_iterations=100,
        ...     final_value=1e-3,
        ...     threshold=1e-6
        ... )
    """

    def __init__(
        self,
        message: str,
        calculation_type: t.Optional[str] = None,
        max_iterations: t.Optional[int] = None,
        final_value: t.Optional[float] = None,
        threshold: t.Optional[float] = None,
        **context: t.Any
    ) -> None:
        """Initialize the convergence error.

        Args:
            message: The error message.
            calculation_type: Type of calculation.
            max_iterations: Maximum iterations allowed.
            final_value: Final unconverged value.
            threshold: Convergence threshold.
            **context: Additional context information.
        """
        if calculation_type is not None:
            context["calculation_type"] = calculation_type
        if max_iterations is not None:
            context["max_iterations"] = max_iterations
        if final_value is not None:
            context["final_value"] = final_value
        if threshold is not None:
            context["threshold"] = threshold
        super().__init__(message, **context)


class ConfigurationError(TerosError):
    """Workflow configuration issues.

    Raised when workflow configuration is invalid or inconsistent, such as:
    - Invalid workflow presets
    - Conflicting configuration options
    - Missing required configuration
    - Invalid builder_inputs structure
    - AiiDA code/computer not found

    Attributes:
        message: The error message.
        config_key: Configuration key that caused the error (optional).
        preset: Workflow preset name (optional).

    Example:
        >>> raise ConfigurationError(
        ...     "Invalid workflow preset",
        ...     preset="invalid_preset",
        ...     config_key="workflow_preset"
        ... )
    """

    def __init__(
        self,
        message: str,
        config_key: t.Optional[str] = None,
        preset: t.Optional[str] = None,
        **context: t.Any
    ) -> None:
        """Initialize the configuration error.

        Args:
            message: The error message.
            config_key: Configuration key that caused the error.
            preset: Workflow preset name.
            **context: Additional context information.
        """
        if config_key is not None:
            context["config_key"] = config_key
        if preset is not None:
            context["preset"] = preset
        super().__init__(message, **context)


class WorkflowError(TerosError):
    """Workflow execution errors.

    Raised when errors occur during workflow execution, such as:
    - Task execution failures
    - Missing required outputs
    - Workflow state inconsistencies
    - Inter-task communication errors
    - Unexpected workflow termination

    Attributes:
        message: The error message.
        task_name: Name of the task that failed (optional).
        process_pk: AiiDA process PK (optional).
        workflow_state: Current workflow state (optional).

    Example:
        >>> raise WorkflowError(
        ...     "Task failed to produce required output",
        ...     task_name="relax_bulk",
        ...     process_pk=12345,
        ...     workflow_state="running"
        ... )
    """

    def __init__(
        self,
        message: str,
        task_name: t.Optional[str] = None,
        process_pk: t.Optional[int] = None,
        workflow_state: t.Optional[str] = None,
        **context: t.Any
    ) -> None:
        """Initialize the workflow error.

        Args:
            message: The error message.
            task_name: Name of the task that failed.
            process_pk: AiiDA process PK.
            workflow_state: Current workflow state.
            **context: Additional context information.
        """
        if task_name is not None:
            context["task_name"] = task_name
        if process_pk is not None:
            context["process_pk"] = process_pk
        if workflow_state is not None:
            context["workflow_state"] = workflow_state
        super().__init__(message, **context)
