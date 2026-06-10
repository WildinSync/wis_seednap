"""Pipeline state management for tracking progress and enabling resumability.

This module provides state tracking functionality that allows the pipeline to:
- Track which steps have been completed
- Store step outputs and metadata
- Resume from the last successful step after failure

Step ordering and dependency validation is not done here; it lives in
config/models/operational.py (PipelineStepsConfig._validate_step_dependencies),
which runs at config load time.
"""

import json
import logging
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field, ValidationError

from seednap.__version__ import __version__ as _SEEDNAP_VERSION
from seednap.errors import SeednapError

logger = logging.getLogger(__name__)


def _serialize_outputs(obj: Any) -> Any:
    """Recursively convert Path objects to strings for JSON serialization.

    Step outputs often contain ``Path`` values (output files); the state is
    persisted as JSON, which has no native path type, so paths are normalized to
    strings before storage. Walks dicts, lists, and tuples; tuples become lists.

    Args:
        obj: Any value to serialize. Paths become strings; dicts/lists/tuples are
            walked recursively; all other values are returned unchanged.

    Returns:
        The same structure with every ``Path`` replaced by its string form and
        every tuple replaced by a list. Non-container, non-Path values pass
        through untouched.
    """
    if isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: _serialize_outputs(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_serialize_outputs(v) for v in obj]
    return obj


class StepStatus(str, Enum):
    """Lifecycle status of a single pipeline step.

    A string-valued enum so the value serializes directly into the state JSON.
    These are the states the orchestrator transitions a step through, and which
    ``--resume`` reads back to decide what to re-run.
    """

    PENDING = "pending"  # Not yet started
    RUNNING = "running"  # Currently executing
    COMPLETED = "completed"  # Successfully completed
    FAILED = "failed"  # Failed with error
    SKIPPED = "skipped"  # Explicitly skipped


class StepState(BaseModel):
    """Recorded state of a single pipeline step (one stage of one run).

    Captures the step's status, timing, error message, output file paths, and
    free-form metadata. This is the per-step unit that makes a run auditable and
    resumable: it records what ran, how long it took, what it produced, and why it
    failed if it did.
    """

    name: str = Field(..., description="Step name (e.g., 'trim', 'dada2')")
    status: StepStatus = Field(default=StepStatus.PENDING, description="Step status")
    started_at: Optional[datetime] = Field(default=None, description="When step started")
    completed_at: Optional[datetime] = Field(
        default=None, description="When step completed"
    )
    duration_seconds: Optional[float] = Field(
        default=None, description="Step duration in seconds"
    )
    error_message: Optional[str] = Field(
        default=None, description="Error message if failed"
    )
    outputs: Dict[str, Any] = Field(
        default_factory=dict, description="Step output files/data"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict, description="Additional step metadata"
    )

    def start(self) -> None:
        """Mark this step as RUNNING and stamp the start time.

        Args:
            None.

        Returns:
            None.
        """
        self.status = StepStatus.RUNNING
        self.started_at = datetime.now()
        logger.info(f"Step '{self.name}' started")

    def complete(self, outputs: Optional[Dict[str, Any]] = None) -> None:
        """
        Mark this step as COMPLETED, stamp the time, and record its outputs.

        Sets the completion time, computes the duration from the recorded start time
        (if any), and stores the step's outputs with Path values normalized to
        strings for JSON serialization.

        Args:
            outputs: Optional dict of output files/data from the step. Path values
                are converted to strings before storage. If falsy, outputs are left
                unchanged.

        Returns:
            None.
        """
        self.status = StepStatus.COMPLETED
        self.completed_at = datetime.now()

        if self.started_at:
            self.duration_seconds = (
                self.completed_at - self.started_at
            ).total_seconds()

        if outputs:
            # Normalize Path values to strings for consistent serialization
            self.outputs = _serialize_outputs(outputs)

        logger.info(
            f"Step '{self.name}' completed in {self.duration_seconds:.1f}s"
            if self.duration_seconds
            else f"Step '{self.name}' completed"
        )

    def fail(self, error: Union[str, Exception]) -> None:
        """
        Mark this step as FAILED, stamp the time, and record the error message.

        Sets the completion time, computes the duration from the recorded start time
        (if any), and stores the error's string form so a resumed run and the audit
        log can show why this step broke.

        Args:
            error: The error message (str) or exception; stored as its string form.

        Returns:
            None.
        """
        self.status = StepStatus.FAILED
        self.completed_at = datetime.now()

        if self.started_at:
            self.duration_seconds = (
                self.completed_at - self.started_at
            ).total_seconds()

        self.error_message = str(error)
        logger.error(f"Step '{self.name}' failed: {self.error_message}")

    def skip(self, reason: Optional[str] = None) -> None:
        """
        Mark this step as SKIPPED, optionally recording why.

        Args:
            reason: Optional human-readable reason for the skip; stored under the
                ``skip_reason`` metadata key when provided.

        Returns:
            None.
        """
        self.status = StepStatus.SKIPPED
        if reason:
            self.metadata["skip_reason"] = reason
        logger.info(f"Step '{self.name}' skipped" + (f": {reason}" if reason else ""))


class PipelineState(BaseModel):
    """Complete execution state of one marker's pipeline run.

    The source of truth for "did this run finish, and what did each step do?".
    Holds the marker, the seednap version that wrote it, config/snapshot paths,
    overall timing, and the per-step StepState map. Serialized to JSON so a run can
    be audited after the fact and resumed from the last good step.
    """

    marker: str = Field(..., description="Marker name")
    seednap_version: Optional[str] = Field(
        default=None,
        description="seednap version that created this run's state "
        "(None = a state file written before version stamping existed)",
    )
    config_path: Optional[Path] = Field(
        default=None, description="Path to config file used"
    )
    config_snapshot_path: Optional[Path] = Field(
        default=None,
        description="Path to the effective merged-config YAML snapshot for this run",
    )
    started_at: datetime = Field(
        default_factory=datetime.now, description="When pipeline started"
    )
    completed_at: Optional[datetime] = Field(
        default=None, description="When pipeline completed"
    )
    steps: Dict[str, StepState] = Field(
        default_factory=dict, description="State of each pipeline step"
    )
    current_step: Optional[str] = Field(
        default=None, description="Currently executing step"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict, description="Pipeline metadata"
    )

    def add_step(self, step_name: str) -> StepState:
        """
        Add a new step to the pipeline state.

        Args:
            step_name: Name of the step

        Returns:
            StepState object for the new step
        """
        if step_name in self.steps:
            logger.warning(f"Step '{step_name}' already exists in state")
            return self.steps[step_name]

        step_state = StepState(name=step_name)
        self.steps[step_name] = step_state
        return step_state

    def get_step(self, step_name: str) -> Optional[StepState]:
        """
        Get state for a specific step.

        Args:
            step_name: Name of the step

        Returns:
            StepState object or None if step doesn't exist
        """
        return self.steps.get(step_name)

    def start_step(self, step_name: str) -> StepState:
        """
        Start a pipeline step.

        Args:
            step_name: Name of the step to start

        Returns:
            StepState object for the step
        """
        if step_name not in self.steps:
            self.add_step(step_name)

        step = self.steps[step_name]
        step.start()
        self.current_step = step_name
        return step

    def complete_step(
        self, step_name: str, outputs: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Mark a named step as completed and clear it as the current step.

        Args:
            step_name: Name of the step to complete.
            outputs: Optional dict of output files/data to record on the step.

        Returns:
            None. No-op if the named step is not present in the state.
        """
        step = self.steps.get(step_name)
        if step:
            step.complete(outputs)
            if self.current_step == step_name:
                self.current_step = None

    def fail_step(self, step_name: str, error: Union[str, Exception]) -> None:
        """
        Mark a named step as failed and clear it as the current step.

        Args:
            step_name: Name of the step to fail.
            error: The error message (str) or exception to record on the step.

        Returns:
            None. No-op if the named step is not present in the state.
        """
        step = self.steps.get(step_name)
        if step:
            step.fail(error)
            if self.current_step == step_name:
                self.current_step = None

    def skip_step(self, step_name: str, reason: Optional[str] = None) -> None:
        """
        Mark a named step as skipped, adding it to the state first if absent.

        Args:
            step_name: Name of the step to skip.
            reason: Optional human-readable reason, recorded on the step's metadata.

        Returns:
            None.
        """
        if step_name not in self.steps:
            self.add_step(step_name)

        step = self.steps[step_name]
        step.skip(reason)

    def is_step_completed(self, step_name: str) -> bool:
        """
        Check if a step has been completed.

        Args:
            step_name: Name of the step

        Returns:
            True if step is completed, False otherwise
        """
        step = self.steps.get(step_name)
        return step.status == StepStatus.COMPLETED if step else False

    def is_step_failed(self, step_name: str) -> bool:
        """
        Check if a step has failed.

        Args:
            step_name: Name of the step

        Returns:
            True if step failed, False otherwise
        """
        step = self.steps.get(step_name)
        return step.status == StepStatus.FAILED if step else False

    def get_completed_steps(self) -> List[str]:
        """
        Get list of completed step names.

        Returns:
            List of step names that are completed
        """
        return [
            name
            for name, step in self.steps.items()
            if step.status == StepStatus.COMPLETED
        ]

    def get_failed_steps(self) -> List[str]:
        """
        Get list of failed step names.

        Returns:
            List of step names that failed
        """
        return [
            name for name, step in self.steps.items() if step.status == StepStatus.FAILED
        ]

    def complete_pipeline(self) -> None:
        """Mark the whole run as finished: stamp completion time, clear current step.

        Note: "completed" here means the orchestrator's step loop ran to the end, not
        that every step succeeded; under continue-on-error some steps may have failed.

        Args:
            None.

        Returns:
            None.
        """
        self.completed_at = datetime.now()
        self.current_step = None
        logger.info("Pipeline completed")

    def get_summary(self) -> Dict[str, Any]:
        """
        Summarize the run's progress for logging and reporting.

        Args:
            None.

        Returns:
            Dictionary with: ``marker`` (str); ``started_at`` / ``completed_at``
            (ISO-8601 strings or None); ``total_duration_seconds`` (float or None);
            the counts ``total_steps``, ``completed``, ``failed``, ``skipped``,
            ``pending`` (ints); ``current_step`` (str or None); and ``steps``, a dict
            keyed by step name, each value holding ``status`` (str),
            ``duration_seconds`` (float or None), and ``error`` (str or None).
        """
        total_steps = len(self.steps)
        completed = len(self.get_completed_steps())
        failed = len(self.get_failed_steps())
        skipped = sum(1 for s in self.steps.values() if s.status == StepStatus.SKIPPED)
        pending = total_steps - completed - failed - skipped

        total_duration = None
        if self.completed_at and self.started_at:
            total_duration = (self.completed_at - self.started_at).total_seconds()

        return {
            "marker": self.marker,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": (
                self.completed_at.isoformat() if self.completed_at else None
            ),
            "total_duration_seconds": total_duration,
            "total_steps": total_steps,
            "completed": completed,
            "failed": failed,
            "skipped": skipped,
            "pending": pending,
            "current_step": self.current_step,
            "steps": {
                name: {
                    "status": step.status.value,
                    "duration_seconds": step.duration_seconds,
                    "error": step.error_message,
                }
                for name, step in self.steps.items()
            },
        }

    def save(self, state_file: Union[str, Path]) -> Path:
        """
        Save pipeline state to JSON file.

        Args:
            state_file: Path to state file

        Returns:
            Path to saved state file
        """
        state_file = Path(state_file)
        state_file.parent.mkdir(parents=True, exist_ok=True)

        # Convert to dict and handle datetime serialization
        state_dict = self.model_dump(mode="python")

        with open(state_file, "w") as f:
            json.dump(state_dict, f, indent=2, default=str)

        logger.debug(f"Saved pipeline state to {state_file}")
        return state_file

    @classmethod
    def load(cls, state_file: Union[str, Path]) -> "PipelineState":
        """
        Load pipeline state from JSON file.

        Args:
            state_file: Path to state file

        Returns:
            PipelineState object

        Raises:
            FileNotFoundError: If state file doesn't exist
            SeednapError: If the state file is unreadable, corrupted, or schema-incompatible
        """
        state_file = Path(state_file)
        if not state_file.exists():
            raise FileNotFoundError(f"State file not found: {state_file}")

        try:
            with open(state_file) as f:
                state_dict = json.load(f)
            state = cls(**state_dict)
        except (json.JSONDecodeError, ValidationError) as exc:
            raise SeednapError(
                f"Could not parse pipeline state file {state_file}: {exc}",
                why=(
                    "The file exists but is not a valid seednap state JSON. It may be "
                    "truncated (a run killed mid-save), corrupted, or written by an "
                    "incompatible seednap version or for a different marker."
                ),
                fix=(
                    "Delete or move it and start a fresh run WITHOUT --resume (a new "
                    "state file is created automatically), or restore a known-good "
                    "copy. The default state file is "
                    "<config.paths.output>/.<marker>_state.json."
                ),
            ) from exc

        logger.info(f"Loaded pipeline state from {state_file}")
        return state

    @classmethod
    def from_config(
        cls, marker: str, config_path: Optional[Path] = None
    ) -> "PipelineState":
        """
        Create new pipeline state from configuration.

        Args:
            marker: Marker name
            config_path: Optional path to config file

        Returns:
            New PipelineState object
        """
        # New runs are stamped with the running version; a None on a loaded state therefore
        # signals a state file written before version stamping existed (handled on resume).
        return cls(
            marker=marker, config_path=config_path, seednap_version=_SEEDNAP_VERSION
        )
