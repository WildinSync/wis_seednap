"""Pipeline state management for tracking progress and enabling resumability.

This module provides state tracking functionality that allows the pipeline to:
- Track which steps have been completed
- Store step outputs and metadata
- Resume from the last successful step after failure
- Validate dependencies between steps
"""

import json
import logging
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)


class StepStatus(str, Enum):
    """Status of a pipeline step."""

    PENDING = "pending"  # Not yet started
    RUNNING = "running"  # Currently executing
    COMPLETED = "completed"  # Successfully completed
    FAILED = "failed"  # Failed with error
    SKIPPED = "skipped"  # Explicitly skipped


class StepState(BaseModel):
    """State of a single pipeline step."""

    name: str = Field(..., description="Step name (e.g., 'trim', 'dada2')")
    status: StepStatus = Field(default=StepStatus.PENDING, description="Step status")
    started_at: Optional[datetime] = Field(None, description="When step started")
    completed_at: Optional[datetime] = Field(None, description="When step completed")
    duration_seconds: Optional[float] = Field(None, description="Step duration in seconds")
    error_message: Optional[str] = Field(None, description="Error message if failed")
    outputs: Dict[str, Any] = Field(
        default_factory=dict, description="Step output files/data"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict, description="Additional step metadata"
    )

    def start(self) -> None:
        """Mark step as started."""
        self.status = StepStatus.RUNNING
        self.started_at = datetime.now()
        logger.info(f"Step '{self.name}' started")

    def complete(self, outputs: Optional[Dict[str, Any]] = None) -> None:
        """
        Mark step as completed.

        Args:
            outputs: Optional dict of output files/data from the step
        """
        self.status = StepStatus.COMPLETED
        self.completed_at = datetime.now()

        if self.started_at:
            self.duration_seconds = (
                self.completed_at - self.started_at
            ).total_seconds()

        if outputs:
            self.outputs = outputs

        logger.info(
            f"Step '{self.name}' completed in {self.duration_seconds:.1f}s"
            if self.duration_seconds
            else f"Step '{self.name}' completed"
        )

    def fail(self, error: Union[str, Exception]) -> None:
        """
        Mark step as failed.

        Args:
            error: Error message or exception
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
        Mark step as skipped.

        Args:
            reason: Optional reason for skipping
        """
        self.status = StepStatus.SKIPPED
        if reason:
            self.metadata["skip_reason"] = reason
        logger.info(f"Step '{self.name}' skipped" + (f": {reason}" if reason else ""))


class PipelineState(BaseModel):
    """Complete pipeline execution state."""

    marker: str = Field(..., description="Marker name")
    config_path: Optional[Path] = Field(None, description="Path to config file used")
    started_at: datetime = Field(
        default_factory=datetime.now, description="When pipeline started"
    )
    completed_at: Optional[datetime] = Field(
        None, description="When pipeline completed"
    )
    steps: Dict[str, StepState] = Field(
        default_factory=dict, description="State of each pipeline step"
    )
    current_step: Optional[str] = Field(None, description="Currently executing step")
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
        Mark a step as completed.

        Args:
            step_name: Name of the step
            outputs: Optional outputs from the step
        """
        step = self.steps.get(step_name)
        if step:
            step.complete(outputs)
            if self.current_step == step_name:
                self.current_step = None

    def fail_step(self, step_name: str, error: Union[str, Exception]) -> None:
        """
        Mark a step as failed.

        Args:
            step_name: Name of the step
            error: Error message or exception
        """
        step = self.steps.get(step_name)
        if step:
            step.fail(error)
            if self.current_step == step_name:
                self.current_step = None

    def skip_step(self, step_name: str, reason: Optional[str] = None) -> None:
        """
        Mark a step as skipped.

        Args:
            step_name: Name of the step
            reason: Optional reason for skipping
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

    def get_pending_steps(self, all_steps: List[str]) -> List[str]:
        """
        Get list of steps that haven't been completed yet.

        Args:
            all_steps: Complete list of pipeline steps

        Returns:
            List of step names that are pending
        """
        completed = set(self.get_completed_steps())
        skipped = set(
            name
            for name, step in self.steps.items()
            if step.status == StepStatus.SKIPPED
        )
        return [step for step in all_steps if step not in completed and step not in skipped]

    def can_resume(self) -> bool:
        """
        Check if pipeline can be resumed.

        Returns:
            True if there are completed steps and no running steps
        """
        has_completed = any(
            step.status == StepStatus.COMPLETED for step in self.steps.values()
        )
        has_running = any(
            step.status == StepStatus.RUNNING for step in self.steps.values()
        )
        return has_completed and not has_running

    def complete_pipeline(self) -> None:
        """Mark entire pipeline as completed."""
        self.completed_at = datetime.now()
        self.current_step = None
        logger.info("Pipeline completed")

    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of pipeline execution.

        Returns:
            Dictionary with pipeline summary statistics
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
            ValueError: If state file is invalid
        """
        state_file = Path(state_file)
        if not state_file.exists():
            raise FileNotFoundError(f"State file not found: {state_file}")

        with open(state_file) as f:
            state_dict = json.load(f)

        logger.info(f"Loaded pipeline state from {state_file}")
        return cls(**state_dict)

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
        return cls(marker=marker, config_path=config_path)
