"""Base class for running R scripts via Rscript subprocess.

Provides shared logic for R availability checks, script execution,
log file writing, and error handling. Used by DADA2, DADA2 taxonomy,
and DECIPHER runners.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

logger = logging.getLogger(__name__)


class RScriptError(Exception):
    """Base exception for R script execution errors."""

    pass


class RScriptRunner:
    """
    Base class for executing R scripts via Rscript.

    Subclasses should set their own error class via `_error_class`
    and call `super().__init__(timeout)` in their constructor.
    """

    _error_class: type = RScriptError

    def __init__(self, timeout: int = 7200):
        """
        Initialize R script runner.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours)
        """
        self.timeout = timeout
        self._check_r_availability()

    def _check_r_availability(self) -> None:
        """
        Check if Rscript is available.

        Raises:
            RScriptError (or subclass): If Rscript is not found
        """
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                check=True,
            )
            logger.debug(f"Found Rscript: {result.stderr.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise self._error_class("Rscript not found. Is R installed?") from e

    def _run_r_script(
        self,
        script_path: Union[str, Path],
        args: List[str],
        log_file: Optional[Union[str, Path]] = None,
    ) -> str:
        """
        Execute R script via Rscript.

        Args:
            script_path: Path to R script file
            args: List of arguments to pass to script
            log_file: Optional path to log file for stdout/stderr

        Returns:
            stdout from Rscript

        Raises:
            RScriptError (or subclass): If Rscript command fails
        """
        script_path = Path(script_path)
        if not script_path.exists():
            raise FileNotFoundError(f"R script not found: {script_path}")

        cmd = ["Rscript", str(script_path)] + [str(arg) for arg in args]
        logger.info(f"Running R script: {script_path.name} {' '.join([str(a) for a in args])}")

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=self.timeout
            )

            # Write to log file if specified
            if log_file:
                log_path = Path(log_file)
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(log_path, "w") as f:
                    f.write(f"Command: {' '.join(cmd)}\n")
                    f.write(f"\n{'='*80}\n")
                    f.write("STDOUT:\n")
                    f.write(result.stdout)
                    f.write(f"\n{'='*80}\n")
                    f.write("STDERR:\n")
                    f.write(result.stderr)

            logger.debug(f"R script completed successfully: {script_path.name}")
            return result.stdout

        except subprocess.CalledProcessError as e:
            error_msg = f"R script failed: {script_path.name}\n{e.stderr}"
            logger.error(error_msg)
            raise self._error_class(error_msg) from e

        except subprocess.TimeoutExpired as e:
            error_msg = f"R script timed out after {self.timeout} seconds: {script_path.name}"
            logger.error(error_msg)
            raise self._error_class(error_msg) from e
