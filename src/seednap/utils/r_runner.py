"""Base class for running R scripts via Rscript subprocess.

Provides shared logic for R availability checks, script execution,
log file writing, and error handling. Used by DADA2, DADA2 taxonomy,
and DECIPHER runners.
"""

import logging
from pathlib import Path
from typing import List, Optional, Union

from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)

# The bundled R scripts live at the repo root in ``scripts/`` (a sibling of ``src/``).
# Anchor them to the installed package so they resolve from any working directory, and
# always use the canonical scripts rather than whatever ``scripts/`` happens to sit under
# the current CWD (e.g. a stale per-project checkout).
SCRIPTS_DIR = Path(__file__).resolve().parents[3] / "scripts"


def r_script_path(name: str) -> Path:
    """Absolute path to a bundled R script (``scripts/<name>``), CWD-independent."""
    return SCRIPTS_DIR / name


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
        run_subprocess(
            ["Rscript", "--version"],
            timeout=10,
            error_class=self._error_class,
        )

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

        return run_subprocess(
            cmd,
            timeout=self.timeout,
            log_file=log_file,
            log_append=False,
            error_class=self._error_class,
        )
