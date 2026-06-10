"""Base class for running R scripts via Rscript subprocess.

Provides shared logic for R availability checks, script execution,
log file writing, and error handling. Used by DADA2, DADA2 taxonomy,
and DECIPHER runners.
"""

import logging
from importlib import resources
from pathlib import Path
from typing import List, Optional, Type, Union

from seednap.utils.subprocess import run_subprocess

logger = logging.getLogger(__name__)

# The bundled R scripts ship inside the package at ``seednap/scripts/``. Anchor them to the
# installed package via importlib.resources (mirroring how darwincore_builder.py resolves its
# bundled CSV templates) so they resolve from any working directory and always use the
# canonical packaged scripts rather than whatever ``scripts/`` happens to sit under the current
# CWD (e.g. a stale per-project checkout).
SCRIPTS_DIR = Path(str(resources.files("seednap").joinpath("scripts")))


def r_script_path(name: str) -> Path:
    """Absolute path to a bundled R script (``seednap/scripts/<name>``), CWD-independent."""
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

    _error_class: Type[Exception] = RScriptError

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
            run_subprocess(
                ["Rscript", "--version"],
                timeout=10,
                error_class=self._error_class,
            )
        except self._error_class as e:
            raise self._error_class(
                f"{e}\n"
                f"  R (Rscript) drives seednap's DADA2 ASV inference and the DADA2-RDP / DECIPHER "
                f"taxonomy steps, so the environment also needs the 'dada2' and 'DECIPHER' R "
                f"packages. Once R is on PATH, verify the packages with: "
                f"Rscript -e 'packageVersion(\"dada2\")'"
            ) from e

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
            raise FileNotFoundError(
                f"Bundled R script not found: {script_path}. seednap resolves its R scripts "
                f"relative to the installed package: they ship inside the package at "
                f"seednap/scripts/ and are bundled in both editable and wheel installs. A missing "
                f"file here means the seednap installation is broken or incomplete -- not a config "
                f"problem. Fix: reinstall from a complete checkout with `pip install -e .` from the "
                f"seednap repo root inside the active conda environment (on the ETH ELE eDNA "
                f"server: /home/shared/edna/envs/seednap), then confirm the scripts are present "
                f"with `ls {SCRIPTS_DIR}` (expect dada2_process.R, taxo_dada2_marker.R, "
                f"taxo_decipher_marker.R)."
            )

        cmd = ["Rscript", str(script_path)] + [str(arg) for arg in args]

        return run_subprocess(
            cmd,
            timeout=self.timeout,
            log_file=log_file,
            log_append=False,
            error_class=self._error_class,
        )
