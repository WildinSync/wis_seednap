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
    """Return the absolute path to a bundled R script, independent of CWD.

    Resolves ``seednap/scripts/<name>`` against the installed package (e.g.
    dada2_process.R, taxo_dada2_marker.R, taxo_decipher_marker.R) so the
    correct packaged script is used regardless of the working directory.

    Args:
        name: File name of the bundled R script (e.g. 'dada2_process.R').

    Returns:
        Absolute Path to the script under the installed package's scripts/
        directory. The path is not checked for existence here; the caller
        (`_run_r_script`) validates that before invoking Rscript.

    Raises:
        None.
    """
    return SCRIPTS_DIR / name


class RScriptError(Exception):
    """Base exception for R script execution errors.

    Raised when Rscript is unavailable or when a bundled R script (DADA2 /
    DECIPHER) exits with an error. Subclassed by the DADA2 and DECIPHER runners
    so callers can distinguish which R step failed.
    """

    pass


class RScriptRunner:
    """
    Base class for executing R scripts via Rscript.

    DADA2 (ASV inference) and DADA2-RDP / DECIPHER (taxonomic assignment) are
    implemented in R; this base wraps invoking them as Rscript subprocesses so
    the Python orchestrator can call them uniformly. Subclasses should set their
    own error class via `_error_class` and call `super().__init__(timeout)` in
    their constructor.

    Attributes:
        _error_class: Exception type raised on failure; overridden per subclass
            so the failing R step is identifiable (defaults to RScriptError).
        timeout: Per-command timeout in seconds, set in __init__.
    """

    _error_class: Type[Exception] = RScriptError

    def __init__(self, timeout: int = 7200):
        """
        Initialize R script runner and verify Rscript is available.

        Args:
            timeout: Command timeout in seconds (default: 7200 = 2 hours).

        Raises:
            RScriptError (or subclass `_error_class`): If Rscript is not on
                PATH (checked eagerly here via `_check_r_availability`).
        """
        self.timeout = timeout
        self._check_r_availability()

    def _check_r_availability(self) -> None:
        """
        Verify that Rscript is installed and on PATH.

        Runs `Rscript --version` with a short timeout; on failure, re-raises
        with guidance that the environment also needs the dada2 and DECIPHER R
        packages. DADA2 and the DADA2-RDP / DECIPHER taxonomy steps cannot run
        without R, so this fails fast at construction rather than mid-pipeline.

        Returns:
            None.

        Raises:
            RScriptError (or subclass `_error_class`): If Rscript is not found.
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
        Execute an R script via Rscript and return its stdout.

        Validates that the bundled script exists, then runs it with the given
        arguments through `run_subprocess` (so the call is logged and any tool
        failure becomes a readable error). The per-step log file is opened in
        write mode, so each invocation overwrites the previous R log.

        Args:
            script_path: Path to the R script file to execute.
            args: Arguments passed to the script (each coerced to str).
            log_file: Optional path to a file that receives the script's
                stdout/stderr (overwritten, not appended).

        Returns:
            The script's standard output as a string.

        Raises:
            FileNotFoundError: If `script_path` does not exist (indicates a
                broken or incomplete seednap installation, not a config error).
            RScriptError (or subclass `_error_class`): If the Rscript command
                exits non-zero, times out, or Rscript is not found.
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
