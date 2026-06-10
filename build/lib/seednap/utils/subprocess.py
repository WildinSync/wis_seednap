"""Shared subprocess execution with logging and error handling."""

import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

logger = logging.getLogger(__name__)


def run_subprocess(
    cmd: List[str],
    *,
    timeout: int = 3600,
    log_file: Optional[Union[str, Path]] = None,
    log_append: bool = True,
    error_class: type = Exception,
) -> str:
    """
    Run a subprocess command with unified logging and error handling.

    Args:
        cmd: Command and arguments to execute.
        timeout: Maximum runtime in seconds (default: 3600).
        log_file: Optional path to write stdout/stderr output.
        log_append: Open log file in append mode (default) or write mode.
        error_class: Exception class to raise on failure.

    Returns:
        Standard output from the command.

    Raises:
        error_class: If the command fails, times out, or is not found.
    """
    cmd_str = " ".join(str(c) for c in cmd)
    logger.info(f"Running: {cmd_str}")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=timeout
        )

        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            mode = "a" if log_append else "w"
            with open(log_path, mode) as f:
                f.write(f"Command: {cmd_str}\n")
                f.write(f"\n{'='*80}\n")
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write(f"\n{'='*80}\n")
                f.write("STDERR:\n")
                f.write(result.stderr)
                f.write(f"\n{'='*80}\n\n")

        logger.debug("Command completed successfully")
        return result.stdout

    except subprocess.CalledProcessError as e:
        stderr = (e.stderr or "").strip() or "(no stderr captured)"
        error_msg = (
            f"External tool '{cmd[0]}' exited with status {e.returncode} (it ran but returned "
            f"an error). The text below is {cmd[0]}'s own output, not a seednap bug; read it "
            f"first.\n"
            f"  command: {cmd_str}\n"
            f"  --- {cmd[0]} stderr ---\n{stderr}\n"
            f"  --- end stderr ---\n"
            f"  Common causes: a malformed or empty input file, a wrong reference/database path, "
            f"or a tool-version mismatch.\n"
            f"  [SDN-TOOL-002] (run `seednap explain SDN-TOOL-002` for more)"
        )
        logger.error(error_msg)
        raise error_class(error_msg) from e

    except subprocess.TimeoutExpired as e:
        error_msg = (
            f"External tool '{cmd[0]}' did not finish within {timeout}s and was killed.\n"
            f"  command: {cmd_str}\n"
            f"  This usually means the dataset is large or the machine is heavily loaded. "
            f"Re-run on a quieter machine, or raise the step timeout."
        )
        logger.error(error_msg)
        raise error_class(error_msg) from e

    except FileNotFoundError as e:
        error_msg = (
            f"Required tool '{cmd[0]}' is not installed or not on PATH: seednap tried to launch "
            f"it and the OS could not find it. Usually the wrong (or no) conda environment is "
            f"active.\n"
            f"  Fix: activate the environment that has seednap's tools, then confirm it is there:\n"
            f"    conda activate /home/shared/edna/envs/seednap   # ETH ELE eDNA server\n"
            f"    {cmd[0]} --version\n"
            f"  seednap's external tools (cutadapt, vsearch, swarm, blastn/makeblastdb, Rscript) "
            f"live in that environment; if '{cmd[0]}' is genuinely missing, install it there.\n"
            f"  [SDN-TOOL-001] (run `seednap explain SDN-TOOL-001` for more)"
        )
        logger.error(error_msg)
        raise error_class(error_msg) from e
