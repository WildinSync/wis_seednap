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
        error_msg = f"Command failed: {cmd_str}\n{e.stderr}"
        logger.error(error_msg)
        raise error_class(error_msg) from e

    except subprocess.TimeoutExpired as e:
        error_msg = f"Command timed out after {timeout}s: {cmd_str}"
        logger.error(error_msg)
        raise error_class(error_msg) from e

    except FileNotFoundError as e:
        error_msg = f"Command not found: {cmd[0]}"
        logger.error(error_msg)
        raise error_class(error_msg) from e
