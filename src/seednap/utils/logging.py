"""Logging configuration and utilities.

Central logging setup for the pipeline. Every CLI command and the orchestrator
route their output through here so a run produces one consistent log stream:
a Rich-formatted console handler (kept at WARNING+ even in quiet mode so the
no-silent-fallback [WARN] messages stay visible) and an optional per-marker
file handler at logs/<marker>_pipeline_run.log. This supports the
reproducibility requirement that every run be reconstructable from its log.
Lives in seednap/utils/ alongside the subprocess and R-runner helpers.
"""

import logging
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.logging import RichHandler


class LogConfig:
    """Logging configuration singleton.

    Re-entrancy contract: __init__ runs its body only once (guarded by
    _initialized), so the stderr Console created on the first construction is
    cached and reused for the lifetime of the process. setup_logging, by
    contrast, may be called more than once (cli.py and orchestrator.py both
    call the module-level setup_logging); each call clears the root logger's
    handlers and re-adds fresh console/file handlers. In short: first call wins
    for the Console, every call resets the root handlers.
    """

    _instance: Optional["LogConfig"] = None
    _initialized: bool = False

    def __new__(cls) -> "LogConfig":
        """Return the single shared LogConfig instance (singleton).

        Args:
            None.

        Returns:
            The one cached LogConfig instance, creating it on first call.

        Raises:
            None.
        """
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """Initialize logging configuration, running its body only once.

        Creates the cached stderr Console and sets default level/log-file
        state on the first construction; subsequent constructions return early
        (guarded by `_initialized`) so the Console is reused for the process.

        Args:
            None.

        Returns:
            None.

        Raises:
            None.
        """
        if LogConfig._initialized:
            return

        self.console = Console(stderr=True)
        self.log_file: Optional[Path] = None
        self.log_level = logging.INFO
        LogConfig._initialized = True

    def setup_logging(
        self,
        level: str = "INFO",
        log_file: Optional[Path] = None,
        format_style: str = "detailed",
        console_output: bool = True,
    ) -> logging.Logger:
        """
        Configure root logging for the application.

        Resets the root logger's handlers and attaches a fresh Rich console
        handler (and a file handler when `log_file` is given). May be called
        more than once; each call re-applies the handlers (last call wins for
        level and file). The console handler is forced to WARNING+ when
        `console_output` is False so safety [WARN] messages remain visible
        under quiet mode while INFO/DEBUG chatter is suppressed.

        Args:
            level: Logging level name (DEBUG, INFO, WARNING, ERROR);
                unrecognised values fall back to INFO.
            log_file: Optional path to a log file; its parent directory is
                created if needed and the file is opened in append mode.
            format_style: Log format style ('simple', 'detailed', or 'json').
                'detailed' adds timestamp/path on the console; 'json' uses a
                JSON-shaped (not escape-safe) file format.
            console_output: If False, the console handler is limited to
                WARNING and above; if True it uses `level`.

        Returns:
            The configured root logger.

        Raises:
            None.
        """
        # Convert string level to logging constant
        numeric_level = getattr(logging, level.upper(), logging.INFO)
        self.log_level = numeric_level
        self.log_file = log_file

        # Get root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(numeric_level)

        # Remove existing handlers
        root_logger.handlers.clear()

        # Console handler with rich formatting. Always attached so safety [WARN]s and errors stay
        # visible on screen; in quiet mode (console_output=False) it is restricted to WARNING+ so
        # only INFO/DEBUG chatter is silenced. This keeps the no-silent-fallback warnings on the
        # console even under --quiet.
        console_handler = RichHandler(
            console=self.console,
            show_time=format_style == "detailed",
            show_path=format_style == "detailed",
            rich_tracebacks=True,
            tracebacks_show_locals=level == "DEBUG",
        )
        console_handler.setLevel(numeric_level if console_output else logging.WARNING)
        root_logger.addHandler(console_handler)

        # File handler with detailed formatting
        if log_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)

            file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
            file_handler.setLevel(numeric_level)

            # Detailed format for file logs
            if format_style == "json":
                # WARNING: this is a JSON-shaped template, not guaranteed-valid
                # JSON. The message is interpolated raw via %(message)s with no
                # escaping, so any log message containing a double-quote,
                # backslash, or newline (e.g. multi-line what/why/fix errors or
                # a subprocess stderr dump) produces malformed JSON. A real JSON
                # formatter is needed before feeding this to a log shipper.
                file_format = (
                    '{"time": "%(asctime)s", "level": "%(levelname)s", '
                    '"logger": "%(name)s", "message": "%(message)s"}'
                )
            else:
                file_format = (
                    "%(asctime)s | %(levelname)-8s | %(name)s:%(lineno)d | %(message)s"
                )

            file_formatter = logging.Formatter(file_format, datefmt="%Y-%m-%d %H:%M:%S")
            file_handler.setFormatter(file_formatter)
            root_logger.addHandler(file_handler)

        return root_logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.

    Thin wrapper over logging.getLogger so modules import logger access from
    one place; the returned logger inherits the root handlers configured by
    setup_logging.

    Args:
        name: Logger name (usually __name__ of the calling module).

    Returns:
        The logging.Logger instance for that name.

    Raises:
        None.
    """
    return logging.getLogger(name)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_style: str = "detailed",
    console_output: bool = True,
) -> logging.Logger:
    """
    Set up logging configuration (module-level convenience function).

    Obtains the shared LogConfig singleton and delegates to its setup_logging.
    Both cli.py and the orchestrator call this; the singleton ensures the
    stderr Console is shared while handlers are refreshed on each call.

    Args:
        level: Logging level name (DEBUG, INFO, WARNING, ERROR);
            unrecognised values fall back to INFO.
        log_file: Optional path to a log file (append mode; parent created).
        format_style: Log format style ('simple', 'detailed', or 'json').
        console_output: Whether INFO/DEBUG logs reach the console; if False,
            the console handler is limited to WARNING and above.

    Returns:
        The configured root logger.

    Raises:
        None.
    """
    config = LogConfig()
    return config.setup_logging(
        level=level, log_file=log_file, format_style=format_style, console_output=console_output
    )


def log_pipeline_step(step_name: str, status: str = "start", logger: Optional[logging.Logger] = None) -> None:
    """
    Emit a one-line log message marking a pipeline step transition.

    A convenience helper for recording when a step (trimming, DADA2, SWARM,
    taxonomy, formatting) starts, completes, or errors, so the run log reads as
    a sequence of step boundaries.

    Args:
        step_name: Name of the pipeline step.
        status: Status of the step. "start" and "complete" log at INFO,
            "error" logs at ERROR. Any other status string is accepted and
            logged generically at INFO ("Pipeline step <name>: <status>");
            no status value is rejected.
        logger: Logger instance to use; the root logger is used if None.

    Returns:
        None.

    Raises:
        None.
    """
    if logger is None:
        logger = logging.getLogger()

    if status == "start":
        logger.info(f"Starting pipeline step: {step_name}")
    elif status == "complete":
        logger.info(f"Completed pipeline step: {step_name}")
    elif status == "error":
        logger.error(f"Error in pipeline step: {step_name}")
    else:
        logger.info(f"Pipeline step {step_name}: {status}")
