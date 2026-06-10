"""Logging configuration and utilities."""

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
        """Ensure only one instance exists."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """Initialize logging configuration (only once)."""
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
        Configure logging for the application.

        Args:
            level: Logging level (DEBUG, INFO, WARNING, ERROR)
            log_file: Optional path to log file
            format_style: Log format style (simple, detailed, json)
            console_output: Whether to output logs to console

        Returns:
            Configured root logger
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

        # Console handler with rich formatting
        if console_output:
            console_handler = RichHandler(
                console=self.console,
                show_time=format_style == "detailed",
                show_path=format_style == "detailed",
                rich_tracebacks=True,
                tracebacks_show_locals=level == "DEBUG",
            )
            console_handler.setLevel(numeric_level)
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

    Args:
        name: Logger name (usually __name__ of the module)

    Returns:
        Logger instance
    """
    return logging.getLogger(name)


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_style: str = "detailed",
    console_output: bool = True,
) -> logging.Logger:
    """
    Set up logging configuration (convenience function).

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional path to log file
        format_style: Log format style (simple, detailed, json)
        console_output: Whether to output logs to console

    Returns:
        Configured root logger
    """
    config = LogConfig()
    return config.setup_logging(
        level=level, log_file=log_file, format_style=format_style, console_output=console_output
    )


def log_pipeline_step(step_name: str, status: str = "start", logger: Optional[logging.Logger] = None) -> None:
    """
    Log a pipeline step.

    Args:
        step_name: Name of the pipeline step
        status: Status of the step. "start" and "complete" log at INFO,
            "error" logs at ERROR. Any other status string is accepted and
            logged generically at INFO ("Pipeline step <name>: <status>");
            no status value is rejected.
        logger: Logger instance (uses root if None)
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
