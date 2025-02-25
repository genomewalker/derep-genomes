"""
Logging utilities for genome dereplication.
"""

import os
import sys
import time
import logging
from typing import Optional

# Default log format
DEFAULT_LOG_FORMAT = "%(levelname)s ::: %(asctime)s ::: %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    log_level: int = logging.INFO,
    log_file: Optional[str] = None,
    log_format: str = DEFAULT_LOG_FORMAT,
    date_format: str = DEFAULT_DATE_FORMAT,
) -> logging.Logger:
    """
    Set up logging configuration.

    Args:
        log_level: Logging level (default: INFO)
        log_file: Path to log file (default: None, log to stdout only)
        log_format: Log message format
        date_format: Date format for log messages

    Returns:
        Configured logger
    """
    # Create logger
    logger = logging.getLogger("derep_genomes")
    logger.setLevel(log_level)

    # Remove existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Create formatter
    formatter = logging.Formatter(log_format, date_format)

    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Add file handler if log file is specified
    if log_file:
        # Create directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)

        # Add file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_logger() -> logging.Logger:
    """
    Get the derep_genomes logger.

    Returns:
        Logger instance
    """
    return logging.getLogger("derep_genomes")


def is_debug_enabled() -> bool:
    """
    Check if debug logging is enabled.

    Returns:
        True if debug logging is enabled, False otherwise
    """
    return get_logger().getEffectiveLevel() == logging.DEBUG


class Timer:
    """Simple timer for measuring execution time."""

    def __init__(self, name: str = ""):
        """
        Initialize timer.

        Args:
            name: Timer name
        """
        self.name = name
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        """Start timer when entering context."""
        self.start_time = time.time()
        return self

    def __exit__(self, *args):
        """Stop timer when exiting context."""
        self.end_time = time.time()
        elapsed = self.end_time - self.start_time

        if self.name:
            get_logger().info(f"{self.name} completed in {elapsed:.2f} seconds")
        else:
            get_logger().info(f"Operation completed in {elapsed:.2f} seconds")

    def elapsed(self) -> float:
        """
        Get elapsed time.

        Returns:
            Elapsed time in seconds
        """
        if self.start_time is None:
            return 0

        end_time = self.end_time or time.time()
        return end_time - self.start_time


def log_memory_usage(prefix: str = "Memory usage"):
    """
    Log current memory usage.

    Args:
        prefix: Prefix for log message
    """
    try:
        import psutil

        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()

        # Convert to MB
        rss_mb = memory_info.rss / (1024 * 1024)
        vms_mb = memory_info.vms / (1024 * 1024)

        get_logger().info(f"{prefix}: RSS={rss_mb:.2f} MB, VMS={vms_mb:.2f} MB")
    except ImportError:
        get_logger().debug("psutil not available, cannot log memory usage")
