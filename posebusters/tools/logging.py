"""Logging utilities."""

from __future__ import annotations

import logging
import os
import sys
from typing import Any

import rdkit
from rdkit import rdBase

# redirect logs to Python logger
rdBase.LogToPythonLogger()


# https://github.com/rdkit/rdkit/discussions/5435
class CaptureLogger(logging.Handler):
    """Helper class that captures Python logger output."""

    def __init__(self, level: str | None = None) -> None:
        """Initialize logger."""
        super().__init__(level=logging.NOTSET if level is None else level)
        self.logs: dict[str, str] = {}
        self.devnull = open(os.devnull, "w")
        rdkit.log_handler.setStream(self.devnull)
        rdkit.logger.addHandler(self)

    def __enter__(self) -> dict[str, str]:
        """Enter context manager."""
        return self.logs

    def __exit__(self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: Any) -> None:
        """Exit context manager."""
        self.release()

    def handle(self, record: logging.LogRecord) -> bool:
        """Handle log record."""
        key = record.levelname
        val = self.format(record)
        self.logs[key] = self.logs.get(key, "") + val
        return False

    def release(self) -> None:
        """Release logger."""
        rdkit.log_handler.setStream(sys.stderr)
        rdkit.logger.removeHandler(self)
        self.devnull.close()

    def emit(self, record: logging.LogRecord) -> None:
        """Emit a log record."""
