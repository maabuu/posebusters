"""Logging utilities."""

import logging
import os
import sys

import rdkit
from rdkit import rdBase

# redirect logs to Python logger
rdBase.LogToPythonLogger()


# https://github.com/rdkit/rdkit/discussions/5435
class CaptureLogger(logging.Handler):
    """Helper class that captures Python logger output."""

    def __init__(self, module=None):
        super().__init__(level=logging.NOTSET)
        self.logs = {}
        self.devnull = open(os.devnull, "w")
        rdkit.log_handler.setStream(self.devnull)
        rdkit.logger.addHandler(self)

    def __enter__(self):
        return self.logs

    def __exit__(self, *args):
        self.release()

    def handle(self, record):
        key = record.levelname
        val = self.format(record)
        self.logs[key] = self.logs.get(key, "") + val
        return False

    def release(self):
        rdkit.log_handler.setStream(sys.stderr)
        rdkit.logger.removeHandler(self)
        self.devnull.close()
        return self.logs
