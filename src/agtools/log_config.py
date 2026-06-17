#!/usr/bin/env python3

import sys
from typing import Optional

from loguru import logger

from agtools import __version__

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"

__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def configure_logger(log_file: Optional[str] = None, use_stderr: bool = False):
    """Configure agtools logging sinks.

    By default, logs are emitted to stdout only. If ``log_file`` is
    provided, logs are written to that file instead. When command output
    is routed to stdout, ``use_stderr`` can be enabled to keep logs on
    stderr.
    """
    logger.remove()
    if log_file:
        logger.add(sink=log_file, level="DEBUG")
    elif use_stderr:
        logger.add(sink=sys.stderr, level="INFO")
    else:
        logger.add(sink=sys.stdout, level="INFO")


# Default configuration for library/programmatic usage: console logs only.
configure_logger()
