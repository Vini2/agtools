#!/usr/bin/env python3

import sys

from loguru import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


# Remove the default logger configuration
logger.remove()

# Console logging (INFO level and above)
logger.add(sink=sys.stdout, level="INFO")

# File logging (DEBUG level and above)
logger.add(sink="agtools.log", level="DEBUG")
