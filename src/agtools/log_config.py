#!/usr/bin/env python3

import sys
from loguru import logger

# Remove the default logger configuration
logger.remove()

# Console logging (INFO level and above)
logger.add(sink=sys.stdout, level="INFO")

# File logging (DEBUG level and above)
logger.add(sink="agtools.log", level="DEBUG")