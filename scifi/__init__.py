#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Commmand line entry point to running the scifi pipeline.
"""

import sys as _sys

from scifi.pipeline import main as __main__


if __name__ == "__main__":
    try:
        _sys.exit(__main__())
    except KeyboardInterrupt:
        _sys.exit(1)
