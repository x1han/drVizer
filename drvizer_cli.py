#!/usr/bin/env python3
"""Compatibility wrapper for running the CLI from a repository checkout."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from drvizer.cli import main


if __name__ == "__main__":
    main()
