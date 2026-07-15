#!/usr/bin/env python3
"""CI-faithful test runner.

Runs the suite the way CI does: with third-party pytest plugin autoload
disabled. CI installs no Qt binding, so the module-level PyQt6 stubs the tests
rely on engage cleanly. Locally, an installed pytest-qt otherwise imports a real
PyQt6 during collection, which defeats those stubs and produces spurious
failures (e.g. the _ClickFilter tests). Disabling autoload keeps local runs
matching CI regardless of what is installed.

Usage:
    python run_tests.py                # run the whole suite
    python run_tests.py -k dialog      # forward any args to pytest
"""
import os
import subprocess
import sys
from pathlib import Path

os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"

ROOT = Path(__file__).resolve().parent
args = sys.argv[1:] or ["tests/"]
raise SystemExit(subprocess.call([sys.executable, "-m", "pytest", *args], cwd=ROOT))
