"""
BayesPrism integration for iobrpy.

This subpackage wraps the (Python) BayesPrism implementation and exposes:

- Core classes:
    * Prism
    * BayesPrism

- Result extraction helpers:
    * get_fraction
    * get_exp

- Tutorial-style pipeline:
    * run_bayesprism  (uses the bundled BP_data reference)

- CLI entry point:
    * cli_main        (called from the top-level `iobrpy` CLI)
"""

# Core classes from the original pybayesprism implementation
from .prism import Prism, BayesPrism

# Helper functions to extract fractions and deconvolved expression
from .extract import get_fraction, get_exp

# Tutorial-style pipeline and its CLI main (bayesprism.py)
from .bayesprism import run_bayesprism, main as cli_main

# (Optional) advanced users may also import process_input / references directly
from . import process_input, references

__all__ = [
    # Core model API
    "Prism",
    "BayesPrism",
    # Result extraction helpers
    "get_fraction",
    "get_exp",
    # Tutorial-style pipeline
    "run_bayesprism",
    # CLI entry point used by the top-level iobrpy.main
    "cli_main",
    # Advanced utilities (not strictly required for normal users)
    "process_input",
    "references",
]