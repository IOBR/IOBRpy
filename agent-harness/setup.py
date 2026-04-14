"""
setup.py - PyPI package configuration for iobrpy-cli.

This package provides a CLI harness for IOBRpy (Immuno-Oncology Biological
Research using Python) with REPL support and JSON output mode.
"""

from pathlib import Path
from setuptools import setup

HERE = Path(__file__).resolve().parent
README_PATH = HERE / "iobrpy_cli" / "iobrpy" / "README.md"
PACKAGES = [
    "iobrpy_cli",
    "iobrpy_cli.iobrpy",
    "iobrpy_cli.iobrpy.core",
    "iobrpy_cli.iobrpy.tests",
    "iobrpy_cli.iobrpy.utils",
]

setup(
    name="iobrpy-cli",
    version="0.1.0",
    description="CLI harness for IOBRpy with REPL support and JSON output",
    long_description=README_PATH.read_text(encoding='utf-8') if README_PATH.exists() else (
        "CLI harness for IOBRpy (Immuno-Oncology Biological Research using Python) "
        "with REPL support, JSON output mode, and project management."
    ),
    long_description_content_type="text/markdown",
    author="Claude Code",
    author_email="noreply@anthropic.com",
    url="https://github.com/IOBR/IOBRpy",
    license="MIT",

    packages=PACKAGES,
    package_dir={"iobrpy_cli": "iobrpy_cli"},
    package_data={
        'iobrpy_cli.iobrpy': ['README.md'],
    },

    # Entry point for CLI
    entry_points={
        'console_scripts': [
            'iobrpy-cli=iobrpy_cli.iobrpy.iobrpy_cli:main',
            'iobrpy-cli-mcp=iobrpy_cli.iobrpy.mcp_server:main',
        ],
    },

    install_requires=[
        "click>=8.0",
        "pandas>=1.5",
        "prompt-toolkit>=3.0",
        "iobrpy>=0.1.8",
    ],

    extras_require={
        'excel': ['openpyxl>=3.0'],
        'full': ['openpyxl>=3.0'],
    },

    python_requires=">=3.11",

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],

    keywords="bioinformatics, immuno-oncology, rna-seq, tme, deconvolution, cli, repl",
)

