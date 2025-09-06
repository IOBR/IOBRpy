# -*- coding: utf-8 -*-
"""
Bootstrap for the 'deside' subcommand with a clean, pinned venv **and**
a safe way to expose only the 'iobrpy' package (no full site-packages leakage).

Key fixes vs. previous versions:
- Do NOT prepend the outer environment's site-packages to PYTHONPATH.
- Instead, create a small shim directory inside the venv and symlink only the
  'iobrpy' package into it, so imports resolve to venv deps (NumPy 1.x etc.).

Env vars:
- IOBRPY_DESIDE_VENV: custom venv path (default: ~/.cache/iobrpy/deside-venv-py{major.minor})
- IOBRPY_DESIDE_REBUILD=1: force delete & rebuild the venv
- PIP_CONSTRAINT: optional pip constraints file to further pin/override
"""

from __future__ import annotations

import os
import sys
import subprocess
import shutil
from pathlib import Path
import venv

# Pins compatible with DeSide==1.3.2 (requires numba<0.57 and umap-learn==0.5.1)
PINNED = [
    "numpy>=1.22,<1.24",
    "llvmlite==0.39.1",
    "numba==0.56.4",
    "umap-learn==0.5.1",
    "scikit-learn==0.24.2",
    "tensorflow==2.11.1",
    "pandas==1.5.3",
    "DeSide==1.3.2",
]

MARKER_FILE = ".iobrpy_deside_ready"
LINKDIR_NAME = "_iobrpy_link"   # lives under the venv; will contain only a symlink to 'iobrpy'


def _default_venv_dir() -> Path:
    ovr = os.environ.get("IOBRPY_DESIDE_VENV")
    if ovr:
        return Path(ovr).expanduser().resolve()
    cache_home = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    pyver = f"py{sys.version_info.major}.{sys.version_info.minor}"
    return (cache_home / "iobrpy" / f"deside-venv-{pyver}").resolve()


def _venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _pip_env(venv_dir: Path) -> dict[str, str]:
    env = os.environ.copy()
    cache_dir = venv_dir / ".pip-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    env.setdefault("PIP_CACHE_DIR", str(cache_dir))
    env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
    env.setdefault("PYTHONNOUSERSITE", "1")
    return env


def _marker_path(venv_dir: Path) -> Path:
    return venv_dir / MARKER_FILE


def _need_install(venv_dir: Path) -> bool:
    return not _marker_path(venv_dir).exists()


def _write_marker(venv_dir: Path) -> None:
    _marker_path(venv_dir).write_text("ok", encoding="utf-8")


def _ensure_venv(venv_dir: Path, rebuild: bool) -> None:
    if rebuild and venv_dir.exists():
        shutil.rmtree(venv_dir, ignore_errors=True)
    if not venv_dir.exists():
        venv_dir.parent.mkdir(parents=True, exist_ok=True)
        venv.EnvBuilder(with_pip=True, clear=False, upgrade=False, symlinks=True).create(str(venv_dir))


def _pip_install(python: Path, venv_dir: Path, packages: list[str]) -> None:
    env = _pip_env(venv_dir)
    subprocess.run([str(python), "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"],
                   check=True, env=env)
    subprocess.run([str(python), "-m", "pip", "install",
                    "--upgrade", "--upgrade-strategy", "eager", "--force-reinstall", *packages],
                   check=True, env=env)


def _verify_numpy_range(python: Path, venv_dir: Path) -> None:
    code = (
        "import numpy,sys; v=numpy.__version__.split('+')[0]; "
        "maj,min=map(int,v.split('.')[:2]); "
        "ok=(maj==1 and 22<=min<24); "
        "print('[iobrpy] Using NumPy', v); sys.exit(0 if ok else 1)"
    )
    env = _pip_env(venv_dir)
    res = subprocess.run([str(python), "-c", code], env=env, capture_output=True, text=True)
    if res.returncode != 0:
        msg = (res.stdout or res.stderr).strip() or "unknown"
        raise SystemExit(f"[iobrpy] FATAL: NumPy must be >=1.22 and <1.24 (got {msg}). Aborting.")


def _setup_iobrpy_shim(venv_dir: Path) -> Path:
    """
    Create a tiny directory in the venv that contains a single entry: 'iobrpy' ->
    the source package directory (…/iobrpy). We'll add this shim directory to PYTHONPATH.
    This avoids importing anything else from the outer site-packages (like NumPy 2.x).
    """
    # __file__ is .../iobrpy/workflow/deside_bootstrap.py
    src_pkg = Path(__file__).resolve().parents[1]  # .../iobrpy
    if not (src_pkg / "__init__.py").exists():
        raise RuntimeError(f"Cannot locate iobrpy package at {src_pkg}")
    link_root = venv_dir / LINKDIR_NAME
    link_root.mkdir(parents=True, exist_ok=True)
    link = link_root / "iobrpy"
    if link.exists() or link.is_symlink():
        try:
            if link.is_dir() and not link.is_symlink():
                shutil.rmtree(link)
            else:
                link.unlink()
        except FileNotFoundError:
            pass
    try:
        link.symlink_to(src_pkg, target_is_directory=True)
    except (OSError, NotImplementedError):
        # Fallback: copytree (slower; larger, but robust where symlink is unavailable)
        shutil.copytree(src_pkg, link)
    return link_root


def run_in_isolated_env(argv: list[str] | None = None) -> None:
    venv_dir = _default_venv_dir()
    rebuild = os.environ.get("IOBRPY_DESIDE_REBUILD", "") == "1"
    _ensure_venv(venv_dir, rebuild=rebuild)

    py = _venv_python(venv_dir)

    first_time = _need_install(venv_dir)
    if first_time or rebuild:
        _pip_install(py, venv_dir, PINNED)
        _verify_numpy_range(py, venv_dir)
        _write_marker(venv_dir)

    # Prepare environment for the worker process
    env = os.environ.copy()
    env["IOBRPY_DESIDE_WORKER"] = "1"
    env.update(_pip_env(venv_dir))

    # Add only the shim directory (containing just 'iobrpy') to PYTHONPATH
    shim = _setup_iobrpy_shim(venv_dir)
    env["PYTHONPATH"] = f"{str(shim)}{os.pathsep}{env.get('PYTHONPATH','')}".rstrip(os.pathsep)

    if argv is None:
        argv = sys.argv[1:]

    cmd = [str(py), "-m", "iobrpy.workflow.deside", *map(str, argv)]
    if first_time:
        print("[iobrpy] deside environment is ready. Launching…", file=sys.stderr)
    subprocess.run(cmd, env=env, check=True)


def main() -> None:
    run_in_isolated_env(argv=None)


if __name__ == "__main__":
    main()
