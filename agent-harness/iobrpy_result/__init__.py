"""Packaged agent assets for IOBRpy result visualization and interpretation."""

from pathlib import Path


ASSET_ROOT = Path(__file__).resolve().parent / "agent_assets"


def packaged_skill_dir() -> Path:
    path = ASSET_ROOT / "iobrpy-result"
    if not path.exists():
        raise FileNotFoundError(f"Bundled iobrpy-result skill not found: {path}")
    return path


def packaged_plugin_dir() -> Path:
    path = ASSET_ROOT / "plugins" / "iobrpy-result"
    if not path.exists():
        raise FileNotFoundError(f"Bundled iobrpy-result plugin not found: {path}")
    return path


def packaged_claude_command_file() -> Path:
    path = ASSET_ROOT / "claude-code-command" / "iobrpy-result.md"
    if not path.exists():
        raise FileNotFoundError(f"Bundled iobrpy-result Claude command not found: {path}")
    return path
