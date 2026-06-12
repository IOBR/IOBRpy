"""
Helpers for installing persistent agent integrations for iobrpy-cli.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import tomllib

from iobrpy_result import (
    packaged_claude_command_file as packaged_result_claude_command_file,
    packaged_plugin_dir as packaged_result_plugin_dir,
    packaged_skill_dir as packaged_result_skill_dir,
)


SUPPORTED_AGENT_CLIENTS = ("codex", "claude-code")
CODEX_SKILL_NAME = "iobrpy"
CODEX_SKILL_SOURCE_DIRNAME = "iobrpy-fastpath"
LEGACY_CODEX_SKILL_NAME = "iobrpy-fastpath"
CODEX_PLUGIN_NAME = "iobrpy"
CODEX_PLUGIN_SOURCE_DIRNAME = "codex-plugin-iobrpy"
CLAUDE_MEMORY_DIRNAME = "iobrpy"
CLAUDE_COMMAND_NAME = "iobrpy"
CLAUDE_COMMAND_SOURCE_DIRNAME = "claude-code-command"
DEFAULT_MCP_SERVER_NAME = "iobrpy"
DEFAULT_CLAUDE_COMMAND = "claude"
RESULT_SKILL_NAME = "iobrpy-result"
RESULT_PLUGIN_NAME = "iobrpy-result"
RESULT_CLAUDE_COMMAND_NAME = "iobrpy-result"
CLAUDE_MEMORY_IMPORT_BLOCK = "\n".join(
    [
        "<!-- iobrpy-fastpath start -->",
        "@~/.claude/iobrpy/CLAUDE.md",
        "<!-- iobrpy-fastpath end -->",
    ]
)
RESULT_CLAUDE_MEMORY_IMPORT_BLOCK = "\n".join(
    [
        "<!-- iobrpy-result start -->",
        "@~/.claude/iobrpy-result/SKILL.md",
        "<!-- iobrpy-result end -->",
    ]
)


@dataclass(frozen=True)
class MCPServerSpec:
    name: str
    command: str
    args: List[str]
    env: Dict[str, str]

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "name": self.name,
            "command": self.command,
            "args": list(self.args),
        }
        if self.env:
            payload["env"] = dict(self.env)
        return payload


@dataclass(frozen=True)
class InstallRecord:
    client: str
    component: str
    status: str
    message: str
    target: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "client": self.client,
            "component": self.component,
            "status": self.status,
            "message": self.message,
        }
        if self.target:
            payload["target"] = self.target
        return payload


@dataclass(frozen=True)
class StatusRecord:
    client: str
    component: str
    status: str
    message: str
    target: Optional[str] = None
    details: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "client": self.client,
            "component": self.component,
            "status": self.status,
            "message": self.message,
        }
        if self.target:
            payload["target"] = self.target
        if self.details is not None:
            payload["details"] = self.details
        return payload


def _codex_home(override: Optional[Path] = None) -> Path:
    if override is not None:
        return override
    env_value = os.environ.get("CODEX_HOME")
    if env_value:
        return Path(env_value).expanduser()
    return Path.home() / ".codex"


def _claude_home(override: Optional[Path] = None) -> Path:
    if override is not None:
        return override
    return Path.home() / ".claude"


def packaged_codex_skill_dir() -> Path:
    skill_dir = Path(__file__).resolve().parent / "agent_assets" / CODEX_SKILL_SOURCE_DIRNAME
    if not skill_dir.exists():
        raise FileNotFoundError(f"Bundled Codex skill not found: {skill_dir}")
    return skill_dir


def packaged_codex_plugin_dir() -> Path:
    plugin_dir = Path(__file__).resolve().parent / "agent_assets" / CODEX_PLUGIN_SOURCE_DIRNAME
    if not plugin_dir.exists():
        raise FileNotFoundError(f"Bundled Codex plugin not found: {plugin_dir}")
    return plugin_dir


def packaged_claude_memory_file() -> Path:
    memory_file = Path(__file__).resolve().parent / "agent_assets" / "claude-code-memory" / "CLAUDE.md"
    if not memory_file.exists():
        raise FileNotFoundError(f"Bundled Claude Code memory file not found: {memory_file}")
    return memory_file


def packaged_claude_command_file() -> Path:
    command_file = Path(__file__).resolve().parent / "agent_assets" / CLAUDE_COMMAND_SOURCE_DIRNAME / f"{CLAUDE_COMMAND_NAME}.md"
    if not command_file.exists():
        raise FileNotFoundError(f"Bundled Claude Code command file not found: {command_file}")
    return command_file


def build_default_mcp_server_spec(server_name: str = DEFAULT_MCP_SERVER_NAME) -> MCPServerSpec:
    return MCPServerSpec(
        name=server_name,
        command=str(Path(sys.executable).resolve()),
        args=["-m", "iobrpy_cli.iobrpy.mcp_server"],
        env={},
    )


def _remove_existing_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
        return
    path.unlink()


def _directory_matches(source: Path, destination: Path) -> bool:
    if not source.is_dir() or not destination.is_dir():
        return False
    source_files = {
        path.relative_to(source)
        for path in source.rglob("*")
        if path.is_file()
    }
    destination_files = {
        path.relative_to(destination)
        for path in destination.rglob("*")
        if path.is_file()
    }
    if source_files != destination_files:
        return False
    return all(
        (source / relative).read_bytes() == (destination / relative).read_bytes()
        for relative in source_files
    )


def _codex_user_root(codex_home: Optional[Path] = None) -> Path:
    return _codex_home(codex_home).parent


def _codex_plugin_path(
    codex_home: Optional[Path] = None,
    plugin_name: str = CODEX_PLUGIN_NAME,
) -> Path:
    return _codex_user_root(codex_home) / "plugins" / plugin_name


def _codex_marketplace_path(codex_home: Optional[Path] = None) -> Path:
    return _codex_user_root(codex_home) / ".agents" / "plugins" / "marketplace.json"


def _expected_codex_marketplace_entry(
    plugin_name: str = CODEX_PLUGIN_NAME,
) -> Dict[str, Any]:
    return {
        "name": plugin_name,
        "source": {
            "source": "local",
            "path": f"./plugins/{plugin_name}",
        },
        "policy": {
            "installation": "AVAILABLE",
            "authentication": "ON_INSTALL",
        },
        "category": "Coding",
    }


def _load_marketplace_payload(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {
            "name": "local-plugins",
            "interface": {
                "displayName": "Local Plugins",
            },
            "plugins": [],
        }

    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Codex marketplace must contain a top-level object: {path}")
    if "plugins" not in payload or not isinstance(payload.get("plugins"), list):
        payload["plugins"] = []
    if "interface" not in payload or not isinstance(payload.get("interface"), dict):
        payload["interface"] = {"displayName": "Local Plugins"}
    if not payload.get("name"):
        payload["name"] = "local-plugins"
    if not payload["interface"].get("displayName"):
        payload["interface"]["displayName"] = "Local Plugins"
    return payload


def _upsert_codex_marketplace_entry(
    payload: Dict[str, Any],
    plugin_name: str = CODEX_PLUGIN_NAME,
) -> tuple[Dict[str, Any], str]:
    expected_entry = _expected_codex_marketplace_entry(plugin_name)
    plugins = payload.get("plugins", [])
    for index, entry in enumerate(plugins):
        if not isinstance(entry, dict):
            continue
        if entry.get("name") != plugin_name:
            continue
        if entry == expected_entry:
            return payload, "already_configured"
        plugins[index] = expected_entry
        payload["plugins"] = plugins
        return payload, "updated"

    plugins.append(expected_entry)
    payload["plugins"] = plugins
    return payload, "installed"


def install_codex_skill(
    *,
    codex_home: Optional[Path] = None,
    force: bool = False,
    dry_run: bool = False,
) -> InstallRecord:
    source = packaged_codex_skill_dir()
    destination = _codex_home(codex_home) / "skills" / CODEX_SKILL_NAME
    if destination.exists() and not force:
        return InstallRecord(
            client="codex",
            component="skill",
            status="already_installed",
            message="Codex skill already exists. Re-run with --force to overwrite it.",
            target=str(destination),
        )

    if dry_run:
        return InstallRecord(
            client="codex",
            component="skill",
            status="dry_run",
            message="Codex skill would be installed.",
            target=str(destination),
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        _remove_existing_path(destination)
        status = "updated"
        message = "Codex skill updated."
    else:
        status = "installed"
        message = "Codex skill installed."
    shutil.copytree(source, destination)
    return InstallRecord(
        client="codex",
        component="skill",
        status=status,
        message=message,
        target=str(destination),
    )


def codex_skill_status(*, codex_home: Optional[Path] = None) -> StatusRecord:
    source = packaged_codex_skill_dir()
    destination = _codex_home(codex_home) / "skills" / CODEX_SKILL_NAME
    legacy_destination = _codex_home(codex_home) / "skills" / LEGACY_CODEX_SKILL_NAME
    if not destination.exists():
        if legacy_destination.exists():
            return StatusRecord(
                client="codex",
                component="skill",
                status="mismatched",
                message="Legacy Codex skill `iobrpy-fastpath` exists, but the `/iobrpy` alias is not installed yet.",
                target=str(destination),
                details={
                    "packaged_source": str(source),
                    "legacy_destination": str(legacy_destination),
                },
            )
        return StatusRecord(
            client="codex",
            component="skill",
            status="not_installed",
            message="Codex skill is not installed.",
            target=str(destination),
            details={"packaged_source": str(source)},
        )

    skill_md = destination / "SKILL.md"
    openai_yaml = destination / "agents" / "openai.yaml"
    matches_packaged = (
        skill_md.exists()
        and openai_yaml.exists()
        and skill_md.read_text(encoding="utf-8") == (source / "SKILL.md").read_text(encoding="utf-8")
        and openai_yaml.read_text(encoding="utf-8") == (source / "agents" / "openai.yaml").read_text(encoding="utf-8")
    )
    return StatusRecord(
        client="codex",
        component="skill",
        status="installed",
        message="Codex skill is installed.",
        target=str(destination),
        details={
            "packaged_source": str(source),
            "matches_packaged": matches_packaged,
            "has_skill_md": skill_md.exists(),
            "has_openai_yaml": openai_yaml.exists(),
        },
    )


def install_codex_plugin(
    *,
    codex_home: Optional[Path] = None,
    force: bool = False,
    dry_run: bool = False,
) -> InstallRecord:
    source = packaged_codex_plugin_dir()
    destination = _codex_plugin_path(codex_home)
    marketplace_path = _codex_marketplace_path(codex_home)

    if dry_run:
        return InstallRecord(
            client="codex",
            component="plugin",
            status="dry_run",
            message="Codex plugin and marketplace entry would be installed.",
            target=str(destination),
        )

    plugin_status = "already_installed"
    if destination.exists():
        if force:
            _remove_existing_path(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination)
            plugin_status = "updated"
        else:
            manifest = destination / ".codex-plugin" / "plugin.json"
            if not manifest.exists():
                return InstallRecord(
                    client="codex",
                    component="plugin",
                    status="failed",
                    message="Existing Codex plugin path is missing `.codex-plugin/plugin.json`. Re-run with --force to replace it.",
                    target=str(destination),
                )
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, destination)
        plugin_status = "installed"

    try:
        payload = _load_marketplace_payload(marketplace_path)
        payload, marketplace_status = _upsert_codex_marketplace_entry(payload)
    except Exception as exc:
        return InstallRecord(
            client="codex",
            component="plugin",
            status="failed",
            message=f"Codex plugin marketplace could not be updated: {exc}",
            target=str(marketplace_path),
        )

    marketplace_path.parent.mkdir(parents=True, exist_ok=True)
    marketplace_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if plugin_status == "already_installed" and marketplace_status == "already_configured":
        status = "already_installed"
        message = "Codex plugin was already installed."
    elif "updated" in {plugin_status, marketplace_status}:
        status = "updated"
        message = "Codex plugin and marketplace entry updated."
    else:
        status = "installed"
        message = "Codex plugin installed and added to the local marketplace."
    return InstallRecord(
        client="codex",
        component="plugin",
        status=status,
        message=message,
        target=str(destination),
    )


def codex_plugin_status(*, codex_home: Optional[Path] = None) -> StatusRecord:
    source = packaged_codex_plugin_dir()
    destination = _codex_plugin_path(codex_home)
    marketplace_path = _codex_marketplace_path(codex_home)
    manifest = destination / ".codex-plugin" / "plugin.json"

    if not destination.exists():
        return StatusRecord(
            client="codex",
            component="plugin",
            status="not_installed",
            message="Codex plugin is not installed.",
            target=str(destination),
            details={
                "packaged_source": str(source),
                "marketplace_path": str(marketplace_path),
            },
        )
    if not manifest.exists():
        return StatusRecord(
            client="codex",
            component="plugin",
            status="mismatched",
            message="Codex plugin directory exists, but `.codex-plugin/plugin.json` is missing.",
            target=str(destination),
            details={
                "packaged_source": str(source),
                "marketplace_path": str(marketplace_path),
            },
        )
    if not marketplace_path.exists():
        return StatusRecord(
            client="codex",
            component="plugin",
            status="mismatched",
            message="Codex plugin files exist, but the local marketplace entry is missing.",
            target=str(destination),
            details={
                "marketplace_path": str(marketplace_path),
            },
        )

    try:
        payload = _load_marketplace_payload(marketplace_path)
    except Exception as exc:
        return StatusRecord(
            client="codex",
            component="plugin",
            status="invalid_config",
            message=f"Codex marketplace could not be parsed: {exc}",
            target=str(marketplace_path),
        )

    expected_entry = _expected_codex_marketplace_entry()
    marketplace_entry = next(
        (
            entry for entry in payload.get("plugins", [])
            if isinstance(entry, dict) and entry.get("name") == CODEX_PLUGIN_NAME
        ),
        None,
    )
    packaged_manifest = (source / ".codex-plugin" / "plugin.json").read_text(encoding="utf-8")
    installed_manifest = manifest.read_text(encoding="utf-8")
    if marketplace_entry == expected_entry and installed_manifest == packaged_manifest:
        return StatusRecord(
            client="codex",
            component="plugin",
            status="installed",
            message="Codex plugin is installed and registered in the local marketplace.",
            target=str(destination),
            details={
                "marketplace_path": str(marketplace_path),
                "marketplace_entry": marketplace_entry,
            },
        )
    return StatusRecord(
        client="codex",
        component="plugin",
        status="mismatched",
        message="Codex plugin exists, but the installed files or marketplace entry do not match the expected iobrpy bundle.",
        target=str(destination),
        details={
            "marketplace_path": str(marketplace_path),
            "marketplace_entry": marketplace_entry,
            "expected_entry": expected_entry,
        },
    )


def install_result_codex_skill(
    *,
    codex_home: Optional[Path] = None,
    force: bool = False,
    dry_run: bool = False,
) -> InstallRecord:
    source = packaged_result_skill_dir()
    destination = _codex_home(codex_home) / "skills" / RESULT_SKILL_NAME
    if destination.exists() and not force:
        status = "already_installed" if _directory_matches(source, destination) else "mismatched"
        return InstallRecord(
            client="codex",
            component="result_skill",
            status=status,
            message=(
                "IOBRpy result skill was already installed."
                if status == "already_installed"
                else "IOBRpy result skill already exists and differs from the packaged version. Re-run with --force to replace it."
            ),
            target=str(destination),
        )

    if dry_run:
        return InstallRecord(
            client="codex",
            component="result_skill",
            status="dry_run",
            message="IOBRpy result skill would be installed.",
            target=str(destination),
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    existed = destination.exists()
    if existed:
        _remove_existing_path(destination)
    shutil.copytree(source, destination)
    return InstallRecord(
        client="codex",
        component="result_skill",
        status="updated" if existed else "installed",
        message="IOBRpy result skill updated." if existed else "IOBRpy result skill installed.",
        target=str(destination),
    )


def result_codex_skill_status(*, codex_home: Optional[Path] = None) -> StatusRecord:
    source = packaged_result_skill_dir()
    destination = _codex_home(codex_home) / "skills" / RESULT_SKILL_NAME
    if not destination.exists():
        return StatusRecord(
            client="codex",
            component="result_skill",
            status="not_installed",
            message="IOBRpy result skill is not installed.",
            target=str(destination),
            details={"packaged_source": str(source)},
        )
    matches = _directory_matches(source, destination)
    return StatusRecord(
        client="codex",
        component="result_skill",
        status="installed" if matches else "mismatched",
        message=(
            "IOBRpy result skill is installed."
            if matches
            else "IOBRpy result skill exists but does not match the packaged version."
        ),
        target=str(destination),
        details={"packaged_source": str(source), "matches_packaged": matches},
    )


def _copy_result_plugin(source: Path, skill_source: Path, destination: Path) -> None:
    shutil.copytree(source, destination)
    installed_skill = destination / "skills" / RESULT_SKILL_NAME
    installed_skill.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(skill_source, installed_skill)


def install_result_codex_plugin(
    *,
    codex_home: Optional[Path] = None,
    force: bool = False,
    dry_run: bool = False,
) -> InstallRecord:
    source = packaged_result_plugin_dir()
    skill_source = packaged_result_skill_dir()
    destination = _codex_plugin_path(codex_home, RESULT_PLUGIN_NAME)
    marketplace_path = _codex_marketplace_path(codex_home)

    if dry_run:
        return InstallRecord(
            client="codex",
            component="result_plugin",
            status="dry_run",
            message="IOBRpy result plugin and marketplace entry would be installed.",
            target=str(destination),
        )

    plugin_status = "already_installed"
    if destination.exists():
        if force:
            _remove_existing_path(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            _copy_result_plugin(source, skill_source, destination)
            plugin_status = "updated"
        else:
            manifest = destination / ".codex-plugin" / "plugin.json"
            installed_skill = destination / "skills" / RESULT_SKILL_NAME
            manifest_matches = (
                manifest.exists()
                and manifest.read_bytes() == (source / ".codex-plugin" / "plugin.json").read_bytes()
            )
            if not manifest_matches or not _directory_matches(skill_source, installed_skill):
                return InstallRecord(
                    client="codex",
                    component="result_plugin",
                    status="failed",
                    message="Existing IOBRpy result plugin is incomplete or stale. Re-run with --force to replace it.",
                    target=str(destination),
                )
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        _copy_result_plugin(source, skill_source, destination)
        plugin_status = "installed"

    try:
        payload = _load_marketplace_payload(marketplace_path)
        payload, marketplace_status = _upsert_codex_marketplace_entry(payload, RESULT_PLUGIN_NAME)
    except Exception as exc:
        return InstallRecord(
            client="codex",
            component="result_plugin",
            status="failed",
            message=f"IOBRpy result plugin marketplace could not be updated: {exc}",
            target=str(marketplace_path),
        )

    marketplace_path.parent.mkdir(parents=True, exist_ok=True)
    marketplace_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if plugin_status == "already_installed" and marketplace_status == "already_configured":
        status = "already_installed"
        message = "IOBRpy result plugin was already installed."
    elif "updated" in {plugin_status, marketplace_status}:
        status = "updated"
        message = "IOBRpy result plugin and marketplace entry updated."
    else:
        status = "installed"
        message = "IOBRpy result plugin installed and added to the local marketplace."
    return InstallRecord(
        client="codex",
        component="result_plugin",
        status=status,
        message=message,
        target=str(destination),
    )


def result_codex_plugin_status(*, codex_home: Optional[Path] = None) -> StatusRecord:
    source = packaged_result_plugin_dir()
    skill_source = packaged_result_skill_dir()
    destination = _codex_plugin_path(codex_home, RESULT_PLUGIN_NAME)
    marketplace_path = _codex_marketplace_path(codex_home)
    manifest = destination / ".codex-plugin" / "plugin.json"
    installed_skill = destination / "skills" / RESULT_SKILL_NAME
    if not destination.exists():
        return StatusRecord(
            client="codex",
            component="result_plugin",
            status="not_installed",
            message="IOBRpy result plugin is not installed.",
            target=str(destination),
            details={"packaged_source": str(source), "marketplace_path": str(marketplace_path)},
        )
    if not manifest.exists() or not installed_skill.exists():
        return StatusRecord(
            client="codex",
            component="result_plugin",
            status="mismatched",
            message="IOBRpy result plugin directory is incomplete.",
            target=str(destination),
        )
    if not marketplace_path.exists():
        return StatusRecord(
            client="codex",
            component="result_plugin",
            status="mismatched",
            message="IOBRpy result plugin files exist, but the local marketplace entry is missing.",
            target=str(destination),
        )
    try:
        payload = _load_marketplace_payload(marketplace_path)
    except Exception as exc:
        return StatusRecord(
            client="codex",
            component="result_plugin",
            status="invalid_config",
            message=f"Codex marketplace could not be parsed: {exc}",
            target=str(marketplace_path),
        )
    expected_entry = _expected_codex_marketplace_entry(RESULT_PLUGIN_NAME)
    marketplace_entry = next(
        (
            entry for entry in payload.get("plugins", [])
            if isinstance(entry, dict) and entry.get("name") == RESULT_PLUGIN_NAME
        ),
        None,
    )
    manifest_matches = manifest.read_bytes() == (source / ".codex-plugin" / "plugin.json").read_bytes()
    skill_matches = _directory_matches(skill_source, installed_skill)
    matches = marketplace_entry == expected_entry and manifest_matches and skill_matches
    return StatusRecord(
        client="codex",
        component="result_plugin",
        status="installed" if matches else "mismatched",
        message=(
            "IOBRpy result plugin is installed and registered in the local marketplace."
            if matches
            else "IOBRpy result plugin exists, but its files or marketplace entry do not match the packaged bundle."
        ),
        target=str(destination),
        details={
            "marketplace_path": str(marketplace_path),
            "marketplace_entry": marketplace_entry,
            "manifest_matches_packaged": manifest_matches,
            "skill_matches_packaged": skill_matches,
        },
    )


def _claude_user_memory_path(claude_home: Optional[Path] = None) -> Path:
    return _claude_home(claude_home) / "CLAUDE.md"


def _claude_managed_memory_path(claude_home: Optional[Path] = None) -> Path:
    return _claude_home(claude_home) / CLAUDE_MEMORY_DIRNAME / "CLAUDE.md"


def _claude_user_command_path(claude_home: Optional[Path] = None) -> Path:
    return _claude_home(claude_home) / "commands" / f"{CLAUDE_COMMAND_NAME}.md"


def _ensure_claude_import_block(text: str) -> tuple[str, str]:
    normalized = text if text.endswith("\n") or not text else text + "\n"
    if CLAUDE_MEMORY_IMPORT_BLOCK in normalized:
        return normalized, "already_configured"

    stripped = normalized.rstrip()
    if stripped:
        stripped += "\n\n"
    stripped += CLAUDE_MEMORY_IMPORT_BLOCK + "\n"
    return stripped, "installed"


def install_claude_code_skill(
    *,
    claude_home: Optional[Path] = None,
    dry_run: bool = False,
) -> InstallRecord:
    packaged = packaged_claude_memory_file()
    user_memory = _claude_user_memory_path(claude_home)
    managed_memory = _claude_managed_memory_path(claude_home)
    existing_main = user_memory.read_text(encoding="utf-8") if user_memory.exists() else ""
    updated_main, import_status = _ensure_claude_import_block(existing_main)
    packaged_contents = packaged.read_text(encoding="utf-8")
    managed_exists = managed_memory.exists()
    managed_matches = managed_exists and managed_memory.read_text(encoding="utf-8") == packaged_contents

    if dry_run:
        return InstallRecord(
            client="claude-code",
            component="skill",
            status="dry_run",
            message="Claude Code user memory would be installed or refreshed.",
            target=str(user_memory),
        )

    managed_memory.parent.mkdir(parents=True, exist_ok=True)
    managed_memory.write_text(packaged_contents, encoding="utf-8")
    user_memory.parent.mkdir(parents=True, exist_ok=True)
    user_memory.write_text(updated_main, encoding="utf-8")

    if managed_matches and import_status == "already_configured":
        status = "already_installed"
        message = "Claude Code user memory was already configured."
    elif managed_exists:
        status = "updated"
        message = "Claude Code user memory updated."
    else:
        status = "installed"
        message = "Claude Code user memory installed."
    return InstallRecord(
        client="claude-code",
        component="skill",
        status=status,
        message=message,
        target=str(user_memory),
    )


def claude_code_skill_status(*, claude_home: Optional[Path] = None) -> StatusRecord:
    packaged = packaged_claude_memory_file()
    user_memory = _claude_user_memory_path(claude_home)
    managed_memory = _claude_managed_memory_path(claude_home)

    if not user_memory.exists() and not managed_memory.exists():
        return StatusRecord(
            client="claude-code",
            component="skill",
            status="not_installed",
            message="Claude Code user memory is not installed.",
            target=str(user_memory),
            details={"managed_memory": str(managed_memory), "packaged_source": str(packaged)},
        )
    if not user_memory.exists():
        return StatusRecord(
            client="claude-code",
            component="skill",
            status="mismatched",
            message="Claude Code managed memory exists, but ~/.claude/CLAUDE.md is missing.",
            target=str(user_memory),
            details={"managed_memory": str(managed_memory), "packaged_source": str(packaged)},
        )

    main_text = user_memory.read_text(encoding="utf-8")
    import_present = CLAUDE_MEMORY_IMPORT_BLOCK in main_text
    if not managed_memory.exists():
        return StatusRecord(
            client="claude-code",
            component="skill",
            status="mismatched",
            message="Claude Code user memory import exists, but the managed iobrpy memory file is missing.",
            target=str(user_memory),
            details={"managed_memory": str(managed_memory), "import_present": import_present},
        )

    managed_matches = managed_memory.read_text(encoding="utf-8") == packaged.read_text(encoding="utf-8")
    if import_present and managed_matches:
        return StatusRecord(
            client="claude-code",
            component="skill",
            status="installed",
            message="Claude Code user memory is installed.",
            target=str(user_memory),
            details={"managed_memory": str(managed_memory), "packaged_source": str(packaged)},
        )
    return StatusRecord(
        client="claude-code",
        component="skill",
        status="mismatched",
        message="Claude Code user memory exists but does not match the expected iobrpy setup.",
        target=str(user_memory),
        details={
            "managed_memory": str(managed_memory),
            "packaged_source": str(packaged),
            "import_present": import_present,
            "managed_matches_packaged": managed_matches,
        },
    )


def install_claude_code_command(
    *,
    claude_home: Optional[Path] = None,
    dry_run: bool = False,
) -> InstallRecord:
    packaged = packaged_claude_command_file()
    destination = _claude_user_command_path(claude_home)
    packaged_contents = packaged.read_text(encoding="utf-8")
    existing_contents = destination.read_text(encoding="utf-8") if destination.exists() else None

    if dry_run:
        return InstallRecord(
            client="claude-code",
            component="slash_command",
            status="dry_run",
            message="Claude Code `/iobrpy` command would be installed.",
            target=str(destination),
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(packaged_contents, encoding="utf-8")

    if existing_contents is None:
        status = "installed"
        message = "Claude Code `/iobrpy` command installed."
    elif existing_contents == packaged_contents:
        status = "already_installed"
        message = "Claude Code `/iobrpy` command was already installed."
    else:
        status = "updated"
        message = "Claude Code `/iobrpy` command updated."
    return InstallRecord(
        client="claude-code",
        component="slash_command",
        status=status,
        message=message,
        target=str(destination),
    )


def claude_code_command_status(*, claude_home: Optional[Path] = None) -> StatusRecord:
    packaged = packaged_claude_command_file()
    destination = _claude_user_command_path(claude_home)
    if not destination.exists():
        return StatusRecord(
            client="claude-code",
            component="slash_command",
            status="not_installed",
            message="Claude Code `/iobrpy` command is not installed.",
            target=str(destination),
            details={"packaged_source": str(packaged)},
        )

    installed_contents = destination.read_text(encoding="utf-8")
    packaged_contents = packaged.read_text(encoding="utf-8")
    if installed_contents == packaged_contents:
        return StatusRecord(
            client="claude-code",
            component="slash_command",
            status="installed",
            message="Claude Code `/iobrpy` command is installed.",
            target=str(destination),
            details={"packaged_source": str(packaged)},
        )
    return StatusRecord(
        client="claude-code",
        component="slash_command",
        status="mismatched",
        message="Claude Code `/iobrpy` command exists but does not match the expected iobrpy command file.",
        target=str(destination),
        details={"packaged_source": str(packaged)},
    )


def _claude_result_skill_path(claude_home: Optional[Path] = None) -> Path:
    return _claude_home(claude_home) / RESULT_SKILL_NAME


def _claude_result_command_path(claude_home: Optional[Path] = None) -> Path:
    return _claude_home(claude_home) / "commands" / f"{RESULT_CLAUDE_COMMAND_NAME}.md"


def _ensure_result_claude_import_block(text: str) -> tuple[str, str]:
    normalized = text if text.endswith("\n") or not text else text + "\n"
    if RESULT_CLAUDE_MEMORY_IMPORT_BLOCK in normalized:
        return normalized, "already_configured"
    stripped = normalized.rstrip()
    if stripped:
        stripped += "\n\n"
    stripped += RESULT_CLAUDE_MEMORY_IMPORT_BLOCK + "\n"
    return stripped, "installed"


def install_claude_code_result_skill(
    *,
    claude_home: Optional[Path] = None,
    dry_run: bool = False,
) -> InstallRecord:
    source = packaged_result_skill_dir()
    user_memory = _claude_user_memory_path(claude_home)
    destination = _claude_result_skill_path(claude_home)
    existing_main = user_memory.read_text(encoding="utf-8") if user_memory.exists() else ""
    updated_main, import_status = _ensure_result_claude_import_block(existing_main)
    destination_matches = _directory_matches(source, destination)

    if dry_run:
        return InstallRecord(
            client="claude-code",
            component="result_skill",
            status="dry_run",
            message="Claude Code IOBRpy result skill would be installed or refreshed.",
            target=str(destination),
        )

    existed = destination.exists()
    if existed:
        _remove_existing_path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination)
    user_memory.parent.mkdir(parents=True, exist_ok=True)
    user_memory.write_text(updated_main, encoding="utf-8")

    if destination_matches and import_status == "already_configured":
        status = "already_installed"
        message = "Claude Code IOBRpy result skill was already configured."
    elif existed:
        status = "updated"
        message = "Claude Code IOBRpy result skill updated."
    else:
        status = "installed"
        message = "Claude Code IOBRpy result skill installed."
    return InstallRecord(
        client="claude-code",
        component="result_skill",
        status=status,
        message=message,
        target=str(destination),
    )


def claude_code_result_skill_status(
    *,
    claude_home: Optional[Path] = None,
) -> StatusRecord:
    source = packaged_result_skill_dir()
    user_memory = _claude_user_memory_path(claude_home)
    destination = _claude_result_skill_path(claude_home)
    if not user_memory.exists() and not destination.exists():
        return StatusRecord(
            client="claude-code",
            component="result_skill",
            status="not_installed",
            message="Claude Code IOBRpy result skill is not installed.",
            target=str(destination),
            details={"packaged_source": str(source)},
        )
    import_present = (
        user_memory.exists()
        and RESULT_CLAUDE_MEMORY_IMPORT_BLOCK in user_memory.read_text(encoding="utf-8")
    )
    skill_matches = _directory_matches(source, destination)
    matches = import_present and skill_matches
    return StatusRecord(
        client="claude-code",
        component="result_skill",
        status="installed" if matches else "mismatched",
        message=(
            "Claude Code IOBRpy result skill is installed."
            if matches
            else "Claude Code IOBRpy result skill or its CLAUDE.md import does not match the packaged setup."
        ),
        target=str(destination),
        details={
            "packaged_source": str(source),
            "import_present": import_present,
            "matches_packaged": skill_matches,
        },
    )


def install_claude_code_result_command(
    *,
    claude_home: Optional[Path] = None,
    dry_run: bool = False,
) -> InstallRecord:
    packaged = packaged_result_claude_command_file()
    destination = _claude_result_command_path(claude_home)
    packaged_contents = packaged.read_text(encoding="utf-8")
    existing_contents = destination.read_text(encoding="utf-8") if destination.exists() else None
    if dry_run:
        return InstallRecord(
            client="claude-code",
            component="result_slash_command",
            status="dry_run",
            message="Claude Code `/iobrpy-result` command would be installed.",
            target=str(destination),
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(packaged_contents, encoding="utf-8")
    if existing_contents is None:
        status = "installed"
        message = "Claude Code `/iobrpy-result` command installed."
    elif existing_contents == packaged_contents:
        status = "already_installed"
        message = "Claude Code `/iobrpy-result` command was already installed."
    else:
        status = "updated"
        message = "Claude Code `/iobrpy-result` command updated."
    return InstallRecord(
        client="claude-code",
        component="result_slash_command",
        status=status,
        message=message,
        target=str(destination),
    )


def claude_code_result_command_status(
    *,
    claude_home: Optional[Path] = None,
) -> StatusRecord:
    packaged = packaged_result_claude_command_file()
    destination = _claude_result_command_path(claude_home)
    if not destination.exists():
        return StatusRecord(
            client="claude-code",
            component="result_slash_command",
            status="not_installed",
            message="Claude Code `/iobrpy-result` command is not installed.",
            target=str(destination),
            details={"packaged_source": str(packaged)},
        )
    matches = destination.read_text(encoding="utf-8") == packaged.read_text(encoding="utf-8")
    return StatusRecord(
        client="claude-code",
        component="result_slash_command",
        status="installed" if matches else "mismatched",
        message=(
            "Claude Code `/iobrpy-result` command is installed."
            if matches
            else "Claude Code `/iobrpy-result` command does not match the packaged version."
        ),
        target=str(destination),
        details={"packaged_source": str(packaged)},
    )


def _toml_literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _codex_section_headers(server_name: str) -> List[str]:
    return [
        f"[mcp_servers.{server_name}]",
        f'[mcp_servers."{server_name}"]',
    ]


def _render_codex_mcp_block(spec: MCPServerSpec) -> str:
    lines = [
        f'[mcp_servers."{spec.name}"]',
        f"command = {_toml_literal(spec.command)}",
        "args = [" + ", ".join(_toml_literal(arg) for arg in spec.args) + "]",
    ]
    if spec.env:
        env_pairs = ", ".join(
            f"{key} = {_toml_literal(value)}"
            for key, value in sorted(spec.env.items())
        )
        lines.append(f"env = {{ {env_pairs} }}")
    return "\n".join(lines) + "\n"


def _read_toml_file(path: Path) -> Dict[str, Any]:
    with path.open("rb") as handle:
        data = tomllib.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"TOML config must contain a table at the top level: {path}")
    return data


def _mcp_entry_matches_spec(entry: Any, spec: MCPServerSpec) -> bool:
    if not isinstance(entry, dict):
        return False

    command = entry.get("command")
    args = entry.get("args")
    env = entry.get("env", {})

    if command != spec.command:
        return False
    if not isinstance(args, list):
        return False
    if list(args) != list(spec.args):
        return False
    if spec.env:
        if not isinstance(env, dict):
            return False
        return dict(env) == dict(spec.env)
    return env in ({}, None)


def _replace_or_append_toml_section(text: str, headers: Sequence[str], new_block: str) -> tuple[str, str]:
    lines = text.splitlines(keepends=True)
    start: Optional[int] = None
    end: Optional[int] = None

    for index, line in enumerate(lines):
        if line.strip() in headers:
            start = index
            break

    if start is None:
        body = text.rstrip()
        if body:
            body += "\n\n"
        body += new_block.rstrip() + "\n"
        return body, "installed"

    for index in range(start + 1, len(lines)):
        stripped = lines[index].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            end = index
            break
    if end is None:
        end = len(lines)

    existing = "".join(lines[start:end]).strip()
    if existing == new_block.strip():
        return text if text.endswith("\n") or not text else text + "\n", "already_configured"

    updated = "".join(lines[:start]) + new_block
    if start > 0 and not updated.endswith("\n"):
        updated += "\n"
    updated += "".join(lines[end:])
    updated = updated.rstrip() + "\n"
    return updated, "updated"


def install_codex_mcp(
    spec: MCPServerSpec,
    *,
    codex_home: Optional[Path] = None,
    dry_run: bool = False,
) -> InstallRecord:
    config_path = _codex_home(codex_home) / "config.toml"
    block = _render_codex_mcp_block(spec)
    existing = config_path.read_text(encoding="utf-8") if config_path.exists() else ""
    updated, status = _replace_or_append_toml_section(existing, _codex_section_headers(spec.name), block)

    if dry_run:
        return InstallRecord(
            client="codex",
            component="mcp",
            status="dry_run",
            message="Codex MCP server would be configured.",
            target=str(config_path),
        )

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(updated, encoding="utf-8")
    message = {
        "installed": "Codex MCP server configured.",
        "updated": "Codex MCP server updated.",
        "already_configured": "Codex MCP server was already configured.",
    }[status]
    return InstallRecord(
        client="codex",
        component="mcp",
        status=status,
        message=message,
        target=str(config_path),
    )


def codex_mcp_status(
    spec: MCPServerSpec,
    *,
    codex_home: Optional[Path] = None,
) -> StatusRecord:
    config_path = _codex_home(codex_home) / "config.toml"
    if not config_path.exists():
        return StatusRecord(
            client="codex",
            component="mcp",
            status="not_configured",
            message="Codex MCP config file does not exist.",
            target=str(config_path),
        )

    try:
        data = _read_toml_file(config_path)
    except Exception as exc:
        return StatusRecord(
            client="codex",
            component="mcp",
            status="invalid_config",
            message=f"Codex config could not be parsed: {exc}",
            target=str(config_path),
        )

    servers = data.get("mcp_servers", {})
    if not isinstance(servers, dict):
        return StatusRecord(
            client="codex",
            component="mcp",
            status="invalid_config",
            message="Codex config has an invalid mcp_servers table.",
            target=str(config_path),
        )

    entry = servers.get(spec.name)
    if entry is None:
        return StatusRecord(
            client="codex",
            component="mcp",
            status="not_configured",
            message="Codex MCP server is not configured.",
            target=str(config_path),
        )
    if _mcp_entry_matches_spec(entry, spec):
        return StatusRecord(
            client="codex",
            component="mcp",
            status="configured",
            message="Codex MCP server matches the expected iobrpy configuration.",
            target=str(config_path),
            details={"configured_entry": entry},
        )
    return StatusRecord(
        client="codex",
        component="mcp",
        status="mismatched",
        message="Codex MCP server exists but does not match the expected iobrpy configuration.",
        target=str(config_path),
        details={"configured_entry": entry, "expected_entry": spec.to_dict()},
    )


def install_claude_code_mcp(
    spec: MCPServerSpec,
    *,
    claude_command: str = DEFAULT_CLAUDE_COMMAND,
    dry_run: bool = False,
    allow_missing_command: bool = False,
) -> InstallRecord:
    executable = shutil.which(claude_command)
    if executable is None:
        status = "skipped" if allow_missing_command else "failed"
        return InstallRecord(
            client="claude-code",
            component="mcp",
            status=status,
            message=(
                "Claude Code CLI was not found on PATH."
                if allow_missing_command
                else "Claude Code CLI was not found on PATH. Install it or pass --client codex only."
            ),
        )

    argv = [executable, "mcp", "add", spec.name, "--scope", "user"]
    for key, value in sorted(spec.env.items()):
        argv.extend(["--env", f"{key}={value}"])
    argv.extend(["--", spec.command, *spec.args])

    if dry_run:
        return InstallRecord(
            client="claude-code",
            component="mcp",
            status="dry_run",
            message="Claude Code MCP server would be configured.",
            target=" ".join(argv),
        )

    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        details = (completed.stderr or completed.stdout or "").strip()
        if "already exists" in details.lower():
            return InstallRecord(
                client="claude-code",
                component="mcp",
                status="already_configured",
                message="Claude Code MCP server is already configured.",
                target=" ".join(argv),
            )
        message = "Claude Code MCP install failed."
        if details:
            message = f"{message} {details}"
        return InstallRecord(
            client="claude-code",
            component="mcp",
            status="failed",
            message=message,
            target=" ".join(argv),
        )

    return InstallRecord(
        client="claude-code",
        component="mcp",
        status="installed",
        message="Claude Code MCP server configured.",
        target=" ".join(argv),
    )


def claude_code_mcp_status(
    spec: MCPServerSpec,
    *,
    claude_command: str = DEFAULT_CLAUDE_COMMAND,
) -> StatusRecord:
    executable = shutil.which(claude_command)
    if executable is None:
        return StatusRecord(
            client="claude-code",
            component="mcp",
            status="cli_not_found",
            message="Claude Code CLI was not found on PATH.",
        )

    argv = [executable, "mcp", "list"]
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        details = (completed.stderr or completed.stdout or "").strip()
        message = "Claude Code CLI could not list MCP servers."
        if details:
            message = f"{message} {details}"
        return StatusRecord(
            client="claude-code",
            component="mcp",
            status="error",
            message=message,
            target=" ".join(argv),
        )

    output = completed.stdout or ""
    pattern = re.compile(rf"(?<![\w-]){re.escape(spec.name)}(?![\w-])", re.IGNORECASE)
    configured = bool(pattern.search(output))
    return StatusRecord(
        client="claude-code",
        component="mcp",
        status="configured" if configured else "not_configured",
        message=(
            "Claude Code MCP server appears in `claude mcp list`."
            if configured
            else "Claude Code CLI is available, but the iobrpy MCP server was not found in `claude mcp list`."
        ),
        target=" ".join(argv),
        details={"list_output": output.strip()},
    )


def _normalize_clients(clients: Sequence[str]) -> List[str]:
    if not clients:
        return ["codex"]
    ordered: List[str] = []
    for client in clients:
        if client == "all":
            for expanded in SUPPORTED_AGENT_CLIENTS:
                if expanded not in ordered:
                    ordered.append(expanded)
            continue
        if client not in ordered:
            ordered.append(client)
    return ordered


def install_agent_bundle(
    clients: Sequence[str],
    *,
    include_skill: bool = True,
    include_mcp: bool = True,
    force: bool = False,
    dry_run: bool = False,
    server_name: str = DEFAULT_MCP_SERVER_NAME,
    codex_home: Optional[Path] = None,
    claude_home: Optional[Path] = None,
    claude_command: str = DEFAULT_CLAUDE_COMMAND,
) -> Dict[str, Any]:
    resolved_clients = _normalize_clients(clients)
    requested_all = "all" in clients
    spec = build_default_mcp_server_spec(server_name)
    records: List[InstallRecord] = []

    for client in resolved_clients:
        if client == "codex":
            if include_skill:
                records.append(
                    install_codex_skill(codex_home=codex_home, force=force, dry_run=dry_run)
                )
                records.append(
                    install_codex_plugin(codex_home=codex_home, force=force, dry_run=dry_run)
                )
                records.append(
                    install_result_codex_skill(codex_home=codex_home, force=force, dry_run=dry_run)
                )
                records.append(
                    install_result_codex_plugin(codex_home=codex_home, force=force, dry_run=dry_run)
                )
            if include_mcp:
                records.append(
                    install_codex_mcp(spec, codex_home=codex_home, dry_run=dry_run)
                )
        elif client == "claude-code":
            if include_skill:
                records.append(
                    install_claude_code_skill(claude_home=claude_home, dry_run=dry_run)
                )
                records.append(
                    install_claude_code_command(claude_home=claude_home, dry_run=dry_run)
                )
                records.append(
                    install_claude_code_result_skill(claude_home=claude_home, dry_run=dry_run)
                )
                records.append(
                    install_claude_code_result_command(claude_home=claude_home, dry_run=dry_run)
                )
            if include_mcp:
                records.append(
                    install_claude_code_mcp(
                        spec,
                        claude_command=claude_command,
                        dry_run=dry_run,
                        allow_missing_command=requested_all,
                    )
                )
        else:
            raise ValueError(f"Unsupported agent client: {client}")

    success = all(record.status not in {"failed", "mismatched"} for record in records)
    next_steps: List[str] = []
    if any(record.client == "codex" and record.component == "skill" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Codex to pick up the `/iobrpy` skill alias.")
    if any(record.client == "codex" and record.component == "plugin" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Codex or reopen the plugin picker so the local `iobrpy` plugin and its namespace commands are visible.")
    if any(record.client == "codex" and record.component == "result_skill" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Codex to pick up the `/iobrpy-result` skill.")
    if any(record.client == "codex" and record.component == "result_plugin" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Codex or reopen the plugin picker so the local `iobrpy-result` plugin is visible.")
    if any(record.client == "codex" and record.component == "mcp" and record.status in {"installed", "updated", "already_configured"} for record in records):
        next_steps.append("Restart Codex after MCP changes so the client reloads its server list.")
    if any(record.client == "claude-code" and record.component == "skill" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Claude Code or start a new session so the updated ~/.claude/CLAUDE.md memory is loaded.")
    if any(record.client == "claude-code" and record.component == "slash_command" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Claude Code or start a new session so `/iobrpy` appears in the slash-command list.")
    if any(record.client == "claude-code" and record.component == "result_skill" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Claude Code or start a new session so the IOBRpy result skill is loaded.")
    if any(record.client == "claude-code" and record.component == "result_slash_command" and record.status in {"installed", "updated", "already_installed"} for record in records):
        next_steps.append("Restart Claude Code or start a new session so `/iobrpy-result` appears in the slash-command list.")
    if any(record.client == "claude-code" and record.component == "mcp" and record.status in {"installed", "already_configured"} for record in records):
        next_steps.append("Run `claude mcp list` or restart Claude Code to confirm the user-scoped server is visible.")

    return {
        "success": success,
        "server": spec.to_dict(),
        "requested_clients": list(clients) if clients else ["codex"],
        "resolved_clients": resolved_clients,
        "components": {
            "skill": include_skill,
            "result_skill": include_skill,
            "result_plugin": include_skill,
            "mcp": include_mcp,
        },
        "results": [record.to_dict() for record in records],
        "next_steps": next_steps,
    }


def agent_status_bundle(
    clients: Sequence[str],
    *,
    include_skill: bool = True,
    include_mcp: bool = True,
    server_name: str = DEFAULT_MCP_SERVER_NAME,
    codex_home: Optional[Path] = None,
    claude_home: Optional[Path] = None,
    claude_command: str = DEFAULT_CLAUDE_COMMAND,
) -> Dict[str, Any]:
    resolved_clients = _normalize_clients(clients or ["all"])
    spec = build_default_mcp_server_spec(server_name)
    records: List[StatusRecord] = []

    for client in resolved_clients:
        if client == "codex":
            if include_skill:
                records.append(codex_skill_status(codex_home=codex_home))
                records.append(codex_plugin_status(codex_home=codex_home))
                records.append(result_codex_skill_status(codex_home=codex_home))
                records.append(result_codex_plugin_status(codex_home=codex_home))
            if include_mcp:
                records.append(codex_mcp_status(spec, codex_home=codex_home))
        elif client == "claude-code":
            if include_skill:
                records.append(claude_code_skill_status(claude_home=claude_home))
                records.append(claude_code_command_status(claude_home=claude_home))
                records.append(claude_code_result_skill_status(claude_home=claude_home))
                records.append(claude_code_result_command_status(claude_home=claude_home))
            if include_mcp:
                records.append(claude_code_mcp_status(spec, claude_command=claude_command))
        else:
            raise ValueError(f"Unsupported agent client: {client}")

    problematic_statuses = {"not_installed", "not_configured", "mismatched", "invalid_config", "cli_not_found", "error"}
    healthy = all(record.status not in problematic_statuses for record in records)
    next_steps: List[str] = []
    if any(record.client == "codex" and record.component == "skill" and record.status == "not_installed" for record in records):
        next_steps.append("Install the Codex skill with `iobrpy-cli agent install --client codex --no-mcp` or rerun the default Codex bootstrap command.")
    if any(record.client == "codex" and record.component == "plugin" and record.status in {"not_installed", "mismatched", "invalid_config"} for record in records):
        next_steps.append("Install or refresh the Codex plugin bundle with `iobrpy-cli agent install --client codex --no-mcp`.")
    if any(record.client == "codex" and record.component in {"result_skill", "result_plugin"} and record.status in {"not_installed", "mismatched", "invalid_config"} for record in records):
        next_steps.append("Install or refresh the IOBRpy result skill and plugin with `iobrpy-cli agent install --client codex --no-mcp --force`.")
    if any(record.client == "codex" and record.component == "mcp" and record.status in {"not_configured", "mismatched", "invalid_config"} for record in records):
        next_steps.append("Install or refresh the Codex MCP config with `iobrpy-cli agent install --client codex --no-skill`.")
    if any(record.client == "claude-code" and record.component == "skill" and record.status in {"not_installed", "mismatched"} for record in records):
        next_steps.append("Install or refresh Claude Code user memory with `iobrpy-cli agent install --client claude-code --no-mcp`.")
    if any(record.client == "claude-code" and record.component == "slash_command" and record.status in {"not_installed", "mismatched"} for record in records):
        next_steps.append("Install or refresh the Claude Code `/iobrpy` command with `iobrpy-cli agent install --client claude-code --no-mcp`.")
    if any(record.client == "claude-code" and record.component in {"result_skill", "result_slash_command"} and record.status in {"not_installed", "mismatched"} for record in records):
        next_steps.append("Install or refresh the Claude Code IOBRpy result integration with `iobrpy-cli agent install --client claude-code --no-mcp`.")
    if any(record.client == "claude-code" and record.component == "mcp" and record.status == "cli_not_found" for record in records):
        next_steps.append("Install the Claude Code CLI first, then run `iobrpy-cli agent install --client claude-code --no-skill`.")
    if any(record.client == "claude-code" and record.component == "mcp" and record.status in {"not_configured", "error"} for record in records):
        next_steps.append("Run `iobrpy-cli agent install --client claude-code --no-skill`, then confirm it with `claude mcp list`.")

    return {
        "success": True,
        "healthy": healthy,
        "server": spec.to_dict(),
        "requested_clients": list(clients) if clients else ["all"],
        "resolved_clients": resolved_clients,
        "components": {
            "skill": include_skill,
            "result_skill": include_skill,
            "result_plugin": include_skill,
            "mcp": include_mcp,
        },
        "results": [record.to_dict() for record in records],
        "next_steps": next_steps,
    }
