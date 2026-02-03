from __future__ import annotations

from pathlib import Path


def resolve_workspace(workspace: str) -> Path:
    ws = Path(workspace).expanduser().resolve()
    ws.mkdir(parents=True, exist_ok=True)
    return ws


def ensure_within_workspace(workspace: Path, target: Path) -> Path:
    ws = workspace.resolve()
    tgt = target.resolve()
    try:
        tgt.relative_to(ws)
    except ValueError as exc:
        raise ValueError(f"Path must be within workspace: {tgt}") from exc
    return tgt


def resolve_run_dir(workspace: Path, run_id: str) -> Path:
    if "/" in run_id or run_id.startswith("."):
        run_dir = Path(run_id).expanduser().resolve()
    else:
        run_dir = workspace / "runs" / run_id
    return ensure_within_workspace(workspace, run_dir)
