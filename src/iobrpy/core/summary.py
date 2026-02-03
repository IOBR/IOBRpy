from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


def build_summary(run_dir: Path) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}
    status_path = run_dir / "status.json"
    if status_path.exists():
        status = json.loads(status_path.read_text(encoding="utf-8"))
        summary["state"] = status.get("state")
        summary["current_step"] = status.get("current_step")
        summary["percent"] = status.get("percent")
    manifest_path = run_dir / "input" / "manifest.tsv"
    if manifest_path.exists():
        sample_count = sum(1 for _ in manifest_path.read_text(encoding="utf-8").splitlines()[1:] if _)
        summary["sample_count"] = sample_count
    summary["run_dir"] = str(run_dir)
    return summary
