from __future__ import annotations

import json
import os
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

from iobrpy.core.inputs import infer_manifest, validate_inputs
from iobrpy.core.planning import build_plan, load_template
from iobrpy.core.report import render_html
from iobrpy.core.summary import build_summary
from iobrpy.core.versions import write_versions
from iobrpy.core.workspace import ensure_within_workspace, resolve_run_dir, resolve_workspace
from iobrpy.runner.local import LocalRunner


def _timestamp_id() -> str:
    return time.strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def validate_inputs_action(workspace: str, inputs: Dict[str, Any]) -> Dict[str, Any]:
    resolve_workspace(workspace)
    result = validate_inputs(inputs.get("fastq_dir"), inputs.get("manifest"))
    return {
        "ok": result.ok,
        "errors": result.errors,
        "warnings": result.warnings,
    }


def create_run(
    workspace: str,
    workflow: str,
    params: Dict[str, Any],
    dry_run: bool = True,
    plan_override: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    ws = resolve_workspace(workspace)
    run_id = _timestamp_id()
    run_dir = ensure_within_workspace(ws, ws / "runs" / run_id)
    for sub in ("input", "logs/steps", "outputs", "artifacts", "tmp"):
        (run_dir / sub).mkdir(parents=True, exist_ok=True)

    inputs = params.get("inputs", {})
    manifest_path = inputs.get("manifest")
    if not manifest_path and inputs.get("fastq_dir"):
        manifest_rows = infer_manifest(Path(inputs["fastq_dir"]).expanduser().resolve())
        manifest_path = run_dir / "input" / "manifest.tsv"
        manifest_path.write_text(
            "sample\tread1\tread2\n"
            + "\n".join(["\t".join(row) for row in manifest_rows])
            + "\n",
            encoding="utf-8",
        )
        inputs["manifest"] = str(manifest_path)

    if plan_override:
        plan = plan_override
    else:
        template = load_template(workflow)
        plan = build_plan(template, params)
    _write_json(run_dir / "plan.json", plan)
    _write_json(run_dir / "params.json", params)
    _write_json(
        run_dir / "status.json",
        {
            "state": "queued",
            "current_step": None,
            "percent": 0,
            "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "dry_run": dry_run,
        },
    )
    write_versions(run_dir / "versions.txt")
    _write_json(run_dir / "metadata.json", {"run_id": run_id, "workspace": str(ws)})
    return {"run_id": run_id, "run_dir": str(run_dir)}


def start_run(run_id: str, apply: bool = False) -> Dict[str, Any]:
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    ws = resolve_workspace(workspace)
    run_dir = resolve_run_dir(ws, run_id)
    status_path = run_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if not apply:
        status.update({"state": "succeeded", "current_step": "dry_run", "percent": 100})
        _write_json(status_path, status)
        return status

    plan = json.loads((run_dir / "plan.json").read_text(encoding="utf-8"))
    params = json.loads((run_dir / "params.json").read_text(encoding="utf-8"))
    runner = LocalRunner(run_dir)
    runner.run(plan, params)
    return json.loads(status_path.read_text(encoding="utf-8"))


def get_status(run_id: str) -> Dict[str, Any]:
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    run_dir = resolve_run_dir(resolve_workspace(workspace), run_id)
    return json.loads((run_dir / "status.json").read_text(encoding="utf-8"))


def tail_log(run_id: str, step: Optional[str] = None, lines: int = 200) -> str:
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    run_dir = resolve_run_dir(resolve_workspace(workspace), run_id)
    log_dir = run_dir / "logs"
    if step:
        path = log_dir / "steps" / f"{step}.log"
    else:
        path = log_dir / "runner.log"
    if not path.exists():
        return ""
    content = path.read_text(encoding="utf-8").splitlines()
    return "\n".join(content[-lines:])


def get_summary(run_id: str) -> Dict[str, Any]:
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    run_dir = resolve_run_dir(resolve_workspace(workspace), run_id)
    summary = build_summary(run_dir)
    summary_path = run_dir / "artifacts" / "summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(summary_path, summary)
    return summary


def render_report(run_id: str, format: str = "html") -> Dict[str, Any]:
    if format != "html":
        raise ValueError("Only html report format is supported in MVP.")
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    run_dir = resolve_run_dir(resolve_workspace(workspace), run_id)
    summary = build_summary(run_dir)
    report_path = run_dir / "artifacts" / "report.html"
    render_html(summary, report_path)
    return {"report_path": str(report_path)}


def cancel_run(run_id: str) -> Dict[str, Any]:
    workspace = os.environ.get("IOBRPY_WORKSPACE", os.getcwd())
    run_dir = resolve_run_dir(resolve_workspace(workspace), run_id)
    status_path = run_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update({"state": "cancelled", "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S")})
    _write_json(status_path, status)
    return status
