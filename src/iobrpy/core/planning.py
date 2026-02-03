from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List


TEMPLATE_DIRNAME = "workflow_templates"


def load_template(name: str) -> Dict[str, Any]:
    template_path = (
        Path(__file__).resolve().parent.parent
        / TEMPLATE_DIRNAME
        / f"{name}.json"
    )
    if not template_path.exists():
        raise FileNotFoundError(f"Workflow template not found: {template_path}")
    return json.loads(template_path.read_text(encoding="utf-8"))


def build_plan(template: Dict[str, Any], params: Dict[str, Any]) -> Dict[str, Any]:
    steps = template.get("steps", [])
    plan_steps: List[Dict[str, Any]] = []
    for idx, step in enumerate(steps, start=1):
        plan_steps.append(
            {
                "id": f"{idx:02d}_{step['id']}",
                "action": step["action"],
                "description": step.get("description", ""),
                "inputs": step.get("inputs", {}),
                "outputs": step.get("outputs", {}),
            }
        )
    return {
        "workflow": template.get("name"),
        "description": template.get("description", ""),
        "steps": plan_steps,
        "params": params,
    }
