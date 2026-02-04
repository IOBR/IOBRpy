from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from iobrpy.ai.backend import LocalPythonBackend
from iobrpy.ai.executor import Executor
from iobrpy.ai.planner import SimplePlanner
from iobrpy.ai.registry import ToolRegistry
from iobrpy.ai.plan import validate_plan_dict


def run_ai(
    prompt: str,
    workspace: Optional[str] = None,
    dry_run: bool = False,
    plan_only: bool = False,
    json_output: bool = False,
    verbose: bool = False,
    allow_unknown: bool = False,
) -> int:
    registry = ToolRegistry.from_main()
    planner = SimplePlanner()
    plan_result = planner.plan(prompt, registry)
    plan_dict = plan_result.to_dict()

    validation = validate_plan_dict(plan_dict)
    if not validation["ok"]:
        message = {"error": "Invalid plan schema", "details": validation["errors"]}
        if json_output:
            print(json.dumps(message, indent=2))
        else:
            print("Invalid plan schema:")
            for err in validation["errors"]:
                print(f"- {err}")
        return 1

    if workspace is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        workspace = f"./iobrpy_ai_runs/{timestamp}"

    executor = Executor(LocalPythonBackend(), registry)
    if plan_dict.get("need_clarification") and not plan_only:
        plan_only = True
    result = executor.run_plan(
        plan=plan_dict,
        workspace=Path(workspace),
        dry_run=dry_run,
        plan_only=plan_only,
        verbose=verbose,
        allow_unknown=allow_unknown,
    )

    output: Dict[str, Any] = {
        "success": result.success,
        "run_dir": str(result.run_dir),
        "need_clarification": plan_dict.get("need_clarification", False),
        "questions": plan_dict.get("questions", []),
    }
    if result.error:
        output["error"] = result.error

    if json_output:
        print(json.dumps(output, indent=2))
    else:
        print(f"AI run directory: {result.run_dir}")
        if plan_dict.get("need_clarification"):
            print("Plan requires clarification:")
            for question in plan_dict.get("questions", []):
                print(f"- {question}")
        if result.error:
            print(f"Error: {result.error}")

    return 0 if result.success else 1
