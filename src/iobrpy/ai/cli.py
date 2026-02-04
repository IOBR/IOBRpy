from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from iobrpy.ai.backend import LocalPythonBackend
from iobrpy.ai.executor import Executor
import json as json_lib
from iobrpy.ai.planner import SimplePlanner, SmartPlanner
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
    interactive: bool = True,
    max_turns: int = 6,
    list_tools: bool = False,
) -> int:
    registry = ToolRegistry.from_main()
    if list_tools:
        for tool in registry.list_tools():
            print(f"{tool.name}: {tool.description}")
        return 0

    planner = SmartPlanner()
    context: Dict[str, object] = {}
    plan_result = planner.plan(prompt, registry, context=context)
    plan_dict = plan_result.to_dict()
    if verbose:
        _print_verbose_plan(plan_result)

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

    if "`" in prompt:
        print(
            "提示：在 bash 中反引号会触发命令替换，请不要用 `...` 包路径；"
            "可改用普通引号 '...' 或直接写路径。"
        )

    turn = 0
    while plan_dict.get("need_clarification") and interactive and turn < max_turns:
        print("缺少必填参数，需补充：")
        for question in plan_dict.get("questions", []):
            print(f"- {question}")
        print("可直接输入参数描述，或输入 JSON（如 {\"args\": {\"outdir\": \"./out\"}}）。")
        user_input = input("> ").strip()
        if not user_input or user_input.lower() in {"quit", "exit"}:
            print("已退出。可重新运行：iobrpy ai \"<你的需求描述>\"")
            return 1
        parsed = _parse_user_input(user_input)
        if parsed.get("args") and isinstance(parsed["args"], dict):
            context.setdefault("args", {})
            context["args"].update(parsed["args"])
        if parsed.get("argv") and isinstance(parsed["argv"], list):
            context.setdefault("argv", [])
            context["argv"].extend(parsed["argv"])
        if parsed.get("prompt"):
            prompt = f"{prompt} {parsed['prompt']}"
        plan_result = planner.plan(prompt, registry, context=context)
        plan_dict = plan_result.to_dict()
        if verbose:
            _print_verbose_plan(plan_result)
        turn += 1

    if plan_dict.get("need_clarification") and not interactive:
        print("需要补充参数，但当前为非交互模式。请提供完整参数后重试。")
        plan_only = True

    executor = Executor(LocalPythonBackend(), registry)
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


def _parse_user_input(user_input: str) -> Dict[str, object]:
    user_input = user_input.strip()
    if user_input.startswith("{") or user_input.startswith("["):
        try:
            parsed = json_lib.loads(user_input)
        except json_lib.JSONDecodeError:
            return {"prompt": user_input}
        if isinstance(parsed, dict):
            return parsed
        if isinstance(parsed, list):
            return {"argv": parsed}
        return {"prompt": user_input}
    return {"prompt": user_input}


def _print_verbose_plan(plan_result: object) -> None:
    steps = getattr(plan_result, "steps", [])
    missing = getattr(plan_result, "required_fields", [])
    if steps:
        step = steps[0]
        print(f"[AI] selected_tool: {step.tool_name}")
        print(f"[AI] extracted_args: {step.args}")
        print(f"[AI] argv: {step.argv}")
    if missing:
        print(f"[AI] missing_required: {missing}")
