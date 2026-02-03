from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from iobrpy.core.actions import (
    create_run,
    get_status,
    get_summary,
    render_report,
    start_run,
    tail_log,
    validate_inputs_action,
)
from iobrpy.core.command_registry import COMMAND_SPECS, find_command
from iobrpy.core.workspace import ensure_within_workspace, resolve_workspace


def _extract_path(message: str) -> Optional[str]:
    match = re.search(r"(/[^\\s]+)", message)
    return match.group(1) if match else None


def _extract_int(message: str) -> Optional[int]:
    match = re.search(r"(\\d+)", message)
    return int(match.group(1)) if match else None


def _default_params() -> Dict[str, Dict[str, object]]:
    return {
        "inputs": {
            "fastq_dir": None,
            "manifest": None,
            "threads": 8,
            "batch_size": 1,
            "species": "hsa",
            "project": None,
        },
        "options": {
            "dry_run": True,
            "share_logs": False,
            "share_paths": False,
            "allow_data": False,
        },
    }


def _parse_kv_pairs(text: str) -> Dict[str, str]:
    pairs: Dict[str, str] = {}
    for key, value in re.findall(r"(\\w+)\\s*=\\s*([^\\s]+)", text):
        pairs[key] = value
    return pairs


def _parse_flags(text: str) -> List[str]:
    return re.findall(r"--\\w+", text)


def _prompt_for_value(prompt: str) -> str:
    while True:
        value = input(f"{prompt}: ").strip()
        if value:
            return value


def _ensure_output_in_workspace(workspace: Path, value: str) -> str:
    output_path = Path(value).expanduser()
    ensure_within_workspace(workspace, output_path.resolve())
    return str(output_path)


def _flag_for_param(param: str, spec_flags: Dict[str, str]) -> str:
    return spec_flags.get(param, f"--{param}")


def _build_cli_plan(command: str, args: List[str]) -> Dict[str, object]:
    return {
        "workflow": f"cli:{command}",
        "description": f"Execute iobrpy {command}",
        "steps": [
            {
                "id": "01_cli_command",
                "action": "cli_command",
                "command": command,
                "args": args,
                "log_name": f"{command}.log",
            }
        ],
        "params": {"command": command, "args": args},
    }


def run_chat(workspace: str) -> None:
    os.environ["IOBRPY_WORKSPACE"] = workspace
    ws_path = resolve_workspace(workspace)
    params = _default_params()
    workflow = "spechla_trust4_qc_report"
    run_id: Optional[str] = None

    print("IOBRpy chat (type 'exit' to quit)")
    while True:
        user_input = input("iobrpy> ").strip()
        if not user_input:
            continue
        if user_input.lower() in {"exit", "quit"}:
            break

        if user_input.lower().startswith("status") and run_id:
            print(get_status(run_id))
            continue
        if user_input.lower().startswith("tail") and run_id:
            print(tail_log(run_id))
            continue
        if user_input.lower().startswith("summary") and run_id:
            print(get_summary(run_id))
            continue
        if user_input.lower().startswith("report") and run_id:
            print(render_report(run_id))
            continue
        if user_input.lower().startswith("confirm") and run_id:
            result = start_run(run_id, apply=True)
            print(result)
            continue

        command = find_command(user_input)
        if not command:
            path = _extract_path(user_input)
            if path:
                params["inputs"]["fastq_dir"] = path
            threads = _extract_int(user_input)
            if threads:
                params["inputs"]["threads"] = threads

            validation = validate_inputs_action(workspace, params["inputs"])
            if not validation["ok"]:
                print("Missing or invalid inputs:")
                for err in validation["errors"]:
                    print(f" - {err}")
                print("You can also specify a command name, e.g., 'runall' or 'tme_profile'.")
                continue

            run_info = create_run(workspace, workflow, params, dry_run=True)
            run_id = run_info["run_id"]
            print(f"Plan created (dry-run). Run ID: {run_id}")
            print("Type 'confirm' to start execution.")
            continue

        spec = COMMAND_SPECS[command]
        kv_pairs = _parse_kv_pairs(user_input)
        flags = _parse_flags(user_input)
        args: List[str] = []
        for key, value in kv_pairs.items():
            args.extend([_flag_for_param(key, spec.param_flags), value])
        for flag in flags:
            if flag not in args:
                args.append(flag)

        for param in spec.output_params:
            flag = _flag_for_param(param, spec.param_flags)
            if flag in args:
                idx = args.index(flag)
                if idx + 1 < len(args):
                    try:
                        args[idx + 1] = _ensure_output_in_workspace(ws_path, args[idx + 1])
                    except ValueError as exc:
                        print(f"{exc}")
                        args[idx + 1] = _prompt_for_value(f"Provide {flag} within workspace")

        for param in spec.required:
            flag = _flag_for_param(param, spec.param_flags)
            if flag not in args:
                while True:
                    value = _prompt_for_value(f"Provide {flag}")
                    try:
                        if param in spec.output_params:
                            value = _ensure_output_in_workspace(ws_path, value)
                        break
                    except ValueError as exc:
                        print(f"{exc}")
                args.extend([flag, value])

        for param in spec.optional:
            if f"--{param}" in args:
                continue
            if param in spec.output_params:
                continue

        plan = _build_cli_plan(command, args)
        run_info = create_run(
            workspace,
            workflow="custom_cli",
            params={"inputs": {}, "options": {"dry_run": True}},
            dry_run=True,
            plan_override=plan,
        )
        run_id = run_info["run_id"]
        print(f"Plan created for command '{command}' (dry-run). Run ID: {run_id}")
        print("Type 'confirm' to start execution.")
