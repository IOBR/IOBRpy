from __future__ import annotations

import json
import os
import time
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from iobrpy.ai.backend import ToolBackend
from iobrpy.ai.registry import ToolRegistry


@dataclass
class ExecutionResult:
    success: bool
    run_dir: Path
    artifacts: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None


class Executor:
    def __init__(self, backend: ToolBackend, registry: ToolRegistry):
        self.backend = backend
        self.registry = registry

    def run_plan(
        self,
        plan: Dict[str, Any],
        workspace: Path,
        dry_run: bool = False,
        plan_only: bool = False,
        verbose: bool = False,
        allow_unknown: bool = False,
    ) -> ExecutionResult:
        run_dir = self._create_run_dir(workspace)
        plan_path = run_dir / "plan.json"
        plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")

        if plan_only:
            return ExecutionResult(success=True, run_dir=run_dir)

        calls_path = run_dir / "calls.jsonl"
        artifacts_path = run_dir / "artifacts.json"
        artifacts: Dict[str, Any] = {}

        steps = plan.get("steps", [])
        for step in steps:
            tool_name = step.get("tool_name")
            tool = self.registry.get(tool_name)
            if tool is None:
                return ExecutionResult(
                    success=False,
                    run_dir=run_dir,
                    error=f"Unknown tool: {tool_name}",
                )

            step_mode = step.get("mode", tool.mode)
            step_args = step.get("args") or {}
            step_argv = step.get("argv") or []

            if step_argv and not (tool.allow_unknown or allow_unknown):
                return ExecutionResult(
                    success=False,
                    run_dir=run_dir,
                    error=(
                        f"Unknown argv passthrough is not allowed for tool '{tool_name}'."
                    ),
                )

            validation_error = self._validate_inputs(tool_name, step_args, step_argv)
            if validation_error:
                return ExecutionResult(
                    success=False,
                    run_dir=run_dir,
                    error=validation_error,
                )

            if dry_run:
                self._write_call(calls_path, {
                    "tool_name": tool_name,
                    "mode": step_mode,
                    "args": step_args,
                    "argv": step_argv,
                    "dry_run": True,
                })
                continue

            if verbose:
                print(f"[AI] Running {tool_name} ({step_mode})")

            start_time = time.time()
            try:
                stdout_buffer = StringIO()
                stderr_buffer = StringIO()
                with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                    output = self.backend.run(
                        tool_name=tool_name,
                        args=step_args if step_mode == "kwargs" else {},
                        argv=step_argv if step_mode == "argv" else step_argv,
                    )
                duration = time.time() - start_time
                step_artifacts = self._record_output(output)
                artifacts[str(step.get("id"))] = step_artifacts
                self._write_call(calls_path, {
                    "tool_name": tool_name,
                    "mode": step_mode,
                    "args": step_args,
                    "argv": step_argv,
                    "start": start_time,
                    "end": time.time(),
                    "duration_s": duration,
                    "return_summary": step_artifacts,
                    "stdout": stdout_buffer.getvalue(),
                    "stderr": stderr_buffer.getvalue(),
                })
            except Exception as exc:  # noqa: BLE001
                duration = time.time() - start_time
                self._write_call(calls_path, {
                    "tool_name": tool_name,
                    "mode": step_mode,
                    "args": step_args,
                    "argv": step_argv,
                    "start": start_time,
                    "end": time.time(),
                    "duration_s": duration,
                    "stdout": stdout_buffer.getvalue(),
                    "stderr": stderr_buffer.getvalue(),
                    "exception": str(exc),
                })
                return ExecutionResult(
                    success=False,
                    run_dir=run_dir,
                    error=(
                        f"Tool '{tool_name}' failed with error: {exc}. "
                        "Check inputs or run with --verbose for details."
                    ),
                )

        artifacts_path.write_text(json.dumps(artifacts, indent=2), encoding="utf-8")
        return ExecutionResult(success=True, run_dir=run_dir, artifacts=artifacts)

    @staticmethod
    def _create_run_dir(workspace: Path) -> Path:
        workspace = Path(workspace)
        workspace.mkdir(parents=True, exist_ok=True)
        run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
        run_dir = workspace / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir

    @staticmethod
    def _write_call(path: Path, record: Dict[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")

    @staticmethod
    def _record_output(output: Any) -> Dict[str, Any]:
        if isinstance(output, dict):
            return output
        if output is None:
            return {"status": "completed"}
        return {"return": str(output)}

    def _validate_inputs(
        self,
        tool_name: str,
        args: Dict[str, Any],
        argv: List[str],
    ) -> Optional[str]:
        def check_threads(value: Any, key: str) -> Optional[str]:
            if value is None:
                return None
            try:
                threads = int(value)
            except (TypeError, ValueError):
                return f"Invalid thread value for '{key}': {value}"
            if threads < 1 or threads > 128:
                return f"Thread count for '{key}' must be between 1 and 128."
            return None

        for key, value in args.items():
            if value is None:
                continue
            if key in {"threads", "num_threads", "num_processes", "parallel_size"}:
                error = check_threads(value, key)
                if error:
                    return error
            if any(term in key for term in ["input", "path"]):
                if any(term in key for term in ["output", "outdir", "path_out", "path2_fastp"]):
                    continue
                if isinstance(value, str) and value:
                    path = Path(value)
                    if not path.exists():
                        return f"Input path does not exist: {value}"
            if any(term in key for term in ["output", "outdir", "path_out", "path2_fastp"]):
                if isinstance(value, str) and value:
                    output_path = Path(value)
                    parent = output_path if output_path.suffix == "" else output_path.parent
                    if not parent.exists():
                        try:
                            parent.mkdir(parents=True, exist_ok=True)
                        except OSError as exc:
                            return f"Cannot create output directory {parent}: {exc}"
                    if not os.access(parent, os.W_OK):
                        return f"Output directory is not writable: {parent}"

        if argv:
            for idx, token in enumerate(argv):
                if token.startswith("--") and idx + 1 < len(argv):
                    value = argv[idx + 1]
                    if token in {"--input", "-i", "--bam", "--path"}:
                        if not Path(value).exists():
                            return f"Input path does not exist: {value}"

        return None
