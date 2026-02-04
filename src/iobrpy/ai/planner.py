from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Optional

from iobrpy.ai.plan import PlanResult, PlanStep
from iobrpy.ai.registry import ToolRegistry


class Planner(ABC):
    @abstractmethod
    def plan(self, prompt: str, registry: ToolRegistry) -> PlanResult:
        raise NotImplementedError


class SimplePlanner(Planner):
    """A lightweight planner that maps prompts to a single tool call."""

    def plan(self, prompt: str, registry: ToolRegistry) -> PlanResult:
        prompt = (prompt or "").strip()
        if not prompt:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=["Please describe the analysis task you want to run."],
                required_fields=["tool_name", "inputs", "outputs"],
                example_payload={
                    "tool_name": "tme_profile",
                    "inputs": {"input": "./data/tpm.tsv"},
                    "outputs": {"output": "./results/tme"},
                },
            )

        tool = self._find_tool(prompt, registry)
        if tool is None:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=[
                    "Which iobrpy command should be used?",
                    "Provide the required input/output parameters for that command.",
                ],
                required_fields=["tool_name", "required_args"],
                example_payload={
                    "tool_name": "bayesprism",
                    "required_args": {
                        "input": "./data/bulk.tsv",
                        "output": "./results/bayesprism",
                    },
                },
            )

        required_args = [
            name
            for name, spec in tool.parameters.items()
            if spec.get("required") and spec.get("action") != "store_true"
        ]

        step = PlanStep(
            id=1,
            tool_name=tool.name,
            mode=tool.mode,
            args={},
            argv=[],
            rationale=f"Matched '{tool.name}' from prompt.",
            expected_outputs=[],
            depends_on=[],
        )

        if tool.mode == "kwargs" and required_args:
            return PlanResult(
                steps=[step],
                need_clarification=True,
                questions=[
                    f"Provide values for required arguments: {', '.join(required_args)}",
                ],
                required_fields=required_args,
                example_payload={
                    "tool_name": tool.name,
                    "args": {arg: f"<{arg}>" for arg in required_args},
                },
            )

        if tool.mode == "argv":
            return PlanResult(
                steps=[step],
                need_clarification=True,
                questions=[
                    "Provide the raw CLI arguments list for this tool.",
                ],
                required_fields=["argv"],
                example_payload={
                    "tool_name": tool.name,
                    "argv": ["--input", "./data/input.bam", "--out", "./results"],
                },
            )

        return PlanResult(steps=[step])

    @staticmethod
    def _find_tool(prompt: str, registry: ToolRegistry) -> Optional[object]:
        prompt_lower = prompt.lower()
        for tool in registry.list_tools():
            if tool.name.lower() in prompt_lower:
                return tool
        return None
