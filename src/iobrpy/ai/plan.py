from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


@dataclass
class PlanStep:
    id: int
    tool_name: str
    mode: str
    args: Dict[str, Any] = field(default_factory=dict)
    argv: List[str] = field(default_factory=list)
    rationale: str = ""
    expected_outputs: List[str] = field(default_factory=list)
    depends_on: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        if not self.args:
            data.pop("args", None)
        if not self.argv:
            data.pop("argv", None)
        return data


@dataclass
class PlanResult:
    steps: List[PlanStep] = field(default_factory=list)
    need_clarification: bool = False
    questions: List[str] = field(default_factory=list)
    required_fields: List[str] = field(default_factory=list)
    example_payload: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "need_clarification": self.need_clarification,
            "questions": self.questions,
            "required_fields": self.required_fields,
            "example_payload": self.example_payload,
            "steps": [step.to_dict() for step in self.steps],
        }


@dataclass
class Plan:
    steps: List[PlanStep]

    def to_dict(self) -> Dict[str, Any]:
        return {"steps": [step.to_dict() for step in self.steps]}


def plan_schema() -> Dict[str, Any]:
    return {
        "type": "object",
        "properties": {
            "need_clarification": {"type": "boolean"},
            "questions": {"type": "array", "items": {"type": "string"}},
            "required_fields": {"type": "array", "items": {"type": "string"}},
            "example_payload": {"type": ["object", "null"]},
            "steps": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["id", "tool_name", "mode", "rationale", "expected_outputs", "depends_on"],
                    "properties": {
                        "id": {"type": "integer"},
                        "tool_name": {"type": "string"},
                        "mode": {"type": "string", "enum": ["kwargs", "argv"]},
                        "args": {"type": "object"},
                        "argv": {"type": "array", "items": {"type": "string"}},
                        "rationale": {"type": "string"},
                        "expected_outputs": {"type": "array", "items": {"type": "string"}},
                        "depends_on": {"type": "array", "items": {"type": "integer"}},
                    },
                },
            },
        },
        "required": ["steps", "need_clarification", "questions", "required_fields"],
    }


def validate_plan_dict(plan: Dict[str, Any]) -> Dict[str, Any]:
    errors: List[str] = []

    if not isinstance(plan, dict):
        return {"ok": False, "errors": ["Plan must be a dictionary."]}

    for field_name in ["steps", "need_clarification", "questions", "required_fields"]:
        if field_name not in plan:
            errors.append(f"Missing required field: {field_name}")

    steps = plan.get("steps", [])
    if not isinstance(steps, list):
        errors.append("Steps must be a list.")
    else:
        for idx, step in enumerate(steps):
            if not isinstance(step, dict):
                errors.append(f"Step {idx} must be a dict.")
                continue
            for required in ["id", "tool_name", "mode", "rationale", "expected_outputs", "depends_on"]:
                if required not in step:
                    errors.append(f"Step {idx} missing '{required}'.")
            if step.get("mode") not in {"kwargs", "argv"}:
                errors.append(f"Step {idx} has invalid mode: {step.get('mode')}")

    return {"ok": not errors, "errors": errors}
