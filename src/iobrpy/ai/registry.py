from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from iobrpy.main import build_parser


@dataclass
class ToolSpec:
    name: str
    description: str
    parameters: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    allow_unknown: bool = False
    mode: str = "kwargs"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "parameters": self.parameters,
            "allow_unknown": self.allow_unknown,
            "mode": self.mode,
        }


class ToolRegistry:
    def __init__(self, tools: List[ToolSpec]):
        self._tools = {tool.name: tool for tool in tools}

    @classmethod
    def from_main(cls) -> "ToolRegistry":
        parser = build_parser()
        return cls.from_parser(parser)

    @classmethod
    def from_parser(cls, parser: argparse.ArgumentParser) -> "ToolRegistry":
        subparsers = None
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                subparsers = action
                break
        if subparsers is None:
            raise ValueError("No subparsers found in parser.")

        tools: List[ToolSpec] = []
        for name, subparser in subparsers.choices.items():
            parameters = cls._extract_parameters(subparser)
            description = subparser.description or subparser.format_help().splitlines()[0]
            allow_unknown = name in {
                "runall",
                "trust4",
                "tme_profile",
                "spechla",
                "bayesprism",
                "extract_hla_read",
                "hla_typing",
            }
            mode = "argv" if name in {"extract_hla_read", "hla_typing", "trust4"} else "kwargs"
            if mode == "argv":
                description = (
                    description
                    + " (This tool expects raw CLI argv-style parameters.)"
                )
            tools.append(
                ToolSpec(
                    name=name,
                    description=description,
                    parameters=parameters,
                    allow_unknown=allow_unknown,
                    mode=mode,
                )
            )

        return cls(tools)

    def list_tools(self) -> List[ToolSpec]:
        return list(self._tools.values())

    def get(self, name: str) -> Optional[ToolSpec]:
        return self._tools.get(name)

    @staticmethod
    def _extract_parameters(parser: argparse.ArgumentParser) -> Dict[str, Dict[str, Any]]:
        parameters: Dict[str, Dict[str, Any]] = {}
        for action in parser._actions:
            if isinstance(action, argparse._HelpAction):
                continue
            if not action.option_strings:
                continue
            dest = action.dest
            param_type = None
            if isinstance(action, argparse._StoreTrueAction) or isinstance(
                action, argparse._StoreFalseAction
            ):
                param_type = "bool"
            elif action.type is not None:
                param_type = getattr(action.type, "__name__", str(action.type))

            parameters[dest] = {
                "flags": action.option_strings,
                "required": bool(getattr(action, "required", False)),
                "action": action.__class__.__name__.lower(),
                "type": param_type,
                "choices": list(action.choices) if action.choices else None,
                "default": action.default,
                "help": action.help,
            }

        return parameters
