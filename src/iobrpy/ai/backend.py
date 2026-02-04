from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

from iobrpy.main import dispatch


class ToolBackend(ABC):
    @abstractmethod
    def run(
        self,
        tool_name: str,
        args: Optional[Dict[str, Any]] = None,
        argv: Optional[List[str]] = None,
    ) -> Any:
        raise NotImplementedError


class LocalPythonBackend(ToolBackend):
    def run(
        self,
        tool_name: str,
        args: Optional[Dict[str, Any]] = None,
        argv: Optional[List[str]] = None,
    ) -> Any:
        return dispatch(command=tool_name, kwargs=args or {}, unknown=argv or [])


class McpBackend(ToolBackend):
    def run(
        self,
        tool_name: str,
        args: Optional[Dict[str, Any]] = None,
        argv: Optional[List[str]] = None,
    ) -> Any:
        raise NotImplementedError("MCP backend is not yet implemented.")
