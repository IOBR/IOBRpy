from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict


class Runner(ABC):
    def __init__(self, run_dir: Path):
        self.run_dir = run_dir

    @abstractmethod
    def run(self, plan: Dict[str, Any], params: Dict[str, Any]) -> None:
        raise NotImplementedError
