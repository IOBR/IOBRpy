from __future__ import annotations

import platform
import sys
from pathlib import Path

from iobrpy.version import VERSION


def write_versions(path: Path) -> None:
    lines = [
        f"iobrpy={VERSION}",
        f"python={sys.version.split()[0]}",
        f"platform={platform.platform()}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
