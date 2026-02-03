from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def render_html(summary: Dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = "\n".join(
        f"<tr><th>{key}</th><td>{value}</td></tr>" for key, value in summary.items()
    )
    html = f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <title>IOBRpy Report</title>
    <style>
      body {{ font-family: Arial, sans-serif; margin: 2rem; }}
      table {{ border-collapse: collapse; width: 100%; }}
      th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
      th {{ background: #f4f4f4; width: 30%; }}
    </style>
  </head>
  <body>
    <h1>IOBRpy Report</h1>
    <p>Summary for this run:</p>
    <table>
      {rows}
    </table>
  </body>
</html>
"""
    output_path.write_text(html, encoding="utf-8")
