from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List

from iobrpy.core.summary import build_summary


class LocalRunner:
    def __init__(self, run_dir: Path):
        self.run_dir = run_dir
        self.logs_dir = run_dir / "logs" / "steps"
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def _update_status(self, state: str, current_step: str | None, percent: int, error: str | None = None) -> None:
        status_path = self.run_dir / "status.json"
        payload = {
            "state": state,
            "current_step": current_step,
            "percent": percent,
            "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "error": error,
        }
        status_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def _run_command(self, command: List[str], log_path: Path) -> None:
        with log_path.open("w", encoding="utf-8") as log_file:
            process = subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT, check=False)
        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(command)}")

    def _run_spechla_batch(self, params: Dict[str, Any]) -> None:
        manifest_path = self.run_dir / "input" / "manifest.tsv"
        if not manifest_path.exists():
            raise FileNotFoundError("Manifest not found for SpecHLA batch.")
        outdir = self.run_dir / "outputs" / "01_spechla"
        outdir.mkdir(parents=True, exist_ok=True)
        rows = manifest_path.read_text(encoding="utf-8").splitlines()[1:]
        for row in rows:
            if not row.strip():
                continue
            sample, read1, read2 = row.split("\t")
            sample_out = outdir / sample
            sample_out.mkdir(parents=True, exist_ok=True)
            cmd = [
                "iobrpy",
                "spechla",
                "--name",
                sample,
                "--read1",
                read1,
                "--read2",
                read2,
                "--outdir",
                str(sample_out),
                "--threads",
                str(params["inputs"].get("threads", 8)),
            ]
            log_path = self.logs_dir / f"spechla_{sample}.log"
            self._run_command(cmd, log_path)

    def _run_trust4(self, params: Dict[str, Any]) -> None:
        inputs = params["inputs"]
        outdir = self.run_dir / "outputs" / "02_trust4"
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "iobrpy",
            "trust4",
            "--fqdir",
            inputs["fastq_dir"],
            "-o",
            str(outdir),
            "-t",
            str(inputs.get("threads", 8)),
        ]
        self._run_command(cmd, self.logs_dir / "trust4.log")

    def _run_fastq_qc(self, params: Dict[str, Any]) -> None:
        inputs = params["inputs"]
        outdir = self.run_dir / "outputs" / "03_qc"
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "iobrpy",
            "fastq_qc",
            "--path1_fastq",
            inputs["fastq_dir"],
            "--path2_fastp",
            str(outdir),
            "--num_threads",
            str(inputs.get("threads", 8)),
            "--batch_size",
            str(inputs.get("batch_size", 1)),
        ]
        self._run_command(cmd, self.logs_dir / "fastq_qc.log")

    def _run_report(self) -> None:
        summary = build_summary(self.run_dir)
        summary_path = self.run_dir / "artifacts" / "summary.json"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    def run(self, plan: Dict[str, Any], params: Dict[str, Any]) -> None:
        steps = plan.get("steps", [])
        total = len(steps)
        self._update_status("running", None, 0)
        for idx, step in enumerate(steps, start=1):
            action = step["action"]
            self._update_status("running", step["id"], int((idx - 1) / total * 100))
            if action == "spechla_batch":
                self._run_spechla_batch(params)
            elif action == "trust4":
                self._run_trust4(params)
            elif action == "fastq_qc":
                self._run_fastq_qc(params)
            elif action == "report":
                self._run_report()
            elif action == "cli_command":
                command = step.get("command")
                args = step.get("args", [])
                if not command:
                    raise ValueError("cli_command step missing 'command'.")
                cmd = ["iobrpy", command] + [str(item) for item in args]
                log_name = step.get("log_name", f"{command}.log")
                self._run_command(cmd, self.logs_dir / log_name)
            else:
                raise ValueError(f"Unsupported action: {action}")
        self._update_status("succeeded", steps[-1]["id"] if steps else None, 100)
