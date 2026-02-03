import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from iobrpy.core.actions import create_run, get_status, tail_log, validate_inputs_action


class RunManagementTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.workspace = Path(self.temp_dir.name)
        os.environ["IOBRPY_WORKSPACE"] = str(self.workspace)
        self.fastq_dir = self.workspace / "fastq"
        self.fastq_dir.mkdir()
        (self.fastq_dir / "sampleA_1.fastq.gz").write_text("R1", encoding="utf-8")
        (self.fastq_dir / "sampleA_2.fastq.gz").write_text("R2", encoding="utf-8")

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_create_run_structure(self):
        params = {
            "inputs": {"fastq_dir": str(self.fastq_dir), "threads": 2, "batch_size": 1},
            "options": {"dry_run": True},
        }
        run_info = create_run(str(self.workspace), "spechla_trust4_qc_report", params, dry_run=True)
        run_dir = Path(run_info["run_dir"])
        self.assertTrue((run_dir / "plan.json").exists())
        self.assertTrue((run_dir / "params.json").exists())
        self.assertTrue((run_dir / "status.json").exists())
        self.assertTrue((run_dir / "versions.txt").exists())
        manifest = (run_dir / "input" / "manifest.tsv").read_text(encoding="utf-8")
        self.assertIn("sampleA", manifest)

    def test_validate_inputs(self):
        result = validate_inputs_action(str(self.workspace), {"fastq_dir": str(self.fastq_dir)})
        self.assertTrue(result["ok"])

    def test_tail_log(self):
        params = {
            "inputs": {"fastq_dir": str(self.fastq_dir), "threads": 2, "batch_size": 1},
            "options": {"dry_run": True},
        }
        run_info = create_run(str(self.workspace), "spechla_trust4_qc_report", params, dry_run=True)
        run_dir = Path(run_info["run_dir"])
        log_dir = run_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "runner.log"
        log_path.write_text("line1\nline2\nline3\n", encoding="utf-8")
        tail = tail_log(run_info["run_id"], lines=2)
        self.assertEqual(tail.strip().splitlines(), ["line2", "line3"])

    def test_plan_override(self):
        plan = {
            "workflow": "cli:test",
            "steps": [{"id": "01_cli", "action": "cli_command", "command": "log2_eset", "args": ["--input", "in.tsv", "--output", "out.tsv"]}],
            "params": {},
        }
        run_info = create_run(
            str(self.workspace),
            "custom_cli",
            {"inputs": {}, "options": {"dry_run": True}},
            dry_run=True,
            plan_override=plan,
        )
        run_dir = Path(run_info["run_dir"])
        plan_on_disk = (run_dir / "plan.json").read_text(encoding="utf-8")
        self.assertIn("cli_command", plan_on_disk)


if __name__ == "__main__":
    unittest.main()
