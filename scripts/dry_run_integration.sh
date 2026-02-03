#!/usr/bin/env bash
set -euo pipefail

WORKSPACE="${1:-/tmp/iobrpy_ws}"
FASTQ_DIR="${WORKSPACE}/fastq"
mkdir -p "${FASTQ_DIR}"
echo "R1" > "${FASTQ_DIR}/sampleA_1.fastq.gz"
echo "R2" > "${FASTQ_DIR}/sampleA_2.fastq.gz"
export WORKSPACE

python - <<'PY'
import os
from iobrpy.core.actions import create_run, start_run

workspace = os.environ.get("WORKSPACE", "/tmp/iobrpy_ws")
params = {
    "inputs": {"fastq_dir": f"{workspace}/fastq", "threads": 2, "batch_size": 1},
    "options": {"dry_run": True},
}
run_info = create_run(workspace, "spechla_trust4_qc_report", params, dry_run=True)
status = start_run(run_info["run_id"], apply=False)
print(status)
PY
