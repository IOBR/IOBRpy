from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple


FASTQ_SUFFIXES = ("_1.fastq.gz", "_2.fastq.gz", "_1.fq.gz", "_2.fq.gz")


@dataclass
class ValidationResult:
    ok: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    manifest_rows: List[Tuple[str, str, str]] = field(default_factory=list)


def _infer_sample_name(filename: str) -> str | None:
    for suffix in FASTQ_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return None


def _infer_pairs(fastq_dir: Path) -> Dict[str, Dict[str, str]]:
    pairs: Dict[str, Dict[str, str]] = {}
    for path in sorted(fastq_dir.iterdir()):
        if not path.is_file():
            continue
        name = _infer_sample_name(path.name)
        if not name:
            continue
        entry = pairs.setdefault(name, {})
        if path.name.endswith(("_1.fastq.gz", "_1.fq.gz")):
            entry["read1"] = str(path)
        elif path.name.endswith(("_2.fastq.gz", "_2.fq.gz")):
            entry["read2"] = str(path)
    return pairs


def infer_manifest(fastq_dir: Path) -> List[Tuple[str, str, str]]:
    pairs = _infer_pairs(fastq_dir)
    rows: List[Tuple[str, str, str]] = []
    for sample, reads in sorted(pairs.items()):
        if "read1" in reads and "read2" in reads:
            rows.append((sample, reads["read1"], reads["read2"]))
    return rows


def validate_inputs(fastq_dir: str | None, manifest_path: str | None) -> ValidationResult:
    errors: List[str] = []
    warnings: List[str] = []
    manifest_rows: List[Tuple[str, str, str]] = []

    if not fastq_dir and not manifest_path:
        errors.append("Either fastq_dir or manifest must be provided.")
        return ValidationResult(False, errors, warnings, manifest_rows)

    if fastq_dir:
        fq_path = Path(fastq_dir).expanduser().resolve()
        if not fq_path.exists():
            errors.append(f"FASTQ directory not found: {fq_path}")
        elif not fq_path.is_dir():
            errors.append(f"FASTQ path is not a directory: {fq_path}")
        else:
            manifest_rows = infer_manifest(fq_path)
            if not manifest_rows:
                warnings.append("No paired FASTQ files inferred from directory.")

    if manifest_path:
        mp = Path(manifest_path).expanduser().resolve()
        if not mp.exists():
            errors.append(f"Manifest not found: {mp}")
        elif not mp.is_file():
            errors.append(f"Manifest path is not a file: {mp}")

    return ValidationResult(len(errors) == 0, errors, warnings, manifest_rows)
