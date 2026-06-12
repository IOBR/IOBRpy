"""
Rendering helpers for agent-facing pipeline map payloads.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .pipeline_map_scan import (
    _is_hla_intermediate_read_path,
    _nonempty_directory_entries,
    _path_has_exact_part,
    _path_has_part_token,
    _path_parts_lower,
    _part_tokens,
)

_EXTERNAL_STATUS_CHECKLIST_ITEM_IDS = {"hla_typing_summary", "tcr_bcr_summary"}
_BAYESPRISM_OPTIONAL_NOTE_EN = "optional standalone; not run by tme_profile/runall"
_BAYESPRISM_OPTIONAL_NOTE_ZH = "独立可选，不由 tme_profile/runall 自动运行"


def _dedupe_matches(matches: Iterable[str], *, limit: int = 5) -> List[str]:
    return sorted(dict.fromkeys(matches))[:limit]


def _render_checklist_marker(checked: bool, status: str) -> str:
    if status == "external":
        return "⚠"
    if status in {"partial", "inferred_existing"}:
        return "△"
    return "✓" if checked else "[ ]"


def _render_checklist_status_labels(status: str) -> Tuple[str, str]:
    labels = {
        "completed": ("Completed", "已完成"),
        "partial": ("Partially completed", "部分完成"),
        "pending": ("Not completed", "未完成"),
        "external": ("External-tool results", "外部工具结果"),
        "inferred_existing": ("Existing results inferred", "已推断存在结果"),
    }
    return labels.get(status, (status.replace("_", " ").title(), status))


def _render_result_source_labels(source_id: str) -> Tuple[str, str]:
    labels = {
        "iobrpy_confirmed_results": ("IOBRpy-confirmed", "IOBRpy 已确认"),
        "agent_inferred_existing_results": ("Existing results inferred", "已推断存在结果"),
        "external_tool_results": ("External-tool results", "外部工具结果"),
        "pending": ("Not yet confirmed", "尚未确认"),
    }
    return labels.get(source_id, (source_id.replace("_", " ").title(), source_id))


def _helper_or_nonresult_file(relpath: str) -> bool:
    lower = relpath.lower()
    suffixes = [suffix.lower() for suffix in Path(relpath).suffixes]
    if lower.endswith(".rdata") or lower.endswith(".rds"):
        return False
    blocked = {
        ".py",
        ".r",
        ".sh",
        ".bash",
        ".zsh",
        ".bat",
        ".cmd",
        ".ps1",
        ".sbatch",
        ".yaml",
        ".yml",
        ".ini",
        ".cfg",
        ".conf",
        ".log",
        ".bam",
        ".bai",
        ".sam",
        ".cram",
        ".crai",
        ".fq",
        ".fastq",
        ".fa",
        ".fasta",
        ".pyc",
    }
    return bool(suffixes and suffixes[-1] in blocked)


def _display_parent_dir(relpath: str) -> str:
    parent = Path(relpath).parent.as_posix()
    if not parent or parent == ".":
        return "./"
    return f"{parent.rstrip('/')}/"


def _display_top_level_dir(relpath: str) -> str:
    parts = Path(relpath).parts
    if len(parts) <= 1:
        return "./"
    return f"{parts[0].rstrip('/')}/"


def _has_late_stage_prefix(relpath: str) -> bool:
    for part in _path_parts_lower(relpath):
        match = re.match(r"^(\d{2})[-_]", part)
        if not match:
            continue
        if int(match.group(1)) < 4:
            continue
        if {"raw", "fastq", "fq"}.intersection(_part_tokens(part)):
            continue
        return True
    return False


def _looks_like_processed_fastq_path(relpath: str) -> bool:
    lower = relpath.lower()
    parts = [part.lower() for part in Path(relpath).parts]
    name = Path(relpath).name.lower()
    if _is_hla_intermediate_read_path(relpath):
        return True
    if _has_late_stage_prefix(relpath):
        return True
    if any(token in name for token in ("trust_", "_toassemble_")):
        return True
    if any(
        part in {"tcr", "bcr", "tcrbcr", "trust4"}
        or re.match(r"^\d{2}[-_]?tcrbcr$", part)
        for part in parts
    ) and any(token in lower for token in ("trust_", "trust4", "toassemble")):
        return True
    if any(
        part in {"03-hr", "hr", "host-remove", "host_remove", "hostremove"}
        or re.match(r"^\d{2}[-_]?hr$", part)
        or re.match(r"^\d{2}[-_]?qc\d*$", part)
        or re.match(r"^qc\d*$", part)
        for part in parts
    ):
        return True
    if any(token in lower for token in ("host_remove", "host-remove", "hostremove")):
        return True
    if any(
        token in lower
        for token in (
            "01-qc",
            "03-fastp",
            "02-salmon",
            "02-star",
            "03-tpm",
            "01-signatures",
            "02-tme",
            "03-lr_cal",
            "04-signatures",
            "05-tme",
            "06-lr_cal",
            "07-tcrbcr",
            "fastp",
            "multiqc",
            "trim",
            "clean",
            "salmon",
            "readspergene",
            "aligned.sortedbycoord.out.bam",
            "trust4",
            "mixcr",
            "immunarch",
            "spechla",
            "extracthlaread",
        )
    ):
        return True
    return False


def _looks_like_raw_fastq_path(relpath: str) -> bool:
    lower = relpath.lower()
    if not any(lower.endswith(ext) for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        return False
    return not _looks_like_processed_fastq_path(relpath)


def _absolute_display_path(base_path: Path, display_path: str) -> str:
    if not display_path:
        return display_path
    normalized = display_path.rstrip("/")
    if normalized in {"", "."}:
        absolute = base_path.resolve()
        return f"{absolute.as_posix()}/"
    absolute = (base_path / normalized).resolve()
    suffix = "/" if display_path.endswith("/") else ""
    return f"{absolute.as_posix()}{suffix}"


def _fallback_evidence_stage_ids(item_id: str) -> List[str]:
    fallback_map = {
        "raw_data": ["raw_fastq"],
        "quality_control": ["fastq_qc"],
        "salmon_quant_merge": ["salmon_quant", "salmon_merge"],
        "star_quant_merge": ["star_quant", "star_merge"],
    }
    return fallback_map.get(item_id, [])


def _has_any_exact_part(relpath: str, parts: Iterable[str]) -> bool:
    return bool(set(_path_parts_lower(relpath)).intersection({part.lower() for part in parts}))


def _looks_like_fastp_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if _has_any_exact_part(relpath, {"badfile", "badfiles", "failed", "failures"}):
        return False
    if name in {"01-qc", "02-fastp", "03-fastp", "fastp", "multiqc"}:
        return True
    if _path_has_exact_part(relpath, "01-qc", "02-fastp", "03-fastp"):
        return True
    if name.endswith(("_fastp.html", "_fastp.json")):
        return True
    if name in {"fastp.html", "fastp.json", "multiqc_report.html", "multiqc_data.json"}:
        return True
    return "multiqc" in lower and "fastp" in lower


def _looks_like_salmon_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if _path_has_part_token(relpath, "star"):
        return False
    if name == "quant.sf":
        return True
    if _path_has_exact_part(relpath, "02-salmon", "salmon"):
        return True
    if _path_has_part_token(relpath, "salmon"):
        return True
    return any(token in lower for token in ("merge_salmon", "salmon_quant", "tximport"))


def _has_result_file_suffix(relpath: str) -> bool:
    suffixes = [suffix.lower() for suffix in Path(relpath.rstrip("/")).suffixes]
    if not suffixes:
        return False
    if suffixes[-1] in {".csv", ".tsv", ".txt", ".sf", ".rdata", ".rds"}:
        return True
    return suffixes[-1] == ".gz" and len(suffixes) >= 2 and suffixes[-2] in {".csv", ".tsv", ".txt"}


def _looks_like_concrete_salmon_result_file(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if relpath.endswith("/") or not _has_result_file_suffix(relpath):
        return False
    if _path_has_part_token(relpath, "star"):
        return False
    if name == "quant.sf":
        return True
    return any(
        token in lower
        for token in (
            "merge_salmon",
            "salmon_quant",
            "tximport",
            "_salmon_tpm",
            "_salmon_count",
            "merged_salmon",
            "salmon_tpm",
            "salmon_count",
        )
    )


def _looks_like_concrete_star_count_file(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if relpath.endswith("/") or not _has_result_file_suffix(relpath):
        return False
    if name in {"biosample_count.txt", "biosample_counts.txt"}:
        return False
    return any(
        token in lower
        for token in (
            "star.count",
            "star_count",
            "merge_star_count",
            "featurecounts",
            "featurecount",
            "count_matrix",
            "counts_matrix",
            "gene_count",
            "read_count",
        )
    )


def _display_salmon_quant_root(relpath: str) -> str:
    path = Path(relpath.rstrip("/"))
    if path.name.lower() != "quant.sf":
        return relpath
    sample_dir = path.parent
    salmon_root = sample_dir.parent
    if str(salmon_root) in {"", "."}:
        return _display_parent_dir(relpath)
    return f"{salmon_root.as_posix().rstrip('/')}/"


def _looks_like_trust4_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if any(token in lower for token in ("mixcr", "immunarch", "immunearch")):
        return False
    return any(
        token in lower
        for token in (
            "trust4",
            "07-tcrbcr",
            "tcrbcr",
            "trust_",
            "_airr.tsv",
            "_airr_align.tsv",
            "_cdr3.out",
            ".trust4.done",
            "trust4_immdata",
            "trust4_immune_indices",
        )
    ) or (
        name.endswith("_report.tsv")
        and any(part in {"tcr", "bcr", "tcrbcr", "trust4"} or re.match(r"^\d{2}[-_]?tcrbcr$", part) for part in _path_parts_lower(relpath))
    )


def _display_trust4_result_root(relpath: str) -> str:
    path = Path(relpath.rstrip("/"))
    parts = path.parts
    lower_parts = [part.lower() for part in parts]
    for idx, part in enumerate(lower_parts):
        tokens = set(_part_tokens(part))
        if (
            part in {"tcr", "bcr", "tcrbcr", "trust4", "07-tcrbcr"}
            or re.match(r"^\d{2}[-_]?tcrbcr$", part)
            or tokens.intersection({"tcrbcr", "trust4"})
        ):
            return f"{Path(*parts[: idx + 1]).as_posix().rstrip('/')}/"
    name = path.name.lower()
    if any(token in name for token in ("trust_", "trust4_", ".trust4.done")) and path.parent != Path("."):
        parent = path.parent
        if parent.parent != Path(".") and parent.name:
            return f"{parent.parent.as_posix().rstrip('/')}/"
        return _display_parent_dir(relpath)
    return relpath


def _looks_like_external_tcr_bcr_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    if _looks_like_trust4_result_path(relpath) or relpath.endswith("/"):
        return False
    return any(
        token in lower
        for token in (
            "mixcr",
            "immunarch",
            "immunearch",
            "clonotype",
            "clonotypes",
            ".clones_",
            "airr",
            "cdr3",
            "vdj",
            "rearrangement",
            "align.report.json",
        )
    )


def _display_tcr_bcr_result_root(relpath: str) -> str:
    path = Path(relpath.rstrip("/"))
    parts = path.parts
    lower_parts = [part.lower() for part in parts]
    for idx, part in enumerate(lower_parts):
        tokens = set(_part_tokens(part))
        if (
            part in {"tcr", "bcr", "tcrbcr", "trust4", "mixcr", "05-mixcr", "07-tcrbcr", "17-tcrbcr", "19-mixcr"}
            or re.match(r"^\d{2}[-_]?(tcrbcr|mixcr)$", part)
            or tokens.intersection({"tcr", "bcr", "tcrbcr", "trust4", "mixcr"})
        ):
            return f"{Path(*parts[: idx + 1]).as_posix().rstrip('/')}/"
    return _display_parent_dir(relpath) if path.suffixes else relpath


def _display_hla_result_root(relpath: str) -> str:
    path = Path(relpath.rstrip("/"))
    parts = path.parts
    lower_parts = [part.lower() for part in parts]
    for idx, part in enumerate(lower_parts):
        tokens = set(_part_tokens(part))
        if (
            part in {"hla", "spechla", "06-hla", "11-spechla", "11-spechla_1", "13-hla"}
            or re.match(r"^\d{2}[-_]?hla$", part)
            or tokens.intersection({"hla", "spechla", "hlatyping"})
        ):
            return f"{Path(*parts[: idx + 1]).as_posix().rstrip('/')}/"
    return _display_parent_dir(relpath) if path.suffixes else relpath


def _looks_like_hla_per_sample_result_path(relpath: str) -> bool:
    name = Path(relpath.rstrip("/")).name.lower()
    return (
        not relpath.endswith("/")
        and name.endswith(("_result.tsv", "_result.csv", "_result.txt"))
        and _path_has_part_token(relpath, "hla", "spechla", "hlatyping")
    )


def _has_star_stage_part(relpath: str) -> bool:
    stage_parts = {
        "02-star",
        "03-star",
        "04-star",
        "star",
        "star-count",
        "star-counts",
        "star_count",
        "star_counts",
        "star-output",
        "star-outputs",
        "star_output",
        "star_outputs",
        "star-results",
        "star_results",
    }
    if _has_any_exact_part(relpath, stage_parts):
        return True
    for part in _path_parts_lower(relpath):
        tokens = set(_part_tokens(part))
        if "star" not in tokens:
            continue
        if re.match(r"^\d{2}[-_]?star$", part):
            return True
        if tokens.intersection({"count", "counts", "output", "outputs", "result", "results"}):
            return True
    return False


def _looks_like_star_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if _has_any_exact_part(relpath, {"ctat_db_dir", "script", "scripts"}):
        return False
    if _has_star_stage_part(relpath):
        return True
    return any(
        token in lower
        for token in (
            "readspergene.out.tab",
            "aligned.sortedbycoord.out.bam",
            "aligned.out.sam",
            "log.final.out",
            "sj.out.tab",
            "merge_star_count",
            "star.count",
        )
    ) or name.startswith("readspergene")


def _filter_item_specific_evidence(item_id: str, evidence_paths: Iterable[str]) -> List[str]:
    filtered: List[str] = []
    for path in evidence_paths:
        lower = path.lower()
        if item_id == "raw_data":
            if _looks_like_raw_fastq_path(path):
                filtered.append(path)
            continue
        if item_id == "quality_control":
            if _looks_like_fastp_result_path(path):
                filtered.append(path)
            continue
        if item_id == "salmon_quant_merge":
            if _looks_like_salmon_result_path(path):
                filtered.append(path)
            continue
        if item_id == "star_quant_merge":
            if _looks_like_star_result_path(path):
                filtered.append(path)
            continue
        if item_id == "immune_deconvolution":
            if any(
                token in lower
                for token in (
                    "cibersort",
                    "estimate",
                    "ips",
                    "mcpcounter",
                    "quantiseq",
                    "epic",
                    "bayesprism",
                    "deconvo_merged",
                    "deconvo-merged",
                )
            ):
                filtered.append(path)
            continue
        if item_id == "ligand_receptor_analysis":
            if any(
                token in lower
                for token in (
                    "lr_cal",
                    "lr_scores",
                    "/03-lr_cal/",
                    "\\03-lr_cal\\",
                    "/06-lr_cal/",
                    "\\06-lr_cal\\",
                )
            ):
                filtered.append(path)
            continue
        if item_id == "feature_scoring":
            if "calculate_sig_score" in lower or any(
                token in lower
                for token in (
                    "/01-signatures/",
                    "\\01-signatures\\",
                    "/04-signatures/",
                    "\\04-signatures\\",
                )
            ):
                filtered.append(path)
            continue
        if item_id == "hla_typing_summary":
            if any(
                token in lower
                for token in (
                    "mixcr",
                    "trust4",
                    "tcr",
                    "bcr",
                    "airr",
                    "clonotype",
                    "clonotypes",
                    "vdj",
                    "cdr3",
                    "immunarch",
                    "immunearch",
                    ".clones_",
                )
            ):
                continue
            filtered.append(path)
            continue
        if item_id == "tcr_bcr_summary":
            if _looks_like_trust4_result_path(path) or any(
                token in lower
                for token in (
                    "tcr",
                    "bcr",
                    "mixcr",
                    "immunarch",
                    "airr",
                    "cdr3",
                    "clonotype",
                    "vdj",
                    "_report.tsv",
                    "immune_indices",
                    "immdata",
                    "repertoire",
                )
            ):
                filtered.append(path)
            continue
        if item_id == "tpm_matrix_ready":
            if any(token in lower for token in ("tpm", "prepare_salmon", "count2tpm", "gene_symbol", "gene_symbols")):
                filtered.append(path)
            continue
        filtered.append(path)
    return filtered


def _prefer_more_specific_dirs(paths: Iterable[str]) -> List[str]:
    normalized = []
    for path in paths:
        cleaned = path.rstrip("/") + "/"
        if cleaned not in normalized:
            normalized.append(cleaned)
    kept: List[str] = []
    for path in normalized:
        if any(other != path and other.startswith(path) for other in normalized):
            continue
        kept.append(path)
    return kept


def _prefer_result_files(evidence_paths: Iterable[str], *, limit: int) -> List[str]:
    stage_dirs = {
        "01-signatures",
        "02-tme",
        "03-lr_cal",
        "04-clustering",
        "04-signatures",
        "05-tme",
        "06-lr_cal",
        "07-tcrbcr",
    }
    file_like = [
        path
        for path in evidence_paths
        if not path.endswith("/")
        and Path(path.rstrip("/")).name.lower() not in stage_dirs
        and bool(Path(path.rstrip("/")).suffixes)
    ]
    return _dedupe_matches(file_like or evidence_paths, limit=limit)


def _bayesprism_display_priority(relpath: str) -> int:
    name = Path(relpath.rstrip("/")).name.lower()
    if name == "theta" or name.startswith("theta."):
        return 0
    if name.startswith("theta_cv"):
        return 1
    if name.startswith("z_"):
        return 2
    return 3


def _prefer_deconvolution_result_files(evidence_paths: Iterable[str], *, limit: int) -> List[str]:
    stage_dirs = {
        "01-signatures",
        "02-tme",
        "03-lr_cal",
        "04-clustering",
        "04-signatures",
        "05-tme",
        "06-lr_cal",
        "07-tcrbcr",
    }
    candidates = [
        path
        for path in evidence_paths
        if not path.endswith("/")
        and Path(path.rstrip("/")).name.lower() not in stage_dirs
        and bool(Path(path.rstrip("/")).suffixes)
    ] or list(evidence_paths)
    unique = list(dict.fromkeys(candidates))

    def _sort_key(path: str) -> Tuple[str, int, str, str]:
        normalized = path.rstrip("/").replace("\\", "/")
        lower = normalized.lower()
        if _path_has_part_token(path, "bayesprism"):
            parent = Path(normalized).parent.as_posix().lower()
            name = Path(normalized).name.lower()
            return (parent, _bayesprism_display_priority(normalized), name, lower)
        return (lower, 0, Path(normalized).name.lower(), lower)

    return sorted(unique, key=_sort_key)[:limit]


def _checklist_evidence_display_paths(
    item_id: str,
    evidence_paths: List[str],
    display_hints: Dict[str, List[str]],
) -> List[str]:
    if item_id == "salmon_quant_merge":
        quant_dirs = [
            _display_salmon_quant_root(path)
            for path in evidence_paths
            if Path(path.rstrip("/")).name.lower() == "quant.sf"
        ]
        concrete_files = [
            path
            for path in evidence_paths
            if _looks_like_concrete_salmon_result_file(path)
            and Path(path.rstrip("/")).name.lower() != "quant.sf"
        ]
        if quant_dirs or concrete_files:
            return _dedupe_matches(quant_dirs + concrete_files, limit=6)

    if item_id in {"raw_data", "quality_control"} and evidence_paths:
        hinted = display_hints.get(item_id, [])
        if hinted:
            return _dedupe_matches(hinted, limit=6)

    if item_id == "raw_data":
        dirs = [
            _display_parent_dir(path)
            for path in evidence_paths
            if _looks_like_raw_fastq_path(path)
        ]
        return _dedupe_matches(dirs or evidence_paths, limit=6)

    if item_id == "quality_control":
        dirs = [_display_parent_dir(path) for path in evidence_paths if _looks_like_fastp_result_path(path)]
        return _dedupe_matches(dirs or evidence_paths, limit=6)

    if item_id == "salmon_quant_merge":
        quant_dirs = [
            _display_salmon_quant_root(path)
            for path in evidence_paths
            if Path(path.rstrip("/")).name.lower() == "quant.sf"
        ]
        concrete_files = [
            path
            for path in evidence_paths
            if _looks_like_concrete_salmon_result_file(path)
            and Path(path.rstrip("/")).name.lower() != "quant.sf"
        ]
        if quant_dirs or concrete_files:
            return _dedupe_matches(quant_dirs + concrete_files, limit=6)
        dirs = [
            _display_top_level_dir(path)
            for path in evidence_paths
            if _looks_like_salmon_result_path(path)
        ]
        return _dedupe_matches(dirs or evidence_paths, limit=6)

    if item_id == "star_quant_merge":
        count_files = [
            path
            for path in evidence_paths
            if _looks_like_concrete_star_count_file(path)
        ]
        dirs = [
            _display_top_level_dir(path)
            for path in evidence_paths
            if _looks_like_star_result_path(path)
        ]
        if count_files:
            return _dedupe_matches(count_files, limit=6)
        return _dedupe_matches(dirs or evidence_paths, limit=6)

    if item_id == "tpm_matrix_ready":
        preferred = [
            path for path in evidence_paths
            if any(token in path.lower() for token in ("merged_salmon_tpm", "tpm_gene_symbols", "prepare_salmon", "count2tpm"))
        ]
        return _dedupe_matches(preferred or evidence_paths, limit=3)

    if item_id == "immune_deconvolution":
        # Keep concrete deconvolution result files visible so the checklist can
        # render full absolute paths, one per line when several outputs exist.
        return _prefer_deconvolution_result_files(evidence_paths, limit=12)

    if item_id in {"feature_scoring", "ligand_receptor_analysis"}:
        return _prefer_result_files(evidence_paths, limit=6)

    if item_id == "tcr_bcr_summary":
        result_roots: List[str] = []
        for path in evidence_paths:
            if _looks_like_trust4_result_path(path):
                result_roots.append(_display_trust4_result_root(path))
                continue
            if _looks_like_external_tcr_bcr_result_path(path):
                result_roots.append(_display_tcr_bcr_result_root(path))
        if result_roots:
            return _dedupe_matches(result_roots, limit=3)
        file_like = _prefer_result_files(evidence_paths, limit=3)
        if file_like:
            return file_like
        dirs = _prefer_more_specific_dirs(_display_parent_dir(path) for path in evidence_paths)
        return _dedupe_matches(dirs or evidence_paths, limit=3)

    if item_id == "hla_typing_summary":
        preferred = [
            path for path in evidence_paths
            if any(token in path.lower() for token in ("hla_result_merged", "cohort_genotype"))
        ]
        if preferred:
            return _dedupe_matches(preferred, limit=2)
        hla_roots = [
            _display_hla_result_root(path)
            for path in evidence_paths
            if (path.endswith("/") or not Path(path.rstrip("/")).suffixes)
            and _path_has_part_token(path, "hla", "spechla", "hlatyping")
            or _looks_like_hla_per_sample_result_path(path)
            or Path(path.rstrip("/")).name.lower().endswith(".task.complete")
        ]
        if hla_roots:
            return _dedupe_matches(hla_roots, limit=3)
        file_like = _prefer_result_files(evidence_paths, limit=3)
        if file_like:
            return file_like
        dirs = [_display_top_level_dir(path) for path in evidence_paths]
        return _dedupe_matches(dirs or evidence_paths, limit=3)

    return _dedupe_matches(evidence_paths, limit=6)


def _scan_checklist_display_hints(entries: Iterable[str]) -> Dict[str, List[str]]:
    raw_dirs: List[str] = []
    qc_dirs: List[str] = []
    salmon_dirs: List[str] = []
    star_dirs: List[str] = []
    normalized_entries = [str(entry).rstrip("/\\") for entry in entries if str(entry).strip()]
    nonempty_dirs = _nonempty_directory_entries(normalized_entries)

    for relpath in entries:
        lower = relpath.lower()
        suffixes = [suffix.lower() for suffix in Path(relpath).suffixes]
        normalized = str(relpath).rstrip("/\\")
        if _looks_like_raw_fastq_path(relpath):
            raw_dirs.append(_display_parent_dir(relpath))
        if _looks_like_fastp_result_path(relpath):
            if not suffixes and normalized not in nonempty_dirs:
                continue
            qc_dirs.append(_display_parent_dir(relpath) if suffixes else _display_top_level_dir(relpath))
        if _looks_like_concrete_salmon_result_file(relpath):
            salmon_dirs.append(_display_top_level_dir(relpath))
        if any(token in lower for token in ("readspergene.out.tab", "aligned.sortedbycoord.out.bam")):
            star_dirs.append(_display_top_level_dir(relpath))
        elif _looks_like_star_result_path(relpath):
            if not suffixes and normalized not in nonempty_dirs:
                continue
            star_dirs.append(_display_top_level_dir(relpath))

    raw_dirs = _dedupe_matches(raw_dirs, limit=6)
    qc_dirs = _dedupe_matches(qc_dirs, limit=6)
    salmon_dirs = _dedupe_matches(salmon_dirs, limit=6)
    star_dirs = _dedupe_matches(star_dirs, limit=6)
    if any(path != "./" for path in salmon_dirs):
        salmon_dirs = [path for path in salmon_dirs if path != "./"]
    if any(path != "./" for path in star_dirs):
        star_dirs = [path for path in star_dirs if path != "./"]
    if any(path != "./" for path in qc_dirs):
        qc_dirs = [path for path in qc_dirs if path != "./"]
    if any(path != "./" for path in raw_dirs):
        raw_dirs = [path for path in raw_dirs if path != "./"]

    return {
        "raw_data": raw_dirs,
        "quality_control": qc_dirs,
        "salmon_quant_merge": salmon_dirs,
        "star_quant_merge": star_dirs,
    }


def _format_map_command(
    path: str,
    *,
    quick: bool = False,
    focus_roots: Iterable[str] = (),
    max_entries: Optional[int] = None,
) -> str:
    parts = ["iobrpy-cli", "map", "--path", path, "--json"]
    if quick:
        parts.append("--quick")
    for focus_root in focus_roots:
        parts.extend(["--focus", focus_root])
    if max_entries is not None:
        parts.extend(["--max-entries", str(max_entries)])
    return " ".join(shlex.quote(part) for part in parts)


def _focus_root_tokens(root: str) -> Set[str]:
    tokens: Set[str] = set()
    for part in Path(root).parts:
        tokens.update(_part_tokens(part))
    return tokens


def _top_level_focus_root_from_path(path: str) -> str:
    normalized = Path(path)
    parts = normalized.parts
    if not parts:
        return ""
    if parts[0] in {".", "/"}:
        return parts[1] if len(parts) > 1 else ""
    return parts[0]


def _recommended_deep_focus_roots(payload: Dict[str, Any]) -> List[str]:
    scan = payload.get("scan", {})
    completed_stages = set(payload.get("completed_stages", []))
    candidates = list(scan.get("focus_roots", []) or [])
    for detection in payload.get("function_detections", []):
        if detection.get("status") == "not_detected":
            continue
        for evidence in detection.get("evidence", []) or []:
            root = _top_level_focus_root_from_path(str(evidence))
            if root:
                candidates.append(root)
    for hint in payload.get("external_analysis_hints", []):
        for evidence in hint.get("evidence", []) or []:
            root = _top_level_focus_root_from_path(str(evidence))
            if root:
                candidates.append(root)
    candidates = _dedupe_matches(candidates, limit=len(candidates))
    if not candidates:
        return []

    preferred: List[str] = []
    fallback: List[str] = []
    bulky_fallback: List[str] = []
    for root in candidates:
        tokens = _focus_root_tokens(root)
        is_profile_like = bool(tokens.intersection({"tme", "profile", "runall", "signature", "signatures", "lr", "cluster", "nmf", "trust4"}))
        is_merge_like = "merge" in tokens
        is_hla_like = bool(tokens.intersection({"hla", "spechla"}))
        is_bulky_upstream = bool(tokens.intersection({"fastq", "fq", "qc", "fastp", "salmon", "star"}))
        is_mixcr_raw = "mixcr" in tokens and "merge" not in tokens

        if is_bulky_upstream:
            bulky_fallback.append(root)
            continue
        if is_mixcr_raw:
            bulky_fallback.append(root)
            continue
        if is_profile_like or is_merge_like or is_hla_like:
            preferred.append(root)
        else:
            fallback.append(root)

    if "tpm_matrix" in completed_stages or "prepare_salmon" in completed_stages or "salmon_merge" in completed_stages:
        fallback = [root for root in fallback if "salmon" not in _focus_root_tokens(root)]
        preferred = [root for root in preferred if "salmon" not in _focus_root_tokens(root)]
    if "tpm_matrix" in completed_stages or "star_merge" in completed_stages:
        fallback = [root for root in fallback if "star" not in _focus_root_tokens(root)]
        preferred = [root for root in preferred if "star" not in _focus_root_tokens(root)]

    preferred_tokens = [_focus_root_tokens(root) for root in preferred]
    if any("merge" in tokens and "mixcr" in tokens for tokens in preferred_tokens):
        fallback = [
            root
            for root in fallback
            if not ("mixcr" in _focus_root_tokens(root) and "merge" not in _focus_root_tokens(root))
        ]
        bulky_fallback = [
            root
            for root in bulky_fallback
            if not ("mixcr" in _focus_root_tokens(root) and "merge" not in _focus_root_tokens(root))
        ]

    if preferred:
        return preferred[:3]

    recommended = fallback[:3]
    if recommended:
        return recommended

    recommended = bulky_fallback[:2]
    return recommended


def _stage_evidence_for_checklist_item(
    item: Dict[str, Any],
    strong_evidence: Dict[str, List[str]],
    weak_evidence: Dict[str, List[str]],
    external_analysis_hints: List[Dict[str, Any]],
    deconvolution_methods: Dict[str, Any],
    display_hints: Dict[str, List[str]],
    *,
    strict_iobrpy_investigation_item_ids: Set[str],
) -> Tuple[List[str], List[str]]:
    evidence: List[str] = []
    for stage_id in item.get("matched_stage_ids", []):
        evidence.extend(strong_evidence.get(stage_id, []))
    if not evidence:
        for stage_id in item.get("matched_stage_ids", []):
            evidence.extend(weak_evidence.get(stage_id, []))
    if not evidence:
        for stage_id in _fallback_evidence_stage_ids(item.get("id", "")):
            evidence.extend(strong_evidence.get(stage_id, []))
            evidence.extend(weak_evidence.get(stage_id, []))
    if item.get("id") == "immune_deconvolution":
        for paths in deconvolution_methods.get("evidence", {}).values():
            evidence.extend(paths)
    if item.get("matched_external_hint_ids"):
        hint_map = {hint.get("id"): hint for hint in external_analysis_hints}
        for hint_id in item.get("matched_external_hint_ids", []):
            evidence.extend(hint_map.get(hint_id, {}).get("evidence", []))
    strict_iobrpy_only = item.get("id") in strict_iobrpy_investigation_item_ids
    bucket_counts = item.get("investigation_bucket_counts", {})
    if not strict_iobrpy_only or bucket_counts.get("iobrpy_confirmed_results"):
        evidence.extend(item.get("investigation_sample_paths", []))
    filtered: List[str] = []
    for path in evidence:
        if not path:
            continue
        if item.get("id") == "raw_data":
            if _looks_like_processed_fastq_path(path):
                continue
            filtered.append(path)
            continue
        if not _helper_or_nonresult_file(path):
            filtered.append(path)
    filtered = _filter_item_specific_evidence(item.get("id", ""), filtered)
    raw_evidence_paths = _dedupe_matches(filtered, limit=12)
    display_paths = _checklist_evidence_display_paths(item.get("id", ""), raw_evidence_paths, display_hints)
    return raw_evidence_paths, display_paths


def _checklist_missing_items(
    item: Dict[str, Any],
    deconvolution_methods: Dict[str, Any],
    language: str = "en",
) -> List[str]:
    if item.get("id") == "immune_deconvolution":
        labels = list(deconvolution_methods.get("missing_labels", []))
        ids = list(deconvolution_methods.get("missing_ids", []))
        optional_ids = set(deconvolution_methods.get("optional_missing_ids", []))
        optional_labels = list(deconvolution_methods.get("optional_missing_labels", []))

        def annotate(label: str, method_id: str = "") -> str:
            is_optional = method_id in optional_ids or label in optional_labels
            if not is_optional:
                return label
            if language == "zh":
                return f"{label}（独立可选，不由 tme_profile/runall 自动运行）"
            return f"{label} ({_BAYESPRISM_OPTIONAL_NOTE_EN})"

        if len(ids) == len(labels):
            return [annotate(label, method_id) for method_id, label in zip(ids, labels)]

        rendered = [annotate(label) for label in labels if label not in optional_labels]
        rendered.extend(annotate(label) for label in optional_labels if label not in labels)
        return rendered
    if item.get("id") == "hla_typing_summary" and item.get("matched_external_hint_ids"):
        key = "pending_stage_labels_zh" if language == "zh" else "pending_stage_labels"
        pending = list(item.get(key, []))
        merged_tokens = ("合并",) if language == "zh" else ("merged",)
        return [label for label in pending if any(token in label.lower() for token in merged_tokens)]
    if item.get("id") == "tpm_matrix_ready" and item.get("checked"):
        return []
    key = "pending_stage_labels_zh" if language == "zh" else "pending_stage_labels"
    return list(item.get(key, []))


def _checklist_wrapper_hint(
    item_id: str,
    function_detections: Dict[str, Dict[str, Any]],
) -> Tuple[str, str]:
    tme_profile = function_detections.get("tme_profile", {})
    status = tme_profile.get("status")
    if item_id not in {"feature_scoring", "immune_deconvolution", "ligand_receptor_analysis"}:
        return ("", "")
    if status == "confirmed_by_content":
        return (
            "These downstream outputs are consistent with a previous iobrpy `tme_profile` run.",
            "这些下游输出与之前运行过的 iobrpy `tme_profile` 一致。",
        )
    if status == "likely_iobrpy_result":
        return (
            "These downstream outputs look consistent with a previous `tme_profile`-like run, but the full wrapper result set is not completely confirmed yet.",
            "这些下游输出看起来与之前的 `tme_profile` 类运行一致，但完整 wrapper 结果尚未完全确认。",
        )
    return ("", "")


def _render_checklist_lines(items: List[Dict[str, Any]], *, language: str) -> List[str]:
    lines: List[str] = []
    for item in items:
        marker = item.get("marker", "[ ]")
        label = item.get("label_zh", item.get("label", "?")) if language == "zh" else item.get("label", "?")
        status_label = item.get("status_label_zh", item.get("status_label", "")) if language == "zh" else item.get("status_label", "")
        lines.append(f"{marker} [{status_label}] {label}")

        evidence = item.get("evidence_display_paths_zh", item.get("evidence_display_paths", [])) if language == "zh" else item.get("evidence_display_paths", [])
        if evidence:
            prefix = "   检测到：" if language == "zh" else "   Detected in: "
            connector = "、" if language == "zh" else ", "
            lines.append(f"{prefix}{connector.join(evidence)}")

        wrapper_hint = item.get("wrapper_hint_zh", item.get("wrapper_hint", "")) if language == "zh" else item.get("wrapper_hint", "")
        if wrapper_hint:
            prefix = "   备注：" if language == "zh" else "   Note: "
            lines.append(f"{prefix}{wrapper_hint}")

        missing = item.get("missing_items_zh", item.get("missing_items", [])) if language == "zh" else item.get("missing_items", [])
        if missing:
            prefix = "   还缺少：" if language == "zh" else "   Still missing: "
            connector = "、" if language == "zh" else ", "
            lines.append(f"{prefix}{connector.join(missing)}")
    return lines


def _table_cell(value: str) -> str:
    cell = (value or "-").replace("|", "\\|").replace("\n", " ").strip()
    return cell or "-"


def _render_checklist_table(items: List[Dict[str, Any]], *, language: str) -> str:
    if language == "zh":
        header = "| 标记 | 状态 | 阶段 | 检测到 | 仍缺少 |"
    else:
        header = "| Marker | Status | Stage | Detected in | Still missing |"
    lines = [header, "| --- | --- | --- | --- | --- |"]
    for item in items:
        marker = item.get("marker", "[ ]")
        label = item.get("label_zh", item.get("label", "?")) if language == "zh" else item.get("label", "?")
        status_label = item.get("status_label_zh", item.get("status_label", "")) if language == "zh" else item.get("status_label", "")
        evidence = item.get("evidence_display_paths_zh", item.get("evidence_display_paths", [])) if language == "zh" else item.get("evidence_display_paths", [])
        missing = item.get("missing_items_zh", item.get("missing_items", [])) if language == "zh" else item.get("missing_items", [])
        detected_text = "<br>".join(evidence) if evidence else "-"
        missing_text = "<br>".join(missing) if missing else "-"
        lines.append(
            "| {0} | {1} | {2} | {3} | {4} |".format(
                _table_cell(marker),
                _table_cell(status_label),
                _table_cell(label),
                _table_cell(detected_text),
                _table_cell(missing_text),
            )
        )
    return "\n".join(lines)


def _render_checklist_lines(items: List[Dict[str, Any]], *, language: str) -> List[str]:
    lines: List[str] = []
    for item in items:
        marker = item.get("marker", "[ ]")
        label = item.get("label_zh", item.get("label", "?")) if language == "zh" else item.get("label", "?")
        status_label = item.get("status_label_zh", item.get("status_label", "")) if language == "zh" else item.get("status_label", "")
        lines.append(f"{marker} [{status_label}] {label}")

        evidence = (
            item.get("evidence_display_paths_zh", item.get("evidence_display_paths", []))
            if language == "zh"
            else item.get("evidence_display_paths", [])
        )
        if evidence:
            lines.append("   检测到：" if language == "zh" else "   Detected in:")
            for entry in evidence:
                lines.append(f"     - {entry}")

        wrapper_hint = (
            item.get("wrapper_hint_zh", item.get("wrapper_hint", ""))
            if language == "zh"
            else item.get("wrapper_hint", "")
        )
        if wrapper_hint:
            lines.append(f"{'   备注：' if language == 'zh' else '   Note: '}{wrapper_hint}")

        missing = (
            item.get("missing_items_zh", item.get("missing_items", []))
            if language == "zh"
            else item.get("missing_items", [])
        )
        if missing:
            lines.append("   还缺少：" if language == "zh" else "   Still missing:")
            for entry in missing:
                lines.append(f"     - {entry}")
    return lines


def _table_list_cell(
    values: List[str],
    *,
    separator: str = " <br/> ",
    continuation_prefix: str = "",
) -> str:
    cleaned = [str(value).strip() for value in values if str(value).strip()]
    if not cleaned:
        return "-"
    if len(cleaned) == 1:
        return cleaned[0]
    if not continuation_prefix:
        return separator.join(cleaned)
    return separator.join([cleaned[0], *[f"{continuation_prefix}{value}" for value in cleaned[1:]]])


def _table_cell(value: str) -> str:
    cell = (value or "-").replace("|", "\\|").replace("<br>", "<br/>").replace("\n", "<br/>").strip()
    return cell or "-"


def _render_checklist_table(items: List[Dict[str, Any]], *, language: str) -> str:
    if language == "zh":
        header = "| 标记 | 状态 | 阶段 | 检测到 | 仍缺少 |"
    else:
        header = "| Marker | Status | Stage | Detected in | Still missing |"
    lines = [header, "| --- | --- | --- | --- | --- |"]
    for item in items:
        item_id = str(item.get("id") or "")
        marker = item.get("marker", "[ ]")
        label = item.get("label_zh", item.get("label", "?")) if language == "zh" else item.get("label", "?")
        status_label = item.get("status_label_zh", item.get("status_label", "")) if language == "zh" else item.get("status_label", "")
        evidence = (
            item.get("evidence_display_paths_zh", item.get("evidence_display_paths", []))
            if language == "zh"
            else item.get("evidence_display_paths", [])
        )
        missing = (
            item.get("missing_items_zh", item.get("missing_items", []))
            if language == "zh"
            else item.get("missing_items", [])
        )
        multiline_separator = " <br/><br/> " if item_id == "immune_deconvolution" else " <br/> "
        lines.append(
            "| {0} | {1} | {2} | {3} | {4} |".format(
                _table_cell(marker),
                _table_cell(status_label),
                _table_cell(label),
                _table_cell(
                    _table_list_cell(
                        list(evidence),
                        separator=multiline_separator,
                        continuation_prefix="• ",
                    )
                ),
                _table_cell(
                    _table_list_cell(
                        list(missing),
                        separator=multiline_separator,
                        continuation_prefix="• ",
                    )
                ),
            )
        )
    return "\n".join(lines)


def _terminal_stage_map(workflow_checklist: List[Dict[str, Any]], language: str = "en") -> str:
    lines = ["工作流程检查清单" if language == "zh" else "IOBRpy workflow checklist"]
    for item in workflow_checklist:
        label = item.get("label_zh") if language == "zh" else item.get("label")
        text = item.get("text_zh") if language == "zh" else item.get("text")
        lines.append(f"{item['marker']} {label}: {text}")
    return "\n".join(lines)


def _enrich_workflow_checklist_payload(
    payload: Dict[str, Any],
    *,
    strict_iobrpy_investigation_item_ids: Set[str],
) -> Dict[str, Any]:
    if not payload.get("success"):
        return payload

    base_path = Path(payload.get("path", "."))
    evidence_map = payload.get("evidence", {})
    strong_evidence = {stage_id: item.get("strong", []) for stage_id, item in evidence_map.items()}
    weak_evidence = {stage_id: item.get("weak", []) for stage_id, item in evidence_map.items()}
    external_analysis_hints = payload.get("external_analysis_hints", [])
    deconvolution_methods = payload.get("deconvolution_methods", {})
    function_detections = {item["id"]: item for item in payload.get("function_detections", [])}
    display_hints = payload.get("checklist_display_hints", {})
    updated_items: List[Dict[str, Any]] = []

    for item in payload.get("workflow_checklist", []):
        updated = dict(item)
        status = updated.get("status", "pending")
        strict_iobrpy_only = updated.get("id") in strict_iobrpy_investigation_item_ids
        allows_external_status = updated.get("id") in _EXTERNAL_STATUS_CHECKLIST_ITEM_IDS
        if status == "external" and not allows_external_status:
            status = "completed" if updated.get("checked") else "pending"
            updated["status"] = status
        raw_evidence_paths, evidence_paths = _stage_evidence_for_checklist_item(
            updated,
            strong_evidence,
            weak_evidence,
            external_analysis_hints,
            deconvolution_methods,
            display_hints,
            strict_iobrpy_investigation_item_ids=strict_iobrpy_investigation_item_ids,
        )
        evidence_display_paths = [_absolute_display_path(base_path, item) for item in evidence_paths]
        missing_en = _checklist_missing_items(updated, deconvolution_methods, language="en")
        missing_zh = _checklist_missing_items(updated, deconvolution_methods, language="zh")
        wrapper_hint_en, wrapper_hint_zh = _checklist_wrapper_hint(updated.get("id", ""), function_detections)

        source_id = "pending"
        bucket_counts = updated.get("investigation_bucket_counts", {})
        if (
            (status == "external" or bucket_counts.get("external_tool_results"))
            and not strict_iobrpy_only
            and allows_external_status
        ):
            source_id = "external_tool_results"
        elif bucket_counts.get("agent_inferred_existing_results") and not strict_iobrpy_only:
            source_id = "agent_inferred_existing_results"
        elif updated.get("checked"):
            source_id = "iobrpy_confirmed_results"

        status_label, status_label_zh = _render_checklist_status_labels(status)
        source_label, source_label_zh = _render_result_source_labels(source_id)
        evidence_text_en = f"Detected in: {', '.join(evidence_display_paths)}." if evidence_display_paths else ""
        evidence_text_zh = f"检测到的结果位置：{'、'.join(evidence_display_paths)}。" if evidence_display_paths else ""
        missing_text_en = f"Still missing: {', '.join(missing_en)}." if missing_en else ""
        missing_text_zh = f"仍缺少：{'、'.join(missing_zh)}。" if missing_zh else ""
        table_separator = " <br/><br/> " if updated.get("id") == "immune_deconvolution" else " <br/> "
        table_detected_value_en = _table_list_cell(
            list(evidence_display_paths),
            separator=table_separator,
            continuation_prefix="• ",
        )
        table_detected_value_zh = table_detected_value_en
        table_missing_value_en = _table_list_cell(
            list(missing_en),
            separator=table_separator,
            continuation_prefix="• ",
        )
        table_missing_value_zh = _table_list_cell(
            list(missing_zh),
            separator=table_separator,
            continuation_prefix="• ",
        )

        text_en = updated.get("text", "").strip()
        text_zh = updated.get("text_zh", "").strip()
        for extra in (wrapper_hint_en, evidence_text_en, missing_text_en):
            if extra and extra not in text_en:
                text_en = f"{text_en} {extra}".strip()
        for extra in (wrapper_hint_zh, evidence_text_zh, missing_text_zh):
            if extra and extra not in text_zh:
                text_zh = f"{text_zh} {extra}".strip()

        updated["marker"] = _render_checklist_marker(bool(updated.get("checked")), status)
        updated["status_label"] = status_label
        updated["status_label_zh"] = status_label_zh
        updated["result_source_id"] = source_id
        updated["result_source_label"] = source_label
        updated["result_source_label_zh"] = source_label_zh
        updated["raw_evidence_paths"] = raw_evidence_paths
        updated["raw_evidence_paths_zh"] = raw_evidence_paths
        updated["evidence_paths"] = evidence_paths
        updated["evidence_paths_zh"] = evidence_paths
        updated["evidence_display_paths"] = evidence_display_paths
        updated["evidence_display_paths_zh"] = evidence_display_paths
        updated["missing_items"] = missing_en
        updated["missing_items_zh"] = missing_zh
        updated["detected_column_value"] = table_detected_value_en
        updated["detected_column_value_zh"] = table_detected_value_zh
        updated["missing_column_value"] = table_missing_value_en
        updated["missing_column_value_zh"] = table_missing_value_zh
        updated["detected_column_source"] = "workflow_checklist[].evidence_display_paths"
        updated["missing_column_source"] = "workflow_checklist[].missing_items"
        updated["wrapper_hint"] = wrapper_hint_en
        updated["wrapper_hint_zh"] = wrapper_hint_zh
        updated["text"] = text_en
        updated["text_zh"] = text_zh
        updated_items.append(updated)

    payload = dict(payload)
    payload["workflow_checklist"] = updated_items
    payload["workflow_checklist_title"] = "Workflow Checklist"
    payload["workflow_checklist_title_zh"] = "工作流检查清单"
    payload["workflow_checklist_lines"] = _render_checklist_lines(updated_items, language="en")
    payload["workflow_checklist_lines_zh"] = _render_checklist_lines(updated_items, language="zh")
    payload["workflow_checklist_table"] = _render_checklist_table(updated_items, language="en")
    payload["workflow_checklist_table_zh"] = _render_checklist_table(updated_items, language="zh")
    payload["workflow_checklist_report"] = "\n".join([payload["workflow_checklist_title"], *payload["workflow_checklist_lines"]])
    payload["workflow_checklist_report_zh"] = "\n".join([payload["workflow_checklist_title_zh"], *payload["workflow_checklist_lines_zh"]])
    payload["default_user_facing_checklist"] = payload["workflow_checklist_table"]
    payload["default_user_facing_checklist_field"] = "workflow_checklist_table"
    payload["default_user_facing_checklist_zh"] = payload["workflow_checklist_table_zh"]
    payload["default_user_facing_checklist_field_zh"] = "workflow_checklist_table_zh"
    payload["terminal_stage_map"] = _terminal_stage_map(updated_items)
    payload["terminal_stage_map_zh"] = _terminal_stage_map(updated_items, language="zh")
    return payload


def _suggest_scan_retry_limits(max_depth: int, max_entries: int) -> Dict[str, int]:
    return {
        "max_depth": max(max_depth + 4, 10),
        "max_entries": 0 if max_entries <= 0 else max(max_entries * 2, 12000),
    }


def _scan_limit_messages(path: Path, scan: Dict[str, Any]) -> Dict[str, Any]:
    return _scan_limit_messages_impl(path, scan)
    warning_bits: List[str] = []
    warning_bits_zh: List[str] = []

    if scan.get("truncated"):
        warning_bits.append(
            f"The directory scan hit the max_entries limit ({scan.get('max_entries')}), so detections may be incomplete."
        )
        warning_bits_zh.append(
            f"目录扫描触发了 max_entries 上限（{scan.get('max_entries')}），因此结果识别可能不完整。"
        )

    depth_limited_dir_count = scan.get("depth_limited_dir_count", 0)
    if depth_limited_dir_count:
        warning_bits.append(
            f"The scan also stopped descending into {depth_limited_dir_count} subdirectories after reaching max_depth={scan.get('max_depth')}."
        )
        warning_bits_zh.append(
            f"扫描在达到 max_depth={scan.get('max_depth')} 后，还停止继续向下进入 {depth_limited_dir_count} 个子目录。"
        )

    if not warning_bits:
        return {
            "scan_retry_recommended": False,
        }

    retry_limits = _suggest_scan_retry_limits(scan.get("max_depth", 0), scan.get("max_entries", 0))
    retry_command_parts = ["iobrpy-cli", "map", "--path", str(path), "--json"]
    for focus_root in scan.get("explicit_focus_roots", []):
        retry_command_parts.extend(["--focus", focus_root])
    retry_command_parts.extend(
        [
            "--max-depth",
            str(retry_limits["max_depth"]),
            "--max-entries",
            str(retry_limits["max_entries"]),
        ]
    )
    retry_command = " ".join(retry_command_parts)

    return {
        "scan_retry_recommended": True,
        "scan_retry_limits": retry_limits,
        "scan_warning": " ".join(warning_bits),
        "scan_warning_zh": "".join(warning_bits_zh),
        "scan_retry_command": retry_command,
        "scan_retry_command_zh": retry_command,
    }


def _scan_limit_messages_impl(path: Path, scan: Dict[str, Any]) -> Dict[str, Any]:
    warning_bits: List[str] = []
    note_bits: List[str] = []
    retry_relevant = False
    entry_limit_hit = bool(scan.get("entry_limit_hit"))
    backlog_limit_hit = bool(scan.get("backlog_limit_hit"))
    entry_limit_enabled = bool(scan.get("entry_limit_enabled", True))

    if entry_limit_hit:
        warning_bits.append(
            f"The directory scan hit the max_entries limit ({scan.get('max_entries')}), so detections may be incomplete."
        )
        retry_relevant = True
    elif backlog_limit_hit:
        warning_bits.append(
            "The directory scan filled its internal entry backlog while walking this path, so detections may still be incomplete."
        )
        retry_relevant = True

    depth_limited_dir_count = scan.get("depth_limited_dir_count", 0)
    if depth_limited_dir_count:
        warning_bits.append(
            f"The scan also stopped descending into {depth_limited_dir_count} subdirectories after reaching max_depth={scan.get('max_depth')}."
        )
        retry_relevant = True

    if scan.get("full_skipped_generic_dir_count") and scan.get("focus_roots"):
        note_bits.append(
            "The deep scan intentionally kept non-focus generic branches shallow while inspecting the selected focus roots more thoroughly."
        )
    elif scan.get("quick_skipped_generic_dir_count"):
        note_bits.append(
            "Quick scan intentionally skipped some generic subdirectories to stay fast, so missing results should be interpreted cautiously."
        )

    if not warning_bits:
        payload = {
            "scan_retry_recommended": False,
        }
        if note_bits:
            payload["scan_note"] = " ".join(note_bits)
            payload["scan_note_zh"] = " ".join(note_bits)
        return payload

    retry_limits = _suggest_scan_retry_limits(scan.get("max_depth", 0), scan.get("max_entries", 0))
    retry_command_parts = ["iobrpy-cli", "map", "--path", str(path), "--json"]
    for focus_root in scan.get("explicit_focus_roots", []):
        retry_command_parts.extend(["--focus", focus_root])
    if retry_relevant:
        retry_command_parts.extend(
            [
                "--max-depth",
                str(retry_limits["max_depth"]),
                "--max-entries",
                str(retry_limits["max_entries"]),
            ]
        )
        retry_command = " ".join(retry_command_parts)
    else:
        retry_command = None

    payload = {
        "scan_retry_recommended": retry_relevant,
        "scan_retry_limits": retry_limits if retry_relevant else {},
        "scan_warning": " ".join(warning_bits),
        "scan_warning_zh": " ".join(warning_bits),
        "scan_retry_command": retry_command,
        "scan_retry_command_zh": retry_command,
        "scan_unbounded_entries": not entry_limit_enabled,
    }
    if note_bits:
        payload["scan_note"] = " ".join(note_bits)
        payload["scan_note_zh"] = " ".join(note_bits)
    return payload


def _agent_rendering_hints_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    auto_retry_performed = bool(payload.get("auto_scan_retry_performed"))
    scan_retry_command = payload.get("scan_retry_command")
    manual_scan_retry_not_recommended = bool(payload.get("manual_scan_retry_not_recommended"))
    path = str(payload.get("path", "<path>"))
    focus_roots = _recommended_deep_focus_roots(payload)
    recommended_deep_map_command = _format_map_command(
        path,
        focus_roots=focus_roots,
        max_entries=0 if focus_roots else None,
    )
    return {
        "confirm_scan_mode_before_first_map": False,
        "default_scan_mode_for_path_requests": "full",
        "always_use_full_scan_for_path_requests": True,
        "do_not_ask_user_to_choose_scan_mode": True,
        "prefer_quick_first_for_path_requests": False,
        "quick_scan_deprecated_for_agent_use": True,
        "wait_for_scan_mode_answer_before_first_map": False,
        "do_not_start_quick_scan_while_awaiting_scan_mode_choice": True,
        "quick_planning_scan_only_after_user_selects_full_scan": False,
        "prefer_interactive_scan_mode_selector": False,
        "interactive_scan_mode_tool": None,
        "interactive_scan_mode_fallback": None,
        "preferred_map_command": "iobrpy-cli map --path <path> --json",
        "preferred_quick_map_command": None,
        "preferred_full_map_command": "iobrpy-cli map --path <path> --json",
        "preferred_deep_scan_planning_command": None,
        "deep_scan_strategy": "direct_full_scan",
        "run_quick_planning_scan_before_full_scan": False,
        "use_focus_roots_from_quick_scan_when_available": False,
        "recommended_deep_focus_roots": focus_roots,
        "path_specific_recommended_deep_map_command": recommended_deep_map_command,
        "path_specific_recommended_deep_map_command_uses_repeated_focus_flags": True,
        "path_specific_recommended_deep_map_command_uses_unbounded_entries": bool(focus_roots),
        "unbounded_entry_limit_value": 0,
        "path_specific_recommended_deep_map_reason_ids": ["focus_roots_available"] if focus_roots else ["no_focus_roots_available"],
        "avoid_fixed_scan_time_promises": True,
        "scan_duration_estimates_are_relative_only": True,
        "scan_duration_caveat_en": "On very large directories or network-mounted storage, full scans can take much longer than a few minutes.",
        "scan_duration_caveat_zh": "在超大目录或网络挂载存储上，完整扫描可能会明显超过几分钟。",
        "deep_scan_command_self_retrying": True,
        "do_not_repeat_identical_deep_scan_commands": True,
        "manual_deep_rescan_only_when": [
            "user_explicitly_requests_rescan",
            "directory_changed_since_last_scan",
            "previous_deep_scan_failed",
        ],
        "scan_mode_prompt_en": "",
        "scan_mode_prompt_zh": "",
        "quick_scan_when": [],
        "full_scan_when": ["all_path_requests"],
        "scan_mode_options": [],
        "quick_followup_conditions": [],
        "adaptive_scan_limits_enabled": True,
        "adaptive_scan_parameter_policy": "let_map_choose_initial_limits",
        "manual_limit_override_only_when": [
            "user_explicitly_requests_custom_limits",
            "map_returns_scan_retry_command_or_scan_retry_limits",
        ],
        "never_lower_limits_on_manual_retry": True,
        "preferred_checklist_table_field": "workflow_checklist_table",
        "preferred_checklist_table_field_zh": "workflow_checklist_table_zh",
        "preferred_checklist_report_field": "workflow_checklist_report",
        "preferred_checklist_report_field_zh": "workflow_checklist_report_zh",
        "preferred_checklist_lines_field": "workflow_checklist_lines",
        "preferred_checklist_lines_field_zh": "workflow_checklist_lines_zh",
        "preferred_rendering_field_order": [
            "default_user_facing_checklist",
            "workflow_checklist_report",
            "workflow_checklist_table",
            "workflow_checklist_lines",
            "workflow_checklist",
        ],
        "preferred_rendering_field_order_zh": [
            "default_user_facing_checklist_zh",
            "workflow_checklist_report_zh",
            "workflow_checklist_table_zh",
            "workflow_checklist_lines_zh",
            "workflow_checklist",
        ],
        "use_compact_default_user_facing_checklist_when_present": True,
        "prefer_prebuilt_checklist_table": True,
        "checklist_table_is_authoritative": True,
        "reuse_prebuilt_checklist_table_verbatim": True,
        "do_not_replace_prebuilt_table_with_custom_summary_table": True,
        "keep_prebuilt_table_layout_even_for_tall_rows": True,
        "do_not_convert_prebuilt_table_to_field_list": True,
        "default_user_facing_checklist_format": "markdown_table",
        "do_not_convert_prebuilt_markdown_table_to_ascii_box_table": True,
        "preserve_absolute_evidence_paths_in_checklist": True,
        "prefer_line_breaks_between_multiple_evidence_paths": True,
        "immune_deconvolution_multiline_spacing": "double_break",
        "do_not_collapse_multiline_evidence_paths_inside_table_cells": True,
        "multi_value_cell_visible_fallback_separator": "•",
        "do_not_shorten_detected_paths_in_tables": True,
        "do_not_use_zh_checklist_fields_in_english_rendering": True,
        "immune_deconvolution_detected_column_source": "workflow_checklist[].evidence_display_paths",
        "immune_deconvolution_missing_column_source_en": "workflow_checklist[].missing_column_value",
        "immune_deconvolution_missing_column_source_zh": "workflow_checklist[].missing_column_value_zh",
        "tpm_matrix_detected_column_source": "workflow_checklist[].evidence_display_paths",
        "do_not_borrow_feature_scoring_outputs_for_tpm_matrix_row": True,
        "if_rebuilding_table_use_evidence_display_paths_not_deconvolution_method_labels": True,
        "use_detected_and_missing_column_value_fields_when_rebuilding_table": True,
        "if_html_line_breaks_are_stripped_keep_visible_cell_separator": True,
        "do_not_rebuild_checklist_from_deconvolution_methods": True,
        "workflow_checklist_missing_items_are_authoritative": True,
        "do_not_use_default_bundle_missing_labels_for_checklist_missing_column": True,
        "immune_deconvolution_missing_column_source": "workflow_checklist[].missing_items_zh",
        "deconvolution_methods_missing_labels_include_bayesprism": True,
        "default_bundle_missing_labels_are_only_for_tme_profile_or_runall_scope": True,
        "bayesprism_is_deconvolution_but_standalone_optional": True,
        "bayesprism_missing_display_en": f"BayesPrism ({_BAYESPRISM_OPTIONAL_NOTE_EN})",
        "bayesprism_missing_display_zh": "BayesPrism（独立可选，不由 tme_profile/runall 自动运行）",
        "prefer_missing_downstream_analysis_suggestions": True,
        "use_missing_downstream_analysis_suggestions_after_checklist": True,
        "do_not_omit_bayesprism_from_recommendations_when_missing": True,
        "use_suggested_command_details_for_parameter_prompts": True,
        "suggested_command_details_include_prompt_en_and_prompt_zh": True,
        "missing_downstream_suggestions_include_reason_en_and_reason_zh": True,
        "show_parameter_placeholders_in_natural_language": True,
        "parameter_hints_source": "src/iobrpy/RAG_MCP/iobrpy_agent_parameter_hints.json",
        "parameter_source_of_truth": "src/iobrpy/RAG_MCP/iobrpy_required_params.json",
        "do_not_drop_pending_rows_or_missing_columns": True,
        "prefer_prebuilt_checklist_report": True,
        "checklist_lines_are_authoritative": True,
        "do_not_shell_parse_mcp_tool_results": True,
        "do_not_read_client_tool_result_files_for_rendering": True,
        "do_not_request_user_approval_for_json_extraction_commands": True,
        "no_extra_confirmation_for_reading_map_payload": True,
        "confirmation_prompts_should_describe_intent_not_raw_command": True,
        "approval_prompt_style": "concise_natural_language_intent",
        "approval_prompt_example_en": "Should I read the checklist and next-step summary from the scan result?",
        "approval_prompt_example_zh": "是否要我读取刚才扫描结果中的检查表和下一步建议？",
        "result_rendering_format": "table",
        "render_full_checklist": True,
        "do_not_regroup_checklist": True,
        "preserve_evidence_and_missing_lines": True,
        "do_not_infer_sibling_results_from_shared_directory": True,
        "only_mark_method_detected_when_supported_by_concrete_file_or_content": True,
        "supported_retry_flags": ["--max-depth", "--max-entries"],
        "unsupported_retry_flags": ["--max-files"],
        "reuse_scan_retry_command_exactly": True,
        "scan_retry_command_available": bool(scan_retry_command),
        "manual_retry_needed": bool(payload.get("scan_retry_recommended")) and not auto_retry_performed and not manual_scan_retry_not_recommended,
        "manual_retry_already_performed": auto_retry_performed,
        "manual_scan_retry_not_recommended": manual_scan_retry_not_recommended,
        "retry_exhausted": bool(payload.get("retry_exhausted")),
        "prefer_integrated_checklist_over_custom_summary": True,
    }


def format_pipeline_map_report(payload: Dict[str, Any]) -> str:
    if not payload.get("success"):
        return f"Pipeline map: failed\n- {payload.get('error', 'Unknown error')}"

    lines = [payload.get("terminal_stage_map", "IOBRpy workflow checklist")]
    lines.append("")
    scan_warning = payload.get("scan_warning")
    if scan_warning:
        lines.append(f"Scan warning: {scan_warning}")
        retry_command = payload.get("scan_retry_command")
        if retry_command:
            lines.append(f"Retry with: {retry_command}")
        lines.append("")
    lines.append(f"Current stage: {payload.get('current_label')}")
    confidence = payload.get("current_stage_confidence")
    if confidence and confidence != "none":
        lines.append(f"Confidence: {confidence}")

    scenario = payload.get("scenario", {})
    if scenario:
        lines.append(f"Scenario: {scenario.get('label')}")
        summary = scenario.get("summary")
        if summary:
            lines.append(f"Summary: {summary}")

    roadmap = payload.get("roadmap", {})
    position_summary = roadmap.get("position_summary")
    if position_summary:
        lines.append(f"Roadmap position: {position_summary}")
    next_summaries = payload.get("next_step_summaries", [])
    if next_summaries:
        lines.append("What can be done next:")
        for item in next_summaries:
            lines.append(f"- {item.get('text')}")

    missing_suggestions = payload.get("missing_downstream_analysis_suggestions", [])
    if missing_suggestions:
        lines.append("Missing downstream analyses to consider:")
        for item in missing_suggestions:
            labels = ", ".join(item.get("missing_analysis_labels", []))
            suffix = f": {labels}" if labels else ""
            lines.append(f"- {item.get('label')}{suffix}")

    recommended_action = payload.get("recommended_action", {})
    if recommended_action:
        lines.append("Recommended next action:")
        label = recommended_action.get("label")
        if label:
            lines.append(f"- {label}")
        reason = recommended_action.get("why")
        if reason:
            lines.append(f"- Why: {reason}")

    external_analysis_hints = payload.get("external_analysis_hints", [])
    if external_analysis_hints:
        lines.append("Other analysis hints:")
        for hint in external_analysis_hints:
            lines.append(f"- {hint.get('label')}: {hint.get('description')}")

    investigation = payload.get("existing_result_investigation", {})
    if investigation.get("enabled") and not payload.get("existing_result_investigation_applied"):
        summary = investigation.get("summary", {})
        lines.append("Existing-result investigation:")
        if investigation.get("auto_scan_retry_performed"):
            lines.append("- Expanded scan limits automatically before investigation.")
        lines.append(f"- Findings: {summary.get('total_findings', 0)}")
        for bucket in investigation.get("result_buckets", []):
            bucket_id = bucket.get("id")
            count = summary.get("bucket_counts", {}).get(bucket_id, 0)
            if count:
                lines.append(f"- {bucket.get('label', bucket_id)}: {count}")

    choice_details = scenario.get("choice_details", [])
    if choice_details:
        lines.append("Choices:")
        for detail in choice_details:
            label = detail.get("label", detail.get("id", "choice"))
            reason = detail.get("when_to_choose")
            if reason:
                lines.append(f"- {label}: {reason}")
            else:
                lines.append(f"- {label}")

    roadmap_targets = [roadmap.get("local_doc"), roadmap.get("url")]
    roadmap_targets = [target for target in roadmap_targets if target]
    if roadmap_targets:
        lines.append("Roadmap references:")
        for target in roadmap_targets:
            lines.append(f"- {target}")

    return "\n".join(lines).strip()
