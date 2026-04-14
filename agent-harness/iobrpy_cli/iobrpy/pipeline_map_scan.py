"""
Scanning and path-matching helpers for pipeline_map.
"""

from __future__ import annotations

import fnmatch
import os
import re
from collections import deque
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


_QUICK_GENERIC_DESCENT_LIMIT = 128
_DEFAULT_MAP_MAX_DEPTH = 8
_DEFAULT_MAP_MAX_ENTRIES = 8000
_BULK_STAGE_CHILD_PROBE_LIMIT = 6
_BULK_STAGE_CHILD_ENTRY_LIMIT = 8
_QUICK_MIN_MAX_DEPTH = 10
_QUICK_MIN_MAX_ENTRIES = 12000
_FULL_GENERIC_DEEP_SCAN_LIMIT = 3
_EXPLICIT_FOCUS_GENERIC_DEEP_SCAN_LIMIT = 1
_ADAPTIVE_SCAN_PROBE_LIMIT = 256
_ADAPTIVE_LARGE_TOP_LEVEL_THRESHOLD = 24
_ADAPTIVE_VERY_LARGE_TOP_LEVEL_THRESHOLD = 96
_ADAPTIVE_STAGE_DIR_THRESHOLD = 4
_ADAPTIVE_STAGE_SCAN_TOKENS = (
    "fastq",
    "fq",
    "qc",
    "fastp",
    "salmon",
    "star",
    "tpm",
    "signature",
    "tme",
    "lr",
    "hla",
    "mixcr",
    "trust4",
    "tcrbcr",
    "tcr",
    "bcr",
    "airr",
    "profile",
    "cluster",
)

_BULK_STAGE_DIR_RULES = {
    "trust4": {"exact_parts": {"07-tcrbcr", "tcrbcr", "tcr", "bcr"}, "semantic_tokens": {"tcrbcr", "trust4"}, "dir_limit_quick": 12, "dir_limit_full": 48, "file_limit_quick": 24, "file_limit_full": 96},
    "raw": {"exact_parts": {"00-fastq", "01-fastq", "02-fastq", "03-fastq"}, "semantic_tokens": {"fastq", "fq"}, "dir_limit_quick": 6, "dir_limit_full": 20, "file_limit_quick": 12, "file_limit_full": 48},
    "qc": {"exact_parts": {"01-qc", "03-fastp"}, "semantic_tokens": {"qc", "fastp"}, "dir_limit_quick": 6, "dir_limit_full": 20, "file_limit_quick": 12, "file_limit_full": 48},
    "salmon": {"exact_parts": {"02-salmon"}, "semantic_tokens": {"salmon"}, "dir_limit_quick": 8, "dir_limit_full": 24, "file_limit_quick": 16, "file_limit_full": 64},
    "star": {"exact_parts": {"02-star"}, "semantic_tokens": {"star"}, "dir_limit_quick": 8, "dir_limit_full": 24, "file_limit_quick": 16, "file_limit_full": 64},
    "hla": {"exact_parts": set(), "semantic_tokens": {"hla", "spechla", "extracthlaread"}, "dir_limit_quick": 8, "dir_limit_full": 24, "file_limit_quick": 16, "file_limit_full": 64},
}


class _EntryMatcher:
    def __init__(self, entries: Iterable[str]):
        normalized_entries = [_normalize_pattern(str(entry)).rstrip("/") for entry in entries if str(entry).strip()]
        self._rows = [
            (entry, Path(entry).name, entry.lower(), Path(entry).name.lower())
            for entry in normalized_entries
        ]
        self._nonempty_dirs = _nonempty_directory_entries(normalized_entries)
        self._match_cache: Dict[Tuple[str, int], List[str]] = {}

    def match_pattern(self, pattern: str, *, limit: int = 3) -> List[str]:
        normalized = _normalize_pattern(pattern)
        cache_key = (normalized, limit)
        cached = self._match_cache.get(cache_key)
        if cached is not None:
            return cached[:]

        normalized_lower = normalized.lower()
        matches: List[str] = []
        for entry, basename, entry_lower, basename_lower in self._rows:
            if (
                fnmatch.fnmatch(entry, normalized)
                or fnmatch.fnmatch(basename, normalized)
                or fnmatch.fnmatch(entry_lower, normalized_lower)
                or fnmatch.fnmatch(basename_lower, normalized_lower)
            ):
                if _is_empty_directory_match(entry, self._nonempty_dirs):
                    continue
                matches.append(entry)
                if len(matches) >= limit:
                    break

        self._match_cache[cache_key] = matches[:]
        return matches

    def collect(self, patterns: Iterable[str], *, limit: int = 3) -> List[str]:
        matches: List[str] = []
        for pattern in patterns:
            for match in self.match_pattern(pattern, limit=limit):
                if match not in matches:
                    matches.append(match)
            if len(matches) >= limit:
                break
        return matches


def _normalize_pattern(pattern: str) -> str:
    return pattern.replace("\\", "/")


def _normalize_relpath(relpath: Path) -> str:
    normalized = relpath.as_posix()
    if normalized.startswith("./"):
        return normalized[2:]
    return normalized


def _nonempty_directory_entries(entries: Iterable[str]) -> Set[str]:
    nonempty: Set[str] = set()
    for entry in entries:
        normalized = _normalize_pattern(str(entry)).rstrip("/")
        if not normalized:
            continue
        parts = Path(normalized).parts
        for idx in range(1, len(parts)):
            nonempty.add(Path(*parts[:idx]).as_posix())
    return nonempty


def _is_empty_directory_match(entry: str, nonempty_dirs: Set[str]) -> bool:
    normalized = _normalize_pattern(str(entry)).rstrip("/")
    if not normalized:
        return False
    if normalized in nonempty_dirs:
        return False
    return not bool(Path(normalized).suffixes)


def _normalize_token(token: str) -> str:
    return "".join(ch for ch in token.lower() if ch.isalnum())


def _path_parts_lower(relpath: str) -> Tuple[str, ...]:
    return tuple(part.lower() for part in Path(relpath).parts)


def _part_tokens(part: str) -> Tuple[str, ...]:
    lower = part.lower()
    tokens = tuple(token for token in re.split(r"[^a-z0-9]+", lower) if token)
    normalized = _normalize_token(lower)
    if normalized and normalized not in tokens:
        return tokens + (normalized,)
    return tokens


def _path_has_part_token(relpath: str, *tokens: str) -> bool:
    wanted = {_normalize_token(token) for token in tokens if token}
    if not wanted:
        return False
    for part in _path_parts_lower(relpath):
        if wanted.intersection(_part_tokens(part)):
            return True
    return False


def _path_has_part_tokens(relpath: str, *tokens: str) -> bool:
    wanted = {_normalize_token(token) for token in tokens if token}
    if not wanted:
        return False
    for part in _path_parts_lower(relpath):
        if wanted.issubset(set(_part_tokens(part))):
            return True
    return False


def _path_has_exact_part(relpath: str, *parts: str) -> bool:
    path_parts = set(_path_parts_lower(relpath))
    return any(part.lower() in path_parts for part in parts)


def _path_has_keyword(relpath: str, keywords: Iterable[str]) -> bool:
    lower = relpath.lower()
    return any(keyword in lower for keyword in keywords)


def _normalized_path_lower(relpath: str) -> str:
    return _normalize_pattern(relpath).lower()


def _has_star_result_context(relpath: str) -> bool:
    lower = _normalized_path_lower(relpath)
    if _path_has_exact_part(relpath, "02-star", "03-star", "04-star", "star"):
        return True
    if _path_has_part_token(relpath, "star"):
        return True
    return any(
        token in lower
        for token in (
            "readspergene.out.tab",
            "aligned.sortedbycoord.out.bam",
            "star.count",
            "star_count",
            "merge_star_count",
            "__stargenome",
            "__starpass1",
        )
    )


def _looks_like_star_count_matrix_name(relpath: str) -> bool:
    lower = _normalized_path_lower(relpath)
    name = Path(relpath).name.lower()
    if not any(lower.endswith(ext) for ext in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")):
        return False
    if name in {"biosample_count.txt", "biosample_counts.txt"}:
        return False
    if any(token in lower for token in ("star.count", "star_count", "readspergene", "featurecounts", "featurecount")):
        return True
    if not _has_star_result_context(relpath):
        return False
    return any(token in name for token in ("count_matrix", "counts_matrix", "gene_count", "read_count"))


def _has_hla_result_context(relpath: str) -> bool:
    return _path_has_part_token(relpath, "hla", "spechla", "extracthlaread", "hlatyping", "hla_typing")


def _is_priority_scan_dir(relpath: str) -> bool:
    lower = _normalized_path_lower(relpath)
    standard_dir_hits = {
        "00-fastq",
        "01-fastq",
        "02-fastq",
        "03-fastq",
        "01-qc",
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
        "tcrbcr",
    }
    return (
        any(token in lower for token in standard_dir_hits)
        or _path_has_part_token(relpath, "runall", "salmon", "star", "hla", "mixcr", "trust4", "tcrbcr", "spechla")
        or _path_has_part_tokens(relpath, "tme", "profile")
    )


def _matches_bulk_stage_file(relpath: str, kind: str, *, for_classification: bool = False) -> bool:
    lower = relpath.lower()
    name = Path(relpath).name.lower()
    if kind == "trust4":
        return (
            "trust4" in lower
            or "trust_" in name
            or name.endswith("_report.tsv")
            or name.endswith("_airr.tsv")
            or name.endswith("_airr_align.tsv")
            or name.endswith("_cdr3.out")
            or name.endswith("_final.out")
            or name.endswith("_raw.out")
            or name.endswith("_annot.fa")
            or name.endswith("_assembled_reads.fa")
            or name.endswith(".trust4.done")
        )
    if kind == "raw":
        return any(lower.endswith(ext) for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"))
    if kind == "qc":
        return "fastp" in lower or "multiqc" in lower or name == "multiqc_report"
    if kind == "salmon":
        salmon_specific = (
            name == "quant.sf"
            or "salmon" in lower
            or "tximport" in lower
            or "merged_salmon" in lower
        )
        if salmon_specific or for_classification:
            return salmon_specific
        return (
            lower.endswith(".rdata")
            or ("count" in name and "mcpcounter" not in name and any(lower.endswith(ext) for ext in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")))
            or ("tpm" in name and any(lower.endswith(ext) for ext in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")))
        )
    if kind == "star":
        star_specific = (
            "readspergene.out.tab" in lower
            or "aligned.sortedbycoord.out.bam" in lower
            or (name.endswith(".task.complete") and _has_star_result_context(relpath))
            or "merge_star_count" in lower
            or "readspergene" in lower
        )
        if star_specific or for_classification:
            return star_specific
        return _looks_like_star_count_matrix_name(relpath)
    if kind == "hla":
        hla_context = _has_hla_result_context(relpath)
        return (
            "hla_result_merged" in lower
            or "hlafinal.type" in lower
            or "hla.result" in lower
            or "hla.results" in lower
            or ("cohort" in lower and "genotype" in lower)
            or ".spechla.done" in lower
            or ".extracthlaread.done" in lower
            or (hla_context and name.endswith(".task.complete"))
            or (hla_context and any(name.endswith(suffix) for suffix in ("_result.tsv", "_result.csv", "_result.txt")))
            or ("extracthlaread" in lower and any(lower.endswith(ext) for ext in ("_1.fq.gz", "_2.fq.gz", "_1.fastq.gz", "_2.fastq.gz")))
        )
    return False


def _probe_bulk_stage_child_kind(dir_path: Path, kind: str) -> bool:
    try:
        with os.scandir(dir_path) as iterator:
            for index, item in enumerate(iterator):
                if index >= _BULK_STAGE_CHILD_ENTRY_LIMIT:
                    break
                name = item.name
                relpath = f"{dir_path.name}/{name}"
                try:
                    is_dir = item.is_dir(follow_symlinks=False)
                except OSError:
                    is_dir = False
                if _matches_bulk_stage_file(relpath, kind, for_classification=True):
                    return True
                if kind == "hla" and is_dir and _path_has_part_token(name, "spechla", "extracthlaread"):
                    return True
    except OSError:
        return False
    return False


def _bulk_stage_dir_kind(
    current_path: Path,
    rel_dir: Path,
    dirnames: Optional[List[Tuple[str, Path, str]]] = None,
    filenames: Optional[List[Tuple[str, Path, str]]] = None,
) -> Optional[str]:
    dirnames = dirnames or []
    filenames = filenames or []
    relpath = rel_dir.as_posix()
    parts = {part.lower() for part in rel_dir.parts if part}
    if relpath in {"", "."} and not parts:
        return None
    for kind, spec in _BULK_STAGE_DIR_RULES.items():
        if parts.intersection(spec["exact_parts"]) or _path_has_part_token(relpath, *spec["semantic_tokens"]):
            return kind
    for kind in _BULK_STAGE_DIR_RULES:
        if any(_matches_bulk_stage_file(normalized_relpath, kind, for_classification=True) for _name, _relpath, normalized_relpath in filenames):
            return kind
        if kind == "hla" and any(_path_has_part_token(name, "spechla", "extracthlaread") for name, _relpath, _normalized_relpath in dirnames):
            return kind
    for kind in _BULK_STAGE_DIR_RULES:
        for index, (name, _relpath, _normalized_relpath) in enumerate(sorted(dirnames, key=lambda item: item[0])):
            if index >= _BULK_STAGE_CHILD_PROBE_LIMIT:
                break
            if _probe_bulk_stage_child_kind(current_path / name, kind):
                return kind
    return None


def _bulk_stage_scan_limits(kind: str, *, quick: bool) -> Tuple[int, int]:
    spec = _BULK_STAGE_DIR_RULES[kind]
    return (
        int(spec["dir_limit_quick"] if quick else spec["dir_limit_full"]),
        int(spec["file_limit_quick"] if quick else spec["file_limit_full"]),
    )


def _is_hla_intermediate_read_path(relpath: str) -> bool:
    return _path_has_part_token(relpath, "extracthlaread", "spechla")


def _scan_file_priority_rank(relpath: str) -> Tuple[int, int, str]:
    lower = relpath.lower()
    name = Path(relpath).name.lower()
    summary_keywords = (
        "hla_result_merged",
        "prepare_salmon",
        "count2tpm",
        "calculate_sig_score",
        "deconvo_merged",
        "lr_cal",
        "tpm_matrix",
        "trust4_immdata",
        "trust4_immune_indices",
        "trust_",
        "airr",
        "cibersort",
        "estimate",
        "mcpcounter",
        "quantiseq",
        "epic",
        "bayesprism",
        "immune_indices",
        "immdata",
    )
    native_summary_names = {
        "quant.sf",
        "lr_cal.csv",
        "lr_cal.tsv",
        "calculate_sig_score.csv",
        "prepare_salmon.csv",
        "count2tpm.csv",
        "deconvo_merged.csv",
        "deconvo_merged.tsv",
        "deconvo_merged.csv.gz",
        "deconvo_merged.tsv.gz",
        "ips_results.csv",
        "ips_results.tsv",
        "ips_results.csv.gz",
        "ips_results.tsv.gz",
        "estimate_results.csv",
        "estimate_results.tsv",
        "estimate_results.csv.gz",
        "estimate_results.tsv.gz",
        "mcpcounter_results.csv",
        "mcpcounter_results.tsv",
        "mcpcounter_results.csv.gz",
        "mcpcounter_results.tsv.gz",
        "quantiseq_results.csv",
        "quantiseq_results.tsv",
        "quantiseq_results.csv.gz",
        "quantiseq_results.tsv.gz",
        "epic_results.csv",
        "epic_results.tsv",
        "epic_results.csv.gz",
        "epic_results.tsv.gz",
        "cibersort_results.csv",
        "cibersort_results.tsv",
        "cibersort_results.csv.gz",
        "cibersort_results.tsv.gz",
        "hla_result_merged.txt",
        "trust4_immdata.csv",
        "trust4_immune_indices.csv",
    }
    if name in native_summary_names or any(keyword in lower for keyword in summary_keywords):
        return (0, len(Path(relpath).parts), relpath)
    if _looks_like_star_count_matrix_name(relpath) or (_has_hla_result_context(relpath) and name.endswith(("_result.tsv", "_result.csv", "_result.txt"))):
        return (0, len(Path(relpath).parts), relpath)
    if any(keyword in lower for keyword in ("tpm", "matrix", "merged")):
        return (1, len(Path(relpath).parts), relpath)
    if any(keyword in lower for keyword in ("01-signatures", "02-tme", "03-lr_cal", "04-signatures", "05-tme", "06-lr_cal", "spec hla", "spechla", "hla.result", "hlafinal.type")):
        return (2, len(Path(relpath).parts), relpath)
    return (3, len(Path(relpath).parts), relpath)


def _probe_directory_entries(path: Path, *, limit: int) -> List[Dict[str, Any]]:
    signature: List[Dict[str, Any]] = []
    try:
        with os.scandir(path) as iterator:
            for index, item in enumerate(iterator):
                if index >= limit:
                    break
                try:
                    is_dir = item.is_dir(follow_symlinks=False)
                except OSError:
                    is_dir = False
                try:
                    mtime_ns = item.stat(follow_symlinks=False).st_mtime_ns
                except OSError:
                    mtime_ns = None
                signature.append(
                    {
                        "name": item.name,
                        "is_dir": is_dir,
                        "mtime_ns": mtime_ns,
                    }
                )
    except OSError:
        return []
    return signature


def _dedupe_normalized_paths(paths: Iterable[str]) -> List[str]:
    deduped: List[str] = []
    seen = set()
    for path in paths:
        normalized = _normalize_relpath(Path(path)).rstrip("/")
        if not normalized or normalized == "." or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def _resolve_focus_roots(root: Path, focus: Optional[Iterable[str]]) -> List[str]:
    if not focus:
        return []

    resolved_focus_roots: List[str] = []
    root_resolved = root.resolve()
    for raw_focus in focus:
        if not raw_focus:
            continue
        focus_path = Path(raw_focus).expanduser()
        if focus_path.is_absolute():
            try:
                relpath = focus_path.resolve().relative_to(root_resolved)
            except ValueError:
                continue
        else:
            relpath = focus_path
        normalized = _normalize_relpath(relpath).rstrip("/")
        if normalized and normalized != ".":
            resolved_focus_roots.append(normalized)

    return _dedupe_normalized_paths(resolved_focus_roots)


def _relpath_is_under_focus(relpath: str, focus_roots: Iterable[str]) -> bool:
    normalized = _normalize_relpath(Path(relpath)).rstrip("/")
    for focus_root in focus_roots:
        if normalized == focus_root or normalized.startswith(f"{focus_root}/"):
            return True
    return False


def _auto_focus_root_candidates(dir_entries: Iterable[Tuple[str, Path, str]]) -> List[str]:
    candidates = [
        normalized_relpath
        for _name, _relpath, normalized_relpath in dir_entries
        if _is_priority_scan_dir(normalized_relpath)
    ]
    return _dedupe_normalized_paths(candidates)


def _should_promote_directory_focus(
    rel_dir: Path,
    dirnames: List[Tuple[str, Path, str]],
    bulk_stage_kind: Optional[str],
) -> bool:
    if not rel_dir.parts:
        return False
    if bulk_stage_kind:
        return True
    nested_priority_count = sum(
        1 for _name, _relpath, normalized_relpath in dirnames if _is_priority_scan_dir(normalized_relpath)
    )
    return nested_priority_count >= 2


def _suggest_initial_scan_limits(
    path: Path,
    *,
    max_depth: int,
    max_entries: int,
    quick: bool,
    focus: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    requested = {
        "max_depth": max_depth,
        "max_entries": max_entries,
        "quick": quick,
        "focus": list(focus or []),
    }
    if not path.exists() or not path.is_dir():
        return {
            "requested": requested,
        "effective": {
            "max_depth": max_depth,
            "max_entries": max_entries,
            "quick": quick,
            "focus": list(focus or []),
        },
            "applied": False,
            "reason_ids": [],
            "probe": {},
        }

    signature = _probe_directory_entries(path.resolve(), limit=_ADAPTIVE_SCAN_PROBE_LIMIT)
    dir_names = [str(item.get("name", "")).lower() for item in signature if item.get("is_dir")]
    file_names = [str(item.get("name", "")).lower() for item in signature if not item.get("is_dir")]
    top_level_dir_count = len(dir_names)
    top_level_file_count = len(file_names)
    stage_like_dir_count = sum(
        1 for name in dir_names if any(token in name for token in _ADAPTIVE_STAGE_SCAN_TOKENS)
    )
    bulk_stage_dir_count = sum(
        1 for name in dir_names if any(token in name for token in ("fastq", "fastp", "qc", "salmon", "star", "hla", "mixcr", "trust4", "tcrbcr", "tcr", "bcr"))
    )
    merged_result_file_count = sum(
        1 for name in file_names if any(token in name for token in ("merged", "tpm", "count", "signature", "deconvo", "lr_cal", "hla_result"))
    )
    probe_limit_reached = len(signature) >= _ADAPTIVE_SCAN_PROBE_LIMIT

    reason_ids: List[str] = []
    if probe_limit_reached:
        reason_ids.append("many_top_level_entries")
    if top_level_dir_count >= _ADAPTIVE_LARGE_TOP_LEVEL_THRESHOLD:
        reason_ids.append("many_top_level_directories")
    if stage_like_dir_count >= _ADAPTIVE_STAGE_DIR_THRESHOLD:
        reason_ids.append("many_stage_like_directories")
    if bulk_stage_dir_count >= 2:
        reason_ids.append("multiple_bulk_stage_directories")
    if merged_result_file_count >= 3:
        reason_ids.append("multiple_summary_result_files")

    effective_depth = max_depth
    effective_entries = max_entries
    unbounded_entries_requested = max_entries <= 0

    if quick:
        if reason_ids and effective_depth <= _QUICK_MIN_MAX_DEPTH:
            effective_depth = max(effective_depth, 12)
        if (
            probe_limit_reached
            or top_level_dir_count >= _ADAPTIVE_VERY_LARGE_TOP_LEVEL_THRESHOLD
            or stage_like_dir_count >= (_ADAPTIVE_STAGE_DIR_THRESHOLD + 2)
        ):
            effective_depth = max(effective_depth, 14)
        if not unbounded_entries_requested and reason_ids and effective_entries <= _QUICK_MIN_MAX_ENTRIES:
            effective_entries = max(effective_entries, 16000)
        if (
            probe_limit_reached
            or top_level_dir_count >= _ADAPTIVE_VERY_LARGE_TOP_LEVEL_THRESHOLD
            or (bulk_stage_dir_count >= 3 and merged_result_file_count >= 2)
        ) and not unbounded_entries_requested:
            effective_entries = max(effective_entries, 24000)
    else:
        if reason_ids and effective_depth <= _DEFAULT_MAP_MAX_DEPTH:
            effective_depth = max(effective_depth, 12)
        if (
            probe_limit_reached
            or top_level_dir_count >= _ADAPTIVE_VERY_LARGE_TOP_LEVEL_THRESHOLD
            or stage_like_dir_count >= (_ADAPTIVE_STAGE_DIR_THRESHOLD + 2)
        ):
            effective_depth = max(effective_depth, 16)
        if not unbounded_entries_requested and reason_ids and effective_entries <= _DEFAULT_MAP_MAX_ENTRIES:
            effective_entries = max(effective_entries, 16000)
        if (
            probe_limit_reached
            or top_level_dir_count >= _ADAPTIVE_VERY_LARGE_TOP_LEVEL_THRESHOLD
            or (bulk_stage_dir_count >= 3 and merged_result_file_count >= 2)
        ) and not unbounded_entries_requested:
            effective_entries = max(effective_entries, 24000)

    return {
        "requested": requested,
        "effective": {
            "max_depth": effective_depth,
            "max_entries": effective_entries,
            "quick": quick,
            "focus": list(focus or []),
        },
        "applied": effective_depth != max_depth or effective_entries != max_entries,
        "reason_ids": reason_ids,
        "probe": {
            "top_level_dir_count": top_level_dir_count,
            "top_level_file_count": top_level_file_count,
            "stage_like_dir_count": stage_like_dir_count,
            "bulk_stage_dir_count": bulk_stage_dir_count,
            "merged_result_file_count": merged_result_file_count,
            "probe_limit_reached": probe_limit_reached,
            "probe_sample_count": len(signature),
        },
    }


def _matches_any_pattern(entry: str, patterns: Iterable[str]) -> bool:
    basename = Path(entry).name
    entry_lower = entry.lower()
    basename_lower = basename.lower()
    for pattern in patterns:
        normalized = _normalize_pattern(pattern)
        normalized_lower = normalized.lower()
        if (
            fnmatch.fnmatch(entry, normalized)
            or fnmatch.fnmatch(basename, normalized)
            or fnmatch.fnmatch(entry_lower, normalized_lower)
            or fnmatch.fnmatch(basename_lower, normalized_lower)
        ):
            return True
    return False


def _build_entry_matcher(entries: Iterable[str]) -> _EntryMatcher:
    return _EntryMatcher(entries)


def _scan_path_entries(
    path: Path,
    *,
    max_depth: int = 8,
    max_entries: int = 8000,
    priority_patterns: Optional[Iterable[str]] = None,
    quick: bool = False,
    focus: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    entry_limit_enabled = max_entries > 0
    entry_limit: Optional[int] = max_entries if entry_limit_enabled else None
    if path.is_file():
        return {
            "entries": [_normalize_relpath(Path(path.name))],
            "entry_count": 1,
            "truncated": False,
            "max_depth": 0,
            "max_entries": max_entries,
            "entry_limit_enabled": entry_limit_enabled,
            "entry_limit_hit": False,
            "backlog_limit_hit": False,
            "unbounded_entry_scan": not entry_limit_enabled,
            "truncated_reason_ids": [],
            "depth_limited_dir_count": 0,
            "depth_limited_dirs": [],
            "strategy": "quick" if quick else "full",
            "quick_skipped_generic_dir_count": 0,
            "full_skipped_generic_dir_count": 0,
            "sampled_bulk_stage_dir_count": 0,
            "sampled_bulk_child_dir_count": 0,
            "sampled_bulk_file_count": 0,
            "focus_mode": "explicit" if focus else "none",
            "focus_roots": _resolve_focus_roots(path.parent if path.parent.exists() else Path("."), focus),
            "explicit_focus_roots": _resolve_focus_roots(path.parent if path.parent.exists() else Path("."), focus),
            "auto_focus_roots": [],
            "promoted_focus_root_count": 0,
            "generic_deep_scan_limit": None,
        }

    entries: List[str] = []
    priority_dir_backlog: List[str] = []
    dir_backlog: List[str] = []
    priority_file_backlog: List[str] = []
    file_backlog: List[str] = []
    truncated = False
    entry_limit_hit = False
    backlog_limit_hit = False
    depth_limited_dir_count = 0
    depth_limited_dirs: List[str] = []
    root = path.resolve()
    priority_queue = deque([(root, Path(), 0)])
    generic_queue: deque[Tuple[Path, Path, int]] = deque()
    file_backlog_limit = entry_limit * 4 if entry_limit is not None else None
    normalized_priority_patterns = list(priority_patterns or [])
    quick_skipped_generic_dir_count = 0
    full_skipped_generic_dir_count = 0
    quick_generic_descents = 0
    sampled_bulk_stage_dir_count = 0
    sampled_bulk_child_dir_count = 0
    sampled_bulk_file_count = 0
    explicit_focus_roots = _resolve_focus_roots(root, focus)
    active_focus_roots = list(explicit_focus_roots)
    auto_focus_roots: List[str] = []
    promoted_focus_root_count = 0
    generic_deep_scan_limit: Optional[int] = None
    focus_mode = "explicit" if explicit_focus_roots else "none"

    while priority_queue or generic_queue:
        if priority_queue:
            current_path, rel_dir, depth = priority_queue.popleft()
        else:
            current_path, rel_dir, depth = generic_queue.popleft()
        try:
            with os.scandir(current_path) as iterator:
                dirnames: List[Tuple[str, Path, str]] = []
                filenames: List[Tuple[str, Path, str]] = []
                for item in iterator:
                    try:
                        is_dir = item.is_dir(follow_symlinks=False)
                    except OSError:
                        continue

                    name = item.name
                    relpath = Path(name) if not rel_dir.parts else rel_dir / name
                    normalized_relpath = _normalize_relpath(relpath)
                    if is_dir:
                        dirnames.append((name, relpath, normalized_relpath))
                        continue
                    filenames.append((name, relpath, normalized_relpath))
        except OSError:
            continue

        bulk_stage_kind = _bulk_stage_dir_kind(current_path, rel_dir, dirnames, filenames)
        sampled_dir_limit: Optional[int] = None
        sampled_file_limit: Optional[int] = None
        if bulk_stage_kind:
            sampled_dir_limit, sampled_file_limit = _bulk_stage_scan_limits(bulk_stage_kind, quick=quick)
            if len(dirnames) > sampled_dir_limit or len(filenames) > sampled_file_limit:
                sampled_bulk_stage_dir_count += 1

        if not rel_dir.parts:
            auto_focus_roots = _auto_focus_root_candidates(dirnames)
            if not active_focus_roots and auto_focus_roots:
                active_focus_roots = list(auto_focus_roots)
                focus_mode = "auto"
            if not quick and active_focus_roots:
                generic_deep_scan_limit = min(
                    max_depth,
                    _EXPLICIT_FOCUS_GENERIC_DEEP_SCAN_LIMIT if explicit_focus_roots else _FULL_GENERIC_DEEP_SCAN_LIMIT,
                )
        elif (
            not quick
            and not _relpath_is_under_focus(rel_dir.as_posix(), active_focus_roots)
            and _should_promote_directory_focus(rel_dir, dirnames, bulk_stage_kind)
        ):
            normalized_rel_dir = _normalize_relpath(rel_dir).rstrip("/")
            if normalized_rel_dir and normalized_rel_dir not in active_focus_roots:
                active_focus_roots.append(normalized_rel_dir)
                promoted_focus_root_count += 1
                if focus_mode == "none":
                    focus_mode = "auto"
            if active_focus_roots:
                generic_deep_scan_limit = min(
                    max_depth,
                    _EXPLICIT_FOCUS_GENERIC_DEEP_SCAN_LIMIT if explicit_focus_roots else _FULL_GENERIC_DEEP_SCAN_LIMIT,
                )

        for index, (name, relpath, normalized_relpath) in enumerate(sorted(dirnames, key=lambda item: item[0])):
            is_focus_dir = _relpath_is_under_focus(normalized_relpath, active_focus_roots)
            is_priority_dir = _is_priority_scan_dir(normalized_relpath) or is_focus_dir
            target_dir_backlog = priority_dir_backlog if is_priority_dir else dir_backlog
            if file_backlog_limit is None or len(target_dir_backlog) < file_backlog_limit:
                target_dir_backlog.append(normalized_relpath)
            else:
                truncated = True
                backlog_limit_hit = True
            if depth < max_depth:
                if sampled_dir_limit is not None and index >= sampled_dir_limit:
                    sampled_bulk_child_dir_count += 1
                    continue
                queue_entry = (current_path / name, relpath, depth + 1)
                if is_priority_dir:
                    priority_queue.append(queue_entry)
                elif quick and depth >= 1:
                    if quick_generic_descents >= _QUICK_GENERIC_DESCENT_LIMIT:
                        quick_skipped_generic_dir_count += 1
                        truncated = True
                    else:
                        generic_queue.append(queue_entry)
                        quick_generic_descents += 1
                elif not quick and generic_deep_scan_limit is not None and depth >= generic_deep_scan_limit:
                    full_skipped_generic_dir_count += 1
                    truncated = True
                else:
                    generic_queue.append(queue_entry)
            else:
                depth_limited_dir_count += 1
                if len(depth_limited_dirs) < 20:
                    depth_limited_dirs.append(normalized_relpath)

        file_matches_kept = 0
        for _name, _relpath, normalized_relpath in sorted(filenames, key=lambda item: _scan_file_priority_rank(item[2])):
            if bulk_stage_kind and not _matches_bulk_stage_file(normalized_relpath, bulk_stage_kind):
                continue
            if sampled_file_limit is not None and file_matches_kept >= sampled_file_limit:
                sampled_bulk_file_count += 1
                continue
            file_matches_kept += 1
            target_backlog = (
                priority_file_backlog
                if normalized_priority_patterns and _matches_any_pattern(normalized_relpath, normalized_priority_patterns)
                else file_backlog
            )
            if file_backlog_limit is not None and len(target_backlog) >= file_backlog_limit:
                truncated = True
                backlog_limit_hit = True
                continue
            target_backlog.append(normalized_relpath)

    priority_dir_backlog.sort()
    priority_file_backlog.sort(key=_scan_file_priority_rank)
    dir_backlog.sort()
    file_backlog.sort(key=_scan_file_priority_rank)

    for bucket in (priority_dir_backlog, priority_file_backlog, dir_backlog, file_backlog):
        for entry in bucket:
            if entry_limit is not None and len(entries) >= entry_limit:
                truncated = True
                entry_limit_hit = True
                break
            entries.append(entry)
        if entry_limit is not None and len(entries) >= entry_limit:
            truncated = True
            entry_limit_hit = True
            break

    truncated_reason_ids: List[str] = []
    if entry_limit_hit:
        truncated_reason_ids.append("entry_limit_hit")
    if backlog_limit_hit:
        truncated_reason_ids.append("backlog_limit_hit")
    if depth_limited_dir_count:
        truncated_reason_ids.append("depth_limit_hit")
    if quick_skipped_generic_dir_count:
        truncated_reason_ids.append("quick_generic_branches_skipped")
    if full_skipped_generic_dir_count:
        truncated_reason_ids.append("generic_branches_shallow_scanned")

    return {
        "entries": entries,
        "entry_count": len(entries),
        "truncated": truncated,
        "max_depth": max_depth,
        "max_entries": max_entries,
        "entry_limit_enabled": entry_limit_enabled,
        "entry_limit_hit": entry_limit_hit,
        "backlog_limit_hit": backlog_limit_hit,
        "unbounded_entry_scan": not entry_limit_enabled,
        "truncated_reason_ids": truncated_reason_ids,
        "depth_limited_dir_count": depth_limited_dir_count,
        "depth_limited_dirs": depth_limited_dirs,
        "strategy": "quick" if quick else "full",
        "quick_skipped_generic_dir_count": quick_skipped_generic_dir_count,
        "full_skipped_generic_dir_count": full_skipped_generic_dir_count,
        "sampled_bulk_stage_dir_count": sampled_bulk_stage_dir_count,
        "sampled_bulk_child_dir_count": sampled_bulk_child_dir_count,
        "sampled_bulk_file_count": sampled_bulk_file_count,
        "focus_mode": focus_mode,
        "focus_roots": _dedupe_normalized_paths(active_focus_roots),
        "explicit_focus_roots": explicit_focus_roots,
        "auto_focus_roots": auto_focus_roots,
        "promoted_focus_root_count": promoted_focus_root_count,
        "generic_deep_scan_limit": generic_deep_scan_limit,
    }


def _match_pattern(
    entries: Iterable[str],
    pattern: str,
    *,
    limit: int = 3,
    matcher: Optional[_EntryMatcher] = None,
) -> List[str]:
    if matcher is not None:
        return matcher.match_pattern(pattern, limit=limit)

    normalized = _normalize_pattern(pattern)
    normalized_lower = normalized.lower()
    matches: List[str] = []
    for entry in entries:
        basename = Path(entry).name
        if (
            fnmatch.fnmatch(entry, normalized)
            or fnmatch.fnmatch(basename, normalized)
            or fnmatch.fnmatch(entry.lower(), normalized_lower)
            or fnmatch.fnmatch(basename.lower(), normalized_lower)
        ):
            matches.append(entry)
            if len(matches) >= limit:
                break
    return matches


def _collect_stage_matches(
    entries: Iterable[str],
    patterns: Iterable[str],
    *,
    matcher: Optional[_EntryMatcher] = None,
) -> List[str]:
    if matcher is not None:
        return matcher.collect(patterns, limit=3)

    normalized_entries = [_normalize_pattern(str(entry)).rstrip("/") for entry in entries if str(entry).strip()]
    nonempty_dirs = _nonempty_directory_entries(normalized_entries)
    matches: List[str] = []
    for pattern in patterns:
        for match in _match_pattern(normalized_entries, pattern, limit=max(3, len(normalized_entries))):
            if _is_empty_directory_match(match, nonempty_dirs):
                continue
            if match not in matches:
                matches.append(match)
                if len(matches) >= 3:
                    break
        if len(matches) >= 3:
            break
    return matches
