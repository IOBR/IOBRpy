#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
runall.py — End-to-end orchestrator for IOBRpy

What this revision does (key points):
- No need to type per-step module names (e.g., fastq_qc, batch_salmon). Long flags are auto-routed.
- 'tme_cluster' step removed. Final directory renamed from '07-LR_cal' to '06-LR_cal'.
- Unified concurrency/batch:
  * Top-level --threads propagates to:
    fastq_qc.num_threads, batch_salmon.num_threads, batch_star_count.num_threads,
    merge_salmon.num_processes, cibersort.threads, calculate_sig_score.parallel_size
  * Top-level --batch_size propagates to:
    fastq_qc.batch_size, batch_salmon.batch_size, batch_star_count.batch_size
  Legacy flags (--num_threads / --parallel_size / --num_processes) are absorbed and normalized into --threads.
- NEW: --remove_version is now OPTIONAL in both modes.
  * salmon: if provided, routed to prepare_salmon
  * star:   if provided, routed to count2tpm
  Defaults do NOT force remove_version anymore.
- Output layout:
    01-qc/                 # fastp outputs
    02-salmon/ or 02-star/ # quantification / align counts
    03-tpm/                # unified TPM matrix (tpm_matrix.csv)
    04-signatures/         # calculate_sig_score outputs
    05-tme/                # deconvolution outputs (+ deconvo_merged.csv)
    06-LR_cal/             # LR_cal outputs (renamed from 07-LR_cal)
"""

import argparse
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import sys

try:
    import pandas as pd
except Exception:
    pd = None  # pandas is optional; used only for merging deconvolution results

# Known step names (kept for backward-compatible "sectioned" command style)
METHOD_SECTIONS = {
    "fastq_qc", "batch_salmon", "merge_salmon", "batch_star_count", "merge_star_count",
    "prepare_salmon", "count2tpm",
    "calculate_sig_score", "sig_score",
    "cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic",
    # "tme_cluster",  # removed
    "LR_cal",
}

# --------------------- Utilities ---------------------

def _ensure_dir(p: Path) -> None:
    """Create a directory if it doesn't exist."""
    p.mkdir(parents=True, exist_ok=True)

def _nonempty(p: Path) -> bool:
    """Return True if the path exists and is a non-empty file or directory."""
    return p.exists() and ((p.is_file() and p.stat().st_size > 0) or (p.is_dir() and any(p.iterdir())))

def _run(cmd: List[str], cwd: Optional[Path] = None, dry: bool = False) -> int:
    """Run a subprocess, streaming output to console."""
    line = " ".join(shlex.quote(x) for x in cmd)
    header = f"[run] {line}"
    print("=" * len(header))
    print(header)
    print("=" * len(header))
    if dry:
        print("[dry-run] skipped execution")
        return 0
    rc = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr, cwd=str(cwd) if cwd else None).returncode
    if rc != 0:
        print(f"[ERROR] {cmd[1]} failed (rc={rc}).")
    else:
        print(f"[ok] {cmd[1]}")
    return rc

def _normalize_flag_token(tok: str) -> str:
    """Normalize a long flag token by replacing '-' with '_' in its name."""
    if tok.startswith("--") and len(tok) > 2:
        return f"--{tok[2:].replace('-', '_')}"
    return tok

def _flag_name(tok: str) -> Optional[str]:
    """Return normalized flag name (lowercase, '_' instead of '-') or None if not a long flag."""
    if tok.startswith("--") and len(tok) > 2:
        return tok[2:].replace("-", "_").lower()
    return None

def _parse_passthrough_blocks(tokens: List[str]) -> Dict[str, List[str]]:
    """
    Parse the legacy "sectioned" style:
      fastq_qc --num_threads 8 ... batch_salmon --index ...
    If the user did not provide module name sections, this returns {}.
    """
    blocks: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    expect_value = False

    for tok in tokens:
        name = tok[2:] if tok.startswith("--") else tok
        # Switch section when encountering a bare step name
        if not tok.startswith("--") and name in METHOD_SECTIONS and not expect_value:
            cur = name
            blocks.setdefault(cur, [])
            continue
        # Collect tokens into the current section
        if cur is not None:
            ntok = _normalize_flag_token(tok)
            blocks[cur].append(ntok)
            expect_value = ntok.startswith("--")
        else:
            expect_value = False
    return blocks

def _find_latest(dirpath: Path, globs: List[str]) -> Optional[Path]:
    """Return the most recently modified file in dirpath matching any pattern, or None."""
    cand: List[Path] = []
    for pat in globs:
        cand += list(dirpath.glob(pat))
    if not cand:
        return None
    cand.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return cand[0]

def _append_passthrough(cmd: List[str], buckets: Dict[str, List[str]], key: str, alt_key: Optional[str] = None) -> None:
    """Append auto/sectioned-routed args for a given step into the subcommand list."""
    for k in filter(None, [key, alt_key]):
        if k in buckets and buckets[k]:
            cmd += buckets[k]

# ------------------ Auto-routing of flags ------------------

# Per-step known long flags (normalized names without '--')
FLAG_BUCKETS: Dict[str, set] = {
    # fastq_qc: threads and batch_size will be injected at top-level
    "fastq_qc": {"se", "length_required", "suffix1"},
    # salmon
    "batch_salmon": {"index", "suffix1", "gtf"},
    "merge_salmon": {"project"},  # num_processes injected from --threads
    # star
    "batch_star_count": {"index", "suffix1"},
    "merge_star_count": {"project"},
    # TPM conversion
    "prepare_salmon": {"return_feature", "remove_version"},
    "count2tpm": {"idtype", "org", "source", "id", "length", "gene_symbol", "check_data", "efflength_csv", "remove_version"},
    # signature scores
    "calculate_sig_score": {"signature", "method", "mini_gene_count", "adjust_eset"},
    # deconvolution
    "cibersort": {"perm", "qn", "absolute", "abs_method"},  # cibersort.threads comes from --threads
    "IPS": set(),
    "estimate": {"platform"},
    "mcpcounter": {"features"},
    "quantiseq": {"arrays", "signame", "tumor", "scale_mrna", "method"},
    "epic": {"reference"},
    # ligand-receptor
    "LR_cal": {"data_type", "id_type", "cancer_type", "verbose"},
}

def _consume_top_level_scalars(unknown: List[str]) -> Tuple[List[str], Optional[int], Optional[int]]:
    """
    Pull legacy concurrency flags out of 'unknown' and promote them to top-level:
      threads: threads / num_threads / parallel_size / num_processes
      batch:   batch_size
    Returns (remaining_unknown, threads_val, batch_val).
    """
    tokens = list(unknown)
    threads_val: Optional[int] = None
    batch_val: Optional[int] = None

    def pop_int(names: List[str]) -> Optional[int]:
        nonlocal tokens
        got: Optional[int] = None
        i = 0
        names_set = {n.lower() for n in names}
        while i < len(tokens):
            name = _flag_name(tokens[i]) or ""
            if name in names_set:
                # If next token is a value, capture and remove both
                if i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
                    try:
                        got = int(tokens[i + 1])
                    except Exception:
                        pass
                    del tokens[i:i + 2]
                    continue
                else:
                    # Flag without value -> remove flag token and continue
                    del tokens[i:i + 1]
                    continue
            i += 1
        return got

    # The last occurrence wins
    v_threads = pop_int(["threads", "num_threads", "parallel_size", "num_processes"])
    if v_threads is not None:
        threads_val = v_threads
    v_batch = pop_int(["batch_size"])
    if v_batch is not None:
        batch_val = v_batch

    return tokens, threads_val, batch_val

def _autobucket(tokens: List[str], mode: str) -> Dict[str, List[str]]:
    """
    Route long flags to steps when the user does NOT provide section names.
    Special-cases:
      - --index targets batch_salmon (salmon) or batch_star_count (star)
      - --project targets merge_salmon (salmon) or merge_star_count (star)
      - --remove_version is routed by mode:
          salmon -> prepare_salmon
          star   -> count2tpm
      - --method is disambiguated between calculate_sig_score and quantiseq by value
    """
    buckets: Dict[str, List[str]] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        name = _flag_name(tok)
        if not name:
            i += 1
            continue

        val: Optional[str] = None
        if i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
            val = tokens[i + 1]

        # Mode-aware routing
        if name == "index":
            target = "batch_salmon" if mode == "salmon" else "batch_star_count"
        elif name == "project":
            target = "merge_salmon" if mode == "salmon" else "merge_star_count"
        elif name == "remove_version":
            target = "prepare_salmon" if mode == "salmon" else "count2tpm"
        else:
            # Generic mapping via FLAG_BUCKETS
            target = None
            for mod, flags in FLAG_BUCKETS.items():
                if name in flags:
                    # Disambiguate 'method' between calc_sig_score and quantiseq
                    if name == "method" and mod in ("calculate_sig_score", "quantiseq") and val:
                        v = str(val).lower()
                        if v in {"integration", "pca", "zscore", "ssgsea"}:
                            target = "calculate_sig_score"
                        elif v in {"lsei", "hampel", "huber", "bisquare"}:
                            target = "quantiseq"
                        else:
                            target = "calculate_sig_score"
                        break
                    target = mod
                    break

        if target:
            buckets.setdefault(target, []).append(_normalize_flag_token(tok))
            if val is not None:
                buckets[target].append(val)
        else:
            print(f"[warn] Unrecognized flag (ignored by router): {tok}{(' ' + val) if val else ''}")

        i += 1 if val is None else 2
    return buckets

# --------------------- Main pipeline ---------------------

def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="iobrpy runall", description="End-to-end orchestrator (salmon/star) with auto routing.")
    parser.add_argument("--mode", choices=["salmon", "star"], required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fastq", required=True, help="Path to raw FASTQ directory (used as fastq_qc --path1_fastq)")
    parser.add_argument("--threads", type=int, default=None, help="Unified concurrency for multiple steps")
    parser.add_argument("--batch_size", type=int, default=None, help="Unified batch size for fastq_qc/salmon/star")
    parser.add_argument("--resume", action="store_true", help="Skip steps if outputs already exist")
    parser.add_argument("--dry_run", action="store_true", help="Print commands without executing")
    ns, unknown = parser.parse_known_args(argv)

    outdir = Path(ns.outdir).resolve()

    # Numbered directories (tme_cluster removed; LR_cal renamed to 06)
    d_fastp    = outdir / "01-qc"
    d_tpm      = outdir / "03-tpm"
    d_sigscore = outdir / "04-signatures"
    d_deconv   = outdir / "05-tme"
    d_lrcal    = outdir / "06-LR_cal"
    d_salmon   = outdir / "02-salmon"
    d_star     = outdir / "02-star"
    # Create common directories (shared across modes)
    for d in [d_fastp, d_tpm, d_sigscore, d_deconv, d_lrcal]:
        _ensure_dir(d)

    # Create the mode-specific directory only
    _ensure_dir(d_salmon if ns.mode == "salmon" else d_star)

    # Legacy "sectioned" style (optional)
    blocks_named = _parse_passthrough_blocks([_normalize_flag_token(t) for t in unknown])
    using_named = any(blocks_named.values())

    # Consume legacy concurrency and batch flags
    unknown, legacy_threads, legacy_batch = _consume_top_level_scalars(unknown)

    # Auto-route long flags if no sections were provided
    blocks_auto = {} if using_named else _autobucket([_normalize_flag_token(t) for t in unknown], ns.mode)

    # Merge routes from both sources
    blocks: Dict[str, List[str]] = {}
    for k in set(list(blocks_named.keys()) + list(blocks_auto.keys())):
        blocks[k] = (blocks_named.get(k) or []) + (blocks_auto.get(k) or [])

    # Final unified values (explicit top-level overrides legacy)
    threads = ns.threads if ns.threads is not None else (legacy_threads if legacy_threads is not None else 8)
    batch_size = ns.batch_size if ns.batch_size is not None else (legacy_batch if legacy_batch is not None else 4)

    # 1) fastq_qc -> 01-qc/
    fastp_done_flag = d_fastp / ".fastq_qc.done"
    if ns.resume and fastp_done_flag.exists() and _nonempty(d_fastp):
        print("[resume] fastq_qc skipped (01-qc/ already has outputs).")
    else:
        cmd = ["iobrpy", "fastq_qc",
               "--path1_fastq", str(Path(ns.fastq).resolve()),
               "--path2_fastp", str(d_fastp),
               "--num_threads", str(threads),
               "--batch_size", str(batch_size)]
        _append_passthrough(cmd, blocks, "fastq_qc", "fastq")
        rc = _run(cmd, dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)
        fastp_done_flag.write_text("done\n", encoding="utf-8")

    # 2/3/4) Quantify & merge -> 02-*/ and 03-tpm/
    if ns.mode == "salmon":
        # 2a) batch_salmon -> 02-salmon/
        if not (ns.resume and _nonempty(d_salmon / ".batch_salmon.done") and _nonempty(d_salmon)):
            cmd = ["iobrpy", "batch_salmon",
                   "--path_fq", str(d_fastp),
                   "--path_out", str(d_salmon),
                   "--num_threads", str(threads),
                   "--batch_size", str(batch_size)]
            _append_passthrough(cmd, blocks, "batch_salmon", "salmon")
            rc = _run(cmd, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_salmon / ".batch_salmon.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] batch_salmon skipped.")

        # 3a) merge_salmon (cwd=02-salmon/)
        if not (ns.resume and _nonempty(d_salmon / ".merge_salmon.done")):
            cmd = ["iobrpy", "merge_salmon",
                   "--path_salmon", str(d_salmon),
                   "--project", "runall",
                   "--num_processes", str(threads)]
            _append_passthrough(cmd, blocks, "merge_salmon", "salmon")
            rc = _run(cmd, cwd=d_salmon, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_salmon / ".merge_salmon.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] merge_salmon skipped.")

        merged_salmon_tpm = _find_latest(d_salmon, ["*_salmon_tpm.tsv", "*_salmon_tpm.tsv.gz"])
        if merged_salmon_tpm is None:
            print("[ERROR] Cannot find merged Salmon TPM in '02-salmon/' (pattern '*_salmon_tpm.tsv*').")
            sys.exit(2)

        # 4a) prepare_salmon -> 03-tpm/tpm_matrix.csv
        tpm_matrix = d_tpm / "tpm_matrix.csv"
        if ns.resume and _nonempty(tpm_matrix):
            print("[resume] prepare_salmon skipped.")
        else:
            cmd = ["iobrpy", "prepare_salmon",
                   "--input", str(merged_salmon_tpm),
                   "--output", str(tpm_matrix)]
            # Apply default return_feature only if user didn't provide one
            ps_args = blocks.get("prepare_salmon") or []
            if not any(a.startswith("--return_feature") for a in ps_args):
                cmd += ["--return_feature", "symbol"]
            _append_passthrough(cmd, blocks, "prepare_salmon")
            # Do NOT force --remove_version; it is passed only if the user provided it.
            rc = _run(cmd, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)

    else:
        # 2b) batch_star_count -> 02-star/
        if not (ns.resume and _nonempty(d_star / ".batch_star_count.done") and _nonempty(d_star)):
            cmd = ["iobrpy", "batch_star_count",
                   "--path_fq", str(d_fastp),
                   "--path_out", str(d_star),
                   "--num_threads", str(threads),
                   "--batch_size", str(batch_size)]
            _append_passthrough(cmd, blocks, "batch_star_count", "star")
            rc = _run(cmd, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_star / ".batch_star_count.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] batch_star_count skipped.")

        # 3b) merge_star_count (cwd=02-star/)
        if not (ns.resume and _nonempty(d_star / ".merge_star_count.done")):
            cmd = ["iobrpy", "merge_star_count",
                   "--path", str(d_star),
                   "--project", "runall"]
            _append_passthrough(cmd, blocks, "merge_star_count", "star")
            rc = _run(cmd, cwd=d_star, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_star / ".merge_star_count.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] merge_star_count skipped.")

        merged_star_counts = _find_latest(d_star, ["*_star_ReadsPerGene.tsv", "*_star_ReadsPerGene.tsv.gz", "*.STAR.count*.gz"])
        if merged_star_counts is None:
            print("[ERROR] Cannot find merged STAR ReadsPerGene in '02-star/' (pattern '*_star_ReadsPerGene.tsv*').")
            sys.exit(2)

        # 4b) count2tpm -> 03-tpm/tpm_matrix.csv
        tpm_matrix = d_tpm / "tpm_matrix.csv"
        if ns.resume and _nonempty(tpm_matrix):
            print("[resume] count2tpm skipped.")
        else:
            cmd = ["iobrpy", "count2tpm",
                   "--input", str(merged_star_counts),
                   "--output", str(tpm_matrix),
                   "--idtype", "ensembl",
                   "--org", "hsa",
                   "--source", "local"]
            # Do NOT force --remove_version; it is passed only if the user provided it.
            _append_passthrough(cmd, blocks, "count2tpm")
            rc = _run(cmd, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)

    # 5) calculate_sig_score -> 04-signatures/
    sig_out = d_sigscore / "calculate_sig_score.csv"
    if ns.resume and _nonempty(sig_out):
        print("[resume] calculate_sig_score skipped.")
    else:
        cmd = ["iobrpy", "calculate_sig_score",
               "--input", str(tpm_matrix),
               "--output", str(sig_out),
               "--parallel_size", str(threads)]
        _append_passthrough(cmd, blocks, "calculate_sig_score", "sig_score")
        rc = _run(cmd, dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)

    # 6) Deconvolution (6 methods) -> 05-tme/
    if not _nonempty(tpm_matrix):
        print("[ERROR] TPM matrix missing. Abort before deconvolution.")
        sys.exit(2)

    produced: List[Path] = []
    for m in ["cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic"]:
        out_file = d_deconv / f"{m}_results.csv"
        if ns.resume and _nonempty(out_file):
            print(f"[resume] {m} skipped.")
            produced.append(out_file)
            continue

        if m == "cibersort":
            cmd = ["iobrpy", "cibersort", "--input", str(tpm_matrix), "--output", str(out_file),
                   "--threads", str(threads)]
        elif m == "IPS":
            cmd = ["iobrpy", "IPS", "--input", str(tpm_matrix), "--output", str(out_file)]
        elif m == "estimate":
            cmd = ["iobrpy", "estimate", "--input", str(tpm_matrix), "--platform", "affymetrix", "--output", str(out_file)]
        elif m == "mcpcounter":
            cmd = ["iobrpy", "mcpcounter", "--input", str(tpm_matrix), "--features", "HUGO_symbols", "--output", str(out_file)]
        elif m == "quantiseq":
            cmd = ["iobrpy", "quantiseq", "--input", str(tpm_matrix), "--output", str(out_file)]
        else:
            cmd = ["iobrpy", "epic", "--input", str(tpm_matrix), "--reference", "TRef", "--output", str(out_file)]

        _append_passthrough(cmd, blocks, m)
        rc = _run(cmd, dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)
        produced.append(out_file)

    # 7) Merge deconvolution results -> 05-tme/deconvo_merged.csv
    merged_wide_path = d_deconv / "deconvo_merged.csv"

    if ns.resume and _nonempty(merged_wide_path):
        print("[resume] merge deconvolution skipped.")
    else:
        if pd is None:
            print("[WARN] pandas is not available; skip merged deconvolution table.")
        else:
            def _read_csv_any(p: Path):
                """Read a CSV with a safe fallback parser."""
                try:
                    return pd.read_csv(p)
                except Exception:
                    return pd.read_csv(p, engine="python")

            def _normalize_id(df: pd.DataFrame) -> pd.DataFrame:
                """
                Standardize the sample identifier column to 'ID' and cleanup.
                - Accept 'ID' or 'Unnamed: 0' as the sample column; otherwise use the first column.
                - Drop redundant 'Unnamed:*' columns.
                - Ensure 'ID' is string type and appears as the first column.
                """
                if "ID" in df.columns:
                    pass
                elif "Unnamed: 0" in df.columns:
                    df = df.rename(columns={"Unnamed: 0": "ID"})
                else:
                    df = df.rename(columns={df.columns[0]: "ID"})

                drop_cols = [c for c in df.columns if c.startswith("Unnamed:") and c != "ID"]
                if drop_cols:
                    df = df.drop(columns=drop_cols)

                df["ID"] = df["ID"].astype(str)
                cols = ["ID"] + [c for c in df.columns if c != "ID"]
                return df[cols]

            # Read all produced method outputs (paths are assumed in `produced`)
            frames = []
            for f in produced:
                df = _read_csv_any(f)
                df = _normalize_id(df)
                # Drop columns that are entirely NaN
                df = df.dropna(axis=1, how="all")
                frames.append(df)

            # Outer-join on 'ID' to build a single wide matrix per sample
            wide = frames[0]
            for df in frames[1:]:
                # If different methods accidentally share identical column names,
                # suffix the incoming duplicates to avoid collisions.
                overlap = (set(wide.columns) & set(df.columns)) - {"ID"}
                if overlap:
                    df = df.rename(columns={c: f"{c}_dup" for c in overlap})
                wide = pd.merge(wide, df, on="ID", how="outer")

            # Sort by ID and write the final wide table
            wide = wide.sort_values("ID").reset_index(drop=True)
            wide.to_csv(merged_wide_path, index=False)
            print(f"[ok] merged deconvolution -> {merged_wide_path}")

    # 8) LR_cal -> 06-LR_cal/
    lr_out = d_lrcal / "lr_cal.csv"
    if ns.resume and _nonempty(lr_out):
        print("[resume] LR_cal skipped.")
    else:
        cmd = ["iobrpy", "LR_cal",
               "--input", str(tpm_matrix),
               "--output", str(lr_out),
               "--data_type", "tpm",
               "--id_type", "ensembl",
               "--cancer_type", "pancan"]
        _append_passthrough(cmd, blocks, "LR_cal")
        rc = _run(cmd, dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)

    print("\n[done] runall finished.")

if __name__ == "__main__":
    main()