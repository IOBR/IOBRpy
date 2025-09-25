#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
runall.py — End-to-end orchestrator for IOBRpy (salmon/star only, no deside).

Pipeline:
  mode=salmon:
    fastq_qc -> batch_salmon -> merge_salmon -> prepare_salmon -> calculate_sig_score
      -> deconv(6) -> merge (CSV) -> tme_cluster -> LR_cal

  mode=star:
    fastq_qc -> batch_star_count -> merge_star_count -> count2tpm -> calculate_sig_score
      -> deconv(6) -> merge (CSV) -> tme_cluster -> LR_cal

Output layout under --outdir (CSV everywhere):
  01-qc/                    # fastq_qc outputs
  02-salmon/ or 02-star/    # batch_* + merge_* outputs (depending on --mode)
  03-tpm/                   # prepare_salmon OR count2tpm -> tpm_matrix.csv
  04-signatures/            # calculate_sig_score.csv
  05-tme/                   # <method>_results.csv + deconvo_merged.csv
  06-tme_cluster/           # tme_cluster.csv
  07-LR_cal/                # lr_cal.csv

tme_cluster behavior in runall:
- Default: use 05-tme/deconvo_merged.csv.
- If 'tme_cluster --pattern <regex>' is provided, use the single file in 05-tme that matches <regex>.
- If you also pass 'tme_cluster --pattern_cols <regex>', it is forwarded to tme_cluster's own '--pattern' to filter columns.
- Example:
    tme_cluster --pattern cibersort --features 1:22 --id ID --max_nc 3 --print_result --scale
    # uses 05-tme/cibersort_results.csv

Notes:
- Long flags in passthrough are normalized: '--id-type' -> '--id_type' (values unchanged).
- Top-level --fastq feeds fastq_qc --path1_fastq.
- No log files; commands stream to console.
"""

from __future__ import annotations
import argparse
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
import sys
import time
import re

try:
    import pandas as pd
except Exception:
    pd = None

# Supported block names that can follow their own arguments
METHOD_SECTIONS = {
    # step-specific
    "fastq_qc", "batch_salmon", "merge_salmon", "batch_star_count", "merge_star_count",
    # generic fallbacks
    "fastq", "salmon", "star", "prepare_salmon", "count2tpm",
    # scoring
    "calculate_sig_score", "sig_score",
    # deconvolution (6 methods; no deside)
    "cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic",
    # finals
    "tme_cluster", "LR_cal",
}

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _run(cmd: List[str], _unused_logf: Path, cwd: Optional[Path] = None, dry: bool = False) -> int:
    """Execute a command and stream output to console (no log files)."""
    line = " ".join(shlex.quote(x) for x in cmd)
    header = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] CMD: {line} (cwd={cwd or Path.cwd()})"
    if dry:
        print(f"[DRY-RUN] {header}")
        return 0
    print(header)
    rc = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr, cwd=str(cwd) if cwd else None).returncode
    if rc != 0:
        print(f"[ERROR] {cmd[1]} failed (rc={rc}).")
    else:
        print(f"[ok] {cmd[1]}")
    return rc

def _normalize_flag_token(tok: str) -> str:
    """Normalize long flags by replacing '-' with '_' (values unchanged)."""
    if tok.startswith("--") and len(tok) > 2:
        return f"--{tok[2:].replace('-', '_')}"
    return tok

def _parse_passthrough_blocks(tokens: List[str]) -> Dict[str, List[str]]:
    """
    Parse unknown argv tokens into blocks keyed by section name.

    Fix: if the previous token was a flag ('--xxx'), the very next non-flag token
    is treated as its VALUE even if it equals a method name (to avoid mis-parsing
    like: 'tme_cluster --pattern cibersort' where 'cibersort' is also a block name).
    """
    blocks: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    expect_value = False

    for tok in tokens:
        if expect_value and not tok.startswith("--"):
            if cur is not None:
                blocks[cur].append(_normalize_flag_token(tok))
            expect_value = False
            continue

        name = tok[2:] if (tok.startswith("--") and len(tok) > 2) else tok

        # Start a new block only when not expecting a value and token is a bare method name
        if not expect_value and not tok.startswith("--") and name in METHOD_SECTIONS:
            cur = name
            blocks.setdefault(cur, [])
            continue

        if cur is not None:
            ntok = _normalize_flag_token(tok)
            blocks[cur].append(ntok)
            expect_value = ntok.startswith("--")
        # else: ignore stray tokens

    return blocks

def _append_passthrough(cmd: List[str], blocks: Dict[str, List[str]], *sections: str) -> None:
    """Append passthrough args for the first existing section among 'sections'."""
    for s in sections:
        extra = blocks.get(s)
        if extra:
            cmd += extra
            return

def _pop_opt(blocks: Dict[str, List[str]], section: str, opt_name: str, has_value: bool = True):
    """
    Pop an option from blocks[section] and return its value (or True if flag-like).
    'opt_name' must be the normalized long-form, e.g., '--pattern' or '--pattern_cols'.
    """
    lst = blocks.get(section)
    if not lst:
        return None
    i = 0
    while i < len(lst):
        if lst[i] == opt_name:
            if has_value:
                val = lst[i + 1] if i + 1 < len(lst) else None
                del lst[i:i + 2]
                return val
            else:
                del lst[i]
                return True
        i += 1
    return None

def _find_latest(dirpath: Path, globs: List[str]) -> Optional[Path]:
    """Find the most recently modified file matching any of the glob patterns."""
    cand: List[Path] = []
    for pat in globs:
        cand += list(dirpath.glob(pat))
    if not cand:
        return None
    cand.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return cand[0]

def _nonempty(p: Path) -> bool:
    """True if path exists and is a non-empty file or a non-empty directory."""
    return p.exists() and (p.is_file() and p.stat().st_size > 0 or (p.is_dir() and any(p.iterdir())))

def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="iobrpy runall", description="End-to-end orchestrator (salmon/star only).")
    parser.add_argument("--mode", choices=["salmon", "star"], required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fastq", required=True, help="Path to raw FASTQ directory for fastq_qc --path1_fastq")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry_run", action="store_true")
    ns, unknown = parser.parse_known_args(argv)

    outdir = Path(ns.outdir).resolve()

    # Numbered directories (RENAMED per your request)
    d_fastp    = outdir / "01-qc"
    d_tpm      = outdir / "03-tpm"
    d_sigscore = outdir / "04-signatures"
    d_deconv   = outdir / "05-tme"
    d_tmeclust = outdir / "06-tme_cluster"
    d_lrcal    = outdir / "07-LR_cal"
    for d in [d_fastp, d_tpm, d_sigscore, d_deconv, d_tmeclust, d_lrcal]:
        _ensure_dir(d)

    # Mode-specific directory (02-salmon OR 02-star)
    d_salmon = d_star = None
    if ns.mode == "salmon":
        d_salmon = outdir / "02-salmon"
        _ensure_dir(d_salmon)
    else:  # star
        d_star = outdir / "02-star"
        _ensure_dir(d_star)

    blocks = _parse_passthrough_blocks(unknown)

    # 1) fastq_qc -> 01-qc/
    fastp_done_flag = d_fastp / ".fastq_qc.done"
    if ns.resume and fastp_done_flag.exists() and _nonempty(d_fastp):
        print("[resume] fastq_qc skipped (01-qc/ already has outputs).")
    else:
        cmd = ["iobrpy", "fastq_qc",
               "--path1_fastq", str(Path(ns.fastq).resolve()),
               "--path2_fastp", str(d_fastp)]
        _append_passthrough(cmd, blocks, "fastq_qc", "fastq")
        rc = _run(cmd, Path(), dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)
        fastp_done_flag.write_text("done\n", encoding="utf-8")

    # 2/3/4) quant + merge + TPM -> 02-*/ + 03-tpm/
    if ns.mode == "salmon":
        # 2a) batch_salmon -> 02-salmon/
        if not (ns.resume and _nonempty((d_salmon / ".batch_salmon.done")) and _nonempty(d_salmon)):
            cmd = ["iobrpy", "batch_salmon",
                   "--path_fq", str(d_fastp),
                   "--path_out", str(d_salmon)]
            _append_passthrough(cmd, blocks, "batch_salmon", "salmon")
            rc = _run(cmd, Path(), dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_salmon / ".batch_salmon.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] batch_salmon skipped.")

        # 3a) merge_salmon (cwd=02-salmon/)
        if not (ns.resume and _nonempty(d_salmon / ".merge_salmon.done")):
            cmd = ["iobrpy", "merge_salmon",
                   "--path_salmon", str(d_salmon),
                   "--project", "runall"]
            _append_passthrough(cmd, blocks, "merge_salmon", "salmon")
            rc = _run(cmd, Path(), cwd=d_salmon, dry=ns.dry_run)
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
                   "--output", str(tpm_matrix),
                   "--return_feature", "symbol",
                   "--remove_version"]
            _append_passthrough(cmd, blocks, "prepare_salmon")
            rc = _run(cmd, Path(), dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)

    else:  # star
        # 2b) batch_star_count -> 02-star/
        if not (ns.resume and _nonempty((d_star / ".batch_star_count.done")) and _nonempty(d_star)):
            cmd = ["iobrpy", "batch_star_count",
                   "--path_fq", str(d_fastp),
                   "--path_out", str(d_star)]
            _append_passthrough(cmd, blocks, "batch_star_count", "star")
            rc = _run(cmd, Path(), dry=ns.dry_run)
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
            rc = _run(cmd, Path(), cwd=d_star, dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)
            (d_star / ".merge_star_count.done").write_text("done\n", encoding="utf-8")
        else:
            print("[resume] merge_star_count skipped.")

        merged_star_counts = _find_latest(d_star, ["*_star_ReadsPerGene.tsv", "*_star_ReadsPerGene.tsv.gz"])
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
                   "--idType", "Ensembl",
                   "--org", "hsa",
                   "--source", "local",
                   "--remove_version"]
            _append_passthrough(cmd, blocks, "count2tpm")
            rc = _run(cmd, Path(), dry=ns.dry_run)
            if rc != 0:
                sys.exit(rc)

    # 5) calculate_sig_score -> 04-signatures/calculate_sig_score.csv
    sig_out = d_sigscore / "calculate_sig_score.csv"
    if ns.resume and _nonempty(sig_out):
        print("[resume] calculate_sig_score skipped.")
    else:
        cmd = ["iobrpy", "calculate_sig_score",
               "--input", str(tpm_matrix),
               "--output", str(sig_out)]
        _append_passthrough(cmd, blocks, "calculate_sig_score", "sig_score")
        rc = _run(cmd, Path(), dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)

    # 6) deconvolution(6) -> 05-tme/<method>_results.csv
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
            cmd = ["iobrpy", "cibersort", "--input", str(tpm_matrix), "--output", str(out_file)]
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
        rc = _run(cmd, Path(), dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)
        produced.append(out_file)

    # 7) Merge deconvolution results -> 05-tme/deconvo_merged.csv (always produced)
    merged_path = d_deconv / "deconvo_merged.csv"
    if ns.resume and _nonempty(merged_path):
        print("[resume] merge deconv skipped.")
    else:
        if pd is None:
            print("[WARN] pandas not available; skip merged deconvolution table (needed by tme_cluster default).")
        else:
            df_all = None
            for f in produced:
                try:
                    df = pd.read_csv(f)  # CSV
                except Exception:
                    df = pd.read_csv(f, engine="python")
                if df.shape[1] < 2:
                    print(f"[WARN] skip malformed file: {f}")
                    continue
                df = df.rename(columns={df.columns[0]: "ID"})
                df_all = df if df_all is None else pd.merge(df_all, df, on="ID", how="outer")
            if df_all is None:
                print("[ERROR] No valid deconvolution results to merge.")
                sys.exit(2)
            df_all.to_csv(merged_path, index=False)
            print(f"[ok] merged deconvolution -> {merged_path}")

    # 8) tme_cluster -> 06-tme_cluster/tme_cluster.csv
    #    Default: use merged CSV. If 'tme_cluster --pattern <regex>' is provided,
    #    use the single matching *_results.csv in 05-tme instead.
    tme_out = d_tmeclust / "tme_cluster.csv"
    if ns.resume and _nonempty(tme_out):
        print("[resume] tme_cluster skipped.")
    else:
        # Extract file regex (runall-level) and optional column regex
        file_pattern = _pop_opt(blocks, 'tme_cluster', '--pattern', has_value=True)
        pattern_cols = _pop_opt(blocks, 'tme_cluster', '--pattern_cols', has_value=True)

        input_for_cluster = merged_path
        if file_pattern:
            all_files = list(d_deconv.glob("*_results.csv"))
            rex = re.compile(file_pattern, re.IGNORECASE)
            matches = [f for f in all_files if rex.search(f.name)]
            if not matches:
                print(f"[ERROR] No deconvolution result matches pattern '{file_pattern}' in {d_deconv}")
                sys.exit(2)
            if len(matches) > 1:
                names = ", ".join(f.name for f in matches)
                print(f"[ERROR] Pattern '{file_pattern}' matched multiple files: {names}. Please be more specific.")
                sys.exit(2)
            input_for_cluster = matches[0]

        cmd = ["iobrpy", "tme_cluster",
               "--input", str(input_for_cluster),
               "--output", str(tme_out)]
        if pattern_cols:
            cmd += ["--pattern", pattern_cols]  # forward as column regex to tme_cluster

        # Pass all other tme_cluster args (e.g., --features, --id, --scale, --max_nc, etc.)
        extra = blocks.get("tme_cluster", [])
        if extra:
            cmd += extra

        rc = _run(cmd, Path(), dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)

    # 9) LR_cal -> 07-LR_cal/lr_cal.csv
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
        rc = _run(cmd, Path(), dry=ns.dry_run)
        if rc != 0:
            sys.exit(rc)

    print("\n[done] runall finished.")