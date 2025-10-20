#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tme_profile.py — Signature scoring + immune deconvolution + LR_cal from a TPM matrix
- Top-level CLI (like main.py): -i/--input, -o/--output, --threads
- Any extra long flags (e.g., --perm 1000, --qn) are auto-routed to the right method
  just like runall.py.
"""

from __future__ import annotations
import argparse, shlex, subprocess, sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pandas as pd
except Exception:
    pd = None  # optional; only needed for merged deconv table

# --------------------- Routing helpers (mirrors runall) ---------------------

METHODS = ["calculate_sig_score", "cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic", "LR_cal"]

# Per-step known long flags (normalized: "--x-y" -> "x_y")
FLAG_BUCKETS: Dict[str, set] = {
    "calculate_sig_score": {"signature", "method", "mini_gene_count", "adjust_eset", "parallel_size"},
    "cibersort": {"perm", "qn", "absolute", "abs_method", "threads"},
    "IPS": set(),
    "estimate": {"platform"},
    "mcpcounter": {"features"},
    "quantiseq": {"arrays", "signame", "tumor", "scale_mrna", "method"},
    "epic": {"reference"},
    "LR_cal": {"data_type", "id_type", "cancer_type", "verbose"},
}

def _normalize_flag_token(tok: str) -> str:
    """Convert --long-flag to normalized form where '-' becomes '_' in the name."""
    if tok.startswith("--") and len(tok) > 2:
        return f"--{tok[2:].replace('-', '_')}"
    return tok

def _flag_name(tok: str) -> Optional[str]:
    """Normalized flag name (without leading dashes), or None if not a long flag."""
    if tok.startswith("--") and len(tok) > 2:
        return tok[2:].replace("-", "_").lower()
    return None

def _parse_named_blocks(tokens: List[str]) -> Dict[str, List[str]]:
    """
    Optional 'sectioned' style: e.g.
      cibersort --perm 1000 quantiseq --method lsei
    Returns a dict: {'cibersort': ['--perm','1000'], 'quantiseq': ['--method','lsei']}
    """
    blocks: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    expect_val = False
    for tok in tokens:
        bare = tok if not tok.startswith("--") else tok[2:]
        if bare in METHODS and not expect_val:
            cur = bare
            blocks.setdefault(cur, [])
            continue
        if cur is not None:
            nt = _normalize_flag_token(tok)
            blocks[cur].append(nt)
            expect_val = nt.startswith("--")
        else:
            expect_val = False
    return blocks

def _autobucket(tokens: List[str]) -> Dict[str, List[str]]:
    """
    Auto-route long flags by looking up known names in FLAG_BUCKETS.
    Example: '--perm 1000' -> 'cibersort'
    """
    buckets: Dict[str, List[str]] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        name = _flag_name(tok)
        if not name:
            i += 1
            continue
        val = tokens[i + 1] if (i + 1 < len(tokens) and not tokens[i + 1].startswith("--")) else None

        target = None
        # Disambiguate 'method' between sig_score and quantiseq by value
        if name == "method" and val is not None:
            v = str(val).lower()
            if v in {"integration", "pca", "zscore", "ssgsea"}:
                target = "calculate_sig_score"
            elif v in {"lsei", "hampel", "huber", "bisquare"}:
                target = "quantiseq"

        if target is None:
            for mod, flags in FLAG_BUCKETS.items():
                if name in flags:
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

def _append(cmd: List[str], buckets: Dict[str, List[str]], key: str) -> None:
    """Append routed args for a given step to the command list."""
    if key in buckets and buckets[key]:
        cmd += buckets[key]

# --------------------- IO helpers ---------------------

def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def _run(cmd: List[str]) -> int:
    line = " ".join(shlex.quote(x) for x in cmd)
    header = f"[run] {line}"
    print("=" * len(header)); print(header); print("=" * len(header))
    rc = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stderr).returncode
    if rc != 0:
        print(f"[ERROR] {cmd[1]} failed (rc={rc}).")
    else:
        print(f"[ok] {cmd[1]}")
    return rc

def _read_csv_any(p: Path):
    try:
        return pd.read_csv(p)
    except Exception:
        return pd.read_csv(p, engine="python")

def _normalize_id(df):
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

# --------------------- Main ---------------------

def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(
        prog="iobrpy tme_profile",
        description="Compute signature scores, immune infiltration, and LR interactions from a TPM matrix."
    )
    p.add_argument("-i", "--input", required=True, help="TPM matrix (genes x samples). CSV/TSV supported.")
    p.add_argument("-o", "--output", required=True, help="Output directory.")
    p.add_argument("--threads", type=int, default=1, help="Threads for ssGSEA and CIBERSORT (default: 1).")
    ns, unknown = p.parse_known_args(argv)

    in_path = Path(ns.input).resolve()
    outdir = Path(ns.output).resolve()
    d_sig = outdir / "01-signatures"
    d_tme = outdir / "02-tme"
    d_lr  = outdir / "03-LR_cal"
    for d in (d_sig, d_tme, d_lr):
        _ensure_dir(d)

    # Route extra flags to method buckets (support both named sections and flat flags)
    tokens = [_normalize_flag_token(t) for t in unknown]
    blocks_named = _parse_named_blocks(tokens)
    using_named = any(blocks_named.values())
    blocks_auto = {} if using_named else _autobucket(tokens)
    buckets: Dict[str, List[str]] = {}
    for k in set(list(blocks_named.keys()) + list(blocks_auto.keys())):
        buckets[k] = (blocks_named.get(k) or []) + (blocks_auto.get(k) or [])

    # 1) calculate_sig_score -> 01-signatures/calculate_sig_score.csv
    sig_out = d_sig / "calculate_sig_score.csv"
    cs_user = buckets.get("calculate_sig_score") or []
    cmd = [
        "iobrpy", "calculate_sig_score",
        "--input", str(in_path),
        "--output", str(sig_out),
        "--parallel_size", str(ns.threads),  # top-level threads; user-supplied parallel_size (if any) will override when appended
    ]
    # Defaults (only if user didn't set them)
    if not any(a.startswith("--signature") for a in cs_user):
        cmd += ["--signature", "all"]
    if not any(a.startswith("--method") for a in cs_user):
        cmd += ["--method", "integration"]
    if not any(a.startswith("--mini_gene_count") for a in cs_user):
        cmd += ["--mini_gene_count", "2"]
    if not any(a.startswith("--adjust_eset") for a in cs_user):
        cmd += ["--adjust_eset"]
    _append(cmd, buckets, "calculate_sig_score")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)

    # 2) Deconvolution (6 methods) -> 02-tme/{method}_results.csv
    produced: List[Path] = []

    # CIBERSORT
    out_file = d_tme / "cibersort_results.csv"
    cmd = ["iobrpy", "cibersort", "--input", str(in_path), "--output", str(out_file),
           "--threads", str(ns.threads)]
    _append(cmd, buckets, "cibersort")  # e.g., --perm 1000, --qn
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # IPS
    out_file = d_tme / "IPS_results.csv"
    cmd = ["iobrpy", "IPS", "--input", str(in_path), "--output", str(out_file)]
    _append(cmd, buckets, "IPS")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # ESTIMATE (default platform unless user sets)
    out_file = d_tme / "estimate_results.csv"
    est_user = buckets.get("estimate") or []
    cmd = ["iobrpy", "estimate", "--input", str(in_path), "--output", str(out_file)]
    if not any(a.startswith("--platform") for a in est_user):
        cmd += ["--platform", "affymetrix"]
    _append(cmd, buckets, "estimate")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # MCPcounter (default features unless user sets)
    out_file = d_tme / "mcpcounter_results.csv"
    mcp_user = buckets.get("mcpcounter") or []
    cmd = ["iobrpy", "mcpcounter", "--input", str(in_path), "--output", str(out_file)]
    if not any(a.startswith("--features") for a in mcp_user):
        cmd += ["--features", "HUGO_symbols"]
    _append(cmd, buckets, "mcpcounter")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # quanTIseq (enable arrays/tumor/scale_mrna unless user set)
    out_file = d_tme / "quantiseq_results.csv"
    q_user = buckets.get("quantiseq") or []
    cmd = ["iobrpy", "quantiseq", "--input", str(in_path), "--output", str(out_file)]
    if not any(a == "--arrays" for a in q_user):
        cmd.append("--arrays")
    if not any(a == "--tumor" for a in q_user):
        cmd.append("--tumor")
    if not any(a.startswith("--scale_mrna") for a in q_user) and not any(a.startswith("--mRNAscale") for a in q_user):
        cmd.append("--scale_mrna")
    _append(cmd, buckets, "quantiseq")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # EPIC (default reference unless user sets)
    out_file = d_tme / "epic_results.csv"
    epic_user = buckets.get("epic") or []
    cmd = ["iobrpy", "epic", "--input", str(in_path), "--output", str(out_file)]
    if not any(a.startswith("--reference") for a in epic_user):
        cmd += ["--reference", "TRef"]
    _append(cmd, buckets, "epic")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)
    produced.append(out_file)

    # 3) Merge deconvolution -> 02-tme/deconvo_merged.csv
    merged_wide_path = d_tme / "deconvo_merged.csv"
    if pd is None:
        print("[WARN] pandas not available; skip merged deconvolution table.")
    else:
        frames = []
        for f in produced:
            df = _read_csv_any(f)
            df = _normalize_id(df)
            df = df.dropna(axis=1, how="all")
            frames.append(df)
        if frames:
            wide = frames[0]
            for df in frames[1:]:
                overlap = (set(wide.columns) & set(df.columns)) - {"ID"}
                if overlap:
                    df = df.rename(columns={c: f"{c}_dup" for c in overlap})
                wide = pd.merge(wide, df, on="ID", how="outer")
            wide = wide.sort_values("ID").reset_index(drop=True)
            wide.to_csv(merged_wide_path, index=False)
            print(f"[ok] merged deconvolution -> {merged_wide_path}")
        else:
            print("[WARN] No deconvolution outputs to merge.")

    # 4) LR_cal (defaults unless user sets)
    lr_out = d_lr / "lr_cal.csv"
    lr_user = buckets.get("LR_cal") or []
    cmd = ["iobrpy", "LR_cal", "--input", str(in_path), "--output", str(lr_out)]
    if not any(a.startswith("--data_type") for a in lr_user):
        cmd += ["--data_type", "tpm"]
    if not any(a.startswith("--id_type") for a in lr_user):
        cmd += ["--id_type", "symbol"]
    if not any(a.startswith("--cancer_type") for a in lr_user):
        cmd += ["--cancer_type", "pancan"]
    if not any(a.startswith("--verbose") for a in lr_user):
        cmd.append("--verbose")
    _append(cmd, buckets, "LR_cal")
    rc = _run(cmd);  assert rc == 0, sys.exit(rc)

    print("\n[done] tme_profile finished.")

if __name__ == "__main__":
    main()