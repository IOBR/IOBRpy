
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch-run Salmon quantification for paired-end FASTQs.

Improvements in this fixed version:
- Robust preflight: print Salmon runtime version and basic index meta info.
- Clear guidance when Salmon fails due to index/runner version mismatch.
- Safe R1->R2 inference from --suffix1 (e.g., "_1.fastq.gz" -> "_2.fastq.gz", "R1.fq.gz" -> "R2.fq.gz").
- Resume-friendly: skip finished samples; no UnboundLocalError on skip path.
- Uses subprocess.run(list) (no shell=True) to avoid quoting/escaping issues.
- Progress bar with tqdm over samples; per-sample logging.
"""
import os
import re
import sys
import glob
import json
import argparse
import subprocess
from multiprocessing import Pool
from functools import partial

from tqdm.auto import tqdm


def _salmon_version_tuple():
    """Return Salmon version as (major, minor, patch) or None if unavailable."""
    try:
        out = subprocess.check_output(["salmon", "--version"], text=True).strip()
        # expected: "salmon 1.10.3"
        m = re.search(r'(\d+)\.(\d+)\.(\d+)', out)
        if m:
            return tuple(map(int, m.groups()))
    except Exception:
        pass
    return None


def _read_index_meta(index_dir):
    """Best-effort read of index meta json; return dict or None."""
    for name in ("versionInfo.json", "meta_info.json", "info.json"):
        p = os.path.join(index_dir, name)
        if os.path.exists(p):
            try:
                with open(p, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return {"__raw__": "unreadable"}
    return None


def _infer_suffix2_from_suffix1(s1: str) -> str:
    """Infer read2 suffix from read1 suffix. Handles common patterns safely."""
    # Prefer specific tokens first
    if "_1" in s1:
        return s1.replace("_1", "_2")
    if "R1" in s1:
        return s1.replace("R1", "R2")
    # Fallback near file stem: ...1.fastq(.gz)? / 1.fq(.gz)?
    s2 = re.sub(r'1(\.f(?:ast)?q(?:\.gz)?)$', r'2\1', s1)
    if s2 != s1:
        return s2
    # Last resort: try delimiter-number pattern
    s2 = re.sub(r'([._-])1(\.f(?:ast)?q(?:\.gz)?)$', r'\g<1>2\2', s1)
    if s2 != s1:
        return s2
    raise ValueError(f"Unable to infer suffix2 from suffix1='{s1}'. Please rename FASTQs or add a clear pattern like '_1'/'_2' or 'R1'/'R2'.")


def _pair_from_r1(path_r1: str, suffix1: str, suffix2: str):
    """Given an R1 path and suffixes, return (sample_id, r1, r2)."""
    if not path_r1.endswith(suffix1):
        raise ValueError(f"File does not end with suffix1: {path_r1} (suffix1={suffix1})")
    base = os.path.basename(path_r1)[:-len(suffix1)]
    dirn = os.path.dirname(path_r1)
    r2 = os.path.join(dirn, base + suffix2)
    return base, path_r1, r2


def _exists_nonempty(p: str) -> bool:
    return os.path.exists(p) and os.path.getsize(p) > 0


def _ensure_dir(d: str):
    os.makedirs(d, exist_ok=True)


def _run_salmon_one(sample_id: str, r1: str, r2: str, index: str, out_root: str, threads: int, gtf: str = None):
    out_dir = os.path.join(out_root, sample_id)
    _ensure_dir(out_dir)

    # Resume logic
    quant_sf = os.path.join(out_dir, "quant.sf")
    done_flag = os.path.join(out_dir, "task.complete")
    if _exists_nonempty(quant_sf) and os.path.exists(done_flag):
        print(f"[skip] {sample_id}: quant.sf exists and task.complete present.")
        print(f"Saved to: {quant_sf}", flush=True)
        return sample_id, True, None

    cmd = [
        "salmon", "quant",
        "-i", index,
        "-l", "ISF",
        "--gcBias",
        "-1", r1,
        "-2", r2,
        "-p", str(threads),
        "-o", out_dir,
        "--validateMappings",
    ]
    if gtf:
        cmd += ["-g", gtf]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        # Friendly guidance for common cause
        msg = (
            "[error] salmon quant failed for sample {sid} (exit {code}).\n"
            "        Command: {cmd}\n"
            "        Common cause: index built by newer Salmon than your runtime,\n"
            "        leading to rapidjson assertion failures or 'invoked improperly'.\n"
            "        Fix: upgrade Salmon to the version used to build the index (or newer),\n"
            "        or rebuild the index using your current Salmon version."
        ).format(sid=sample_id, code=e.returncode, cmd=" ".join(cmd))
        return sample_id, False, msg

    # Mark done
    try:
        with open(done_flag, "w") as f:
            f.write("ok\n")
    except Exception:
        pass
    print(f"Saved to: {quant_sf}", flush=True)
    return sample_id, True, None


def process_sample(args_tuple):
    """Pool wrapper that returns (sample_id, success, error_message_or_None)."""
    return _run_salmon_one(*args_tuple)


def main():
    parser = argparse.ArgumentParser(description="Batch Salmon quant for paired-end FASTQs")
    parser.add_argument("--index", required=True, help="Path to Salmon index directory")
    parser.add_argument("--path_fq", required=True, help="Directory containing FASTQ files")
    parser.add_argument("--path_out", required=True, help="Output directory for per-sample results")
    parser.add_argument("--suffix1", default="_1.fastq.gz", help="R1 file suffix (default: _1.fastq.gz)")
    parser.add_argument("--batch_size", type=int, default=1, help="Number of concurrent samples (processes)")
    parser.add_argument("--num_threads", type=int, default=8, help="Threads per Salmon process")
    parser.add_argument("--gtf", default=None, help="Optional GTF file for -g to produce gene-level quant")
    args = parser.parse_args()

    # Preflight info
    sv = _salmon_version_tuple()
    if sv:
        print(f"[preflight] salmon version: {'.'.join(map(str, sv))}")
    else:
        print("[preflight] salmon version: <unknown> (salmon not found on PATH?)")

    meta = _read_index_meta(args.index)
    if isinstance(meta, dict):
        # Print a short summary of meta keys; don't over-parse
        keys = list(meta.keys())[:6]
        print(f"[preflight] index meta keys: {keys if keys else '<none>'}")
    else:
        print("[preflight] index meta: <not found>")

    # Resolve suffix2 from suffix1
    try:
        suffix2 = _infer_suffix2_from_suffix1(args.suffix1)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

    # Discover R1s and pair
    pattern = os.path.join(args.path_fq, f"*{args.suffix1}")
    r1_files = sorted(glob.glob(pattern))
    if not r1_files:
        print(f"[error] No FASTQ files matched: {pattern}", file=sys.stderr)
        sys.exit(2)

    pairs = []
    for r1 in r1_files:
        sample_id, r1p, r2p = _pair_from_r1(r1, args.suffix1, suffix2)
        if not os.path.exists(r2p):
            print(f"[warn] Missing R2 for {sample_id}: {r2p}. Skipping this sample.", file=sys.stderr)
            continue
        pairs.append((sample_id, r1p, r2p, args.index, args.path_out, args.num_threads, args.gtf))

    if not pairs:
        print("[error] No valid R1/R2 pairs to process.", file=sys.stderr)
        sys.exit(2)

    # Ensure output root exists
    _ensure_dir(args.path_out)

    # Run
    print(f"[plan] {len(pairs)} samples; batch_size={args.batch_size}; threads_per_job={args.num_threads}")
    failures = []
    with Pool(processes=max(1, int(args.batch_size))) as pool:
        for sid, ok, err in tqdm(pool.imap_unordered(process_sample, pairs),
                                 total=len(pairs), desc="Samples", unit="sample"):
            if not ok:
                failures.append((sid, err))
                # Fail fast? comment next 2 lines to continue others
                # for continuing, don't sys.exit here; we'll report after loop
                # break

    if failures:
        print("\n[summary] Some samples failed:", file=sys.stderr)
        for sid, err in failures:
            print(f"--- {sid} ---", file=sys.stderr)
            print(err, file=sys.stderr)
        # Non-zero exit to signal error overall
        sys.exit(1)
    else:
        print("\n[summary] All samples finished successfully.")

if __name__ == "__main__":
    main()
