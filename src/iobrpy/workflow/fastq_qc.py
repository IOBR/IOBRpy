#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import random
import argparse
from multiprocessing import Pool

# try tqdm for a nice progress bar; fall back gracefully if not installed
try:
    from tqdm.auto import tqdm  # type: ignore
except Exception:  # pragma: no cover
    def tqdm(iterable, **kwargs):
        return iterable

from iobrpy.utils.print_colorful_message import print_colorful_message

def process_sample(file, path1_fastq, path2_fastp, num_threads, suffix1, se, length_required):
    """
    Processes a single FASTQ file using fastp.
    NOTE: Parameters & fastp options are unchanged.
    Returns a dict with sample id, status, and output file paths.
    """
    suffix2 = suffix1.replace("1", "2")
    outputs = []

    if not file.endswith(suffix1):
        return {"sample": file, "status": "ignored", "outputs": []}

    forward_file = os.path.join(path1_fastq, file)
    sample_id = file[:-len(suffix1)]
    output_forward = os.path.join(path2_fastp, file)
    task_file = os.path.join(path2_fastp, f"{sample_id}.task.complete")

    # anticipate paired file path for printing, even if se is True
    reverse_file = forward_file[:-len(suffix1)] + suffix2
    output_reverse = output_forward[:-len(suffix1)] + suffix2

    if os.path.exists(output_forward) and os.path.exists(task_file) and (se or os.path.exists(output_reverse)):
        # Already processed; report outputs for summary
        outputs.append(output_forward)
        if not se:
            outputs.append(output_reverse)
        return {"sample": sample_id, "status": "skipped", "outputs": outputs}

    print_colorful_message(f"Processing: {sample_id}", "green")

    try:
        if se:
            command = [
                "fastp", "-i", forward_file, "-o", output_forward,
                "--thread", str(num_threads),
                "--disable_length_filtering",
                "--disable_quality_filtering",
                "--n_base_limit", "6", "--compression", "6",
                "--html", f"{path2_fastp}/{sample_id}_fastp.html",
                "--json", f"{path2_fastp}/{sample_id}_fastp.json"
            ]
        else:
            command = [
                "fastp", "-i", forward_file, "-o", output_forward,
                "-I", reverse_file, "-O", output_reverse,
                "--thread", str(num_threads),
                "--disable_length_filtering",
                "--disable_quality_filtering",
                "--n_base_limit", "6", "--compression", "6",
                "--html", f"{path2_fastp}/{sample_id}_fastp.html",
                "--json", f"{path2_fastp}/{sample_id}_fastp.json"
            ]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(task_file, 'w') as f:
            f.write("Processing complete for " + sample_id)
        print_colorful_message(f"Task complete for {sample_id}", "cyan")

        outputs.append(output_forward)
        if not se:
            outputs.append(output_reverse)
        return {"sample": sample_id, "status": "processed", "outputs": outputs}

    except subprocess.CalledProcessError as e:
        print(f"Error processing {sample_id}: {e.stderr.decode()}")
        return {"sample": sample_id, "status": "error", "outputs": []}

def _worker(args_tuple):
    """Unpack arguments for Pool.imap_unordered -> process_sample."""
    return process_sample(*args_tuple)

def _run_multiqc(path2_fastp: str):
    """
    Generate a MultiQC report that summarizes all fastp JSON results.
    Output directory: {path2_fastp}/multiqc_report
    Output file:      multiqc_fastp_report.html (+ multiqc_data/)
    """
    # collect fastp JSONs; if none, skip quietly
    json_reports = [f for f in os.listdir(path2_fastp) if f.endswith("_fastp.json")]
    if not json_reports:
        print("No fastp JSON files found; skipping MultiQC.")
        return None

    out_dir = os.path.join(path2_fastp, "multiqc_report")
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        "multiqc",
        "--module", "fastp",          # restrict to fastp module
        "--force",                    # overwrite existing outputs
        "--outdir", out_dir,
        "--filename", "multiqc_fastp_report",
        path2_fastp                   # scan the fastp output directory
    ]
    try:
        completed = subprocess.run(
            cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        report_html = os.path.join(out_dir, "multiqc_fastp_report.html")
        data_dir = os.path.join(out_dir, "multiqc_data")
        print("\nMultiQC report saved to:")
        print(report_html)
        return report_html
    except FileNotFoundError:
        print("MultiQC not found. Please install it, e.g.: conda install -c bioconda multiqc")
    except subprocess.CalledProcessError as e:
        print(f"MultiQC failed: {e.stderr.decode()}")
    return None

def step1_fastq_qc(path1_fastq, path2_fastp, num_threads=8, suffix1="_1.fastq.gz", batch_size=5, se=False, length_required=50):
    """
    Preprocess FASTQ files using fastp in parallel.
    Args:
        path1_fastq (str): Path to raw FASTQ files.
        path2_fastp (str): Path to preprocessed FASTQ files.
        num_threads (int): Number of threads for fastp.
        suffix1 (str): Suffix for the forward reads files.
        batch_size (int): Number of samples to process simultaneously.
        se (bool): Flag to indicate if the sequencing data is single-end.
        length_required (int): Minimum length of reads to keep after processing.
    """
    # Keep your original banner at the beginning (unchanged)

    print("### FASTQ files quality control using fastp ###")

    os.makedirs(path2_fastp, exist_ok=True)
    fastq_files = [f for f in os.listdir(path1_fastq) if f.endswith(suffix1)]
    random.shuffle(fastq_files)

    tasks = [(file, path1_fastq, path2_fastp, num_threads, suffix1, se, length_required)
             for file in fastq_files]

    results = []
    if tasks:
        with Pool(processes=batch_size) as pool:
            for res in tqdm(pool.imap_unordered(_worker, tasks),
                            total=len(tasks),
                            desc="Running fastp",
                            unit="sample"):
                results.append(res)

    # ------- FINAL: print output file paths first -------
    all_outputs = []
    for r in results:
        all_outputs.extend(r.get("outputs", []))

    # Unique & sorted for readability
    unique_outputs = sorted(set(all_outputs))

    print("\nSaved/Existing output files:")
    if unique_outputs:
        for p in unique_outputs:
            print(p)
    else:
        print("(No outputs produced or files already present.)")

    _run_multiqc(path2_fastp)

    # ------- Then print the IOBRpy banner specified by user -------
    print("   ")
    print_colorful_message("#########################################################", "blue")
    print_colorful_message(" IOBRpy: Immuno-Oncology Biological Research using Python ", "cyan")
    print_colorful_message(" If you encounter any issues, please report them at ", "cyan")
    print_colorful_message(" https://github.com/IOBR/IOBRpy/issues ", "cyan")
    print_colorful_message("#########################################################", "blue")
    print(" Author: Haonan Huang, Dongqiang Zeng")
    print(" Email: interlaken@smu.edu.cn ")
    print_colorful_message("#########################################################", "blue")
    print("   ")

def main():
    parser = argparse.ArgumentParser(description="Preprocess FASTQ files using fastp")
    parser.add_argument("--path1_fastq", type=str, required=True, help="Path to raw FASTQ files")
    parser.add_argument("--path2_fastp", type=str, required=True, help="Path to preprocessed FASTQ files")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads for fastp")
    parser.add_argument("--suffix1", type=str, default="_1.fastq.gz", help="Suffix of the forward reads file")
    parser.add_argument("--batch_size", type=int, default=5, help="Number of samples to process simultaneously")
    parser.add_argument("--se", action='store_true', help="Indicate if the sequencing data is single-end")
    parser.add_argument("--length_required", type=int, default=50, help="Minimum length of reads to keep after processing")
    args = parser.parse_args()

    step1_fastq_qc(args.path1_fastq, args.path2_fastp, args.num_threads, args.suffix1, args.batch_size, args.se, args.length_required)

if __name__ == "__main__":
    main()
