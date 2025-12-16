---
title: "TRUST4 tutorial"
page-layout: article
toc: true
---

## Overview
`iobrpy trust4` wraps the TRUST4 TCR/BCR reconstruction CLI so you can run single samples or batch runs from BAMs/FASTQs and automatically summarize immune repertoires.

## Inputs
- **Single BAM**: pass a path to `-b <sample.bam>`.
- **Batch BAMs**: point `-b` to a directory; every `*.bam` will be processed.
- **Paired FASTQs**: use `--fqdir <dir>`; files must end with `_1.fastq.gz` / `_2.fastq.gz` per sample.

## Basic usage
```bash
# Single BAM
iobrpy trust4 -b /path/to/sample.bam -o /path/to/outdir

# Batch over a BAM directory
iobrpy trust4 -b /data/bam_dir -o /data/trust4_results -t 8

# Batch over paired FASTQs
# Results for each sample are grouped under the output root
iobrpy trust4 --fqdir /data/fastqs -o /data/trust4_results -t 8
```

Key behavior:
- In batch modes, each sample is placed in its own subfolder under `-o` with `TRUST_<sample>` prefixes.
- After TRUST4 finishes, iobrpy scans the output tree for `*_report.tsv` files and produces `trust4_immdata.csv` and `trust4_immune_indices.csv` for downstream analysis.

## Tips
- Set `TRUST4_BIN` if the TRUST4 executable is not on `PATH`.
- Use `-t/--threads` to control TRUST4 threading.
- Combine with `iobrpy runall --mode star` to reuse STAR BAMs directly with `iobrpy trust4 -b <star_output_dir>`.
