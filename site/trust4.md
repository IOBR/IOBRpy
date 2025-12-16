---
title: "TRUST4 tutorial"
page-layout: article
toc: true
---

## Overview
`iobrpy trust4` wraps the TRUST4 TCR/BCR reconstruction CLI so you can run single samples or batch runs from BAMs/FASTQs and automatically summarize immune repertoires.

- **What it computes / why it matters**: Based on TRUST4 assemblies, iobrpy derives TCR/BCR clone counts, clone frequencies, and diversity statistics (e.g., Shannon, Gini) to quantify clonal expansion and immune activity. These summaries feed downstream comparisons of immune infiltration, therapy response, or survival analyses.

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

### Metrics in `trust4_immune_indices.csv`
When all samples finish, the output root contains `trust4_immune_indices.csv`, which aggregates per-sample repertoire metrics:

- `Nreads`: total reads used for clone analysis.
- `Nclones`: number of non-zero clones detected.
- `Length_CDR3`: read-weighted mean CDR3 length.
- `Shannon_Index`: Shannon diversity index measuring clone diversity.
- `Evenness`: normalized Shannon evenness across clones.
- `Top_clone` / `Second_top_clone`: frequencies of the most and second-most abundant clones.
- `Rare_clone` / `Second_Rare_clone`: frequencies of the rarest and second-rarest clones.
- `Gini`: Gini coefficient for clone frequency concentration.
- `Gini_Simpson`: Gini-Simpson index reflecting diversity and unevenness.

These metrics summarize clonal expansion, richness, and evenness across samples, making it easier to compare immune repertoire shifts across cohorts or treatment arms.

## Tips
- Set `TRUST4_BIN` if the TRUST4 executable is not on `PATH`.
- Use `-t/--threads` to control TRUST4 threading.
- Combine with `iobrpy runall --mode star` to reuse STAR BAMs directly with `iobrpy trust4 -b <star_output_dir>`.
