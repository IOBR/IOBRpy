---
title: "TRUST4 & TCR-BCR"
page-layout: article
toc: true
---

## Overview
- TRUST4 is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from fluid and solid tissues, including tumors. TRUST4 performs de novo assembly on V, J, C genes including the hypervariable complementarity-determining region 3 (CDR3) and reports consensus contigs of BCR/TCR sequences. TRUST4 then realigns the contigs to IMGT reference gene sequences to identify the corresponding gene and CDR3 details.
- IOBRpy wraps the TRUST4 TCR/BCR reconstruction CLI, enabling single-sample or batch runs from BAMs/FASTQs with automatic immune repertoire summarization.
- Based on TRUST4 assemblies, iobrpy derives TCR/BCR clone counts, clone frequencies, and diversity statistics (e.g., Shannon, Gini) to quantify clonal expansion and immune activity. These summaries feed downstream comparisons of immune infiltration, therapy response, or survival analyses.

## Inputs
- **Single BAM**: pass a path to `-b <sample.bam>`.
- **Batch BAMs**: point `-b` to a directory; every `*.bam` will be processed.
- **Single Sample Paired FASTQ**: point `-1/-2` to paired-end read files
- **Batch Paired FASTQs**: use `--fqdir <dir>`; files must end with `_1.fastq.gz` / `_2.fastq.gz` or `_1.fq.gz` / `_2.fq.gz` per sample.

## Basic usage
```bash
# Single BAM
iobrpy trust4 -b /path/to/sample.bam -o /path/to/outdir -t 8

# Batch over a BAM directory
iobrpy trust4 -b /data/bam_dir -o /data/trust4_results -t 8

# Single sample with paired FASTQ
iobrpy trust4 -1 /path/to/sample_1.fastq.gz -2 /path/to/sample_2.fastq.gz -o /path/to/outdir -t 8
iobrpy trust4 -1 /path/to/sample_1.fq.gz -2 /path/to/sample_2.fq.gz -o /path/to/outdir -t 8

# Batch over paired FASTQs
iobrpy trust4 --fqdir /data/fastqs -o /data/trust4_results -t 8
```

**Notes:**
- Use `-t` to control TRUST4 threading.
- After TRUST4 finishes, iobrpy scans the output tree for `*_report.tsv` files and produces `trust4_immdata.csv` and `trust4_immune_indices.csv` for downstream analysis.

### Metrics in `trust4_immune_indices.csv`
When all samples finish, the output root contains `trust4_immune_indices.csv`, which aggregates per-sample repertoire metrics:

- **What it computes**: These metrics summarize clonal expansion, richness, and evenness across samples, making it easier to compare immune repertoire shifts across cohorts or treatment arms.

- `Nreads`: total reads used for clone analysis.
- `Nclones`: number of non-zero clones detected.
- `Length_CDR3`: read-weighted mean CDR3 length.
- `Shannon_Index`: Shannon diversity index measuring clone diversity.
- `Evenness`: normalized Shannon evenness across clones.
- `Top_clone` / `Second_top_clone`: frequencies of the most and second-most abundant clones.
- `Rare_clone` / `Second_Rare_clone`: frequencies of the rarest and second-rarest clones.
- `Gini`: Gini coefficient for clone frequency concentration.
- `Gini_Simpson`: Gini-Simpson index reflecting diversity and unevenness.