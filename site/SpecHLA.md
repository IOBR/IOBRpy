---
title: "SpecHLA & HLA typing"
page-layout: article
toc: true
---

## Single-sample SpecHLA (`iobrpy spechla`)
Run SpecHLA in RNA-seq mode when you already have paired FASTQs.

- **What it computes / why it matters**: SpecHLA performs high-resolution HLA genotyping (HLA-A/B/C, DR/DQ/DP, and other loci), producing allele calls that support immune-escape assessment, transplant compatibility review, immunotherapy prognosis, and neoantigen prediction.

```bash
iobrpy spechla \
  -n SAMPLE1 \
  -1 /data/SAMPLE1_1.fq.gz \
  -2 /data/SAMPLE1_2.fq.gz \
  -o /data/spechla_results \
  -j 16
```
- `-n/--name`: sample identifier.
- `-1/--read1` and `-2/--read2`: FASTQ.gz mates.
- `-o/--outdir`: result directory (created if missing).
- `-j/--threads`: threads for SpecHLA (default 8).

## Extract HLA reads from BAM/CRAM (`iobrpy extract_hla_read`)
Use SpecHLA's ExtractHLAread helper to pull HLA reads and convert them to FASTQ.

- **What it computes / why it matters**: ExtractHLAread isolates HLA-mapped reads from aligned BAM/CRAM files and converts them to FASTQ to provide cleaner inputs for HLA typing, improving downstream SpecHLA accuracy and enabling QC-friendly batch pipelines.

```bash
iobrpy extract_hla_read \
  -s SAMPLE1 \
  -b /data/bams/SAMPLE1.bam \
  -r hg38 \
  -o /data/ExtractHLAread/SAMPLE1
```
- `-s/--sample`: sample name written into outputs.
- `-b/--bam`: sorted + indexed BAM/CRAM.
- `-r/--ref`: `hg19` or `hg38` reference for ExtractHLAread.
- `-o/--outdir`: where the extracted FASTQs and logs go.
- Add `--no-auto-install` to skip automatic dependency checks.

## Batch HLA typing (`iobrpy hla_typing`)
End-to-end wrapper that runs **ExtractHLAread** for every BAM in a directory, launches **SpecHLA** on the extracted FASTQs, and merges per-sample calls.

- **What it computes / why it matters**: The batch workflow outputs per-sample HLA genotypes and a merged `hla.results.txt`, simplifying cohort-wide comparisons, clinical association tests, and downstream immune analyses such as neoantigen or donor matching workflows.

```bash
iobrpy hla_typing \
  -b /data/bams \
  -r hg38 \
  -o /data/hla_typing_results \
  -j 16
```
Workflow details:
- Input directory `-b/--bam-dir` is scanned for `*.bam`; each sample gets `ExtractHLAread/<sample>` and `SpecHLA/<sample>` under the output root.
- Resume-friendly markers (`*.ExtractHLAread.done`, `*.SpecHLA.done`) prevent reruns on completed samples.
- A merged `hla.results.txt` is written alongside the per-sample folders for quick review.
