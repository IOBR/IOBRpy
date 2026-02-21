---
title: "SpecHLA & HLA typing"
page-layout: article
toc: true
---

## Automatic dependency installation

Both `iobrpy extract_hla_read` and `iobrpy spechla` support **lazy (on-demand) dependency installation**.  
When you run the command, iobrpy will first check whether required dependencies are available; if something is missing, it will **automatically install** via `conda install` and then continue the workflow.

## Single-sample HLA typing

### Extract HLA reads from BAM/CRAM

Use `iobrpy extract_hla_read` to extract HLA reads and convert them to FASTQ.

- **What it computes**: ExtractHLAread isolates HLA-mapped reads from aligned BAM/CRAM files and converts them to FASTQ to provide cleaner inputs for HLA typing, improving downstream SpecHLA accuracy and enabling QC-friendly batch pipelines.

```bash
iobrpy extract_hla_read \
  -s SAMPLE1 \
  -b /data/bam/SAMPLE1.bam \
  -r hg38 \
  -o /data/ExtractHLAread/SAMPLE1
```
- `-s`: desired sample name.
- `-b`: sorted and indexed bam or cram.
- `-r`: `hg19` or `hg38` reference.
- `-o`: folder to save extracted reads.

### Perform HLA typing using HLA reads

Run `iobrpy spechla` to perform HLA typing on HLA reads.

- **What it computes**: SpecHLA performs high-resolution HLA genotyping (HLA-A/B/C, DR/DQ/DP, and other loci), producing allele calls that support immune-escape assessment, transplant compatibility review, immunotherapy prognosis, and neoantigen prediction.

```bash
iobrpy spechla \
  -n SAMPLE1 \
  -1 /data/ExtractHLAread/SAMPLE1/SAMPLE1_1.fq.gz \
  -2 /data/ExtractHLAread/SAMPLE1/SAMPLE1_2.fq.gz \
  -o /data/spechla_results \
  -j 8
```
- `-n`: sample ID.
- `-1` and `-2`: The first and the second fastq file of paired-end data.
- `-o`: result directory (created if missing).
- `-j`: threads for SpecHLA (default 8).

## Batch HLA typing (`iobrpy hla_typing`)
End-to-end wrapper that runs **ExtractHLAread** for every BAM in a directory, launches **SpecHLA** on the extracted HLA reads, and merges per-sample results.

- **What it computes**: The batch workflow outputs per-sample HLA genotypes and a merged `hla_result_merged.txt`.

```bash
iobrpy hla_typing \
  -b /data/bam \
  -r hg38 \
  -o /data/hla_typing_results \
  -j 8
```

Workflow details:
- Input directory `-b` is scanned for `*.bam`; each sample gets `ExtractHLAread/<sample>` and `SpecHLA/<sample>` under the output root.
- Resume-friendly markers (`*.ExtractHLAread.done`, `*.SpecHLA.done`) for each sample to prevent reruns on completed samples.
- A merged `hla_result_merged.txt` is written alongside the per-sample folders for quick review.