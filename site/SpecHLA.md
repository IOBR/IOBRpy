---
title: "SpecHLA & HLA typing"
page-layout: article
toc: true
---

## Single-sample SpecHLA (`iobrpy spechla`)
Run SpecHLA in RNA-seq mode when you already have paired FASTQs.

- **计算了什么 / 有什么用**：SpecHLA对样本进行HLA等位基因分型，输出高分辨率的HLA-A/B/C、DR/DQ/DP等基因型，可用于免疫逃逸/移植风险评估、免疫治疗预后分析及新抗原预测的基础数据。

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

- **计算了什么 / 有什么用**：从对齐好的BAM/CRAM中提取HLA相关reads并转换成FASTQ，保证HLA分型输入更纯净，提高后续SpecHLA结果的准确性；适合集成到批量HLA分型流程或质量复核。

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

- **计算了什么 / 有什么用**：批量自动产出每个样本的HLA分型结果并汇总成 `hla.results.txt`，便于跨队列比较、与临床特征关联，或直接作为免疫相关分析（如neoantigen预测、免疫匹配）的输入。

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
