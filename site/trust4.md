---
title: "TRUST4 tutorial"
page-layout: article
toc: true
---

## Overview
`iobrpy trust4` wraps the TRUST4 TCR/BCR reconstruction CLI so you can run single samples or batch runs from BAMs/FASTQs and automatically summarize immune repertoires。

- **计算了什么 / 有什么用**：基于TRUST4的组装结果，iobrpy会提取TCR/BCR克隆信息并计算克隆数、克隆频率、克隆多样性等指标，用于评估免疫克隆扩增和微环境免疫活性；输出可直接用于后续的免疫浸润比较或与临床结局关联分析。

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

### trust4_immune_indices.csv包含的指标
当所有样本的TRUST4运行结束后，会在输出根目录生成 `trust4_immune_indices.csv`，按样本汇总了以下免疫多样性参数：

- `Nreads`：用于克隆分析的总reads数；
- `Nclones`：检测到的非零克隆数量；
- `Length_CDR3`：CDR3序列长度的reads加权平均值；
- `Shannon_Index`：香农多样性指数，衡量克隆多样性；
- `Evenness`：克隆均匀度（Shannon指数归一化）；
- `Top_clone` / `Second_top_clone`：最丰富及次丰富克隆的频率；
- `Rare_clone` / `Second_Rare_clone`：最稀有及次稀有克隆的频率；
- `Gini`：克隆分布的基尼系数，越高表示集中度越强；
- `Gini_Simpson`：Gini-Simpson指数，度量多样性和不均匀度。

## Tips
- Set `TRUST4_BIN` if the TRUST4 executable is not on `PATH`.
- Use `-t/--threads` to control TRUST4 threading.
- Combine with `iobrpy runall --mode star` to reuse STAR BAMs directly with `iobrpy trust4 -b <star_output_dir>`.
