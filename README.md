<p align="center">
  <a href="#pypi"><img alt="PyPI" src="https://img.shields.io/badge/PyPI-pip%20install-blue"></a>
  <a href="#conda"><img alt="Conda" src="https://img.shields.io/badge/Conda-bioconda-brightgreen"></a>
  <a href="#docker"><img alt="Docker" src="https://img.shields.io/badge/Docker-pull-2496ED"></a>
</p>

# IOBRpy

**IOBRpy** is a **command-line toolkit** for bulk RNA-seq tumor microenvironment (TME) analysis. It wires together FASTQ QC, quantification (Salmon or STAR), matrix assembly, signature scoring, immune deconvolution, clustering, and ligand–receptor scoring.

![IOBRpy logo](./IOBRpy.png)

---

## Documentation

A complete documentation for IOBRpy can be found at https://iobr.github.io/IOBRpy/.

---

## Installation

### Quick install

```bash
# Method 1 : PyPI
pip install iobrpy

# Method 2 : Conda (bioconda via conda-forge + bioconda)
conda install -c conda-forge -c bioconda iobrpy=0.1.4

# Method 3 : Docker
docker pull hhn123123/iobrpy:latest
```

### PyPI

<details><summary><strong>Show full PyPI steps</strong></summary>

```bash
# Creating a virtual environment is recommended
conda create -n iobrpy python=3.9 -y
conda activate iobrpy
```
```bash
# Update pip
python -m pip install --upgrade pip
```
```bash
# Install iobrpy
pip install iobrpy
```
```bash
# Install fastp, salmon, STAR and MultiQC
# Recommended: use mamba for faster solves (if available)
mamba install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  trust4

# If you don't have mamba, use conda instead
conda install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  trust4
```
</details>

### Conda

> **Prerequisite (Conda):** Please install Miniconda or Anaconda first. We recommend [Miniconda](https://docs.anaconda.com/miniconda/).

<details><summary><strong>Show full Conda steps</strong></summary>

```bash
# Creating a virtual environment is recommended
conda create -n iobrpy python=3.9 -y
conda activate iobrpy
```
```bash
# Install iobrpy 0.1.4 (from bioconda via conda-forge + bioconda)
# Recommended: use mamba for faster solves (if available)
mamba install -y -c conda-forge -c bioconda iobrpy=0.1.4

# If you don't have mamba, use conda instead
conda install -y -c conda-forge -c bioconda iobrpy=0.1.4
```
</details>

### Docker

> **Docker Hub website:** [Docker Hub](https://hub.docker.com/)

<details><summary><strong>Show Docker pull</strong></summary>

```bash
# Option 1: Pull the latest image from Docker Hub
docker pull hhn123123/iobrpy:latest
```
```bash
# Option 2: Offline install (from GitHub Release)
# 1) Download iobrpy.tar.gz from https://github.com/IOBR/IOBRpy/releases/tag/v1.0.0
# 2) Change to the directory where the archive is saved and load the image
cd /path/to/iobrpy.tar.gz
docker load -i iobrpy.tar.gz
```
</details>

---

## Features

**End-to-End Pipeline Runner**
- `runall` — A single command that wires the full Salmon or STAR pipeline end-to-end and writes the standardized layout:
  The pipeline creates the following directories, in order: `01-qc/`, `02-salmon/` or `02-star/`, `03-tpm/`, `04-signatures/`, `05-tme/`, and `06-LR_cal/`.

**All-in-one TME profiling**
- `tme_profile` - A single command that inputs a TPM (genes×samples) matrix, performs signature scoring, runs six immune deconvolution methods, merges their outputs, and computes ligand–receptor scores, using the functions `calculate_sig_score`, `cibersort`, `IPS`, `estimate`, `mcpcounter`, `quantiseq`, `epic`, and `LR_cal`.

**Preprocessing**
- `fastq_qc` — Parallel FASTQ QC/trimming via **fastp**, with per-sample HTML/JSON and an optional **MultiQC** summary report under `01-qc/multiqc_report/`. Resume-friendly and prints output paths first.

**Salmon submodule (quantification, merge, and TPM)**
- `batch_salmon` — Batch **salmon quant** on paired-end FASTQs; safe R1/R2 inference; per-sample `quant.sf`; progress and preflight checks (salmon version, index meta).  
- `merge_salmon` — Recursively collect per-sample `quant.sf` and produce two matrices: **TPM** and **NumReads**.  
- `prepare_salmon` — Clean up Salmon outputs into a TPM matrix; strip version suffixes; keep `symbol`/`ENSG`/`ENST` identifiers.

**STAR submodule (alignment, counts, and TPM)**
- `batch_star_count` — Batch **STAR** alignment with `--quantMode GeneCounts`, sorted BAM + `_ReadsPerGene.out.tab`; resume-friendly summary.  
- `merge_star_count` — Merge multiple `_ReadsPerGene.out.tab` into one wide count matrix.  
- `count2tpm` — Convert counts to TPM (supports Ensembl/Entrez/Symbol/MGI; optional effective length CSV).

**Expression Annotation & Mouse to Human Mapping & log2(x+1) (Optional)**
- `anno_eset` — Harmonize/annotate an expression matrix (choose symbol/probe columns; deduplicate; aggregation method).
- `mouse2human_eset` — Convert mouse gene symbols to human gene symbols. Supports two modes: **matrix mode** (rows = genes) or **table mode** (input contains a symbol column). 
- `log2_eset` — Apply log2(x+1) to a **genes × samples** expression matrix.

**Pathway / signature scoring**
- `calculate_sig_score` — Sample‑level signature scores via `pca`, `zscore`, `ssgsea`, or `integration`. 
  Supports the following signature **groups** (space‑ or comma‑separated), or `all` to merge them:
  - `go_bp`, `go_cc`, `go_mf`
  - `signature_collection`, `signature_tme`, `signature_sc`, `signature_tumor`, `signature_metabolism`
  - `kegg`, `hallmark`, `reactome`

**Immune deconvolution and scoring**
- `cibersort` — CIBERSORT wrapper/implementation with permutations, quantile normalization, absolute mode.
- `quantiseq` — quanTIseq deconvolution with `lsei` or robust norms (`hampel`, `huber`, `bisquare`); tumor‑gene filtering; mRNA scaling.
- `epic` — EPIC cell fractions using `TRef`/`BRef` references.
- `estimate` — ESTIMATE immune/stromal/tumor purity scores.
- `mcpcounter` — MCPcounter infiltration scores.
- `IPS` — Immunophenoscore (AZ/SC/CP/EC + total).
- `deside` — Deep learning–based deconvolution (requires pre‑downloaded model; supports pathway‑masked mode via KEGG/Reactome GMTs).

**Clustering / decomposition**
- `tme_cluster` — k‑means with **automatic k** via KL index (Hartigan–Wong), feature selection and standardization.
- `nmf` — NMF‑based clustering (auto‑selects k; excludes k=2) with PCA plot and top features.

**Ligand–receptor**
- `LR_cal` — Ligand–receptor interaction scoring using cancer‑type specific networks.
---

## Input Requirements
- **FASTQ layout**: paired-end by default. Filenames end with `*_1.fastq.gz` / `*_2.fastq.gz` (configurable via `--suffix1`).
- **Expression matrix orientation**: **genes × samples** by default.
- **Output file delimiters**: automatically inferred from the file extension; .csv and .tsv/.txt are recommended.

---

## Command‑line usage

### From FASTQ to TME - `runall`

#### How `runall` passes options
`runall` defines a small set of top-level options (e.g., `--mode/--outdir/--fastq/--threads/--batch_size`). Any unrecognized options are forwarded to the corresponding sub-steps. This keeps `runall` flexible as sub-commands evolve.

Below are **two fully wired workflows** handled by `iobrpy runall`.  

#### Salmon mode
```bash
iobrpy runall \
  --mode salmon \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/salmon/index" \
  --project MyProj
```
#### STAR mode
```bash
iobrpy runall \
  --mode star \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/star/index" \
  --project MyProj
```

---

### Option legend for the `runall` examples

#### Common options
| Flag | Purpose |
|---|---|
| `--mode {salmon / star}` | Select backend (Salmon quant vs. STAR align+count) |
| `--outdir <DIR>` | Root output directory (creates the standardized layout) |
| `--fastq <DIR>` | Raw FASTQ dir |
| `--index <DIR>` | Salmon : path to **Salmon index**; STAR : path to **STAR index** |
| `--project <STR>` | Prefix for merged outputs |
| `--threads <INT>` / `--batch_size <INT>` | Global concurrency/batching |

---

### Expected layout
```
# Salmon mode：
/path/to/outdir
|-- 01-qc
|   |-- <sample>_1.fastq.gz
|   |-- <sample>_2.fastq.gz
|   |-- <sample>_fastp.html
|   |-- <sample>_fastp.json
|   |-- <sample>.task.complete
|   `-- multiqc_report
|       `-- multiqc_fastp_report.html
|-- 02-salmon
|   |-- <sample>
|   |   `-- quant.sf
|   |-- MyProj_salmon_count.tsv.gz
|   `-- MyProj_salmon_tpm.tsv.gz
|-- 03-tpm
|   |-- prepare_salmon.csv
|   `-- tpm_matrix.csv
|-- 04-signatures
|   `-- calculate_sig_score.csv
|-- 05-tme
|   |-- cibersort_results.csv
|   |-- epic_results.csv
|   |-- quantiseq_results.csv
|   |-- IPS_results.csv
|   |-- estimate_results.csv
|   |-- mcpcounter_results.csv
|   `-- deconvo_merged.csv
`-- 06-LR_cal
    `-- lr_cal.csv
# STAR mode：
/path/to/outdir
|-- 01-qc
|   |-- <sample>_1.fastq.gz
|   |-- <sample>_2.fastq.gz
|   |-- <sample>_fastp.html
|   |-- <sample>_fastp.json
|   |-- <sample>.task.complete
|   `-- multiqc_report
|       `-- multiqc_fastp_report.html
|-- 02-star
|   |-- <sample>/
|   |-- <sample>__STARgenome/
|   |-- <sample>__STARpass1/
|   |-- <sample>_STARtmp/
|   |-- <sample>_Aligned.sortedByCoord.out.bam
|   |-- <sample>_Log.final.out
|   |-- <sample>_Log.out
|   |-- <sample>_Log.progress.out
|   |-- <sample>_ReadsPerGene.out.tab
|   |-- <sample>_SJ.out.tab
|   |-- <sample>.task.complete
|   |-- .batch_star_count.done
|   |-- .merge_star_count.done
|   `-- MyProj.STAR.count.tsv.gz
|-- 03-tpm
|   |-- count2tpm.csv
|   `-- tpm_matrix.csv
|-- 04-signatures
|   `-- calculate_sig_score.csv
|-- 05-tme
|   |-- cibersort_results.csv
|   |-- epic_results.csv
|   |-- quantiseq_results.csv
|   |-- IPS_results.csv
|   |-- estimate_results.csv
|   |-- mcpcounter_results.csv
|   `-- deconvo_merged.csv
`-- 06-LR_cal
    `-- lr_cal.csv
```

---

### Output Reference

#### Standard layout (produced by `iobrpy runall`)
- `01-qc/` — fastp outputs; a resume flag `.fastq_qc.done` is written when the step completes.
- `02-salmon/` **or** `02-star/` — quantification/alignment + merged matrices; resume flags like `.batch_salmon.done`, `.merge_salmon.done`, or `.merge_star_count.done`.
- `03-tpm/` — unified TPM matrix `tpm_matrix.csv`. For Salmon mode it comes from `prepare_salmon`; for STAR mode it comes from `count2tpm`.
- `04-signatures/` — signature scoring results (file: `calculate_sig_score.csv`).
- `05-tme/` — deconvolution outputs from multiple methods + `deconvo_merged.csv`.
- `06-LR_cal/` — ligand–receptor results `lr_cal.csv`.

#### Salmon mode (`02-salmon/`)
- Per-sample Salmon folders containing `quant.sf` (from `batch_salmon`). A `.batch_salmon.done` flag is written after completion.
- Merged matrices (from `merge_salmon`):
  - `<PROJECT>_salmon_tpm.tsv[.gz]`
  - `<PROJECT>_salmon_count.tsv[.gz]`  
  A `.merge_salmon.done` flag is written after completion.
- `03-tpm/prepare_salmon.csv` — cleaned genes × samples TPM matrix produced by `prepare_salmon` (default `--return_feature symbol` unless overridden).
- `03-tpm/tpm_matrix.csv` — **log2(x+1)** matrix produced by `log2_eset` from `prepare_salmon.csv`.

#### STAR mode (`02-star/`)
- Per-sample STAR outputs (BAM, logs, `*_ReadsPerGene.out.tab`, etc.).
- Merged counts (from `merge_star_count`):
  - `<PROJECT>.STAR.count.tsv.gz` . A `.merge_star_count.done` flag is written after completion.
- `03-tpm/count2tpm.csv` — TPM matrix produced by `count2tpm` from the merged STAR ReadPerGene/count matrix.
- `03-tpm/tpm_matrix.csv` — **log2(x+1)** matrix produced by `log2_eset` from `count2tpm.csv`.

#### Signatures (`04-signatures/`)
- `calculate_sig_score.csv` — per-sample pathway/signature scores. Columns correspond to the selected signature set and method (`integration`, `pca`, `zscore`, or `ssgsea`). 

#### Deconvolution (`05-tme/`)
Each method writes a single table named `<method>_results.csv`:

- `cibersort_results.csv` — columns suffixed with `_CIBERSORT`. Note whether `--perm` and `--QN` were used.
- `quantiseq_results.csv` — quanTIseq fractions. Document the chosen `--method {lsei|hampel|huber|bisquare}` and flags like `--arrays`, `--tumor`, `--scale_mrna`, `--signame`.
- `epic_results.csv` — EPIC fractions; record the reference profile used (`--reference {TRef|BRef|both}`).
- `estimate_results.csv` — ESTIMATE immune/stromal/purity scores; columns suffixed `_estimate`.
- `mcpcounter_results.csv` — MCPcounter scores; columns suffixed `_MCPcounter`.
- `IPS_results.csv` — IPS sub-scores and total score.

**Merged table**
- `deconvo_merged.csv` — produced by `runall` after all deconvolution methods finish; normalizes the sample index to a column named `ID` and outer-joins by sample ID across methods.

#### Ligand–receptor (`06-LR_cal/`)
- `lr_cal.csv` — ligand–receptor scoring table from `LR_cal`. Record the `--data_type {count|tpm}` and the `--id_type` you used.

---

## Contact / Support
- Issues: https://github.com/IOBR/IOBRpy/issues
- Maintainers: [ Haonan Huang ] (email = 2905611068@qq.com); [ Dongqiang Zeng ] (email = interlaken@smu.edu.cn)