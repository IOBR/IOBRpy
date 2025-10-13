---
title: From FASTQ to TME - runall
layout: default
nav_order: 4
---

# **From FASTQ to TME - `runall`**

## How `runall` passes options
`runall` defines a small set of top-level options (e.g., `--mode/--outdir/--fastq/--threads/--batch_size`). Any unrecognized options are forwarded to the corresponding sub-steps. This keeps `runall` flexible as sub-commands evolve.

Below are **two fully wired workflows** handled by `iobrpy runall`.  

## Salmon mode
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
## STAR mode
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

## Option legend for the `runall` examples

### Common options

- `--mode {salmon|star}` — Select backend (Salmon quant vs. STAR align+count)
- `--outdir <DIR>` — Root output directory (creates the standardized layout)
- `--fastq <DIR>` — Raw FASTQ dir, forwarded to `fastq_qc --path1_fastq`
- `--threads <INT>` / `--batch_size <INT>` — Global concurrency / batching
- `--resume` — Skip steps whose outputs already exist
- `--dry_run` — Print planned commands without executing

### Salmon-only

- `--index <DIR>` — Salmon index for `batch_salmon`
- `--project <STR>` — Prefix for merged outputs in `merge_salmon`
- `--return_feature {symbol|ENSG|ENST}` — Output gene ID type in `prepare_salmon`
- `--remove_version` — Strip version suffix in `prepare_salmon`

### STAR-only

- `--index <DIR>` — STAR genomeDir for `batch_star_count`
- `--project <STR>` — Prefix for merged counts in `merge_star_count`
- `--idtype {ensembl|entrez|symbol|mgi}` — Gene ID type for `count2tpm`
- `--org {hsa|mmus}` — Organism for `count2tpm`
- `--remove_version` — Strip version suffix before `count2tpm`

### Signature scoring

- `--method {integration|pca|zscore|ssgsea}` — Scoring method for `calculate_sig_score`
- `--signature <set>` — Which signature set to use (`all`, etc.)
- `--mini_gene_count <INT>` — Minimum genes per signature
- `--adjust_eset` — Extra filtering after log transform

### Deconvolution

- `--perm <INT>` / `--QN {true|false}` — CIBERSORT permutations / quantile normalization
- `--platform <STR>` — ESTIMATE platform
- `--features HUGO_symbols` — MCPcounter feature type
- `--arrays` `--tumor` `--scale_mrna` — quanTIseq options
- `--reference {TRef|BRef|both}` — EPIC reference profile

### Ligand–receptor

- `--data_type {tpm|count}` — Input matrix type for `LR_cal`
- `--id_type {symbol|ensembl|...}` — Gene ID type for `LR_cal`
- `--verbose` — Verbose logging

---

## Expected layout
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

## Output Reference

### Standard layout (produced by `iobrpy runall`)
- `01-qc/` — fastp outputs; a resume flag `.fastq_qc.done` is written when the step completes.
- `02-salmon/` **or** `02-star/` — quantification/alignment + merged matrices; resume flags like `.batch_salmon.done`, `.merge_salmon.done`, or `.merge_star_count.done`.
- `03-tpm/` — unified TPM matrix `tpm_matrix.csv`. For Salmon mode it comes from `prepare_salmon`; for STAR mode it comes from `count2tpm`.
- `04-signatures/` — signature scoring results (file: `calculate_sig_score.csv`).
- `05-tme/` — deconvolution outputs from multiple methods + `deconvo_merged.csv`.
- `06-LR_cal/` — ligand–receptor results `lr_cal.csv`.

### Salmon mode (`02-salmon/`)
- Per-sample Salmon folders containing `quant.sf` (from `batch_salmon`). A `.batch_salmon.done` flag is written after completion.
- Merged matrices (from `merge_salmon`):
  - `<PROJECT>_salmon_tpm.tsv[.gz]`
  - `<PROJECT>_salmon_count.tsv[.gz]`  
  A `.merge_salmon.done` flag is written after completion.
- `03-tpm/prepare_salmon.csv` — cleaned genes × samples TPM matrix produced by `prepare_salmon` (default `--return_feature symbol` unless overridden).
- `03-tpm/tpm_matrix.csv` — **log2(x+1)** matrix produced by `log2_eset` from `prepare_salmon.csv`.

### STAR mode (`02-star/`)
- Per-sample STAR outputs (BAM, logs, `*_ReadsPerGene.out.tab`, etc.).
- Merged counts (from `merge_star_count`):
  - `<PROJECT>.STAR.count.tsv.gz` . A `.merge_star_count.done` flag is written after completion.
- `03-tpm/count2tpm.csv` — TPM matrix produced by `count2tpm` from the merged STAR ReadPerGene/count matrix.
- `03-tpm/tpm_matrix.csv` — **log2(x+1)** matrix produced by `log2_eset` from `count2tpm.csv`.

### Signatures (`04-signatures/`)
- `calculate_sig_score.csv` — per-sample pathway/signature scores. Columns correspond to the selected signature set and method (`integration`, `pca`, `zscore`, or `ssgsea`). 

### Deconvolution (`05-tme/`)
Each method writes a single table named `<method>_results.csv`:

- `cibersort_results.csv` — columns suffixed with `_CIBERSORT`. Note whether `--perm` and `--QN` were used.
- `quantiseq_results.csv` — quanTIseq fractions. Document the chosen `--method {lsei|hampel|huber|bisquare}` and flags like `--arrays`, `--tumor`, `--scale_mrna`, `--signame`.
- `epic_results.csv` — EPIC fractions; record the reference profile used (`--reference {TRef|BRef|both}`).
- `estimate_results.csv` — ESTIMATE immune/stromal/purity scores; columns suffixed `_estimate`.
- `mcpcounter_results.csv` — MCPcounter scores; columns suffixed `_MCPcounter`.
- `IPS_results.csv` — IPS sub-scores and total score.

**Merged table**
- `deconvo_merged.csv` — produced by `runall` after all deconvolution methods finish; normalizes the sample index to a column named `ID` and outer-joins by sample ID across methods.

### Ligand–receptor (`06-LR_cal/`)
- `lr_cal.csv` — ligand–receptor scoring table from `LR_cal`. Record the `--data_type {count|tpm}` and the `--id_type` you used.