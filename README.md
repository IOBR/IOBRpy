# iobrpy

A Python **command‑line toolkit** for bulk RNA‑seq analysis of the tumor microenvironment (TME): data prep → signature scoring → immune deconvolution → clustering → ligand–receptor scoring.

---

## Features

**Data preparation**
- `prepare_salmon` — Clean up Salmon outputs into a TPM matrix; strip version suffixes; keep `symbol`/`ENSG`/`ENST` identifiers.
- `count2tpm` — Convert read counts to TPM (supports Ensembl/Entrez/Symbol/MGI; biomart/local annotation; effective length CSV).
- `anno_eset` — Harmonize/annotate an expression matrix (choose symbol/probe columns; deduplicate; aggregation method).

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

## Installation

```bash
# Option A: Conda
conda env create -f environment.yml
conda activate iobrpy

# Option B: pip
pip install -r requirements.txt
# Dev install (recommended while editing code)
pip install -e .
```

---

## Command‑line usage

### Global
```bash
iobrpy -h
iobrpy <command> --help
# Example: show help for count2tpm
iobrpy count2tpm --help
```

### General I/O conventions
- **Input orientation**: genes × samples by default.
- **Separators**: auto‑detected from file extension (`.csv` vs `.tsv`/`.txt`); you can override via command options where available.
- **Outputs**: CSV/TSV/TXT

### Typical end‑to‑end workflow

1) **Prepare an expression matrix**
```bash
# a) From Salmon outputs → TPM
iobrpy prepare_salmon -i salmon_tpm.tsv.gz -o TPM_matrix.csv --return_feature symbol --remove_version

# b) From raw gene counts → TPM
iobrpy count2tpm -i counts.tsv.gz -o TPM_matrix.csv --idType Ensembl --org hsa --source local
# (Optionally provide transcript effective lengths)
#   --effLength_csv efflen.csv --id id --length eff_length --gene_symbol symbol
```

2) **(Optional) Annotate / de‑duplicate**
```bash
iobrpy anno_eset -i TPM_matrix.csv -o TPM_anno.csv --annotation anno_hug133plus2 --symbol symbol --probe id --method mean 
# You can also use: --annotation-file my_anno.csv --annotation-key gene_id
```

3) **Signature scoring**
```bash
iobrpy calculate_sig_score -i TPM_anno.csv -o sig_scores.csv --signature signature_collection --method pca --mini_gene_count 2 --parallel_size 1
# Accepts space‑separated or comma‑separated groups; use "all" for a full merge.
```

4) **Immune deconvolution (choose one or many)**
```bash
# CIBERSORT
iobrpy cibersort -i TPM_anno.csv -o cibersort.csv --perm 100 --QN True --absolute Flase --abs_method sig.score --threads 1

# quanTIseq (method: lsei / robust norms)
iobrpy quantiseq -i TPM_anno.csv -o quantiseq.csv --signame TIL10 --method lsei --tumor --arrays --scale_mrna

# EPIC
iobrpy epic -i TPM_anno.csv -o epic.csv --reference TRef

# ESTIMATE
iobrpy estimate -i TPM_anno.csv -o estimate.csv --platform affymetrix

# MCPcounter
iobrpy mcpcounter -i TPM_anno.csv -o mcpcounter.csv --features HUGO_symbols

# IPS
iobrpy IPS -i TPM_anno.csv -o IPS.csv

# DeSide
iobrpy deside --model_dir path/to/your/DeSide_model -i TPM_anno.csv -o deside.csv --result_dir path/to/your/plot/folder --exp_type TPM --method_adding_pathway add_to_end --scaling_by_constant --transpose --print_info
```

5) **TME clustering / NMF clustering**
```bash
# KL index auto‑select k (k‑means)
iobrpy tme_cluster -i cibersort.csv -o tme_cluster.csv --features 1:22 --id "ID" --min_nc 2 --max_nc 5 --print_result --scale

# NMF clustering (auto k, excludes k=2)
iobrpy nmf -i cibersort.csv -o path/to/your/result/folder --kmin 2 --kmax 10 --features 1:22 --max-iter 10000
```

6) **Ligand–receptor scoring (optional)**
```bash
iobrpy LR_cal -i TPM_anno.csv -o LR_score.csv --data_type tpm --id_type "symbol" --cancer_type pancan --verbose
```

---

## Commands & common options

### Data preparation
- **prepare_salmon**
  - `-i/--input` Salmon combined TSV/TSV.GZ
  - `-o/--output` cleaned TPM matrix
  - `-r/--return_feature {symbol,ENSG,ENST}`
  - `--remove_version`
- **count2tpm**
  - `-i/--input`, `-o/--output`
  - `--idType {Ensembl,entrez,symbol,mgi}`, `--org {hsa,mmus}`, `--source {local,biomart}`
  - `--effLength_csv`, `--id`, `--length`, `--gene_symbol`, `--check_data`
- **anno_eset**
  - `-i/--input`, `-o/--output`
  - `--annotation anno_grch38|anno_gc_vm32` (or `--annotation-file` + `--annotation-key`)
  - `--symbol`, `--probe`, `--method mean|sd|sum`

### Signature scoring
- **calculate_sig_score**
  - `--signature` groups: `go_bp`, `go_cc`, `go_mf`, `signature_collection`, `signature_tme`, `signature_sc`, `signature_tumor`, `signature_metabolism`, `kegg`, `hallmark`, `  reactome`, or `all`
  - `--method pca|zscore|ssgsea|integration`
  - `--mini_gene_count`, `--adjust_eset`, `--parallel_size`

### Deconvolution / scoring
- **cibersort**: `--perm`, `--QN true|false`, `--absolute true|false`, `--abs_method sig.score|no.sumto1`, `--threads`
- **quantiseq**: `--arrays`, `--signame TIL10`, `--tumor`, `--scale_mrna`, `--method lsei|hampel|huber|bisquare`, `--rmgenes default|none|<list>`
- **epic**: `--reference TRef|BRef|both`
- **estimate**: `-p/--platform affymetrix|agilent|illumina`
- **mcpcounter**: `-f/--features affy133P2_probesets|HUGO_symbols|ENTREZ_ID|ENSEMBL_ID`
- **IPS**: input/output only (matrix → scores)
- **deside** (deep learning–based deconvolution)
  - **Required**: `-m/--model_dir <dir>`, `-i/--input <expr.csv/tsv>`, `-o/--output`
  - **Input format**: `--exp_type {TPM|log_space|linear}`  
  - **Pathway options**: `--gmt <one or more .gmt>`, `--method_adding_pathway {add_to_end|convert}`
  - **Scaling / transforms**: `--scaling_by_constant`, `--scaling_by_sample`, `--one_minus_alpha`
  - **Outputs & logs**: `--print_info`, `--add_cell_type`, `-r/--result_dir <dir>` (save result plots)
  - **Matrix orientation**: `--transpose` (use if your file is samples×genes instead of genes×samples)


### Clustering / decomposition
- **tme_cluster**: `--features 1:K` or regex `--pattern`, `--id`, `--scale/--no-scale`, `--min_nc`, `--max_nc`, `--max_iter`,`--print_result`,`--input_sep`,`--output_sep`
- **nmf**: `--kmin`, `--kmax`, `--features`, `--log1p`, `--normalize`, `--shift`, `--random-state`, `--max-iter`

### Ligand–receptor
- **LR_cal**: `--data_type count|tpm`, `--id_type`, `--cancer_type`, `--verbose`

---

## Troubleshooting

- **Wrong input orientation**  
  Deconvolution commands expect **genes × samples**. For `deside`, `--transpose` can be helpful depending on your file.

- **Mixed separators / encoding**  
  Prefer `.csv`, `.txt` or `.tsv` consistently. Auto‑detection works in most subcommands but you can override with explicit flags where provided.

- **DeSide model missing**
  The `deside` subcommand requires pretrained model files. If you get errors like `FileNotFoundError: DeSide_model not found` , download the official model archive from:
  https://figshare.com/articles/dataset/DeSide_model/25117862/1?file=44330255

---

## Citation & acknowledgments

This toolkit implements or wraps well‑known methods (CIBERSORT, quanTIseq, EPIC, ESTIMATE, MCPcounter, DeSide, etc.). For academic use, please cite the corresponding original papers in addition to this package.

---

## License

Add your license of choice (MIT/BSD/GPL, etc.) here.
