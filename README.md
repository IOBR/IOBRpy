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
# creating a virtual environment is recommended
conda create -n iobrpy python=3.9
conda activate iobrpy
# update pip
python3 -m pip install --upgrade pip
# install deside
pip install iobrpy
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
iobrpy prepare_salmon \
  -i salmon_tpm.tsv.gz \
  -o TPM_matrix.csv \
  --return_feature symbol \
  --remove_version

# b) From raw gene counts → TPM
iobrpy count2tpm \
  -i counts.tsv.gz \
  -o TPM_matrix.csv \
  --idType Ensembl \
  --org hsa \
  --source local
# (Optionally provide transcript effective lengths)
#   --effLength_csv efflen.csv --id id --length eff_length --gene_symbol symbol
```

2) **(Optional) Annotate / de‑duplicate**
```bash
iobrpy anno_eset \
  -i TPM_matrix.csv \
  -o TPM_anno.csv \
  --annotation anno_hug133plus2 \
  --symbol symbol \
  --probe id \
  --method mean  
# You can also use: --annotation-file my_anno.csv --annotation-key gene_id
```

3) **Signature scoring**
```bash
iobrpy calculate_sig_score \
  -i TPM_anno.csv \
  -o sig_scores.csv \
  --signature signature_collection \
  --method pca \
  --mini_gene_count 2 \
  --parallel_size 1
# Accepts space‑separated or comma‑separated groups; use "all" for a full merge.
```

4) **Immune deconvolution (choose one or many)**
```bash
# CIBERSORT
iobrpy cibersort \
  -i TPM_anno.csv \
  -o cibersort.csv \
  --perm 100 \
  --QN True \
  --absolute False \
  --abs_method sig.score \
  --threads 1

# quanTIseq (method: lsei / robust norms)
iobrpy quantiseq \
  -i TPM_anno.csv \
  -o quantiseq.csv \
  --signame TIL10 \
  --method lsei \
  --tumor \
  --arrays \
  --scale_mrna

# EPIC
iobrpy epic \
  -i TPM_anno.csv \
  -o epic.csv \
  --reference TRef

# ESTIMATE
iobrpy estimate \
  -i TPM_anno.csv \
  -o estimate.csv \
  --platform affymetrix

# MCPcounter
iobrpy mcpcounter \
  -i TPM_anno.csv \
  -o mcpcounter.csv \
  --features HUGO_symbols

# IPS
iobrpy IPS \
  -i TPM_anno.csv \
  -o IPS.csv

# DeSide
iobrpy deside \
  --model_dir path/to/your/DeSide_model \
  -i TPM_anno.csv \
  -o deside.csv \
  -r path/to/your/plot/folder \
  --exp_type TPM \
  --method_adding_pathway add_to_end \
  --scaling_by_constant \
  --transpose \
  --print_info
```

5) **TME clustering / NMF clustering**
```bash
# KL index auto‑select k (k‑means)
iobrpy tme_cluster \
  -i cibersort.csv \
  -o tme_cluster.csv \
  --features 1:22 \
  --id ID \
  --min_nc 2 \
  --max_nc 5 \
  --print_result \
  --scale

# NMF clustering (auto k, excludes k=2)
iobrpy nmf \
  -i cibersort.csv \
  -o path/to/your/result/folder \
  --kmin 2 \
  --kmax 10 \
  --features 1:22 \
  --max-iter 10000 \
  --skip_k_2
```

6) **Ligand–receptor scoring (optional)**
```bash
iobrpy LR_cal \
  -i TPM_anno.csv \
  -o LR_score.csv \
  --data_type tpm \
  --id_type symbol \
  --cancer_type pancan \
  --verbose
```

---

## Commands & common options

### Data preparation
- **prepare_salmon**
  - `-i/--input <TSV|TSV.GZ>` (required): Salmon-combined gene TPM table
  - `-o/--output <CSV/TSV>` (required): cleaned TPM matrix (genes × samples)
  - `-r/--return_feature {ENST|ENSG|symbol}` (default: `symbol`): which identifier to keep
  - `--remove_version`: strip version suffix from gene IDs (e.g., `ENSG000001.12 → ENSG000001`)

- **count2tpm**
  - `-i/--input <CSV/TSV[.gz]>` (required): raw count matrix (genes × samples)
  - `-o/--output <CSV/TSV>` (required): output TPM matrix
  - `--effLength_csv <CSV>`: optional effective-length file with columns `id`, `eff_length`, `symbol`
  - `--idType {Ensembl|entrez|symbol|mgi}` (default: `Ensembl`)
  - `--org {hsa|mmus}` (default: `hsa`)
  - `--source {local|biomart}` (default: `local`)
  - `--id <str>` (default: `id`): ID column name in `--effLength_csv`
  - `--length <str>` (default: `eff_length`): length column
  - `--gene_symbol <str>` (default: `symbol`): gene symbol column
  - `--check_data`: check & drop missing/invalid entries before conversion

- **anno_eset**
  - `-i/--input <CSV/TSV/TXT>` (required)
  - `-o/--output <CSV/TSV/TXT>` (required)
  - `--annotation {anno_hug133plus2|anno_rnaseq|anno_illumina|anno_grch38}` (required unless using external file)
  - `--annotation-file <pkl/csv/tsv/xlsx>`: external annotation (overrides built-in)
  - `--annotation-key <str>`: key to pick a table if external `.pkl` stores a dict of DataFrames
  - `--symbol <str>` (default: `symbol`): column used as gene symbol
  - `--probe  <str>` (default: `id`): column used as probe/feature ID
  - `--method {mean|sd|sum}` (default: `mean`): duplicate-ID aggregation

### Signature scoring
- **calculate_sig_score**
  - `-i/--input <CSV/TSV/TXT>` (required), `-o/--output <CSV/TSV/TXT>` (required)
  - `--signature <one or more groups>` (required; space- or comma-separated; `all` uses every group)  
    Groups: `go_bp`, `go_cc`, `go_mf`, `signature_collection`, `signature_tme`, `signature_sc`, `signature_tumor`, `signature_metabolism`, `kegg`, `hallmark`, `reactome`
  - `--method {pca|zscore|ssgsea|integration}` (default: `pca`)
  - `--mini_gene_count <int>` (default: `3`)
  - `--adjust_eset`: apply extra filtering after log2 transform
  - `--parallel_size <int>` (default: `1`; threads for `ssgsea`)

### Deconvolution / scoring
- **cibersort**
  - `-i/--input <CSV/TSV>` (required), `-o/--output <CSV/TSV>` (required)
  - `--perm <int>` (default: `100`)
  - `--QN <True|False>` (default: `True`): quantile normalization
  - `--absolute <True|False>` (default: `False`): absolute mode
  - `--abs_method {sig.score|no.sumto1}` (default: `sig.score`)
  - `--threads <int>` (default: `1`)  
  *Output: columns are suffixed with `_CIBERSORT`, index name is `ID`, separator inferred from output extension.*

- **quantiseq**
  - `-i/--input <CSV/TSV>` (required; genes × samples), `-o/--output <TSV>` (required)
  - `--arrays`: perform quantile normalization for arrays
  - `--signame <str>` (default: `TIL10`)
  - `--tumor`: remove genes highly expressed in tumors
  - `--scale_mrna`: enable mRNA scaling (otherwise raw signature proportions)
  - `--method {lsei|hampel|huber|bisquare}` (default: `lsei`)
  - `--rmgenes <str>` (default: `unassigned`; allowed: `default`, `none`, or comma-separated list)

- **epic**
  - `-i/--input <CSV/TSV>` (required; genes × samples)
  - `-o/--output <CSV/TSV>` (required)
  - `--reference {TRef|BRef|both}` (default: `TRef`)

- **estimate**
  - `-i/--input <CSV/TSV/TXT>` (required; genes × samples)
  - `-p/--platform {affymetrix|agilent|illumina}` (default: `affymetrix`)
  - `-o/--output <CSV/TSV/TXT>` (required)  
  *Output is transposed; columns are suffixed with `_estimate`; index label is `ID`; separator inferred from extension.*

- **mcpcounter**
  - `-i/--input <TSV>` (required; genes × samples)
  - `-f/--features {affy133P2_probesets|HUGO_symbols|ENTREZ_ID|ENSEMBL_ID}` (required)
  - `-o/--output <CSV/TSV>` (required)  
  *Output: columns normalized (spaces → `_`) and suffixed with `_MCPcounter`; index label `ID`; separator inferred from extension.*

- **IPS**
  - `-i/--input <matrix>` (required), `-o/--output <file>` (required)  
  *No extra flags (expression matrix → IPS sub-scores + total).*

- **deside** (deep learning–based deconvolution)
  - `-m/--model_dir <dir>` (required): path to the pre-downloaded DeSide model directory
  - `-i/--input <CSV/TSV>` (required): rows = genes, columns = samples
  - `-o/--output <CSV>` (required)
  - `--exp_type {TPM|log_space|linear}` (default: `TPM`)  
    - `TPM`: already log2 processed  
    - `log_space`: `log2(TPM+1)`  
    - `linear`: linear space (TPM/counts)
  - `--gmt <file1.gmt file2.gmt ...>`: optional one or more GMT files for pathway masking
  - `--method_adding_pathway {add_to_end|convert}` (default: `add_to_end`)
  - `--scaling_by_constant`, `--scaling_by_sample`, `--one_minus_alpha`: optional scaling/transforms
  - `--print_info`: verbose logs
  - `--add_cell_type`: append predicted cell-type labels
  - `--transpose`: use if your file is *samples × genes*
  - `-r/--result_dir <dir>`: optional directory to save result plots/logs


### Clustering / decomposition
- **tme_cluster**
  - `-i/--input <CSV/TSV/TXT>` (required): input table for clustering.
    - Expected shape: first column = sample ID (use `--id` if not first), remaining columns = features.
  - `-o/--output <CSV/TSV/TXT>` (required): output file for clustering results.
  - `--features <spec>`: select feature columns by 1-based inclusive range, e.g. `1:22` (intended for CIBERSORT outputs; **exclude** the sample ID column when counting).
  - `--pattern <regex>`: alternatively select features by a regex on column names (e.g. `^CD8|^NK`).  
    *Tip: use one of `--features` or `--pattern`.*
  - `--id <str>` (default: first column): column name containing sample IDs.
  - `--scale` / `--no-scale`: toggle z-score scaling of features (help text: default = **True**).
  - `--min_nc <int>` (default: `2`): minimum number of clusters to try.
  - `--max_nc <int>` (default: `6`): maximum number of clusters to try.
  - `--max_iter <int>` (default: `10`): maximum iterations for k-means.
  - `--tol <float>` (default: `1e-4`): convergence tolerance for centroid updates.
  - `--print_result`: print intermediate KL scores and cluster counts.
  - `--input_sep <str>` (default: auto): input delimiter (e.g. `,` or `\t`); auto-detected if unset.
  - `--output_sep <str>` (default: auto): output delimiter; inferred from filename if unset.

- **nmf**
  - `-i/--input <CSV/TSV>` (required): matrix to factorize; first column should be sample names (index).
  - `-o/--output <DIR>` (required): directory to save results.
  - `--kmin <int>` (default: `2`): minimum `k` (inclusive).
  - `--kmax <int>` (default: `8`): maximum `k` (inclusive).
  - `--features <spec>`: 1-based inclusive selection of feature columns (e.g. `2-10` or `1:5`), typically cell-type columns.
  - `--log1p`: apply `log1p` to the input (useful for counts).
  - `--normalize`: L1 row normalization (each sample sums to 1).
  - `--shift <float>` (default: `None`): if data contain negatives, add a constant to make all values non-negative.
  - `--random-state <int>` (default: `42`): random seed for NMF.
  - `--max-iter <int>` (default: `1000`): NMF max iterations.
  - `--skip_k_2`: skip evaluating `k = 2` when searching for the best `k`.

### Ligand–receptor
- **LR_cal**
  - `-i/--input <CSV/TSV>` (required): expression matrix (genes × samples).
  - `-o/--output <CSV/TSV>` (required): file to save LR scores.
  - `--data_type {count|tpm}` (default: `tpm`): type of the input matrix.
  - `--id_type <str>` (default: `ensembl`): gene ID type expected by the LR backend.
  - `--cancer_type <str>` (default: `pancan`): cancer-type network to use.
  - `--verbose`: verbose logging.

---

## Troubleshooting

- **Wrong input orientation**  
  Deconvolution commands expect **genes × samples**. For `deside`, `--transpose` can be helpful depending on your file.

- **Mixed separators / encoding**  
  Prefer `.csv` , `.txt` or `.tsv` consistently. Auto‑detection works in most subcommands but you can override with explicit flags where provided.

- **DeSide model missing**
  The `deside` subcommand requires pretrained model files. If you get errors like `FileNotFoundError: DeSide_model not found` , download the official model archive from:
  https://figshare.com/articles/dataset/DeSide_model/25117862/1?file=44330255

---

## Citation & acknowledgments

This toolkit implements or wraps well‑known methods (CIBERSORT, quanTIseq, EPIC, ESTIMATE, MCPcounter, DeSide, etc.). For academic use, please cite the corresponding original papers in addition to this package.

---

## License

MIT License

Copyright (c) 2024 Dongqiang Zeng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.