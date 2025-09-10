# iobrpy

A Python **command‑line toolkit** for bulk RNA‑seq analysis of the tumor microenvironment (TME): data prep → signature scoring → immune deconvolution → clustering → ligand–receptor scoring.

---

## Features

**Data preparation**
- `prepare_salmon` — Clean up Salmon outputs into a TPM matrix; strip version suffixes; keep `symbol`/`ENSG`/`ENST` identifiers.
- `count2tpm` — Convert read counts to TPM (supports Ensembl/Entrez/Symbol/MGI; biomart/local annotation; effective length CSV).
- `anno_eset` — Harmonize/annotate an expression matrix (choose symbol/probe columns; deduplicate; aggregation method).
- `mouse2human_eset` — Convert mouse gene symbols to human gene symbols. Supports two modes: **matrix mode** (rows = genes) or **table mode** (input contains a symbol column). 

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

### Typical end‑to‑end workflow — output file structure examples

1) **Prepare an expression matrix**
```bash
# a) From Salmon outputs → TPM
iobrpy prepare_salmon \
  -i salmon_tpm.tsv.gz \
  -o TPM_matrix.csv \
  --return_feature symbol \
  --remove_version

Gene        TS99       TC89       TC68       TC40       813738     1929563
5S_rRNA     0.000      0.000      0.000      0.000      0.000      0.000
5_8S_rRNA   0.000      0.000      0.000      0.000      0.000      0.000
7SK         0.000      0.000      954.687    1488.249   3691.321   5399.889
A1BG        0.479      1.717      1.844      0.382      1.676      1.126
A1BG-AS1    0.149      0.348      0.755      0.000      0.314      0.400

# b) From raw gene counts → TPM
iobrpy count2tpm \
  -i counts.tsv.gz \
  -o TPM_matrix.csv \
  --idType Ensembl \
  --org hsa \
  --source local
# (Optionally provide transcript effective lengths)
#   --effLength_csv efflen.csv --id id --length eff_length --gene_symbol symbol

| Name        | SAMPLE-2e394f45066d\_20180921 | SAMPLE-88dc3e3cd88e\_20180921 | SAMPLE-b80d019c9afa\_20180921 | SAMPLE-586259880b46\_20180926 | SAMPLE-e95813c8875d\_20180921 | SAMPLE-7bd449ae436b\_20180921 |
| 5S\_rRNA    |                         5.326 |                         2.314 |                         2.377 |                         3.439 |                         6.993 |                         3.630 |
| 5\_8S\_rRNA |                         0.000 |                         0.000 |                         0.000 |                         0.000 |                         0.000 |                         0.000 |
| 7SK         |                         8.006 |                        13.969 |                        11.398 |                         5.504 |                         8.510 |                         6.418 |
| A1BG        |                         3.876 |                         2.576 |                         2.874 |                         2.533 |                         2.034 |                         2.828 |
| A1BG-AS1    |                         5.512 |                         4.440 |                         7.725 |                         4.610 |                         6.292 |                         5.336 |
```

2) (Optional) Mouse → Human symbol mapping
```bash
# Matrix mode: rows are mouse gene symbols, columns are samples
iobrpy mouse2human_eset \
  -i mouse_matrix.tsv \
  -o human_matrix.tsv \
  --is_matrix \
  --verbose

# Table mode: input has a symbol column (e.g., SYMBOL), will de-duplicate then map
iobrpy mouse2human_eset \
  -i mouse_table.csv \
  -o human_matrix.csv \
  --column_of_symbol SYMBOL \
  --verbose

Gene        Sample1    Sample2    Sample3    Sample4    Sample5    Sample6
SCMH1       0.905412   0.993271   0.826294   0.535761   0.515038   0.733388
NARF        0.116423   0.944370   0.847920   0.441993   0.736983   0.467756
CD52        0.988616   0.784523   0.303614   0.886433   0.608639   0.351713
CAV2        0.063843   0.993835   0.891718   0.702293   0.703912   0.248690
HOXB6       0.716829   0.555838   0.638682   0.971783   0.868208   0.802464

```

3) **(Optional) Annotate / de‑duplicate**
```bash
iobrpy anno_eset \
  -i TPM_matrix.csv \
  -o TPM_anno.csv \
  --annotation anno_hug133plus2 \
  --symbol symbol \
  --probe id \
  --method mean  
# You can also use: --annotation-file my_anno.csv --annotation-key gene_id

Gene        GSM1523727   GSM1523728   GSM1523729   GSM1523744   GSM1523745   GSM1523746
SH3KBP1     4.3279743    4.316195     4.3514247    4.2957463    4.2566543    4.2168822
RPL41       4.2461486    4.2468076    4.2579398    4.2955956    4.2426114    4.3464246
EEF1A1      4.2937622    4.291038     4.2621994    4.2718415    4.1992331    4.2639275
HUWE1       4.2255821    4.2111235    4.1993775    4.2192063    4.2214823    4.2046394
LOC1019288  4.2193027    4.2196698    4.2132521    4.1819267    4.2345738    4.2104611

```

4) **Signature scoring**
```bash
iobrpy calculate_sig_score \
  -i TPM_anno.csv \
  -o sig_scores.csv \
  --signature signature_collection \
  --method pca \
  --mini_gene_count 2 \
  --parallel_size 1
# Accepts space‑separated or comma‑separated groups; use "all" for a full merge.
ID          CD_8_T_effector_PCA   DDR_PCA    APM_PCA    Immune_Checkpoint_PCA   CellCycle_Reg_PCA   Pan_F_TBRs_PCA
GSM1523727  -3.003007             0.112244   1.046749   -3.287490               1.226469            -3.836552
GSM1523728  0.631973              1.138303   1.999972   0.405965                1.431343            0.164805
GSM1523729  -2.568384             -1.490780  -0.940420  -2.087635               0.579742            -1.208286
GSM1523744  -0.834788             4.558424   -0.274724  -0.873015               1.400215            -2.880584
GSM1523745  -1.358852             4.754705   -2.215926  -1.086041               1.342590            -1.054318

```

5) **Immune deconvolution (choose one or many)**
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

ID          B_cells_naive_CIBERSORT  B_cells_memory_CIBERSORT  Plasma_cells_CIBERSORT  T_cells_CD8_CIBERSORT  T_cells_CD4_naive_CIBERSORT  T_cells_CD4_memory_resting_CIBERSORT
GSM1523727  0.025261644              0.00067545                0.174139691             0.060873405             0                           0.143873862
GSM1523728  0.007497053              0.022985466               0.079320853             0.052005437             0                           0.137097071
GSM1523729  0.005356156              0.010721794               0.114171733             0                       0                           0.191541779
GSM1523744  0                        0.064645073               0.089539616             0.024437887             0                           0.147821928
GSM1523745  0                        0.014678117               0.121834835             0                       0                           0.176046775

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

6) **TME clustering / NMF clustering**
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

7) **Ligand–receptor scoring (optional)**
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

- **mouse2human\_eset**
  - `-i/--input <CSV|TSV|TXT[.gz]>` (required): input expression **matrix** or **table**
  - `-o/--output <CSV|TSV|TXT[.gz]>` (required): converted matrix indexed by **human** symbols (genes × samples)
  - `--is_matrix`: treat input as a matrix (rows = **mouse** gene symbols, columns = samples); if omitted, runs in **table mode**
  - `--column_of_symbol <str>` (required in table mode): column name that contains **mouse** gene symbols
  - `--sep <,|\t>`: override input separator; if omitted, inferred by extension.
  - `--out_sep <,|\t>`: override output separator; if omitted, inferred by **output** path extension
  - `--verbose`: print shapes and basic run info

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

- **Python version for DeSide**
  The `deside` subcommand runs **ONLY on Python 3.9**. Other versions (3.8/3.10/3.11/…) are **not supported** .When invoked via the `iobrpy CLI`, it **automatically creates/uses an isolated virtual environment with pinned dependencies** so it doesn’t leak packages from your outer env. You can override the venv location with IOBRPY_DESIDE_VENV or force a clean rebuild with IOBRPY_DESIDE_REBUILD=1; the CLI wires iobrpy into that venv through a small shim and then launches the worker. 

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