# IOBRpy

**IOBRpy** is a **command-line toolkit** for bulk RNA-seq tumor microenvironment (TME) analysis. It wires together FASTQ QC, quantification (Salmon or STAR), matrix assembly, signature scoring, immune deconvolution, clustering, and ligand–receptor scoring.

![IOBRpy logo](./IOBRpy.png)

---

## Features

**End-to-End Pipeline Runner**
- `runall` — A single command that wires the full Salmon or STAR pipeline end-to-end and writes the standardized layout:
  The pipeline creates the following directories, in order: `01-qc/`, `02-salmon/` or `02-star/`, `03-tpm/`, `04-signatures/`, `05-tme/`, and `06-LR_cal/`.

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

**Expression Annotation & Mouse to Human Mapping(Optional)**
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

### `runall` — From FASTQ to TME

Below are **two fully wired workflows** handled by `iobrpy runall`.  

#### Salmon mode
```bash
iobrpy runall \
  --mode salmon \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 16 \
  --batch_size 4 \
  --index "/path/to/salmon/index" \
  --project MyProj \
  --return_feature symbol \
  --remove_version \
  --method integration \
  --signature all \
  --mini_gene_count 2 \
  --adjust_eset \
  --perm 1000 \
  --QN true \
  --platform affymetrix \
  --features HUGO_symbols \
  --arrays \
  --tumor \
  --scale_mrna \
  --reference TRef \
  --data_type tpm \
  --id_type "symbol" \
  --verbose
```
#### STAR mode
```bash
iobrpy runall \
  --mode star \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 16 \
  --batch_size 1 \
  --index "/path/to/star/index" \
  --project MyProj \
  --idtype ensembl \
  --org hsa \
  --remove_version \
  --method integration \
  --signature all \
  --mini_gene_count 2 \
  --adjust_eset \
  --perm 100 \
  --QN true \
  --platform affymetrix \
  --features HUGO_symbols \
  --arrays \
  --tumor \
  --scale_mrna \
  --reference TRef \
  --data_type tpm \
  --id_type "symbol" \
  --verbose
```
### Option legend for the `runall` examples

#### Common options
| Flag | Purpose |
|---|---|
| `--mode {salmon|star}` | Select backend (Salmon quant vs. STAR align+count) |
| `--outdir <DIR>` | Root output directory (creates the standardized layout) |
| `--fastq <DIR>` | Raw FASTQ dir, forwarded to `fastq_qc --path1_fastq` |
| `--threads <INT>` / `--batch_size <INT>` | Global concurrency/batching |
| `--resume` | Skip steps whose outputs already exist |
| `--dry_run` | Print planned commands without executing |

#### Salmon-only
| Flag | Purpose |
|---|---|
| `--index <DIR>` | Salmon index for `batch_salmon` |
| `--project <STR>` | Prefix for merged outputs in `merge_salmon` |
| `--return_feature {symbol / ENSG / ENST}` | Output gene ID type in `prepare_salmon` |
| `--remove_version` | Strip version suffix in `prepare_salmon` |

#### STAR-only
| Flag | Purpose |
|---|---|
| `--index <DIR>` | STAR genomeDir for `batch_star_count` |
| `--project <STR>` | Prefix for merged counts in `merge_star_count` |
| `--idtype {ensembl / entrez / symbol / mgi}` | Gene ID type for `count2tpm` |
| `--org {hsa / mmus}` | Organism for `count2tpm` |
| `--remove_version` | Strip version suffix before `count2tpm` |

#### Signature scoring
| Flag | Purpose |
|---|---|
| `--method {integration / pca / zscore / ssgsea}` | Scoring method for `calculate_sig_score` |
| `--signature <set>` | Which signature set to use (`all`, etc.) |
| `--mini_gene_count <INT>` | Min genes per signature |
| `--adjust_eset` | Extra filtering after log transform |

#### Deconvolution
| Flag | Purpose |
|---|---|
| `--perm <INT>` / `--QN {true / false}` | CIBERSORT permutations / quantile normalization |
| `--platform <str>` | ESTIMATE platform |
| `--features HUGO_symbols` | MCPcounter features |
| `--arrays` `--tumor` `--scale_mrna` | quanTIseq options |
| `--reference {TRef / BRef / both}` | EPIC reference profile |

#### Ligand–receptor
| Flag | Purpose |
|---|---|
| `--data_type {tpm / count}` | Input matrix type for `LR_cal` |
| `--id_type {symbol / ensembl / ...}` | Gene ID type for `LR_cal` |
| `--verbose` | Verbose logging |

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

### Typical end‑to‑end workflow — output file structure examples

1) **FASTQ Quality Control**
```bash
iobrpy fastq_qc \
  --path1_fastq "/path/to/fastq" \
  --path2_fastp "/path/to/fastp" \
  --num_threads 16 \
  --batch_size 4
```
```
/path/to/fastp/
  <sample>_1.fastq.gz
  <sample>_2.fastq.gz
  <sample>_fastp.html
  <sample>_fastp.json
  <sample>.task.complete
  multiqc_report/multiqc_fastp_report.html
```
2) **Prepare TPM**
```bash
# From FASTQ_QC to Salmon
iobrpy batch_salmon \
  --index "/path/to/salmon/index" \
  --path_fq "/path/to/fastp" \
  --path_out "/path/to/salmon" \
  --num_threads 16 \
  --batch_size 4
```
```
/path/to/salmon/
  <sample>/quant.sf
```
```bash
iobrpy merge_salmon \
  --project MyProj \
  --path_salmon "/path/to/salmon" \
  --num_processes 16
```
```
/path/to/salmon/
  MyProj_salmon_count.tsv.gz
  MyProj_salmon_tpm.tsv.gz
```
```bash
# From Salmon to TPM
iobrpy prepare_salmon \
  -i MyProj_salmon_tpm.tsv.gz \
  -o TPM_matrix.csv \
  --return_feature symbol \
  --remove_version
```
```
Gene        TS99       TC89       TC68       TC40       813738     1929563
5S_rRNA     0.000      0.000      0.000      0.000      0.000      0.000
5_8S_rRNA   0.000      0.000      0.000      0.000      0.000      0.000
7SK         0.000      0.000      954.687    1488.249   3691.321   5399.889
A1BG        0.479      1.717      1.844      0.382      1.676      1.126
A1BG-AS1    0.149      0.348      0.755      0.000      0.314      0.400
```
```bash
# From FASTQ_QC to STAR
iobrpy batch_star_count \
  --index "/path/to/star/index" \
  --path_fq "/path/to/fastp" \
  --path_out "/path/to/star" \
  --num_threads 16 \
  --batch_size 1
```
```
/path/to/star/
  <sample>/
  <sample>__STARgenome/
  <sample>__STARpass1/
  <sample>_STARtmp/
  <sample>_Aligned.sortedByCoord.out.bam
  <sample>_Log.final.out
  <sample>_Log.out
  <sample>_Log.progress.out
  <sample>_ReadsPerGene.out.tab
  <sample>_SJ.out.tab
  <sample>.task.complete
  .batch_star_count.done
  .merge_star_count.done
```
```bash
iobrpy merge_star_count \
  --project MyProj \
  --path "/path/to/star"
```
```
/path/to/star/
  MyProj.STAR.count.tsv.gz
```
```bash
# b) From STAR to TPM
iobrpy count2tpm \
  -i MyProj.STAR.count.tsv.gz \
  -o TPM_matrix.csv \
  --idtype ensembl \
  --org hsa \
  --remove_version
# (Optionally provide transcript effective lengths)
#   --effLength_csv efflen.csv --id id --length eff_length --gene_symbol symbol
```
```
Name       SAMPLE-2e394f45066d_20180921  SAMPLE-88dc3e3cd88e_20180921  SAMPLE-b80d019c9afa_20180921  SAMPLE-586259880b46_20180926  SAMPLE-e95813c8875d_20180921  SAMPLE-7bd449ae436b_20180921
5S_rRNA    5.326                         2.314                         2.377                         3.439                         6.993                         3.630
5_8S_rRNA  0.000                         0.000                         0.000                         0.000                         0.000                         0.000
7SK        8.006                         13.969                        11.398                        5.504                         8.510                         6.418
A1BG       3.876                         2.576                         2.874                         2.533                         2.034                         2.828
A1BG-AS1   5.512                         4.440                         7.725                         4.610                         6.292                         5.336

```

3) **(Optional) Mouse to Human symbol mapping**
```bash
# Matrix mode: rows are mouse gene symbols, columns are samples
iobrpy mouse2human_eset \
  -i mouse_matrix.tsv \
  -o human_matrix.tsv \
  --is_matrix \
  --verbose
```
```bash
# Table mode: input has a symbol column (e.g., SYMBOL), will de-duplicate then map
iobrpy mouse2human_eset \
  -i mouse_table.csv \
  -o human_matrix.csv \
  --column_of_symbol SYMBOL \
  --verbose
```
```
Gene        Sample1    Sample2    Sample3    Sample4    Sample5    Sample6
SCMH1       0.905412   0.993271   0.826294   0.535761   0.515038   0.733388
NARF        0.116423   0.944370   0.847920   0.441993   0.736983   0.467756
CD52        0.988616   0.784523   0.303614   0.886433   0.608639   0.351713
CAV2        0.063843   0.993835   0.891718   0.702293   0.703912   0.248690
HOXB6       0.716829   0.555838   0.638682   0.971783   0.868208   0.802464

```

4) **(Optional) Annotate / de‑duplicate**
```bash
iobrpy anno_eset \
  -i TPM_matrix.csv \
  -o TPM_anno.csv \
  --annotation anno_grch38 \
  --symbol symbol \
  --probe id \
  --method mean \
  --remove_version
```
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
```
Gene        GSM1523727   GSM1523728   GSM1523729   GSM1523744   GSM1523745   GSM1523746
SH3KBP1     4.3279743    4.316195     4.3514247    4.2957463    4.2566543    4.2168822
RPL41       4.2461486    4.2468076    4.2579398    4.2955956    4.2426114    4.3464246
EEF1A1      4.2937622    4.291038     4.2621994    4.2718415    4.1992331    4.2639275
HUWE1       4.2255821    4.2111235    4.1993775    4.2192063    4.2214823    4.2046394
LOC1019288  4.2193027    4.2196698    4.2132521    4.1819267    4.2345738    4.2104611

```

5) **Signature scoring**
```bash
iobrpy calculate_sig_score \
  -i TPM_anno.csv \
  -o sig_scores.csv \
  --signature signature_collection \
  --method pca \
  --mini_gene_count 2 \
  --parallel_size 1 \
  --adjust_eset
# Accepts space‑separated or comma‑separated groups; use "all" for a full merge.
```
```
ID          CD_8_T_effector_PCA   DDR_PCA    APM_PCA    Immune_Checkpoint_PCA   CellCycle_Reg_PCA   Pan_F_TBRs_PCA
GSM1523727  -3.003007             0.112244   1.046749   -3.287490               1.226469            -3.836552
GSM1523728  0.631973              1.138303   1.999972   0.405965                1.431343            0.164805
GSM1523729  -2.568384             -1.490780  -0.940420  -2.087635               0.579742            -1.208286
GSM1523744  -0.834788             4.558424   -0.274724  -0.873015               1.400215            -2.880584
GSM1523745  -1.358852             4.754705   -2.215926  -1.086041               1.342590            -1.054318

```

6) **Immune deconvolution (choose one or many)**
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
```
```
ID          B_cells_naive_CIBERSORT  B_cells_memory_CIBERSORT  Plasma_cells_CIBERSORT  T_cells_CD8_CIBERSORT  T_cells_CD4_naive_CIBERSORT  T_cells_CD4_memory_resting_CIBERSORT
GSM1523727  0.025261644              0.00067545                0.174139691             0.060873405             0                           0.143873862
GSM1523728  0.007497053              0.022985466               0.079320853             0.052005437             0                           0.137097071
GSM1523729  0.005356156              0.010721794               0.114171733             0                       0                           0.191541779
GSM1523744  0                        0.064645073               0.089539616             0.024437887             0                           0.147821928
GSM1523745  0                        0.014678117               0.121834835             0                       0                           0.176046775
```
```bash
# quanTIseq (method: lsei / robust norms)
iobrpy quantiseq \
  -i TPM_anno.csv \
  -o quantiseq.csv \
  --signame TIL10 \
  --method lsei \
  --tumor \
  --arrays \
  --scale_mrna
```
```
ID          B_cells_quantiseq   Macrophages_M1_quantiseq   Macrophages_M2_quantiseq   Monocytes_quantiseq   Neutrophils_quantiseq   NK_cells_quantiseq
GSM1523727  0.098243385         0.050936602                0.059696474                0                      0.208837962            0.057777168
GSM1523728  0.096665146         0.079422458                0.060696168                0                      0.247916520            0.057952322
GSM1523729  0.102140568         0.044950190                0.075727597                0                      0.230014524            0.060158368
GSM1523744  0.095363945         0.072341346                0.058039861                0                      0.213903654            0.059082891
GSM1523745  0.099119729         0.066757223                0.061254450                0                      0.236191857            0.056277179
```
```bash
# EPIC
iobrpy epic \
  -i TPM_anno.csv \
  -o epic.csv \
  --reference TRef
```
```
ID          Bcells_EPIC           CAFs_EPIC           CD4_Tcells_EPIC      CD8_Tcells_EPIC      Endothelial_EPIC      Macrophages_EPIC
GSM1523727  0.029043394           0.008960087         0.145125027          0.075330211          0.087619386           0.005567638
GSM1523728  0.029268307           0.010942391         0.159158789          0.074554506          0.095359587           0.007104695
GSM1523729  0.030334561           0.010648890         0.148159994          0.074191268          0.094116333           0.006359346
GSM1523744  0.027351486           0.010870086         0.144756807          0.070363208          0.085913230           0.006341159
GSM1523745  0.027688157           0.011024014         0.148947183          0.072791879          0.092757138           0.006766186
```
```bash
# ESTIMATE
iobrpy estimate \
  -i TPM_anno.csv \
  -o estimate.csv \
  --platform affymetrix
```
```
ID          StromalSignature_estimate   ImmuneSignature_estimate   ESTIMATEScore_estimate   TumorPurity_estimate
GSM1523727  -1250.182509                267.9107094                -982.2718                0.895696565
GSM1523728  197.4176128                 1333.936386                1531.353999              0.675043839
GSM1523729  -110.7937025                821.7451865                710.951484               0.758787601
GSM1523744  -118.685488                 662.3002928                543.6148048              0.774555972
GSM1523745  323.7935623                 1015.007089                1338.800651              0.695624427
```
```bash
# MCPcounter
iobrpy mcpcounter \
  -i TPM_anno.csv \
  -o mcpcounter.csv \
  --features HUGO_symbols
```
```
ID          T_cells_MCPcounter   CD8_T_cells_MCPcounter   Cytotoxic_lymphocytes_MCPcounter   B_lineage_MCPcounter   NK_cells_MCPcounter   Monocytic_lineage_MCPcounter
GSM1523727  1.4729234            1.1096225                1.3252089                          1.7530587              1.3129832             1.9197157
GSM1523728  1.5288218            1.0466424                1.5997275                          1.8069543              1.3283454             2.2191597
GSM1523729  1.4688324            1.0731858                1.3722626                          1.8967154              1.3185674             2.0802533
GSM1523744  1.4561831            1.0241529                1.440144                           1.7485736              1.3176502             2.2423225
GSM1523745  1.5078415            1.0987011                1.4883308                          1.7068269              1.3165186             2.27452
```
```bash
# IPS
iobrpy IPS \
  -i TPM_anno.csv \
  -o IPS.csv
```
```
ID          MHC_IPS    EC_IPS     SC_IPS     CP_IPS     AZ_IPS     IPS_IPS
GSM1523727  2.252749   0.403792   -0.19162   0.219981   2.684902   9
GSM1523728  2.373568   0.608176   -0.578189  -0.234406  2.16915    7
GSM1523729  2.101158   0.479571   -0.321637  0.099342   2.358434   8
GSM1523744  2.120172   0.535005   -0.332785  0.013166   2.335558   8
GSM1523745  1.911082   0.558811   -0.479384  0.087989   2.078497   7
```
```bash
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
```
                  Plasma_B_cells_deside  Non_plasma_B_cells_deside  CD4_T_deside  CD8_T_effector_deside  CD8_T_\(GZMK_high\)_deside  Double_neg_like_T_deside
TCGA-55-8508-01A  0.138                  0.014                      0.019         0.003                  0.001                       0
TCGA-67-3771-01A  0.05                   0.005                      0.016         0.002                  0.017                       0.001
TCGA-55-A4DG-01A  0.042                  0.049                      0.014         0.001                  0.035                       0.005
TCGA-91-7771-01A  0.032                  0.014                      0.032         0.006                  0.023                       0.01
TCGA-91-6849-01A  0.07                   0.011                      0.007         0.001                  0.014                       0
```

7) **TME clustering / NMF clustering**
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
```
```
ID          cluster   B_cells_naive_CIBERSORT   B_cells_memory_CIBERSORT   Plasma_cells_CIBERSORT   T_cells_CD8_CIBERSORT   T_cells_CD4_naive_CIBERSORT
GSM1523727  TME1      -0.218307125              -0.588626398               0.824242243              1.136773711             -0.142069534
GSM1523728  TME3      -0.531705309              0.093328188                -0.892611283             1.086091448             -0.142069534
GSM1523729  TME1      -0.359692153              -0.432511044               -0.481593953             -0.685959226            -0.142069534
GSM1523744  TME3      -0.531705309              0.952517071                -0.873856851             0.370938418             -0.142069534
GSM1523745  TME2      -0.531705309              -0.798612476               -0.132728742             -0.685959226            -0.142069534
```
```bash
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
```
sample      cluster   B_cells_naive_CIBERSORT  B_cells_memory_CIBERSORT  Plasma_cells_CIBERSORT  T_cells_CD8_CIBERSORT  T_cells_CD4_naive_CIBERSORT
GSM1523727  cluster2  0.006101201              0.013615524               0.149377703             0.049747382            0
GSM1523728  cluster3  0                        0.033869265               0.076470323             0.048364124            0
GSM1523729  cluster1  0.003348733              0.018252079               0.09392446              0                      0
GSM1523744  cluster2  0                        0.059386784               0.077266743             0.028845636            0
GSM1523745  cluster3  0                        0.007379033               0.108739264             0                      0

cluster   top_1                                 top_2                         top_3                                 top_4                             top_5                                   top_6
cluster1  T_cells_CD4_memory_resting_CIBERSORT  Plasma_cells_CIBERSORT        Macrophages_M2_CIBERSORT              T_cells_gamma_delta_CIBERSORT     Mast_cells_resting_CIBERSORT            T_cells_follicular_helper_CIBERSORT
cluster2  Macrophages_M2_CIBERSORT              Macrophages_M1_CIBERSORT      T_cells_follicular_helper_CIBERSORT   Plasma_cells_CIBERSORT            T_cells_CD4_memory_activated_CIBERSORT  Neutrophils_CIBERSORT
cluster3  T_cells_CD4_memory_resting_CIBERSORT  Neutrophils_CIBERSORT         Macrophages_M0_CIBERSORT              Macrophages_M2_CIBERSORT          Plasma_cells_CIBERSORT                  Mast_cells_activated_CIBERSORT

```

8) **Ligand–receptor scoring (optional)**
```bash
iobrpy LR_cal \
  -i TPM_anno.csv \
  -o LR_score.csv \
  --data_type tpm \
  --id_type symbol \
  --cancer_type pancan \
  --verbose
```
```
ID          A2M_APP_CALR_LRPAP1_PSAP_SERPING1_LRP1   ADAM10_AXL    ADAM10_EFNA1_EPHA3   ADAM12_ITGA9   ADAM12_ITGB1_SDC4   ADAM12_SDC4
GSM1523727  1.547225629                              1.566540118   1.017616452          1.476739407     1.492157038        1.492157038
GSM1523728  1.477988945                              1.757804434   1.408624847          1.492926847     1.492926847        1.492926847
GSM1523729  1.504309415                              1.730361606   1.5367173            1.473255496     1.473255496        1.473255496
GSM1523744  1.514383163                              1.73870604    1.308314516          1.469082453     1.492761796        1.492761796
GSM1523745  1.478643424                              1.76013689    1.552305282          1.449499815     1.449499815        1.449499815

```

---

## Commands & common options

### `runall` — From FASTQ to TME
- **runall**
  - `--mode {salmon|star}` (required)
  - `--outdir <DIR>` (required): root output directory
  - `--fastq <DIR>` (required): forwarded to `fastq_qc --path1_fastq`
  - `--threads <INT>` (per-block): CPU/concurrency control set via block-level flags (e.g., fastq_qc --num_threads, batch_salmon --num_threads, batch_star_count --num_threads, merge_salmon --num_processes, cibersort --threads, calculate_sig_score --parallel_size).
  - `--batch_size <INT>` (per-block): batching size set via block-level flags (e.g., fastq_qc --batch_size, batch_salmon --batch_size, batch_star_count --batch_size).
  - `--resume`: skip steps if outputs already exist
  - `--dry_run`: print planned commands without executing

### From FASTQ through FASTQ Quality Control and Salmon/STAR to TPM
- **fastq_qc**
  - `--path1_fastq <DIR>` (required): raw FASTQ directory
  - `--path2_fastp <DIR>` (required): output directory for fastp results (`01-qc/`)
  - `--num_threads <int>` (default: `8`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`): forward read suffix
  - `--batch_size <int>` (default: `5`)
  - `--se`: single-end mode
  - `--length_required <int>` (default: `50`)
  - Notes: Writes per-sample `*_fastp.html/json`; if **multiqc** is present, also writes `01-qc/multiqc_report/multiqc_fastp_report.html`.  
  *(Implementation: automatic MultiQC invocation and output path)*

#### Salmon mode
- **batch_salmon**
  - `--index <DIR>` (required): salmon index
  - `--path_fq <DIR>` (required): directory of FASTQs (after `fastq_qc`)
  - `--path_out <DIR>` (required): output root (e.g., `02-salmon/`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`)
  - `--batch_size <int>` (default: `1`): concurrent samples (processes)
  - `--num_threads <int>` (default: `8`): threads per salmon
  - `--gtf <FILE>`: optional GTF for `-g` gene-level quant
  - Behavior: safe R1 to R2 inference; per-sample `task.complete`; progress; preflight prints salmon version & index meta keys.

- **merge_salmon**
  - `--path_salmon <DIR>` (required): root containing per-sample salmon outputs (searched recursively)
  - `--project <STR>` (required): prefix for outputs
  - `--num_processes <int>`: I/O threads (default: CPU count)
  - Output: `<project>_salmon_tpm.tsv.gz`, `<project>_salmon_count.tsv.gz` under `--path_salmon` with progress and head preview.

- **prepare_salmon**
  - `-i/--input <TSV|TSV.GZ>` (required): Salmon-combined gene TPM table
  - `-o/--output <CSV/TSV>` (required): cleaned TPM matrix (genes × samples)
  - `-r/--return_feature {ENST|ENSG|symbol}` (default: `symbol`): which identifier to keep
  - `--remove_version`: strip version suffix from gene IDs (e.g., `ENSG000001.12 to ENSG000001`)

#### STAR mode
- **batch_star_count**
  - `--index <DIR>` (required): STAR genomeDir
  - `--path_fq <DIR>` (required): directory of FASTQs (after `fastq_qc`)
  - `--path_out <DIR>` (required): outputs (e.g., `02-star/`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`)
  - `--batch_size <int>` (default: `1`)
  - `--num_threads <int>` (default: `8`)
  - Notes: generates sorted BAM and `_ReadsPerGene.out.tab` per sample and a summary of paths.

- **merge_star_count**
  - `--path <DIR>` (required): directory containing multiple `*_ReadsPerGene.out.tab`
  - `--project <STR>` (required): output prefix
  - Output: `<project>.STAR.count.tsv.gz` (gzipped TSV with gene IDs as rows and samples as columns)

- **count2tpm**
  - `-i/--input <CSV/TSV[.gz]>` (required): raw count matrix (genes × samples)
  - `-o/--output <CSV/TSV>` (required): output TPM matrix
  - `--effLength_csv <CSV>`: optional effective-length file with columns `id`, `eff_length`, `symbol`
  - `--idtype {ensembl|entrez|symbol|mgi}` (default: `ensembl`)
  - `--org {hsa|mmus}` (default: `hsa`)
  - `--id <str>` (default: `id`): ID column name in `--effLength_csv`
  - `--length <str>` (default: `eff_length`): length column
  - `--gene_symbol <str>` (default: `symbol`): gene symbol column
  - `--check_data`: check & drop missing/invalid entries before conversion
  - `--remove_version`: strip version suffix from gene IDs

### (Optional) Mouse to Human symbol mapping
- **mouse2human\_eset**
  - `-i/--input <CSV|TSV|TXT[.gz]>` (required): input expression **matrix** or **table**
  - `-o/--output <CSV|TSV|TXT[.gz]>` (required): converted matrix indexed by **human** symbols (genes × samples)
  - `--is_matrix`: treat input as a matrix (rows = **mouse** gene symbols, columns = samples); if omitted, runs in **table mode**
  - `--column_of_symbol <str>` (required in table mode): column name that contains **mouse** gene symbols
  - `--sep <,|\t>`: override input separator; if omitted, inferred by extension.
  - `--out_sep <,|\t>`: override output separator; if omitted, inferred by **output** path extension
  - `--verbose`: print shapes and basic run info

### (Optional) Annotate / de‑duplicate
- **anno_eset**
  - `-i/--input <CSV/TSV/TXT>` (required)
  - `-o/--output <CSV/TSV/TXT>` (required)
  - `--annotation {anno_hug133plus2|anno_rnaseq|anno_illumina|anno_grch38}` (required unless using external file)
  - `--annotation-file <pkl/csv/tsv/xlsx>`: external annotation (overrides built-in)
  - `--annotation-key <str>`: key to pick a table if external `.pkl` stores a dict of DataFrames
  - `--symbol <str>` (default: `symbol`): column used as gene symbol
  - `--probe  <str>` (default: `id`): column used as probe/feature ID
  - `--method {mean|sd|sum}` (default: `mean`): duplicate-ID aggregation
  - `--remove_version`: strip version suffix from gene IDs

### Signature scoring
- **calculate_sig_score**
  - `-i/--input <CSV/TSV/TXT>` (required), `-o/--output <CSV/TSV/TXT>` (required)
  - `--signature <one or more groups>` (required; space- or comma-separated; `all` uses every group)  
    Groups: `go_bp`, `go_cc`, `go_mf`, `signature_collection`, `signature_tme`, `signature_sc`, `signature_tumor`, `signature_metabolism`, `kegg`, `hallmark`, `reactome`
  - `--method {pca|zscore|ssgsea|integration}` (default: `pca`)
  - `--mini_gene_count <int>` (default: `3`)
  - `--adjust_eset`: apply extra filtering after log2 transform
  - `--parallel_size <int>` (default: `1`; threads for scoring (`PCA`/`zscore`/`ssGSEA`))

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
  *Output: suffixed with `_MCPcounter`; index label `ID`; separator inferred from extension.*

- **IPS**
  - `-i/--input <matrix>` (required), `-o/--output <file>` (required)  
  *No extra flags (the expression matrix yields IPS sub-scores and a total score).*

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
  - `--id_type <str>` (default: `ensembl`): gene ID type expected by the LR backend.Choices: `ensembl`, `entrez`, `symbol`, `mgi`.
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