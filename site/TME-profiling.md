---
title: "TME profiling"
page-layout: article
toc: true
---

## All-in-one TME profiling

### Minimal usage for `tme_profile`

- `tme_profile` runs the whole TME profiling stack **from a TPM matrix** in one command.  It wraps and orchestrates:
  - **Signature scoring** → `calculate_sig_score`
  - **Immune deconvolution (six methods)** → `cibersort`, `IPS`, `estimate`, `mcpcounter`, `quantiseq`, `epic`
  - **Ligand–receptor scoring** → `LR_cal`
  - It also **merges** the deconvolution outputs into a single table

> **Inputs**: genes × samples **TPM** matrix
> **Outputs**: standardized subfolders for signatures, TME deconvolution, and Ligand–receptor scoring  
> **Heads-up:** `tme_profile` does **not** include `deside`, `BayesPrism` and any clustering (`tme_cluster`, `nmf`).  

```bash
iobrpy tme_profile \
  -i TPM.csv \
  -o /path/to/outdir \
  --threads 1
```
```
# Expected layout

/path/to/outdir
|-- 01-signatures
|   `-- calculate_sig_score.csv
|-- 02-tme
|   |-- cibersort_results.csv
|   |-- epic_results.csv
|   |-- quantiseq_results.csv
|   |-- IPS_results.csv
|   |-- estimate_results.csv
|   |-- mcpcounter_results.csv
|   `-- deconvo_merged.csv
`-- 03-LR_cal
    `-- lr_cal.csv
```

### Signature scoring
```bash
iobrpy calculate_sig_score \
  -i TPM.csv \
  -o sig_scores.csv \
  --signature all \
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

### Immune deconvolution

#### CIBERSORT
```bash
# CIBERSORT
iobrpy cibersort \
  -i TPM.csv \
  -o cibersort.csv \
  --perm 100 \
  --QN True \
  --threads 1
```
```bash
# CIBERSORT absolute mode
iobrpy cibersort \
  -i TPM.csv \
  -o cibersort.csv \
  --perm 100 \
  --QN True \
  --absolute True \
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

#### BayesPrism
```bash
# Run with single-cell data provided by IOBRpy
iobrpy bayesprism \
  -i TPM.csv \
  -o results/bayesprism \
  --threads 8
```
```bash
# Run with your single-cell reference and labels
iobrpy bayesprism \
  -i TPM.csv \
  -o results/bayesprism \
  --sc_dat sc_dat.csv \    # single-cell count matrix (cells × genes)
  --cell_state_labels cell_state_labels.csv \    # A single-column CSV file, with cells ordered identically to --sc_dat, representing each cell’s state
  --cell_type_labels cell_type_labels.csv \    # A single-column CSV file, with cells ordered identically to --sc_dat, representing each cell’s type
  --key tumor_cell_type \    # One of the cell types in cell_type_label must be designated as tumor
  --threads 8
```
```
# The output includes theta.csv, theta_cv.csv, and Z_tumor.csv; here we show only a partial preview of the format of theta.csv.
# theta.csv: estimated cell-type proportions.
# theta_cv.csv: coefficient of variation for the proportions.
# Z_tumor.csv: tumor expression estimates (as a matrix).

ID          Malignant_cells_BayesPrism  T_lymphocytes_BayesPrism  Myeloid_cells_BayesPrism  B_lymphocytes_BayesPrism  Fibroblasts_BayesPrism  NK_cells_BayesPrism
GSM1523727  0.344942                    0.018638                  0.061527                  0.120722                  0.123203                0.008854
GSM1523728  0.304567                    0.020351                  0.077356                  0.119172                  0.145886                0.010876
GSM1523729  0.302902                    0.020212                  0.068916                  0.121607                  0.146904                0.008710
GSM1523744  0.316309                    0.019605                  0.069270                  0.125118                  0.143838                0.009714
GSM1523745  0.296863                    0.020212                  0.074517                  0.124529                  0.150886                0.008780
GSM1523746  0.307967                    0.020205                  0.066109                  0.125893                  0.137038                0.008482
```

#### quanTIseq
```bash
# quanTIseq (method: lsei / robust norms)
iobrpy quantiseq \
  -i TPM.csv \
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

#### EPIC
```bash
# EPIC
iobrpy epic \
  -i TPM.csv \
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

#### ESTIMATE
```bash
# ESTIMATE
iobrpy estimate \
  -i TPM.csv \
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

#### MCPcounter
```bash
# MCPcounter
iobrpy mcpcounter \
  -i TPM.csv \
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

#### IPS
```bash
# IPS
iobrpy IPS \
  -i TPM.csv \
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

#### DeSide
```bash
# DeSide
iobrpy deside \
  --model_dir path/to/your/DeSide_model \
  -i TPM.csv \
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

### Ligand–receptor scoring
```bash
iobrpy LR_cal \
  -i TPM.csv \
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

## TME clustering / NMF clustering

### TME_cluster

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

### NMF

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