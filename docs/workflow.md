---
title: Typical end-to-end workflow
layout: default
nav_order: 2
---

# Typical end‑to‑end workflow — output file structure examples

## 1) FASTQ Quality Control
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

## 2) Prepare TPM
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

## 3) (Optional) Mouse to Human symbol mapping
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

## 4) (Optional) Annotate / de‑duplicate
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

## 5) (Optional) Log2 transform
```bash
iobrpy log2_eset \
  -i expr.csv \
  -o expr.log2.csv
```
```
Name      SRR35344563_GSM8516765_Normal4   SRR35344561_GSM8516763_Normal2   SRR35344562_GSM8516764_Normal3   SRR35344560_GSM8516762_Normal1
A1BG      2.229246496                      0.636390662                      2.140913236                      1.420200061
A1BG-AS1  4.206586844                      3.591817651                      0.614426747                      6.842377234
A1CF      0.128261135                      0                                0.414914625                      0.205743238
A2M       1.999453226                      0.679106252                      2.816410018                      2.898826563
```

## 6) Signature scoring
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

## 7) Immune deconvolution (choose one or many)
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

## 8) TME clustering / NMF clustering
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

## 9) Ligand–receptor scoring
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