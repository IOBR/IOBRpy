---
title: "From FASTQ to TPM"
page-layout: article
toc: true
---

## FASTQ Quality Control
```bash
iobrpy fastq_qc \
  --path1_fastq "/path/to/fastq" \
  --path2_fastp "/path/to/fastp" \
  --num_threads 8 \
  --batch_size 1
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

---

## Prepare TPM

### batch_salmon
```bash
# From FASTQ_QC to Salmon
iobrpy batch_salmon \
  --index "/path/to/salmon/index" \
  --path_fq "/path/to/fastp" \
  --path_out "/path/to/salmon" \
  --num_threads 8 \
  --batch_size 1
```
```
/path/to/salmon/
  <sample>/quant.sf
```

### merge_salmon
```bash
iobrpy merge_salmon \
  --project MyProj \
  --path_salmon "/path/to/salmon" \
  --num_processes 8
```
```
/path/to/salmon/
  MyProj_salmon_count.tsv.gz
  MyProj_salmon_tpm.tsv.gz
```

### prepare_salmon
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

### batch_star_count
```bash
# From FASTQ_QC to STAR
iobrpy batch_star_count \
  --index "/path/to/star/index" \
  --path_fq "/path/to/fastp" \
  --path_out "/path/to/star" \
  --num_threads 8 \
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

### merge_star_count
```bash
iobrpy merge_star_count \
  --project MyProj \
  --path "/path/to/star"
```
```
/path/to/star/
  MyProj.STAR.count.tsv.gz
```

### count2tpm
```bash
# From STAR to TPM
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

---

## (Optional) Mouse to Human symbol mapping
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

---

## (Optional) Annotate / de‑duplicate
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

---

## (Optional) Log2 transform
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