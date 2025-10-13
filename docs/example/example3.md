---
title: From TPM to TME
layout: default
parent: Example
nav_order: 3
---

# **From TPM to TME**

This page takes a **TPM** matrix and runs downstream **TME** analyses: signature scoring, immune deconvolution (multiple methods), clustering, and ligand–receptor scoring. 

---

## Inputs

- **TPM matrix**: `TPM_matrix.csv`
- (Optional) **log2 transform**: if desired, apply:
  
```bash
iobrpy log2_eset \
  -i TPM_matrix.csv \
  -o TPM_matrix.log2.csv
```

---

## 1) Signature scoring

Compute pathway/signature scores from TPM.

```bash
iobrpy calculate_sig_score \
  -i TPM_matrix.csv \
  -o sig_scores.csv \
  --signature all \
  --method integration \
  --mini_gene_count 2 \
  --parallel_size 1 \
  --adjust_eset
```

---

## 2) Immune deconvolution

Choose one or several methods below; each writes one result file.

### CIBERSORT
```bash
iobrpy cibersort \
  -i TPM_matrix.csv \
  -o cibersort.csv \
  --perm 100 \
  --QN True \
  --absolute False \
  --abs_method sig.score \
  --threads 1
```

### quanTIseq
```bash
iobrpy quantiseq \
  -i TPM_matrix.csv \
  -o quantiseq.csv \
  --signame TIL10 \
  --method lsei \
  --tumor \
  --arrays \
  --scale_mrna
```

### EPIC
```bash
iobrpy epic \
  -i TPM_matrix.csv \
  -o epic.csv \
  --reference TRef
```

### ESTIMATE
```bash
iobrpy estimate \
  -i TPM_matrix.csv \
  -o estimate.csv \
  --platform affymetrix
```

### MCPcounter
```bash
iobrpy mcpcounter \
  -i TPM_matrix.csv \
  -o mcpcounter.csv \
  --features HUGO_symbols
```

### IPS
```bash
iobrpy IPS \
  -i TPM_matrix.csv \
  -o IPS.csv
```

### DeSide
```bash
iobrpy deside \
  --model_dir path/to/your/DeSide_model \
  -i TPM_matrix.csv \
  -o deside.csv \
  --result_dir path/to/your/plot/folder \
  --exp_type TPM \
  --scaling_by_constant \
  --transpose \
  --print_info
```

---

## 3) TME clustering

You can cluster samples by cell fractions or signature scores.

### k-means with KL index auto-k (recommended)
```bash
iobrpy tme_cluster \
  -i cibersort.csv \
  -o tme_cluster.csv \
  --features 1:22 \
  --id "ID" \
  --min_nc 2 \
  --max_nc 5 \
  --print_result \
  --scale
```

### NMF clustering (auto-k, excluding k=2)
```bash
iobrpy nmf \
  -i cibersort.csv \
  -o path/to/your/result/folder \
  --kmax 10 \
  --features 1:22 \
  --skip_k_2
```

---

## 4) Ligand–receptor scoring

Compute bulk ligand–receptor interaction scores from TPM:

```bash
iobrpy LR_cal \
  -i TPM_matrix.csv \
  -o LR_score.csv \
  --data_type "tpm" \
  --id_type "symbol" \
  --cancer_type pancan \
  --verbose
```