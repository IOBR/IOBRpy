---
title: TPM conversion
layout: default
parent: Example
nav_order: 2
---

# TPM conversion

This page shows four common entry points to a TPM matrix and the final `log2(x+1)` transform you should apply after each path.

> Quick rule of thumb  
> - **Raw counts → TPM**: use `count2tpm`.  
> - **Salmon quant → TPM**: use `prepare_salmon`.  
> - **Gene-expression tables (e.g., arrays) → gene-level matrix**: use `anno_eset` to map/aggregate to symbols.  
> - **Mouse → Human**: use `mouse2human_eset` to map symbols.
> - After any of the above, run `log2_eset`.

---

## A) From **count matrix** to **TPM**

```bash
# 1) counts → TPM
iobrpy count2tpm \
  -i MyProj.STAR.count.tsv.gz \
  -o TPM_matrix.csv \
  --idtype ensembl \
  --org hsa \
  --remove_version
# (Optional) Add effective transcript lengths if available:
#   --effLength_csv efflen.csv --id id --length eff_length --gene_symbol symbol
```

```bash
# 2) TPM → log2(x+1)
iobrpy log2_eset \
  -i TPM_matrix.csv \
  -o TPM_matrix.log2.csv
```

---

## B) From **Salmon matrix** to **TPM**

```bash
# 1) Salmon TPM (gene/transcript) → cleaned gene-level TPM
iobrpy prepare_salmon \
  -i MyProj_salmon_tpm.tsv.gz \
  -o TPM_matrix.csv \
  --return_feature symbol \
  --remove_version
```

```bash
# 2) TPM → log2(x+1)
iobrpy log2_eset \
  -i TPM_matrix.csv \
  -o TPM_matrix.log2.csv
```

---

## C) From **gene-expression matrix** to **gene-level matrix** with annotation (`anno_eset`)

Use when your input is an expression table that needs **ID mapping / de-duplication** (e.g., microarray probes → symbols, or TPM tables with mixed identifiers).  

```bash
# Map/aggregate to symbols using a built-in annotation set
iobrpy anno_eset \
  -i expression_matrix.csv \
  -o expression_anno.csv \
  --annotation anno_grch38 \
  --symbol symbol \
  --probe id \
  --method mean \
  --remove_version
# Alternative platform example:
# iobrpy anno_eset -i expression_matrix.csv -o expression_anno.csv \
#   --annotation anno_hug133plus2 --symbol symbol --probe id --method mean
```

```bash
# if your input was already TPM-like, finish with log2(x+1)
iobrpy log2_eset \
  -i expression_anno.csv \
  -o expression_anno.log2.csv
```

---

## D) **Mouse → Human** gene conversion (`mouse2human_eset`)

Two common modes:

```bash
# Matrix mode: rows = mouse gene symbols, columns = samples
iobrpy mouse2human_eset \
  -i mouse_matrix.tsv \
  -o human_matrix.tsv \
  --is_matrix \
  --verbose
```

```bash
# Table mode: has a symbol column (e.g., SYMBOL); will de-duplicate then map
iobrpy mouse2human_eset \
  -i mouse_table.csv \
  -o human_matrix.csv \
  --column_of_symbol SYMBOL \
  --verbose
```

```bash
# log2(x+1) after mapping
iobrpy log2_eset \
  -i human_matrix.tsv \
  -o human_matrix.log2.tsv
```
