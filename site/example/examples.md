---
title: "Examples"
toc: true
number-sections: true
---

## From downloading data to TME

### Data download

#### Prepare the SRR list
- Retrieve the SRR accessions for **PRJNA1161405** from NCBI SRA: <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1161405>.
- Save the accessions (one per line) into `PRJNA1161405.txt` and upload it to: `path/to/PRJNA1161405/`.

#### High-speed download with `prefetch`
> Requires **SRA Toolkit** installed and on your `PATH`.

```bash
# (Optional) load/activate your environment
# module load sra-tools            # or: conda activate sra-tools

cd path/to/PRJNA1161405
# Download all SRR runs listed in PRJNA1161405.txt into the current directory
prefetch -O ./ --option-file PRJNA1161405.txt
```

#### Convert `.sra` to FASTQ with `fasterq-dump`
This loop finds each run directory produced by `prefetch` and converts the `.sra` file to paired FASTQ files.

```bash
folder="path/to/PRJNA1161405/"
cd path/to/PRJNA1161405

for dir in "${folder}"SRR*; do
  if [[ -d "${dir}" ]]; then
    dir_name="$(basename "${dir}")"
    input_file="${dir}/${dir_name}.sra"
    # -3: skip technical reads, -p: show progress, -e 64: threads, -O . : output to current dir
    fasterq-dump -3 "${input_file}" -p -e 64 -O .
  fi
done
```

#### Multi-thread compression with `pigz`
Compress all `.fastq` files in the folder using 8 threads.

```bash
cd path/to/PRJNA1161405
for file in SRR*.fastq; do
  if [ -f "$file" ]; then
    pigz "$file" -p 8
  fi
done
```

#### (Optional) Direct downloads from ENA FTP with `curl`
If you prefer pulling FASTQ files directly from ENA:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Normal samples
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/063/SRR35344563/SRR35344563_1.fastq.gz -o SRR35344563_GSM8516765_Normal4_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/063/SRR35344563/SRR35344563_2.fastq.gz -o SRR35344563_GSM8516765_Normal4_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/061/SRR35344561/SRR35344561_1.fastq.gz -o SRR35344561_GSM8516763_Normal2_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/061/SRR35344561/SRR35344561_2.fastq.gz -o SRR35344561_GSM8516763_Normal2_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/060/SRR35344560/SRR35344560_1.fastq.gz -o SRR35344560_GSM8516762_Normal1_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/060/SRR35344560/SRR35344560_2.fastq.gz -o SRR35344560_GSM8516762_Normal1_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/062/SRR35344562/SRR35344562_1.fastq.gz -o SRR35344562_GSM8516764_Normal3_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/062/SRR35344562/SRR35344562_2.fastq.gz -o SRR35344562_GSM8516764_Normal3_Homo_sapiens_RNA-Seq_2.fastq.gz

# HCC samples
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/068/SRR35344568/SRR35344568_1.fastq.gz -o SRR35344568_GSM8516770_HCC3_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/068/SRR35344568/SRR35344568_2.fastq.gz -o SRR35344568_GSM8516770_HCC3_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/069/SRR35344569/SRR35344569_1.fastq.gz -o SRR35344569_GSM8516771_HCC4_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/069/SRR35344569/SRR35344569_2.fastq.gz -o SRR35344569_GSM8516771_HCC4_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/070/SRR35344570/SRR35344570_1.fastq.gz -o SRR35344570_GSM8516772_HCC5_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/070/SRR35344570/SRR35344570_2.fastq.gz -o SRR35344570_GSM8516772_HCC5_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/071/SRR35344571/SRR35344571_1.fastq.gz -o SRR35344571_GSM8516773_HCC6_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/071/SRR35344571/SRR35344571_2.fastq.gz -o SRR35344571_GSM8516773_HCC6_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/072/SRR35344572/SRR35344572_1.fastq.gz -o SRR35344572_GSM8516774_HCC7_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/072/SRR35344572/SRR35344572_2.fastq.gz -o SRR35344572_GSM8516774_HCC7_Homo_sapiens_RNA-Seq_2.fastq.gz

# CLD samples
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/073/SRR35344573/SRR35344573_1.fastq.gz -o SRR35344573_GSM8516775_CLD1_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/073/SRR35344573/SRR35344573_2.fastq.gz -o SRR35344573_GSM8516775_CLD1_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/074/SRR35344574/SRR35344574_1.fastq.gz -o SRR35344574_GSM8516776_CLD2_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/074/SRR35344574/SRR35344574_2.fastq.gz -o SRR35344574_GSM8516776_CLD2_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/075/SRR35344575/SRR35344575_1.fastq.gz -o SRR35344575_GSM8516777_CLD3_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/075/SRR35344575/SRR35344575_2.fastq.gz -o SRR35344575_GSM8516777_CLD3_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/076/SRR35344576/SRR35344576_1.fastq.gz -o SRR35344576_GSM8516778_CLD4_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/076/SRR35344576/SRR35344576_2.fastq.gz -o SRR35344576_GSM8516778_CLD4_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/077/SRR35344577/SRR35344577_1.fastq.gz -o SRR35344577_GSM8516779_CLD5_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR353/077/SRR35344577/SRR35344577_2.fastq.gz -o SRR35344577_GSM8516779_CLD5_Homo_sapiens_RNA-Seq_2.fastq.gz
```

---

### From FASTQ to TME - runall

#### Salmon mode
```bash
iobrpy runall \
  --mode salmon \
  --outdir "/path/to/salmon_dir" \
  --fastq "path/to/PRJNA1161405" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/salmon/index" \
  --project SRR
```

#### STAR mode
```bash
iobrpy runall \
  --mode star \
  --outdir "/path/to/star_dir" \
  --fastq "path/to/PRJNA1161405" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/star/index" \
  --project SRR
```

## TPM conversion

This page shows four common entry points to a TPM matrix and the final `log2(x+1)` transform you should apply after each path.

> Quick rule of thumb  
> - **Raw counts → TPM**: use `count2tpm`.  
> - **Salmon quant → TPM**: use `prepare_salmon`.  
> - **Gene-expression tables (e.g., arrays) → gene-level matrix**: use `anno_eset` to map/aggregate to symbols.  
> - **Mouse → Human**: use `mouse2human_eset` to map symbols.
> - After any of the above, run `log2_eset`.

---

### From **count matrix** to **TPM**

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

### From **Salmon matrix** to **TPM**

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

### From **gene-expression matrix** to **gene-level matrix** with annotation (`anno_eset`)

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

### **Mouse → Human** gene conversion (`mouse2human_eset`)

Two common modes:

```bash
# Matrix mode: rows = mouse gene symbols, columns = samples
iobrpy mouse2human_eset \
  -i mouse_matrix.csv \
  -o human_matrix.csv \
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
  -i human_matrix.csv \
  -o human_matrix.log2.csv
```

## From TPM to TME

This page takes a **TPM** matrix and runs downstream **TME** analyses: signature scoring, immune deconvolution (multiple methods), clustering, and ligand–receptor scoring. 

---

### Inputs

- **TPM matrix**: `TPM_matrix.csv`
- (Optional) **log2 transform**: if desired, apply:
  
```bash
iobrpy log2_eset \
  -i TPM_matrix.csv \
  -o TPM_matrix.log2.csv
```

---

### All-in-one TME profiling - `tme_profile`

- `tme_profile` wraps the following functions into one command:
  - **Signature scoring** → `calculate_sig_score`
  - **Immune deconvolution (six methods)** → `cibersort`, `IPS`, `estimate`, `mcpcounter`, `quantiseq`, `epic`
  - **Ligand–receptor scoring** → `LR_cal`
  - It also **merges** the deconvolution outputs into a single table

> **Not included:** `deside`, `bayesprism` and any clustering (`tme_cluster`, `nmf`).  
> **Tip:** You can either run each function **step-by-step** (see the sections below for individual commands and options), **or** use `tme_profile` to execute the full chain in one go.

#### Minimal usage
```bash
iobrpy tme_profile \
  -i TPM_matrix.csv \
  -o out/tme \
  --threads 1
```

---

### Signature scoring

Per-sample signature scores. Columns correspond to the selected signature set and method (`integration`, `pca`, `zscore`, or `ssgsea`).

```bash
iobrpy calculate_sig_score \
  -i TPM_matrix.csv \
  -o sig_scores.csv \
  --signature all \
  --method integration \
  --mini_gene_count 2 \
  --parallel_size 1 \
  --adjust_eset
# Accepts space‑separated or comma‑separated groups; use "all" for a full merge.
```

---

### Immune deconvolution

Choose one or several methods below; each writes one result file.

#### CIBERSORT
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

#### BayesPrism
```bash
iobrpy bayesprism \
  -i TPM_matrix.csv \
  -o results/bayesprism \
  --threads 8
```

#### quanTIseq
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

#### EPIC
```bash
iobrpy epic \
  -i TPM_matrix.csv \
  -o epic.csv \
  --reference TRef
```

#### ESTIMATE
```bash
iobrpy estimate \
  -i TPM_matrix.csv \
  -o estimate.csv \
  --platform affymetrix
```

#### MCPcounter
```bash
iobrpy mcpcounter \
  -i TPM_matrix.csv \
  -o mcpcounter.csv \
  --features HUGO_symbols
```

#### IPS
```bash
iobrpy IPS \
  -i TPM_matrix.csv \
  -o IPS.csv
```

#### DeSide
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

### TME clustering

You can cluster samples by cell fractions or signature scores.

#### tme_cluster: k-means with KL index auto-k (recommended)
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

#### NMF clustering (auto-k, excluding k=2)
```bash
iobrpy nmf \
  -i cibersort.csv \
  -o path/to/your/result/folder \
  --kmax 10 \
  --features 1:22 \
  --skip_k_2
```

---

### Ligand–receptor scoring

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

---

## HLA typing

Run `iobrpy hla_typing` to perform batch HLA typing on the samples.

```bash
iobrpy hla_typing \
  -b /path/to/star_dir/02-star \    # Generated by the STAR mode of iobrpy runall
  -r hg38 \
  -o /data/hla_typing_results \
  -j 8
```

---

## TCR-BCR

Run `iobrpy trust4` to analyze TCR and BCR sequences

```bash
# Batch over a BAM directory
iobrpy trust4 \
  -b /path/to/star_dir/02-star \    # Generated by the STAR mode of iobrpy runall
  -o /data/trust4_results \
  -t 8
```
```bash
# Batch over paired FASTQs
iobrpy trust4 \
  --fqdir /path/to/salmon_dir/01-qc \    # If you use the Salmon mode of iobrpy runall
  -o /data/trust4_results \
  -t 8
```