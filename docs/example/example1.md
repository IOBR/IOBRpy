---
title: From downloading data to TME
layout: default
parent: Example
nav_order: 1
---

# **From downloading data to TME**

## Part 1 — Data download

### Prepare the SRR list
- Retrieve the SRR accessions for **PRJNA1161405** from NCBI SRA: <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1161405>.
- Save the accessions (one per line) into `PRJNA1161405.txt` and upload it to: `path/to/PRJNA1161405/`.

### High-speed download with `prefetch`
> Requires **SRA Toolkit** installed and on your `PATH`.

```bash
# (Optional) load/activate your environment
# module load sra-tools            # or: conda activate sra-tools

cd path/to/PRJNA1161405
# Download all SRR runs listed in PRJNA1161405.txt into the current directory
prefetch -O ./ --option-file PRJNA1161405.txt
```

### Convert `.sra` to FASTQ with `fasterq-dump`
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

### Multi-thread compression with `pigz`
Compress all `.fastq` files in the folder using 8 threads.

```bash
cd path/to/PRJNA1161405
for file in SRR*.fastq; do
  if [ -f "$file" ]; then
    pigz "$file" -p 8
  fi
done
```

### (Optional) Direct downloads from ENA FTP with `curl`
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

## Part 2 — From FASTQ to TME - runall

### Salmon mode
```bash
iobrpy runall \
  --mode salmon \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/salmon/index" \
  --project SRR
```

### STAR mode
```bash
iobrpy runall \
  --mode star \
  --outdir "/path/to/outdir" \
  --fastq "/path/to/fastq" \
  --threads 8 \
  --batch_size 1 \
  --index "/path/to/star/index" \
  --project SRR
```