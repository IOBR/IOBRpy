---
title: Build Salmon & STAR Index
layout: default
nav_order: 3
---

# Build Salmon & STAR Index

## Prerequisites

- Conda (or mamba) installed
- Tools:
  ```bash
  # create/activate an environment (example)
  conda create -n rnaseq -y
  conda activate rnaseq

  # install tools (choose bioconda channel)
  conda install -c bioconda -c conda-forge salmon star -y
  ```

- Choose a base directory for index:
  ```bash
  export BASE=/path/to/index/dir
  mkdir -p "$BASE"
  ```

---

## Build Salmon index

```bash
# Move to the base directory for references
cd "$BASE"

# Download GENCODE v44 annotation & transcript FASTA
# (Use -c to resume if the download is interrupted)
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz

# Decompress the downloads
gunzip -f gencode.v44.annotation.gtf.gz
gunzip -f gencode.v44.transcripts.fa.gz

# Paths
fa="$BASE/gencode.v44.transcripts.fa"

# Output directory for Salmon index
mkdir -p "$BASE/salmon"

# Build Salmon index
# -t: transcript FASTA
# -i: index output directory
# -k: k-mer size (31 is common for human)
# -p: threads
salmon index -t "$fa" -i "$BASE/salmon/gencode44" -k 31 -p 8
```

**Tip:** Adjust `-p` to match your CPU cores; `-k` can be tuned depending on read length and species.

---

## 2) Build STAR index

```bash
# Move to the base directory for references
cd "$BASE"

# Download the primary assembly and GTF (GENCODE v44)
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Decompress the downloads
gunzip -f GRCh38.primary_assembly.genome.fa.gz
gunzip -f gencode.v44.annotation.gtf.gz

# Paths
STAR_INDEX="$BASE/star"
genome_fa="$BASE/GRCh38.primary_assembly.genome.fa"
gtf="$BASE/gencode.v44.annotation.gtf"

# Create output directory for STAR index
mkdir -p "$STAR_INDEX"

# Build STAR genome index
STAR --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX" \
     --genomeFastaFiles "$genome_fa" \
     --sjdbGTFfile "$gtf" \
     --runThreadN 16 \
     --sjdbOverhang 100
```

---

## Notes

- **Release/assembly:** This guide uses **GENCODE v44** (GRCh38). If you need a different release or assembly (e.g., T2T-CHM13), update the URLs and filenames accordingly.
- **Threads:** Increase `--runThreadN`/`-p` to speed up indexing if you have more cores.
- **Storage:** STAR indexes are large (tens of GB). Ensure you have enough disk space.
- **Read length:** For STAR, set `--sjdbOverhang = read_length - 1` for best splice junction annotations.