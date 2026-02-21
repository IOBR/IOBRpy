---
title: "Installation"
page-layout: article
toc: true
---

::: {.panel-tabset .tabset-pills}

## Conda {.unnumbered}

> **Prerequisite (Conda):** Please install Miniconda or Anaconda first. We recommend [Miniconda](https://docs.anaconda.com/miniconda/).

```bash
# Creating a virtual environment is recommended
conda create -n iobrpy python=3.11 -y
conda activate iobrpy
```
```bash
# Install iobrpy 0.1.6 (from bioconda via conda-forge + bioconda)
# Recommended: use mamba for faster solves (if available)
mamba install -y -c conda-forge -c bioconda iobrpy=0.1.7

# If you don't have mamba, use conda instead
conda install -y -c conda-forge -c bioconda iobrpy=0.1.7
```

## PyPI {.unnumbered}

```bash
# Creating a virtual environment is recommended
conda create -n iobrpy python=3.11 -y
conda activate iobrpy
```
```bash
# Update pip
python -m pip install --upgrade pip
```
```bash
# Install iobrpy
pip install iobrpy
```
```bash
# Install fastp, salmon, STAR and MultiQC
# Recommended: use mamba for faster solves (if available)
mamba install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  trust4

# If you don't have mamba, use conda instead
conda install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  trust4
```

## Docker {.unnumbered}

> **Docker Hub website:** [Docker Hub](https://hub.docker.com/)

### Download

```bash
# Option 1: Pull the latest image from Docker Hub
docker pull hhn123123/iobrpy:latest
```
```bash
# Option 2: Offline install (from GitHub Release)
# 1) Download iobrpy.tar.gz from https://github.com/IOBR/IOBRpy/releases/tag/v1.0.0
# 2) Change to the directory where the archive is saved and load the image
cd /path/to/iobrpy.tar.gz
docker load -i iobrpy.tar.gz
```

### Usage

#### Non-DeSide

```bash
# Example: run 'iobrpy runall' with STAR mode
docker run --rm -it \
  --security-opt apparmor=unconfined \
  -e PYTHONUNBUFFERED=1 \
  -v /work:/work \
  iobrpy:latest \
  iobrpy runall \
    --mode star \
    --outdir /work/path/to/outdir \
    --fastq /work/path/to/fastq \
    --index /work/path/to/star/index \
    --project MyProj \
    --threads 8 \
    --batch_size 1
```

#### DeSide (recommended to run via PyPI or Conda)

```bash
# Use a relative ext4-backed directory for the DeSide venv
docker run --rm -it \
  --security-opt apparmor=unconfined \
  -e PYTHONUNBUFFERED=1 \
  -e IOBRPY_NO_DESIDE_BOOTSTRAP=0 \
  -e IOBRPY_DESIDE_VENV=/ext4/.iobr-venv/deside-venv-py3.9 \
  -v /work:/work \
  -v /ext4:/ext4 \
  --entrypoint /bin/sh iobrpy:latest -lc '
    apt-get update &&
    apt-get install -y --no-install-recommends g++ make &&
    exec /bin/sh /usr/local/bin/docker-entrypoint.sh \
      iobrpy deside \
        --model_dir "/work/path/to/DeSide_model" \
        --input "/work/path/to/TPM_matirx.csv" \
        --output "/work/path/to/deside_results.csv" \
        --exp_type TPM \
        --result_dir "/work/path/to/deside" \
        --scaling_by_constant \
        --transpose \
        --print_info
  '
```