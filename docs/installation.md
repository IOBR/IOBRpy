---
title: Installation
layout: default
nav_order: 2
---

# **Installation**

```bash
# Creating a virtual environment is recommended
conda create -n iobrpy python=3.9 -y
conda activate iobrpy

# Update pip
python -m pip install --upgrade pip

# Install iobrpy
pip install iobrpy

#Install fastp, salmon, STAR and MultiQC
# Recommended: use mamba for faster solves (if available)
# Channels order matters: conda-forge first, then bioconda
mamba install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  multiqc

# If you don't have mamba, use conda instead
# (slower dependency solving; otherwise equivalent)
conda install -y -c conda-forge -c bioconda \
  fastp \
  salmon \
  star \
  multiqc

# (Optional) Verify tools are available
fastp --version
salmon --version
STAR --version
multiqc --version
```