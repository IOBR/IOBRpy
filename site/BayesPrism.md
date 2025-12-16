---
title: "BayesPrism tutorial"
page-layout: article
toc: true
---

## Overview
`iobrpy bayesprism` runs the Python BayesPrism implementation bundled with iobrpy. It performs bulk RNA-seq deconvolution using a default single-cell reference or your own custom matrices.

## Required inputs
- `-i/--input`: bulk expression matrix (genes × samples), CSV/TSV (optionally gzipped).
- `-o/--output`: directory to receive `theta.csv`, `theta_cv.csv`, and `Z_tumor.csv`.

## Optional inputs
- `--threads`: CPU cores for BayesPrism (default 1).
- `--sc_dat`: custom single-cell count matrix (genes × cells) if you want to override the bundled `BP_data/sc_dat.csv`.
- `--cell_state_labels`: optional labels file matching the custom single-cell matrix.
- `--cell_type_labels`: optional cell type labels file matching the custom single-cell matrix.

## Example commands
```bash
# Minimal run with bundled single-cell reference
iobrpy bayesprism \
  -i data/bulk_tpm.tsv \
  -o results/bayesprism \
  --threads 8

# Run with custom single-cell reference and labels
iobrpy bayesprism \
  -i data/bulk_tpm.tsv \
  -o results/bayesprism_custom \
  --sc_dat refs/custom_sc_counts.csv \
  --cell_state_labels refs/custom_state_labels.csv \
  --cell_type_labels refs/custom_type_labels.csv \
  --threads 8
```

## Outputs
- `theta.csv`: estimated cell-type proportions.
- `theta_cv.csv`: coefficient of variation for the proportions.
- `Z_tumor.csv`: tumor expression estimates (as a matrix).
