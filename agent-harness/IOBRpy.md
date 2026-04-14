# IOBRpy - CLI Harness Software-Specific SOP

## Overview

**IOBRpy** (Immuno-Oncology Biological Research using Python) is a bulk RNA-seq immuno-oncology toolkit. The `iobrpy-cli` harness wraps the native IOBRpy CLI with:
- project management,
- grouped analysis commands,
- JSON output for agents,
- REPL support,
- and clearer entry points for TME, HLA, and immune repertoire workflows.

At a high level, the software supports four major workflow families:
- **FASTQ to TME**: QC, quantification, TPM generation, signature scoring, deconvolution, and ligand-receptor analysis.
- **TPM to TME**: downstream TME analysis starting from an existing expression matrix.
- **HLA typing**: either per-sample `extract-hla-read -> spechla` or batch `hla-typing`.
- **Immune repertoire analysis**: TRUST4-based TCR/BCR reconstruction.

## Domain Analysis

### Scientific Domain
- **Field**: Bioinformatics / cancer immunology / tumor microenvironment analysis.
- **Primary data type**: bulk RNA-seq.
- **Primary goals**:
  - quantify expression from FASTQ,
  - derive TPM-like matrices,
  - score pathways and immune signatures,
  - estimate TME composition,
  - compute ligand-receptor interaction scores,
  - type HLA alleles from RNA-seq-compatible inputs,
  - reconstruct TCR/BCR repertoires.

### Canonical Workflows

#### 1. FASTQ to TME via `quantify runall`

This is the main end-to-end workflow when starting from raw sequencing reads.

- **Common first stage**: `fastq_qc` performs trimming and QC on FASTQ inputs.
- **Salmon branch**:
  - `batch_salmon`
  - `merge_salmon`
  - `prepare_salmon`
- **STAR branch**:
  - `batch_star_count`
  - `merge_star_count`
  - `count2tpm`
- **Downstream TME stage**:
  - signature scoring,
  - selected deconvolution methods,
  - optional ligand-receptor analysis,
  - optional TRUST4 execution.

This distinction between `prepare_salmon` and `count2tpm` is important:
- `prepare_salmon` is the Salmon-specific TPM preparation path.
- `count2tpm` is the STAR-count-specific TPM conversion path.

#### 2. TPM to TME via `analyze tme-profile`

This is the main downstream workflow when the user already has a TPM matrix.

The default bundle includes:
- `calculate_sig_score`,
- `cibersort`,
- `IPS`,
- `estimate`,
- `mcpcounter`,
- `quantiseq`,
- `epic`,
- `lr-analysis`.

It is the best high-level entry point for users who want "TME profiling from expression matrix" without manually running each method.

#### 3. HLA typing workflows

There are two distinct HLA stories that should be represented clearly:

- **Per-sample HLA workflow**:
  - `immune extract-hla-read`
  - then `immune spechla`
  - This is the explicit "extract reads from BAM/CRAM, convert to FASTQ, then type with SpecHLA" path.

- **Batch HLA workflow**:
  - `immune hla-typing`
  - This is the directory-level wrapper that performs the full extraction + SpecHLA flow across multiple samples and produces merged outputs.

#### 4. Immune repertoire workflow

- `immune trust4` performs TCR/BCR repertoire reconstruction.
- The underlying analyzer supports both single-input and directory/batch style execution depending on whether the input is a BAM file, BAM directory, or FASTQ directory.

### Data Models

#### Expression Matrix
- **Format**: CSV/TSV/TXT depending on command.
- **Orientation**: genes x samples by default.
- **Common gene ID types**: Ensembl, Entrez, gene symbol, MGI.
- **Normalization context**:
  - Salmon-derived data are prepared through `prepare_salmon`.
  - STAR-derived counts are converted through `count2tpm`.
  - downstream commands generally expect TPM-like matrices.

#### Signature and Score Outputs
- **Signature scores**: sample x signature matrix.
- **`calculate_sig_score` methods**:
  - `pca`: one score per signature from the leading principal component.
  - `zscore`: mean-style signature summary after preprocessing.
  - `ssgsea`: sample-wise enrichment scoring.
  - `integration`: combined PCA, z-score, and ssGSEA style outputs.
- **Dedicated score commands**:
  - `estimate`: stromal / immune / ESTIMATE-related outputs.
  - `ips`: immunophenoscore-related outputs.

#### Deconvolution Outputs
- **Typical shape**: sample x cell-state or cell-type matrix.
- **Default TME-profile / runall methods**:
  - `cibersort`
  - `ips`
  - `estimate`
  - `mcpcounter`
  - `quantiseq`
  - `epic`
- **Standalone advanced method**:
  - `bayesprism`, which uses a single-cell reference and is exposed as a separate analysis path rather than part of the default six-method bundle.

#### HLA Outputs
- `extract-hla-read` produces extracted HLA-focused FASTQ outputs.
- `spechla` produces per-sample HLA typing results.
- `hla-typing` produces batch-oriented results, including merged HLA summaries.

#### TRUST4 Outputs
- repertoire reports,
- immune summary tables,
- and derived immune indices where available.

## Command Groups

### 1. Project and Session Management

The harness adds a stateful layer on top of IOBRpy:
- `project create`
- `project list`
- `project status`
- `project show`
- `project delete`
- `repl`

These commands are not biological analyses themselves. They organize analysis state, paths, and outputs for agent-friendly usage.

### 2. TME Analysis (`analyze`)

#### `analyze tme-profile`
Best entry point for downstream TME analysis from a TPM matrix.

#### `analyze signature-score`
Exposes `calculate_sig_score` directly.

**Supported scoring methods**:
- `pca`
- `zscore`
- `ssgsea`
- `integration`

**Supported signature groups** include pathway and immune collections such as:
- `go_bp`
- `go_cc`
- `go_mf`
- `signature_collection`
- `signature_tme`
- `signature_sc`
- `signature_tumor`
- `signature_metabolism`
- `kegg`
- `hallmark`
- `reactome`
- `all`

#### Dedicated TME scoring commands
- `analyze estimate`
- `analyze ips`

These are score-centric TME outputs rather than classical cell-fraction deconvolution.

#### `analyze deconvolute`
General entry point for method-specific deconvolution.

**Supported methods**:
- `cibersort`
- `ips`
- `estimate`
- `mcpcounter`
- `quantiseq`
- `epic`
- `bayesprism`

#### Deconvolution method notes
- **CIBERSORT**: immune cell fraction deconvolution.
- **quanTIseq**: TIL-oriented deconvolution for bulk RNA-seq.
- **EPIC**: immune/stromal/tumor-context deconvolution with configurable references.
- **MCPcounter**: marker-based abundance scores for immune and stromal populations.
- **BayesPrism**: Bayesian deconvolution with a single-cell reference and uncertainty-aware outputs.

#### Other analysis commands
- `analyze lr-analysis` - ligand-receptor analysis.
- `analyze combine` - combine multiple deconvolution outputs into a summary table.

### 3. HLA and Immune Repertoire (`immune`)

#### `immune extract-hla-read`
Starts from BAM/CRAM. Extracts HLA-related reads and converts them to FASTQ. This is the prerequisite step for SpecHLA when the starting point is aligned reads.

#### `immune spechla`
Runs per-sample SpecHLA on FASTQ inputs.

#### `immune hla-typing`
Runs the batch HLA workflow over a BAM directory. Conceptually this is the "extract + SpecHLA" pipeline at directory scale.

#### `immune trust4`
Runs TCR/BCR repertoire reconstruction.

### 4. Quantification (`quantify`)

#### `quantify qc`
FASTQ QC and trimming only.

#### `quantify runall`
Main end-to-end FASTQ workflow.

**Salmon mode**:
- QC -> Salmon quantification -> Salmon merge -> `prepare_salmon` -> downstream TME analysis.

**STAR mode**:
- QC -> STAR counting -> STAR merge -> `count2tpm` -> downstream TME analysis.

**Downstream stage**:
- signature scoring,
- selected deconvolution methods,
- optional ligand-receptor analysis,
- optional TRUST4.

#### Supporting quantification utilities
- `quantify clean`
- `quantify validate`
- `quantify scan-fastq`
- `version --external`

### 5. Single-Step Workflow Utilities (`workflow`)

These commands expose native steps independently for users who want manual control.

**Expression matrix preparation**:
- `workflow prepare-salmon`
- `workflow count2tpm`
- `workflow anno-eset`
- `workflow log2-eset`
- `workflow remove-version`
- `workflow mouse2human`
- `workflow mouse2human-eset`

**Quantification steps**:
- `workflow fastq-qc`
- `workflow batch-salmon`
- `workflow merge-salmon`
- `workflow batch-star-count`
- `workflow merge-star-count`

**Downstream analysis steps**:
- `workflow calculate-sig-score`
- `workflow lr-cal`
- `workflow tme-cluster`
- `workflow nmf-cluster`

### 6. AI Assistant

The underlying repository also contains an AI orchestration layer in `src/iobrpy/ai/` for tool registration, plan generation, and execution. This is separate from the `iobrpy-cli` harness but part of the broader codebase.

## State Model

### Project State
A project represents a managed analysis workspace with:
- input data files,
- processing parameters,
- output directories,
- execution status,
- and reusable context for REPL / JSON workflows.

### Session State
A session tracks:
- current working directory or active project,
- last-used parameters,
- result locations,
- and command history.

### Output Formats

#### Expression Matrix
- CSV/TSV with genes as rows and samples as columns.
- Identifier column contains gene IDs.
- Values are TPM-like expression or derived transformed values depending on the step.

#### Deconvolution Results
- CSV/TSV with samples as rows and cell states, cell types, or scores as columns.
- Values are proportions or abundance-like scores depending on method.

#### Signature Scores
- CSV with samples as rows and signatures as columns.
- Values depend on the scoring method (`pca`, `zscore`, `ssgsea`, `integration`).

#### LR Scores
- CSV with ligand-receptor scores derived from bulk expression input.

## Flag Auto-Routing

Both native `runall` / `tme_profile` style workflows and the harness expose auto-routing patterns:
- unrecognized long flags can be forwarded to the appropriate sub-step,
- legacy concurrency flags can be normalized,
- and top-level workflow commands can remain stable even as subcommands evolve.

This behavior is important when reading the code because orchestration layers may accept options that are not explicitly declared on every wrapper function.

## External Dependencies

### Required System Tools
- **fastp** - FASTQ QC/trimming.
- **salmon** - transcript quantification.
- **star** - genome alignment / gene counting.
- **trust4** - TCR/BCR analysis.

### Additional HLA-related Runtime Dependencies
For `extract-hla-read`, the workflow expects tools such as:
- `samtools`
- `bamutil`
- `libdeflate`
- `htslib`

### Python Dependencies
- `numpy`, `pandas`, `scipy`, `scikit-learn`
- `joblib`, `gseapy`, `tqdm`
- plus method-specific scientific dependencies used by the native IOBRpy modules

## Resource Files

Reference data in `src/iobrpy/resources/` include:
- `calculate_data.pkl` - signature collections used by `calculate_sig_score`.
- `lm22.txt` - CIBERSORT reference matrix.
- `epic_TRef_BRef.pkl` - EPIC references.
- `quantiseq_data.pkl`, `mcp_data.pkl`, `estimate_data.pkl` - method-specific reference data.
- pathway GMTs such as KEGG and Reactome collections.

