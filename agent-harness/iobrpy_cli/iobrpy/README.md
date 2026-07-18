# iobrpy-cli

A stateful CLI harness for **IOBRpy** (Immuno-Oncology Biological Research using Python) with REPL support and JSON output mode for agent consumption.

## MCP Server

This package also provides a real MCP stdio server:

```bash
iobrpy-cli-mcp
```

The MCP server exposes the native `iobrpy` subcommands as tools by deriving the
tool catalog from the real parser in `iobrpy.main`. This keeps command names and
parameter semantics aligned with upstream `iobrpy` instead of maintaining a
second, drift-prone command surface.

Intentional exclusions:
- `ai`

All other native `iobrpy` subcommands are exposed as MCP tools, and their
agent-facing parameter sets are constrained by
`src/iobrpy/RAG_MCP/iobrpy_required_params.json`.

## Persistent Agent Setup

Use the built-in installer when you want persistent agent integrations instead of repeating `iobrpy-cli` in each prompt:

```bash
# Install the bundled `/iobrpy` and `/iobrpy-result` skills/plugins plus Codex MCP
iobrpy-cli agent install --client codex

# Install both Claude Code commands/memories plus Claude Code MCP
iobrpy-cli agent install --client claude-code

# Configure every supported client in one pass
iobrpy-cli agent install --client all

# Inspect the current setup without writing files
iobrpy-cli agent status
```

These commands default to a short human-readable summary. Add `--json` if you want the full structured payload for automation.

What the installer does:
- `codex`: installs the `iobrpy` and `iobrpy-result` skills, their local plugins and marketplace entries, and an MCP server entry in `~/.codex/config.toml`.
- `claude-code`: installs `/iobrpy` and `/iobrpy-result`, refreshes their managed guidance, injects import blocks into `~/.claude/CLAUDE.md`, and runs `claude mcp add <name> --scope user -- <python> -m iobrpy_cli.iobrpy.mcp_server`.

Use `--dry-run` to preview changes and `--force` when you want to overwrite the packaged Codex skill.

## Overview

IOBRpy is a command-line toolkit for bulk RNA-seq immuno-oncology analysis. This harness is now **native-first**: for every non-`ai` IOBRpy command, the agent-facing parameter surface in `iobrpy-cli` is constrained by `src/iobrpy/RAG_MCP/iobrpy_required_params.json`.

Preferred pattern:

```bash
iobrpy-cli <native_command> [native_args...]
```

Examples:

```bash
iobrpy-cli runall --mode salmon --fastq ./fastq --outdir ./out --index /path/to/index
iobrpy-cli tme_profile -i TPM_matrix.csv -o ./outdir --threads 1
iobrpy-cli trust4 -o ./trust4_out -t 8 --fqdir ./fastq
iobrpy-cli runall --help
```

Per-command help is semantic-aware: for commands such as `runall` and `tme_profile`, `--help` is driven by `src/iobrpy/RAG_MCP/iobrpy_required_params.json`, so the parameter set and parameter descriptions stay close to the JSON used by the RAG MCP layer.

Legacy namespaces such as `analyze`, `quantify`, `workflow`, `immune`, and `hla_tcr` are still accepted as compatibility shims, but they should not be the first choice for agents.

## Capability Map

### TME Analysis

- `tme_profile` runs the common downstream bulk-TME workflow from an existing TPM matrix: signature scoring, six default TME methods, and ligand-receptor analysis.
- `calculate_sig_score` exposes four methods: `pca`, `zscore`, `ssgsea`, and `integration`.
- `estimate` and `IPS` expose dedicated score-based outputs: ESTIMATE stromal/immune context and IPS immunophenoscore outputs.
- `cibersort`, `mcpcounter`, `quantiseq`, `epic`, and `bayesprism` stay available as standalone native analysis commands.
- Compatibility aliases such as `analyze tme-profile` and `analyze signature-score` still forward to the matching native command, but the top-level native form is preferred.

### Deconvolution Methods

- `cibersort`: immune cell fraction deconvolution for bulk expression matrices.
- `quantiseq`: TIL-focused bulk RNA-seq deconvolution with optional array normalization, tumor-gene removal, and mRNA scaling.
- `epic`: immune, stromal, and tumor-context deconvolution using configurable EPIC references.
- `mcpcounter`: marker-based abundance scores for immune and stromal populations.
- `bayesprism`: Bayesian deconvolution using a single-cell reference, exposed as a standalone analysis method rather than part of the default `tme-profile` bundle.

### HLA and Immune Repertoire Workflows

- `extract_hla_read` extracts HLA-related reads from BAM/CRAM and converts them to FASTQ. This is the preparation step for SpecHLA when starting from aligned reads.
- `spechla` runs per-sample SpecHLA directly from extracted FASTQ pairs.
- `hla_typing` runs the batch HLA workflow over a BAM directory and produces merged HLA outputs. Conceptually this is the "extract + SpecHLA" pipeline at directory scale.
- `trust4` runs TCR/BCR repertoire reconstruction.

### End-to-End Quantification

- `runall` starts from raw FASTQ and performs FASTQ QC first.
- In **Salmon mode**, the pipeline flows through `batch_salmon -> merge_salmon -> prepare_salmon`, then enters downstream TME analysis.
- In **STAR mode**, the pipeline flows through `batch_star_count -> merge_star_count -> count2tpm`, then enters downstream TME analysis.
- The downstream TME stage covers signature scoring, configurable deconvolution methods, optional ligand-receptor analysis, and optional TRUST4 execution.
- `prepare_salmon` and `count2tpm` are also exposed separately when users want to control TPM generation step by step.

## Installation

### From Source

```bash
cd /path/to/iobrpy-cli/agent-harness
pip install -e .
```

### Dependencies

The package requires:
- Python 3.11+
- click >= 8.0
- pandas >= 1.5
- prompt-toolkit >= 3.0
- iobrpy >= 0.2.0

### Optional Dependencies

For Excel export support:
```bash
pip install openpyxl
```

Or install with all extras:
```bash
pip install -e ".[full]"
```

## Quick Start

### Project-Based Workflow

```bash
# Create a new project
iobrpy-cli project create my_tme_analysis \
  --root ./iobrpy_projects \
  --description "TME analysis for cancer samples" \
  --input-type tpm

# List all projects
iobrpy-cli project list --root ./iobrpy_projects

# Show project status
iobrpy-cli project status my_tme_analysis --root ./iobrpy_projects

# Show project details
iobrpy-cli project show my_tme_analysis --root ./iobrpy_projects
```

### Running TME Analysis

```bash
# Run complete TME profile from TPM matrix
iobrpy-cli tme_profile \
  -i input_data/tpm_matrix.csv \
  -o ./results \
  --threads 4
```

### Running End-to-End Pipeline

```bash
# Run complete pipeline from FASTQ to TME results (Salmon mode)
iobrpy-cli runall \
  --mode salmon \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/salmon_index \
  --threads 8

# Run complete pipeline from FASTQ to TME results (STAR mode)
iobrpy-cli runall \
  --mode star \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/star_index \
  --threads 8

# Run pipeline with resume support
iobrpy-cli runall \
  --mode salmon \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/salmon_index \
  --resume

# Dry run to preview pipeline steps
iobrpy-cli runall \
  --mode salmon \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/salmon_index \
  --dry_run

# Pipeline with custom scoring method forwarded to calculate_sig_score
iobrpy-cli runall \
  --mode salmon \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/salmon_index \
  --signature hallmark kegg \
  --method integration

# Clean pipeline intermediate files
iobrpy-cli quantify clean ./pipeline_output

# Validate pipeline output
iobrpy-cli quantify validate ./pipeline_output --mode salmon

# Scan FASTQ directory
iobrpy-cli quantify scan-fastq /path/to/fastq_files --suffix1 _R1 --suffix2 _R2

# Show version information
iobrpy-cli version --external

# Combine deconvolution results
iobrpy-cli analyze combine ./results/deconvolution --output combined_results.csv

# Pipeline with specific signature groups
iobrpy-cli runall \
  --mode salmon \
  --fastq /path/to/fastq_files \
  --outdir ./pipeline_output \
  --index /path/to/salmon_index \
  --signature hallmark kegg
```

```bash
# Run complete TME profile from TPM matrix
iobrpy-cli tme_profile \
  -i input_data/tpm_matrix.csv \
  -o ./results \
  --threads 4

# Calculate signature scores only
iobrpy-cli calculate_sig_score \
  -i input_data/tpm_matrix.csv \
  -o signature_scores.csv \
  --signature all \
  --method integration

# Run a specific deconvolution method
iobrpy-cli cibersort \
  -i input_data/tpm_matrix.csv \
  -o deconv_results.csv \
  --threads 4

# Run ligand-receptor analysis
iobrpy-cli LR_cal \
  -i input_data/tpm_matrix.csv \
  -o lr_scores.csv \
  --cancer_type pancan
```

### Export Results

```bash
# Export as JSON
iobrpy-cli export json results/tme/cibersort_results.csv \
  --output results/cibersort.json

# Create summary report
iobrpy-cli export summary ./results \
  --output summary.json
```

### REPL Mode

Start an interactive REPL session:

```bash
iobrpy-cli repl --root ./iobrpy_projects
```

In REPL mode, you can:
- Use standard IOBRpy commands directly
- Set active project: `use my_project`
- List projects: `projects`
- Show status: `status`
- View history: `history`
- Exit: `exit`

## JSON Output Mode

For programmatic/agent consumption, use the `--json` flag:

```bash
# List projects in JSON format
iobrpy-cli project list --root ./iobrpy_projects --json

# Get project status as JSON
iobrpy-cli project status my_project --root ./iobrpy_projects --json

# Agent-friendly workflow discovery in JSON
iobrpy-cli recommend --path input.csv --task "analyze tumor microenvironment" --json
```

JSON responses always include:
- `success`: Boolean indicating if operation succeeded
- `message`: Human-readable message
- `data`: Operation-specific result data

For `recommend --json`, the `recommendation` object includes `confirm_parameters` so agents can confirm the JSON-defined parameters before execution.

## Agent Helper Commands

The harness also exposes a few shell-first helper commands so agents can
discover and use IOBRpy quickly without scanning the repository source tree
first:

```bash
iobrpy-cli doctor --json
iobrpy-cli commands --json
iobrpy-cli recommend --path /path/to/input --task "parse tumor microenvironment" --json
```

These commands are meant to be the preferred fast path for agent workflows. When `recommend --json` returns `confirm_parameters`, agents should review those before launching long workflows.
They complement the MCP server instead of replacing it.

## Command Reference

### Project Commands

| Command | Description |
|----------|-------------|
| `project create <name>` | Create a new project |
| `project list` | List all projects |
| `project status <name>` | Show project status |
| `project show <name>` | Show project details |
| `project delete <name>` | Delete a project |

### Native Analysis Commands

| Command | Description |
|----------|-------------|
| `tme_profile -i <input> -o <dir>` | Run complete TME profile analysis |
| `calculate_sig_score -i <input> -o <file>` | Calculate signature scores |
| `cibersort -i <input> -o <file>` | Run CIBERSORT deconvolution |
| `quantiseq -i <input> -o <file>` | Run quanTIseq deconvolution |
| `epic -i <input> -o <file>` | Run EPIC deconvolution |
| `bayesprism -i <input> -o <dir>` | Run BayesPrism deconvolution |
| `LR_cal -i <input> -o <file>` | Run ligand-receptor analysis |
| `trust4 -o <dir> -t <threads> (--bam <bam>|--fqdir <dir>|-1 <r1> -2 <r2>)` | Run TRUST4 TCR/BCR repertoire reconstruction |
| `spechla <native_args>` | Run SpecHLA RNA-seq exon-level HLA typing |
| `extract_hla_read <native_args>` | Extract HLA reads from BAM/CRAM files for HLA typing |
| `hla_typing -b <bam_dir> -r <hg19|hg38> -o <outdir>` | Run batch HLA typing workflow |

### Export Commands

| Command | Description |
|----------|-------------|
| `export json <file>` | Export file as JSON |
| `export summary <dir>` | Create summary report |

### Utility Commands

| Command | Description |
|----------|-------------|
| `info <path>` | Show file/directory information |
| `ls <path>` | List directory contents |
| `version [--external]` | Show version information for IOBRpy and external tools |

### Quantification Commands

| Command | Description |
|----------|-------------|
| `runall --mode <mode> --fastq <dir> --outdir <dir> --index <dir>` | Run complete end-to-end pipeline from FASTQ to TME results |
| `fastq_qc --path1_fastq <dir> --path2_fastp <dir>` | Run FASTQ quality control |
| `quantify clean <output-dir>` | Clean up intermediate files from pipeline run |
| `quantify validate <output-dir>` | Validate pipeline output directory for completeness |
| `quantify scan-fastq <fastq-dir>` | Scan FASTQ directory and report structure |

### Workflow Commands

| Command | Description |
|----------|-------------|
| `prepare_salmon -i <input> -o <file>` | Prepare Salmon data matrix |
| `count2tpm -i <input> -o <file>` | Convert count matrix to TPM |
| `anno_eset -i <input> -o <file>` | Annotate expression set |
| `calculate_sig_score -i <input> -o <file>` | Calculate signature scores |
| `batch_salmon <native_args>` | Batch Salmon quantification |
| `merge_salmon <native_args>` | Merge Salmon results |
| `batch_star_count <native_args>` | Batch STAR quantification |
| `merge_star_count <native_args>` | Merge STAR results |
| `fastq_qc <native_args>` | Run FASTQ quality control |
| `log2_eset -i <input> -o <file>` | Apply log2 transformation |
| `tme_cluster -i <input> -o <file>` | Run TME clustering |
| `nmf -i <input> -o <dir>` | Run NMF clustering |
| `mouse2human_eset -i <input> -o <file>` | Convert mouse to human symbols |
| `LR_cal -i <input> -o <file>` | Calculate ligand-receptor scores |

## Project Structure

A project has the following structure:

```
my_project/
鈹溾攢鈹€ project.json      # Project configuration
鈹溾攢鈹€ input/           # Input data files
鈹溾攢鈹€ output/          # Analysis results
鈹斺攢鈹€ logs/            # Log files
```

## Deconvolution Methods

The following deconvolution methods are supported:

| Method | Description | Cell Types |
|--------|-------------|-------------|
| CIBERSORT | 22 immune cell types | B cells, T cells, NK cells, Macrophages, etc. |
| quanTIseq | 10 TIL cell types | B cells, CD4+ T, CD8+ T, etc. |
| EPIC | Endothelial, fibroblast, immune | B cells, T cells, Endothelial, etc. |
| MCPcounter | Immune/stromal scores | T cells, B cells, NK cells, etc. |
| ESTIMATE | Stromal/immune scores | Stromal, Immune, Tumor Purity |
| IPS | Immunophenoscore | CTL, Effector cells, etc. |

## Signature Groups

Available signature groups for scoring:

- `go_bp` - Gene Ontology Biological Process
- `go_cc` - Gene Ontology Cellular Component
- `go_mf` - Gene Ontology Molecular Function
- `signature_collection` - General signatures
- `signature_tme` - Tumor microenvironment signatures
- `signature_sc` - Single-cell signatures
- `signature_tumor` - Tumor-specific signatures
- `signature_metabolism` - Metabolism signatures
- `kegg` - KEGG pathways
- `hallmark` - MSigDB hallmark gene sets
- `reactome` - Reactome pathways
- `all` - All signature groups

## Python API

You can also use the harness as a Python library:

```python
from pathlib import Path
from iobrpy_cli.iobrpy.core import (
    ProjectManager,
    TMEAnalyzer,
    QuantificationWorkflow,
    QuantificationMode,
    DeconvolutionMethod,
)

# Create a project
manager = ProjectManager(Path("./projects"))
project = manager.create_project("my_analysis", input_type="tpm")

# Run TME analysis
analyzer = TMEAnalyzer()
results = analyzer.run_tme_profile(
    input_matrix=Path("data/tpm_matrix.csv"),
    output_dir=Path("results"),
    threads=4,
)

# Check results
print(f"Signatures: {results['signatures'].success}")
print(f"Deconvolution: {len(results['deconvolution'])} methods")
print(f"LR analysis: {results['lr'].success}")

# Run complete end-to-end pipeline from FASTQ
workflow = QuantificationWorkflow()
pipeline_results = workflow.runall(
    fastq_dir=Path("fastq_files"),
    output_dir=Path("pipeline_output"),
    mode=QuantificationMode.SALMON,
    index=Path("salmon_index"),
    threads=8,
)
```

## Examples

### Example 1: Complete Workflow

```bash
# 1. Create project
iobrpy-cli project create cancer_study \
  --root ./iobrpy_projects \
  --input-type tpm

# 2. Run TME profile
iobrpy-cli tme_profile \
  -i data/cancer_tpm.csv \
  -o ./iobrpy_projects/cancer_study/output/tme_profile \
  --threads 8

# 3. Export summary
iobrpy-cli export summary \
  ./iobrpy_projects/cancer_study/output/tme_profile \
  --output summary_report.json
```

### Example 2: Analysis Pipeline

```bash
# Calculate signatures
iobrpy-cli calculate_sig_score \
  -i data/tpm_matrix.csv \
  -o signatures.csv \
  --signature hallmark kegg \
  --method ssgsea

# Run CIBERSORT
iobrpy-cli cibersort \
  -i data/tpm_matrix.csv \
  -o deconv/cibersort.csv \
  --threads 8

# Run LR analysis
iobrpy-cli LR_cal \
  -i data/tpm_matrix.csv \
  -o lr_scores.csv \
  --cancer_type pancan
```

### Example 4: Individual Workflow Steps

```bash
# Convert counts to TPM
iobrpy-cli count2tpm \
  -i counts_matrix.csv \
  -o tpm_matrix.csv \
  --idtype ensembl --org hsa

# Annotate expression matrix
iobrpy-cli anno_eset \
  -i tpm_matrix.csv \
  -o annotated_matrix.csv \
  --annotation anno_grch38

# Calculate signature scores standalone
iobrpy-cli calculate_sig_score \
  -i annotated_matrix.csv \
  -o signature_scores.csv \
  --signature hallmark kegg go_bp \
  --method integration

# Run TME clustering
iobrpy-cli tme_cluster \
  -i deconvolution_results.csv \
  -o clustering_results.csv \
  --features 1:8 --scale

# Run NMF clustering
iobrpy-cli nmf \
  -i signature_scores.csv \
  -o nmf_results \
  --kmin 2 --kmax 6 --normalize

# Convert mouse to human symbols
iobrpy-cli mouse2human_eset \
  -i mouse_matrix.csv \
  -o human_matrix.csv \
  --is-matrix

# Calculate ligand-receptor scores
iobrpy-cli LR_cal \
  -i human_matrix.csv \
  -o lr_scores.csv \
  --cancer_type brca
```

### Example 6: HLA Typing Workflow

```bash
# Extract HLA reads from BAM file
iobrpy-cli extract_hla_read \
  --sample NA12878 \
  --bam alignment/NA12878.bam \
  --ref hg38 \
  --outdir hla_extracted

# Extract HLA reads with hg19 reference
iobrpy-cli extract_hla_read \
  --sample SAMPLE001 \
  --bam alignment/sample.bam \
  --ref hg19 \
  --outdir hla_extracted \
  --no-auto-install
```

### Example 5: Batch Processing Workflow

```bash
# 1. Run FASTQ QC on raw data
iobrpy-cli fastq_qc \
  --path1_fastq raw_fastq/ \
  --path2_fastp cleaned_fastq/ \
  --num_threads 8 \
  --batch_size 4

# 2. Batch Salmon quantification
iobrpy-cli batch_salmon \
  --index salmon_index/ \
  --path_fq cleaned_fastq/ \
  --path_out salmon_quant/ \
  --suffix1 "_clean_1.fastq.gz" \
  --batch_size 4 \
  --num_threads 8

# 3. Merge Salmon results
iobrpy-cli merge_salmon \
  --path_salmon salmon_quant/ \
  --project my_study

# 4. Prepare final matrix
iobrpy-cli prepare_salmon \
  -i salmon_quant/my_study_TPM_matrix.csv \
  -o final_tpm.csv \
  -r symbol
```

### Example 3: JSON Output for Scripts

```bash
# Create project and get JSON response
PROJECT_ID=$(iobrpy-cli project create my_project --json | jq -r '.data.path')

# Discover the right workflow
iobrpy-cli recommend --path data.csv --task "analyze tumor microenvironment" --json | jq

# Export summary
iobrpy-cli export summary ./out --json | jq '.data.sections'
```

## Development

### Running Tests

```bash
# Run unit tests
pytest tests/test_core.py -v

# Run E2E tests (requires package installation)
pytest tests/test_full_e2e.py -v

# Run all tests with coverage
pytest --cov=iobrpy_cli.iobrpy -v
```

### Code Structure

```
iobrpy_cli/iobrpy/
鈹溾攢鈹€ __init__.py
鈹溾攢鈹€ iobrpy_cli.py      # Main CLI entry point
鈹溾攢鈹€ core/              # Core modules
鈹?  鈹溾攢鈹€ __init__.py
鈹?  鈹溾攢鈹€ project.py      # Project management
鈹?  鈹溾攢鈹€ session.py      # Session management
鈹?  鈹溾攢鈹€ export.py       # Export functionality
鈹?  鈹溾攢鈹€ tme.py         # TME analysis
鈹?  鈹斺攢鈹€ quantification.py  # Quantification workflows including runall pipeline
鈹溾攢鈹€ utils/             # Utilities
鈹?  鈹溾攢鈹€ __init__.py
鈹?  鈹溾攢鈹€ io.py          # I/O utilities
鈹?  鈹斺攢鈹€ json_output.py # JSON output formatting
鈹斺攢鈹€ tests/             # Tests
    鈹溾攢鈹€ __init__.py
    鈹溾攢鈹€ test_core.py    # Unit tests
    鈹溾攢鈹€ test_full_e2e.py # E2E tests
    鈹斺攢鈹€ TEST.md        # Test documentation
```

## License

MIT License - See LICENSE file for details.

## Support

For issues related to:
- This CLI harness: Report in the issue tracker
- IOBRpy functionality: Report at https://github.com/IOBR/IOBRpy/issues

