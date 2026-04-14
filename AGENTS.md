# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project Overview

**IOBRpy** is a command-line toolkit for bulk RNA-seq immuno-oncology analysis.

Its main user journeys are:
- **FASTQ to TME** via `runall`, covering QC, quantification, TPM generation, signature scoring, six default TME methods, ligand-receptor analysis, and TRUST4 outputs.
- **TPM to TME** via `tme_profile` and the individual downstream analysis commands.
- **Immune-focused workflows** for HLA typing and TCR/BCR repertoire analysis.

Mode-specific TPM generation is important in this repository:
- **Salmon mode** uses `batch_salmon -> merge_salmon -> prepare_salmon`.
- **STAR mode** uses `batch_star_count -> merge_star_count -> count2tpm`.

## Agent Fast Path

When the task is about using IOBRpy on data, prefer the CLI entrypoints before reading source files.

Recommended first steps:
- Run `iobrpy-cli doctor --json` to verify the Python environment and discover the preferred entrypoints.
- For agent-driven path parsing, do not ask the user to choose quick versus full scanning. Use the default deep/full path directly: `iobrpy-cli map --path <path> --json`. When using MCP, call `map_path` once with its default `response_detail=summary` unless technical detail is explicitly needed.
- Run `iobrpy-cli recommend --path <path> --task "<task>" --json` to classify the input, get the most likely workflow, and inspect `important_optional_parameters` before execution.
- Run `iobrpy-cli commands --json` to inspect the supported native `iobrpy` command surface; parameter sets come from `src/iobrpy/RAG_MCP/iobrpy_required_params.json`.
- Run `iobrpy-cli <native_command> --help` when you need command semantics; this help is intentionally constrained by `src/iobrpy/RAG_MCP/iobrpy_required_params.json`.
- Prefer top-level native invocations such as `iobrpy-cli runall ...` or `iobrpy-cli tme_profile ...` instead of the older `analyze/quantify/workflow/immune` wrapper namespaces.
- When users want persistent agent setup across directories, use `iobrpy-cli agent install --client codex|claude-code|all` instead of telling them to copy skill or MCP files manually.
- Use `iobrpy-cli agent status` when users want to audit whether Codex or Claude Code has already been wired up.

Agent behavior rules:
- At the start of each new IOBRpy-agent request, render this banner once before scanning or planning, matching the CLI banner in `src/iobrpy/main.py`:
  ```text
  #########################################################
   IOBRpy: Immuno-Oncology Biological Research using Python
   If you encounter any issues, please report them at
   https://github.com/IOBR/IOBRpy/issues
  #########################################################
   Author: Haonan Huang, Dongqiang Zeng
   Email: interlaken@smu.edu.cn
  #########################################################
  ```
  Do not repeat it before every tool call or native command.
- Prefer `iobrpy-cli <native_command> ...` or `iobrpy-cli-mcp` before scanning the repository source.
- Do not substitute the R `IOBR` package when Python `iobrpy` is available.
- Treat `src/iobrpy/RAG_MCP/iobrpy_required_params.json` as the source of truth for non-`deside`/`ai` command parameters in `iobrpy-cli`.
- Before recommending, writing, or executing a native command invocation, validate the command and flags against `src/iobrpy/RAG_MCP/iobrpy_required_params.json` by using `iobrpy-cli commands --json` or `iobrpy-cli <native_command> --help`. Do not invent flags from general workflow concepts. In particular, `tme_profile` has no `--mode`; its JSON-defined parameters are `--input`, `--output`, and optional `--threads`. Reserve `--mode` for commands that declare it, such as `runall`.
- For path-driven requests, run the deep/full `map` scan directly before `recommend --json` so you can tell the user whether the data is raw, partially processed, or ready for downstream analysis. Do not ask a scan-mode question first.
- When `map --json` shows existing intermediate outputs, do not default to rerunning `runall`; summarize the current stage and ask whether the user wants to continue downstream, rerun the current stage, or rerun the full pipeline.
- Reuse `map --json` fields such as `workflow_checklist`, `scenario.summary`, `scenario.communication_goals`, `scenario.required_user_decisions`, `current_label`, `current_stage_confidence`, `completed_stage_summaries`, and `next_step_summaries` when explaining the next choices to the user.
- For path-driven dataset parsing requests, always present the grouped `workflow_checklist` first from the JSON in a table before summarizing the current stage or the rerun/continue choices, even when the user did not explicitly ask “what is already done”.
- Render the full grouped checklist in order, including unchecked items, rather than filtering it down to only completed or currently relevant entries.
- Do not use a client-native scan-mode selector or ask the user to choose between quick and full scans. The default path scan is already the deep/full scan.
- Do not promise fixed minute ranges for scan duration. Warn only when relevant that full scans on very large or network-mounted directories can take much longer than expected.
- `iobrpy-cli map` does not reuse prior scan payloads; rerun it when the user wants a fresh scan.
- On the first scan, do not invent manual `--max-depth` or `--max-entries` values. Let `iobrpy-cli map` choose adaptive initial limits unless the user explicitly asks for custom limits.
- Render that checklist in user-facing text with clear markers such as `✓` for checked items and `[ ]` for unchecked items. When the user is speaking Chinese, prefer localized fields such as `label_zh`, `text_zh`, `details_zh`, and `current_label_zh` when they are present. Use the grouped checklist entries directly instead of expanding them back into many command-level stages unless the user explicitly asks for technical detail.
- In Chinese conversations, use localized checklist fields such as `workflow_checklist_title_zh`, `status_label_zh`, `result_source_label_zh`, `evidence_paths`, and `missing_items_zh` when they are present. In English conversations, do not pull `_zh` checklist fields into English table cells or summaries.
- When `workflow_checklist_table` or `workflow_checklist_table_zh` is present, prefer that prebuilt table first and reuse it verbatim.
- Do not replace that prebuilt checklist table with a shorter custom table. Keep unchecked rows and the `Still missing` / `仍缺少` column visible.
- When `default_user_facing_checklist_format` is `markdown_table`, render `default_user_facing_checklist` exactly as provided. Do not convert it into an ASCII box table, key-value blocks, or another hand-built layout.
- For checklist rows that contain multiple evidence paths, preserve each path verbatim and keep the provided in-cell line separation. Do not collapse multiple absolute paths onto one visual line, and do not concatenate adjacent paths back into one long string.
- For the `immune_deconvolution` checklist row specifically, prefer the full evidence paths from the prebuilt checklist cell or `workflow_checklist[].evidence_display_paths`. Do not replace that `检测到` / `Detected in` cell with method-name summaries such as `CIBERSORT`, `EPIC`, or `BayesPrism`.
- If a client cannot render the provided prebuilt table and absolutely must rebuild it, use `workflow_checklist[].evidence_display_paths` as the `检测到` / `Detected in` column source. Do not rebuild that column from `deconvolution_methods.detected_labels`.
- If `workflow_checklist[].detected_column_value` or `workflow_checklist[].missing_column_value` is present, use those pre-rendered cell values directly when rebuilding the checklist table.
- If a client strips HTML line breaks inside table cells, preserve a visible separator such as `•` between adjacent evidence paths or missing items instead of concatenating them into one string.
- For the `immune_deconvolution` checklist row, prefer the prebuilt checklist table first. If you must rebuild that row, use `workflow_checklist[].missing_column_value` or `workflow_checklist[].missing_items` in English, and use `workflow_checklist[].missing_column_value_zh` or `workflow_checklist[].missing_items_zh` only in Chinese. BayesPrism is an immune deconvolution method, but it is a standalone optional branch outside the default `tme_profile` / `runall` six-method bundle; keep the language-matched standalone-optional note instead of mixing Chinese text into English output.
- For the `tpm_matrix_ready` checklist row, use the prebuilt checklist table or `workflow_checklist[].detected_column_value` / `workflow_checklist[].evidence_display_paths` directly. Do not borrow downstream signature-scoring files such as `calculate_sig_score.csv` into the TPM row.
- Do not use `deconvolution_methods.default_bundle_missing_labels` to render the checklist `Still missing` / `仍缺少` column. Those labels only describe what `tme_profile` or `runall` can backfill automatically, and intentionally exclude standalone BayesPrism.
- When `workflow_checklist_report` or `workflow_checklist_report_zh` is present, prefer that prebuilt report first. Treat it as the default user-facing checklist block and reuse it verbatim before adding any short summary.
- When `workflow_checklist_lines` or `workflow_checklist_lines_zh` are present, prefer those pre-rendered lines directly instead of re-parsing the JSON payload into another ad hoc checklist formatter.
- When those pre-rendered checklist lines are present, render them verbatim and in order, including their evidence and missing-item sub-lines. Do not regroup them under custom headings such as `数据预处理阶段`, `表达矩阵准备`, or `免疫微环境分析`.
- When citing detected directories or files from `map`, prefer the full absolute paths provided by the checklist/report instead of shortening them back to relative directory names.
- Do not truncate the checklist after the earlier sections. Keep later items such as HLA typing and TCR/BCR analysis in the rendered order as well.
- Do not replace the checklist with custom taxonomy-style sections, hand-picked bullet lists, or ad hoc tables when a prebuilt checklist table or report is present.
- Before the first `map` command on a user-provided path, do not ask a scan-mode question. Run `iobrpy-cli map --path <path> --json` or the MCP `map_path` helper directly, and let the scanner choose adaptive limits.
- Once a `map_path` or `iobrpy-cli map --json` scan has succeeded, do not call `analyze_path_task` merely to render `workflow_checklist`, current stage, evidence paths, missing items, or next-step summaries. Answer directly from the map payload. Use `analyze_path_task` only when you need richer command planning that the map payload does not already provide.
- Ignore legacy quick-scan hints in older payloads. Current agent behavior is to use the deep/full map scan directly and avoid exposing quick-planning details to the user.
- For MCP helper calls, keep the default `response_detail=summary` for normal user-facing analysis. Use `response_detail=structured` only when the user explicitly asks for technical/debug details or the map summary is insufficient for a real command-planning decision.
- Do not say "the scan result is large" and then call `analyze_path_task` to extract fields. A large or folded MCP display is a rendering issue, not a reason to call another large helper.
- Do not run the same focused deep scan command repeatedly. A successful `iobrpy-cli map --path <path> --json --focus ...` already includes one internal automatic retry when needed; only launch another deep scan if the previous one failed, the directory changed, or the user explicitly asked to rescan.
- If a manual retry becomes necessary, never lower `--max-depth` or `--max-entries` compared with the current scan. Reuse `scan_retry_command` exactly when it is present; otherwise use `scan_retry_limits` from the payload instead of guessing your own values.
- After `iobrpy-cli map --json` returns, do not search memories or read Claude/Codex temporary task-output files such as `/tmp/claude.../tasks/*.output` just to render `workflow_checklist`. Do not shell out to `python3 -c`, `jq`, or similar one-off parsers for that purpose unless the user explicitly asks for raw JSON manipulation.
- Do not ask the user to approve helper shell commands whose only purpose is extracting fields from a previous `map` or MCP JSON result. In normal summary responses, render `default_user_facing_checklist` directly and use compact `missing_downstream_analysis_suggestions` from the tool result. Use structured fields such as `workflow_checklist` or `suggested_command_details` only when `response_detail=structured` explicitly returned them. In Chinese conversations, use the corresponding `_zh` fields when present.
- For MCP `map_path` summary responses, prefer `default_user_facing_checklist` first; it may already be localized and avoids duplicating the same checklist across several fields.
- If a real action still needs user confirmation, translate the intent into a concise natural-language question rather than exposing raw commands, temp-file paths, or parser snippets. For example, ask `Should I summarize the checklist and next-step suggestions from the scan result?` in English or `是否要我读取刚才扫描结果中的检查表和下一步建议？` in Chinese instead of showing a `python3 -c` command or a Claude/Codex `tool-results` path.
- If an MCP helper response is folded or saved to a client `tool-results` file because it is too large, do not read that file just to recover checklist fields and do not run `python3 -c`, `jq`, `cat`, or similar parser commands. Fall back to the most recent visible `map_path` payload fields, or rerun a compact `map_path` scan with the correct language value if genuinely necessary.
- Treat `workflow_checklist` as the single integrated rendering source. If checklist items already carry `investigation_applied`, `result_source_label`, or evidence/missing fields, do not print a second standalone section for `IOBRpy-confirmed results`, `Agent-inferred existing results`, or `External-tool results`.
- Before concluding that a stage is missing, check `scan.truncated`, `scan.depth_limited_dir_count`, `scan_warning`, and `scan_retry_recommended`. If the scan may be incomplete, automatically rerun `iobrpy-cli map --path <path> --json` once with the higher limits from `scan_retry_command` or `scan_retry_limits` before answering. Do not ask the user for permission for this retry. After the retry, use the retried payload as the main basis for the answer and mention briefly that you expanded the scan limits automatically. If the retried payload still carries `scan_warning`, keep the uncertainty explicit instead of treating the missing stages as definitely absent.
- If `auto_scan_retry_performed` is already true in the payload, do not launch another manual retry unless the user explicitly asks you to rescan again.
- Never invent retry flags. The only supported retry options are `--max-depth` and `--max-entries`; `--max-files` is not a valid `iobrpy-cli map` option. If `scan_retry_command` is present, run that command exactly instead of reconstructing your own variant.
- If the retried `map` payload still appears to miss user-reported existing outputs, or if `agent_fallback.recommended` is true because the directory looks nonstandard or externally processed, rerun `iobrpy-cli map --path <path> --json --investigate-existing` before trusting the negative checklist entries.
- Prefer the CLI-native `existing_result_investigation` payload over ad hoc agent-only searching when that option is available.
- In that fallback mode, reuse `existing_result_investigation`, `agent_fallback.investigation_targets`, `function_detections`, `content_verified_functions`, `likely_iobrpy_functions`, `reusable_result_functions`, and `external_analysis_hints` to explain what the CLI-native investigation found.
- Keep the three provenance buckets (`IOBRpy-confirmed results`, `Agent-inferred existing results`, and `External-tool results`) for internal reasoning, but do not emit them as separate user-facing sections when `workflow_checklist` already integrates that provenance.
- Do not infer additional downstream methods or sibling outputs just because files share a directory or because `tme_profile` / `runall` was detected. Only call a specific method completed when a concrete file path or content-verified evidence supports it.
- When `map.choice_details` includes `rerun_full_pipeline`, mention that option explicitly instead of collapsing the choices into a shorter custom list.
- Treat `external_analysis_hints` such as generic HLA or external TCR/BCR result directories as external analysis hints, not as proof that the corresponding iobrpy TRUST4 or HLA workflow was run.
- When `recommend --json` returns `confirm_parameters`, confirm those with the user before execution.
- When a command has resource-sensitive parameters such as `threads`, `num_threads`, `parallel_size`, `num_processes`, `batch_size`, `-t`, or `-j`, do not blindly keep the native default or hard-code `8`. If the user did not specify the value, first check CPU cores, current load, and available memory on the same execution host, then ask a concise natural-language question such as `我看了一下当前机器资源，建议用 N 个线程继续，是否可以？`.
- Read `agent-manifest.json` when you need a compact machine-readable summary of the preferred workflow entrypoints.
- Only inspect source files when implementing code changes, debugging a native command, or when the CLI output is insufficient.

## Installation and Development

### Installation
```bash
# Install in editable mode for development
pip install -e .

# External tool dependencies (via conda/mamba)
mamba install -y -c conda-forge -c bioconda \
  fastp salmon star trust4
```

### Building
```bash
# Build package
python -m build

# Build source distribution only
python -m build --sdist
```

### Testing
```bash
# Run tests
pytest

# Run specific test file
pytest tests/test_ai_plan.py

# Run with coverage
pytest --cov=iobrpy
```

## Architecture

### Entry Point
- `src/iobrpy/main.py` - CLI dispatcher using argparse with subcommands. Each subcommand maps to a function in a workflow module.

### Module Structure

**Core workflow modules** (`src/iobrpy/workflow/`):
- `runall.py` - End-to-end orchestrator from FASTQ to TME outputs. Salmon mode runs `fastq_qc -> batch_salmon -> merge_salmon -> prepare_salmon -> log2_eset`; STAR mode runs `fastq_qc -> batch_star_count -> merge_star_count -> count2tpm -> log2_eset`; both branches then continue into signature scoring, six TME methods, deconvolution merge, `LR_cal`, and `trust4`.
- `tme_profile.py` - All-in-one downstream TME profiling from an existing TPM matrix: `calculate_sig_score`, `cibersort`, `IPS`, `estimate`, `mcpcounter`, `quantiseq`, `epic`, merged deconvolution outputs, and `LR_cal`.
- `fastq_qc.py` - FASTQ QC/trimming via fastp with MultiQC summary.
- `batch_salmon.py`, `merge_salmon.py`, `prepare_salmon.py` - Salmon quantification pipeline.
- `batch_star_count.py`, `merge_star_count.py`, `count2tpm.py` - STAR alignment + count-to-TPM pipeline.
- Deconvolution / TME scoring methods: `cibersort.py`, `quantiseq.py`, `epic.py`, `estimate.py`, `mcpcounter.py`, `IPS.py`, `deside_bootstrap.py`. Note that the main repository still contains `deside_bootstrap.py`, but the `iobrpy-cli` harness currently exposes the non-DeSide methods plus standalone `bayesprism`.
- `calculate_sig_score.py` - Pathway/signature scoring with `pca`, `zscore`, `ssgsea`, and `integration`.
- `LR_cal.py` - Ligand-receptor scoring.
- `tme_cluster.py`, `nmf.py` - Clustering methods.
- `trust4.py` - TCR/BCR repertoire analysis wrapper.
- `hla_typing.py` - Batch HLA typing workflow wrapper.

**Specialized submodules**:
- `src/iobrpy/SpecHLA/` - HLA typing via SpecHLA and HLA-read extraction. This supports both per-sample `extract_hla_read -> spechla` style workflows and batch `hla_typing` execution.
- `src/iobrpy/bayesprism/` - BayesPrism deconvolution implementation.

**AI orchestration layer** (`src/iobrpy/ai/`):
- `registry.py` - Tool registration system.
- `planner.py` - Planning logic (SimplePlanner, SmartPlanner).
- `executor.py` - Plan execution.
- `backend.py` - Tool backends (LocalPythonBackend, McpBackend).
- `cli.py` - AI command-line interface.

**Resources** (`src/iobrpy/resources/`):
- Contains `.pkl` files with reference data for deconvolution methods.
- GMT files for pathway analysis (KEGG, Reactome).
- Pre-computed annotation datasets.

### Key Design Patterns

**Flag routing in `runall` and `tme_profile`**:
These orchestrators use auto-bucketing. Unrecognized long flags are automatically routed to the appropriate sub-step based on `FLAG_BUCKETS` dictionaries. This keeps the orchestrators flexible as sub-commands evolve.

**Resume-friendly execution**:
Pipeline steps write `.done` flag files such as `.batch_salmon.done` and `.fastq_qc.done` to enable resume capability.

**Standard output layout**:
Both `runall` modes produce a consistent directory structure:
- `01-qc/` - fastp outputs.
- `02-salmon/` or `02-star/` - quantification.
- `03-tpm/` - unified TPM matrix.
- `04-signatures/` - signature scores.
- `05-tme/` - deconvolution results.
- `06-LR_cal/` - ligand-receptor scores.
- `07-TCRBCR/` - TRUST4 TCR/BCR outputs.

### Analysis Semantics

**Signature scoring methods**:
- `pca` - derives one signature score from the leading principal component of signature genes.
- `zscore` - summarizes each signature with per-sample expression averaging after preprocessing.
- `ssgsea` - computes sample-wise enrichment scores.
- `integration` - combines PCA, z-score, and ssGSEA style outputs into one integrated table.

**Dedicated TME scoring commands**:
- `estimate` - stromal score, immune score, and ESTIMATE-derived purity context.
- `IPS` - immunophenoscore-related outputs.

**Default downstream TME methods in `runall` / `tme_profile`**:
- `cibersort` - immune cell fraction deconvolution.
- `quantiseq` - TIL-oriented bulk deconvolution.
- `epic` - immune/stromal/tumor-context deconvolution with configurable references.
- `mcpcounter` - abundance-style immune and stromal marker scores.
- `estimate` and `IPS` - score-based TME outputs that are run alongside the deconvolution methods.

**Standalone but not part of the default six-method `runall` / `tme_profile` bundle**:
- `bayesprism` - Bayesian deconvolution using a single-cell reference and uncertainty-aware outputs.

## Command Patterns

### Full pipeline (runall)
```bash
# Salmon mode
iobrpy runall --mode salmon --outdir ./out --fastq ./fastq \
  --index /path/to/salmon/index --threads <resource_checked_threads> --project MyProj

# STAR mode
iobrpy runall --mode star --outdir ./out --fastq ./fastq \
  --index /path/to/star/index --threads <resource_checked_threads> --project MyProj
```

### TME profiling (from existing TPM)
```bash
iobrpy tme_profile -i TPM_matrix.csv -o ./outdir --threads <resource_checked_threads>
```

### Individual workflow steps
Each workflow step can be run independently. See `iobrpy --help` and individual subcommand help for options.

## Data Conventions

- **Expression matrix orientation**: genes x samples by default.
- **File delimiters**: auto-inferred from extension (`.csv`, `.tsv`, `.txt`).
- **Gene ID types**: Supports Ensembl, Entrez, Symbol, MGI (varies by command).
- **TPM vs counts**: Most downstream TME methods expect TPM-like input; salmon workflows reach TPM through `prepare_salmon`, while STAR workflows reach TPM through `count2tpm`.

## Resource Files

Key `.pkl` files in `resources/`:
- `anno_eset.pkl` - Annotation data for `anno_eset`.
- `calculate_data.pkl` - Signature gene sets for `calculate_sig_score`.
- `epic_TRef_BRef.pkl` - EPIC reference profiles.
- `lm22.txt` - CIBERSORT LM22 signature matrix.
- `quantiseq_data.pkl`, `mcp_data.pkl`, `estimate_data.pkl` - Reference data for respective methods.

