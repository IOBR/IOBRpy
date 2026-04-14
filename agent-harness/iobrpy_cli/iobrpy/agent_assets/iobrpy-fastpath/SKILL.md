---
name: iobrpy
description: Route bulk RNA-seq, TME, 微环境, immune infiltration, FASTQ, TPM, HLA, and TCR/BCR analysis requests to iobrpy-cli or registered iobrpy MCP tools before scanning source code. Use when the task is to analyze data with IOBRpy, choose the right workflow for an input path, or run native workflows such as runall, tme_profile, hla_typing, or trust4 without depending on the current working directory.
---

# Iobrpy Fastpath

## Startup Banner

At the start of each new IOBRpy-agent request, render this banner once before scanning or planning:

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

Do not repeat this banner before every tool call or native command.

## Workflow

1. Prefer registered iobrpy MCP tools when they are available. For a user-provided path, start with `map_path`, not `analyze_path_task`. Pass `language: "en"` for English user requests and `language: "zh"` for Chinese user requests so compact checklist and confirmation fields match the conversation language. Use `recommend_workflow` only when you still need a narrower command recommendation that the current helper payload did not already provide. Use `list_native_commands` when the native command surface is unclear. Otherwise use `iobrpy-cli`.
2. Start with `iobrpy-cli doctor --json` only when the environment or command availability is genuinely uncertain. Do not call it by default for an ordinary path scan.
3. When the user provides a path or a dataset task, do not ask them to choose quick versus full scanning. Run `map_path` directly, or use `iobrpy-cli map --path <path> --json`; the default scan path is the deep/full path for agent use.
3a. Do not use `request_user_input`, clickable choice UI, or a text prompt for scan mode. Legacy quick mode is deprecated for agent use and should not be exposed to the user.
4. Use the map output to present `workflow_checklist` first as a table for every path-driven dataset parsing request, even when the user did not explicitly ask what has already been completed.
5. Render the full grouped `workflow_checklist` in order, including unchecked items, instead of filtering it down to only completed or currently relevant entries.
5a. Do not ask a scan-mode question before the first `map` command on a user-provided path.
5b. Do not show a quick/full scan option table. The default path scan is already the deep/full scan.
5bb. Do not promise fixed minute ranges for scan duration. Warn only when relevant that full scans on very large or network-mounted directories can take much longer than expected.
5c. `iobrpy-cli map` does not reuse prior scan payloads; rerun it when the user wants a fresh scan.
5d. On the first scan, do not invent manual `--max-depth` or `--max-entries` values. Let `iobrpy-cli map` choose adaptive initial limits unless the user explicitly asks for custom limits.
6. Render that checklist directly from `workflow_checklist` using user-facing markers such as `✓` for checked items and `[ ]` for unchecked items. When the user is speaking Chinese, prefer localized fields such as `label_zh`, `text_zh`, `details_zh`, and `current_label_zh` when they are present. Use the grouped checklist entries directly instead of expanding them back into many command-level stages.
7. When `workflow_checklist_title_zh`, `status_label_zh`, `result_source_label_zh`, `evidence_paths`, or `missing_items_zh` are present, use those localized checklist fields directly. In Chinese conversations, keep checklist section titles in Chinese rather than mixing in English titles.
7.0a. When `workflow_checklist_table` or `workflow_checklist_table_zh` is present, prefer that prebuilt table first and reuse it verbatim.
7.0aa. Do not replace that prebuilt checklist table with a shorter custom table. Keep unchecked rows and the `Still missing` / `仍缺少` column visible.
7.0aaa. When `default_user_facing_checklist_format` is `markdown_table`, render `default_user_facing_checklist` exactly as provided. Do not convert it into an ASCII box table, key-value blocks, or another hand-built layout.
7.0aab. For checklist rows that contain multiple evidence paths, preserve each path verbatim and keep the provided in-cell line separation. Do not shorten absolute paths back to basenames or summarize them as method labels such as `CIBERSORT` or `EPIC`.
7.0aac. For the `immune_deconvolution` row specifically, prefer the full evidence paths from the prebuilt checklist cell. Do not replace that cell with method-name summaries, and do not squeeze multiple file paths onto one visual line.
7.0aad. If `workflow_checklist[].detected_column_value` or `workflow_checklist[].missing_column_value` is present, use those pre-rendered cell values directly when rebuilding the checklist table.
7.0aae. If a client strips HTML line breaks inside table cells, preserve a visible separator such as `•` between adjacent evidence paths or missing items instead of concatenating them into one string.
7.0ab. For the `immune_deconvolution` checklist row, treat `workflow_checklist[].missing_items_zh` or the prebuilt checklist table as authoritative. BayesPrism is an immune deconvolution method, but it is a standalone optional branch outside the default `tme_profile` / `runall` six-method bundle; when it is missing, keep it in the `仍缺少` column as `BayesPrism（独立可选，不由 tme_profile/runall 自动运行）`.
7.0ac. Do not use `deconvolution_methods.default_bundle_missing_labels` to render the checklist `Still missing` / `仍缺少` column. Those labels only describe what `tme_profile` or `runall` can backfill automatically, and intentionally exclude standalone BayesPrism.
7.0af. In English conversations, do not pull `missing_items_zh`, `missing_column_value_zh`, or other `_zh` checklist fields into English table cells. For the `immune_deconvolution` row in English, use `missing_column_value` or `missing_items` so BayesPrism keeps the English standalone-optional note.
7.0. When `workflow_checklist_report` or `workflow_checklist_report_zh` is present, prefer that prebuilt report first and reuse it verbatim before adding any short summary.
7a. When `workflow_checklist_lines` or `workflow_checklist_lines_zh` are present, prefer those pre-rendered checklist lines directly instead of re-parsing the JSON payload into another formatter.
7aa. When those pre-rendered checklist lines are present, render them verbatim and in order, including their evidence and missing-item sub-lines. Do not regroup them under custom headings such as `数据预处理阶段`, `表达矩阵准备`, or `免疫微环境分析`.
7aab. When citing detected directories or files from `map`, prefer the full absolute paths provided by the checklist/report instead of shortening them back to relative directory names.
7aaa. Do not truncate the checklist after the earlier sections. Keep later items such as HLA typing and TCR/BCR analysis in the rendered order as well.
7ab. Do not replace the checklist with custom taxonomy-style sections, hand-picked bullet lists, or ad hoc tables when a prebuilt checklist table or report is present.
7ac. Before the first `map` command on a user-provided path, do not ask a scan-mode question. Run `iobrpy-cli map --path <path> --json` or MCP `map_path` directly, and let the scanner choose adaptive limits.
7ac-pre. Once a `map_path` or `iobrpy-cli map --json` scan has succeeded, do not call `analyze_path_task` merely to render `workflow_checklist`, current stage, evidence paths, missing items, or next-step summaries. Answer directly from the map payload. Use `analyze_path_task` only when you need richer command planning that the map payload does not already provide.
7ac0. When `agent_rendering_hints.path_specific_recommended_deep_map_command_uses_unbounded_entries` is true, keep its `--max-entries 0` segment intact. That means the focused deep scan should not stop at the usual entry cap.
7aca. Ignore legacy quick-scan hints in older payloads. Current agent behavior is to use the deep/full map scan directly and avoid exposing quick-planning details to the user.
7acb. If `agent_rendering_hints.do_not_ask_user_to_choose_scan_mode` is present, treat it as authoritative.
7acb1. For MCP helper calls, keep the default `response_detail=summary` for normal user-facing analysis. Use `response_detail=structured` only when the user explicitly asks for technical/debug details or the map summary is insufficient for a real command-planning decision.
7acb2. Do not say "the scan result is large" and then call `analyze_path_task` to extract fields. A large or folded MCP display is a rendering issue, not a reason to call another large helper.
7acc. Do not run the same focused deep scan command repeatedly. A successful `iobrpy-cli map --path <path> --json --focus ...` already includes one internal automatic retry when needed; only launch another deep scan if the previous one failed, the directory changed, or the user explicitly asked to rescan.
7ad. Do not expose quick-scan fallback logic to users. If a scan looks incomplete, follow `scan_retry_command` or `scan_retry_limits` automatically.
7ae. If a previous quick scan payload exists from an older client, prefer a fresh full/deep map payload before making missing-stage claims.
7af. If a manual retry becomes necessary, never lower `--max-depth` or `--max-entries` compared with the current scan. Reuse `scan_retry_command` exactly when it is present; otherwise use `scan_retry_limits` from the payload instead of guessing your own values.
7b. After `iobrpy-cli map --json` returns, do not search memories or read temporary task-output files such as `/tmp/claude.../tasks/*.output` just to render `workflow_checklist`. Do not shell out to `python3 -c`, `jq`, or similar one-off parsers for that purpose unless the user explicitly asks for raw JSON manipulation.
7ba. Do not ask the user to approve helper shell commands whose only purpose is extracting fields from a previous `map` or MCP JSON result. In normal summary responses, render `default_user_facing_checklist` directly and use compact `missing_downstream_analysis_suggestions` from the tool result. Use structured fields such as `workflow_checklist` or `suggested_command_details` only when `response_detail=structured` explicitly returned them. In Chinese conversations, use the corresponding `_zh` fields when present.
7ba1. For MCP `map_path` summary responses, prefer `default_user_facing_checklist` first; it may already be localized and avoids duplicating the same checklist across several fields.
7bb. If a real action still needs user confirmation, translate the intent into a concise natural-language question rather than exposing raw commands, temp-file paths, or parser snippets. For example, ask `Should I summarize the checklist and next-step suggestions from the scan result?` in English or `是否要我读取刚才扫描结果中的检查表和下一步建议？` in Chinese instead of showing a `python3 -c` command or a Claude/Codex `tool-results` path.
7bc. If an MCP helper response is folded or saved to a client `tool-results` file because it is too large, do not read that file just to recover checklist fields and do not run `python3 -c`, `jq`, `cat`, or similar parser commands. Fall back to the most recent visible `map_path` payload fields, or rerun a compact `map_path` scan with the correct `language` value if genuinely necessary.
8. Treat `workflow_checklist` as the single integrated rendering source. If checklist items already carry `investigation_applied`, `result_source_label`, or evidence/missing fields, do not print a second standalone section for `IOBRpy-confirmed results`, `Agent-inferred existing results`, or `External-tool results`.
8a. Do not infer additional downstream methods or sibling outputs just because files share a directory or because `tme_profile` / `runall` was detected. Only call a specific method completed when a concrete file path or content-verified result supports it.
9. Before concluding that a stage is missing, check `scan.truncated`, `scan.depth_limited_dir_count`, `scan_warning`, and `scan_retry_recommended`. If the scan may be incomplete, automatically rerun `iobrpy-cli map --path <path> --json` once with the higher limits from `scan_retry_command` or `scan_retry_limits` before answering. Do not ask the user for permission for that retry. Use the retried payload as the main basis for the answer and mention briefly that you expanded the scan limits automatically. If the retried payload still contains `scan_warning`, keep the uncertainty explicit instead of treating the missing stages as definitely absent.
9a. If `auto_scan_retry_performed` is already true in the payload, do not launch another manual retry unless the user explicitly asks you to rescan again.
9b. Never invent retry flags. The only supported retry options are `--max-depth` and `--max-entries`; `--max-files` is not a valid `iobrpy-cli map` option. If `scan_retry_command` is present, run that command exactly instead of reconstructing your own variant.
10. If the retried `map` payload still appears to miss user-reported existing outputs, or if `agent_fallback.recommended` is true because the directory looks nonstandard or externally processed, rerun `iobrpy-cli map --path <path> --json --investigate-existing` before trusting the negative checklist entries.
11. Prefer the CLI-native `existing_result_investigation` payload over ad hoc agent-only searching when that option is available.
12. In that fallback mode, use `existing_result_investigation`, `agent_fallback.investigation_targets`, `function_detections`, `content_verified_functions`, `likely_iobrpy_functions`, `reusable_result_functions`, and `external_analysis_hints` to explain what the CLI-native investigation found.
13. Report fallback findings inside the integrated checklist rather than as a second checklist-shaped section whenever the checklist already carries source/evidence fields.
14. After the checklist, explain the current stage, confidence, which outputs already exist, and whether it is better to continue downstream, rerun the current stage, or rerun the full pipeline. Reuse `workflow_checklist`, `current_stage_confidence`, `scenario.summary`, `scenario.communication_goals`, `scenario.required_user_decisions`, and `roadmap.position_summary` when available.
14aa. When `analyze_path_task` returns `current_state_reasoning_surface`, `ranked_next_step_options`, and `next_step_reasoning_surface`, use them as the main explanation surface. Treat those as structured inputs for your own summary rather than canned text to repeat verbatim.
14a. Do not hard-limit yourself to three suggestions. Present as many materially relevant next actions as the current state supports.
14a1. When `missing_downstream_analysis_suggestions` is present, use it as the recommendation source after the checklist. It is expected to include every missing downstream analysis, including standalone BayesPrism; do not collapse it to only the default `tme_profile` bundle.
14a2. When `suggested_command_details` contains parameter hints or placeholders, show concise natural-language prompts in the user's language, such as `path_salmon_index="替换你的 Salmon index 的路径"` for Chinese or `path_salmon_index="Replace with your Salmon index path."` for English. These hints are display guidance layered on top of `src/iobrpy/RAG_MCP/iobrpy_required_params.json`; do not invent flags outside that JSON.
14b. When you recommend a command or next action, explain why it fits the current state, what input data it needs, and what result it will produce.
14c. Assume the user may not know the workflow details. Prefer beginner-friendly explanations over terse labels such as "continue downstream" without context.
15. When `map.choice_details` includes a choice, reflect that choice back to the user explicitly instead of replacing it with a shorter custom list. If `rerun_full_pipeline` is present, mention it explicitly.
16. Treat `external_analysis_hints` such as generic HLA or external TCR/BCR result directories as external analysis hints, not as proof that iobrpy TRUST4 or iobrpy HLA workflows were run.
17. When speaking to users, prefer `workflow_checklist`, `current_label`, `next_step_summaries`, and `supplementary_stage_summaries`. Do not re-expand `workflow_checklist` into a second long list of `completed_stage_summaries` unless the user explicitly asks for more technical detail. Do not present raw function names, stage ids, or technical labels unless the user explicitly asks for command-level detail.
17a. When `analyze_path_task` returns `recommended_commands`, `alternative_commands`, or `latent_specialized_commands`, use them to decide what to explain next. Mention matrix-preparation commands such as `anno_eset`, `mouse2human_eset`, and `log2_eset` only when the task actually needs them.
17aa. For clustering requests, prefer `tme_cluster` before `nmf` unless the user explicitly asks for NMF or latent-factor discovery.
17ab. When a clustering card includes `existing_input_candidates`, `execution_guardrails`, or `pre_execution_questions`, use them before execution. Prefer a detected CIBERSORT matrix first when it is present, do not silently auto-switch to MCPcounter, and do not auto-use a very wide signature-score matrix without asking the user which subset they want.
17ac. For `tme_cluster` and `nmf`, remember that `--features` counts columns after excluding the first sample-ID column. Tutorial-style ranges therefore often start at `1`, such as `1:22`.
17ad. Before executing `tme_cluster` or `nmf`, make sure the output parent directory exists.
17b. When `analyze_path_task.agent_guidance.directory_recognition.layer_2_llm_reasoning_surface.should_expand_reasoning` is true, use its `reason_details`, `evidence_samples`, `candidate_function_hypotheses`, and `suggested_user_questions` to resolve ambiguity and widen function judgment instead of assuming the checklist and current stage already cover every relevant iobrpy capability.
17c. Prefer `current_state_reasoning_surface.summary_signals`, `current_state_reasoning_surface.default_posture_reason_ids`, `ranked_next_step_options[*].ranking_reasons`, and `next_step_reasoning_surface.ranking_criteria` over prewritten one-line assessments. The host LLM should compose the final wording itself.
17d. When a downstream matrix command card includes `upstream_recovery_paths`, keep them in mind before execution and especially after failure. If the current TPM-like matrix is unusable, prefer the best detected rollback path from Salmon or STAR outputs instead of giving up immediately.
17e. When a native command fails and the payload includes `failure_guidance.upstream_recovery_paths`, use the preferred recovery path before concluding that the user needs a different biology-specific input. For downstream functions, regenerated TPM should go through `log2_eset` before retrying the failed command.
17f. For BayesPrism specifically, do not jump straight to demanding a custom single-cell reference if the current TPM-like matrix failed once and a standardized TPM regeneration path from existing Salmon or STAR outputs is still available.
18. After `map`, run `iobrpy-cli recommend --path <path> --task "<task>" --json` only when you still need a lower-level native-command recommendation. Use `analyze_path_task` only after `map_path` when the scan result still needs richer planning. Do not call `recommend_workflow` reflexively right after `analyze_path_task`, because that would repeat reasoning the high-level helper already produced.
19. Inspect `recommendation.confirm_parameters`, `missing_parameters`, and `important_optional_parameters`. Confirm the parameters explicitly called out there when they materially affect biology, reference choice, or compute cost.
19a. When `concise_confirmation_plan` or `concise_confirmation_prompts` are present, prefer them over raw `confirm_parameters`. Ask only the most important 1-2 questions first, and phrase them as short natural-language numbered choices such as `1. ... 2. ... 3. ...` when possible.
19b. When a command has resource-sensitive parameters such as `threads`, `num_threads`, `parallel_size`, `num_processes`, `batch_size`, `-t`, or `-j`, do not blindly keep the native default or hard-code `8`. If the user did not specify the value, first check CPU cores, current load, and available memory on the same execution host, then ask a concise natural-language question such as `I checked the current CPU and memory load; can I continue with N threads?` in English or `我看了一下当前机器资源，建议用 N 个线程继续，是否可以？` in Chinese.
20. Prefer top-level native commands such as `iobrpy-cli runall`, `iobrpy-cli tme_profile`, `iobrpy-cli hla_typing`, and `iobrpy-cli trust4`.
21. Before recommending, writing, or executing a native command invocation, validate the command and flags against `src/iobrpy/RAG_MCP/iobrpy_required_params.json` by using `iobrpy-cli commands --json` or `iobrpy-cli <native_command> --help`. Do not invent flags from general workflow concepts. In particular, `tme_profile` has no `--mode`; its JSON-defined parameters are `--input`, `--output`, and optional `--threads`. Reserve `--mode` for commands that declare it, such as `runall`.
22. Inspect source code only when implementing code changes, debugging native command behavior, or when the CLI output is insufficient.
23. Do not substitute the R `IOBR` package when Python `iobrpy` is available.

## Checklist Cell Rendering

- If `workflow_checklist[].detected_column_value` or `workflow_checklist[].missing_column_value` is present, use those pre-rendered cell values directly when rebuilding the checklist table.
- If a client strips HTML line breaks inside table cells, preserve a visible separator such as `•` between adjacent evidence paths or missing items instead of concatenating them into one string.
- For the `immune_deconvolution` row, prefer `detected_column_value` or `workflow_checklist[].evidence_display_paths` over method-name summaries.

- For the `tpm_matrix_ready` row, prefer `detected_column_value` or `workflow_checklist[].evidence_display_paths` directly and do not add downstream signature-scoring files such as `calculate_sig_score.csv`.

## Input Routing

- FASTQ directory: prefer `runall`
- TPM-like expression matrix: prefer `tme_profile`
- BAM or CRAM directory for HLA typing: prefer `hla_typing`
- TRUST4 or repertoire request: prefer `trust4`
- Standalone deconvolution request from an expression matrix: use the matching native command, but stay inside `iobrpy-cli` or the registered MCP tools

## Guardrails

- Keep the work inside `iobrpy-cli` or iobrpy MCP tools unless the user explicitly asks for source inspection or code changes.
- Use `iobrpy-cli commands --json` or native `--help` to inspect the supported parameter surface when the exact workflow is unclear.
- For path-driven tasks, do not default to rerunning `runall` when `map --json` already shows completed intermediate outputs.
- If `recommend --json` returns `confirm_parameters`, ask a concise follow-up before execution instead of silently guessing those values.
