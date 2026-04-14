---
description: Route a bulk RNA-seq, TME, HLA, or TCR/BCR request through iobrpy-cli and iobrpy MCP tools.
---

# IOBRpy

Treat `"$ARGUMENTS"` as the active IOBRpy request. If no arguments were provided, ask the user for a path or a concrete analysis goal before doing anything else.

## Startup Banner

At the start of each new `/iobrpy` request, render this banner once before scanning or planning:

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

1. Prefer iobrpy MCP helper tools first. For a user-provided path, start with `map_path`, not `analyze_path_task`. Pass `language: "en"` for English user requests and `language: "zh"` for Chinese user requests so compact checklist and confirmation fields match the conversation language.
2. Start with `iobrpy-cli doctor --json` only when the environment or command availability is genuinely unclear. Do not call it by default for an ordinary path scan.
3. For a user-provided path, do not ask the user to choose quick versus full scanning. Run `map_path` directly; the default scan path is the deep/full path for agent use.
4. Use `analyze_path_task` only after `map_path` when the scan result still needs richer planning. Do not call `recommend_workflow` immediately after `analyze_path_task` unless the user changed goals or you still need a narrower low-level recommendation.
5. Do not use a scan-mode selector, clickable choice UI, or text prompt for quick/full mode. Legacy quick mode is deprecated for agent use and should not be exposed to the user.
6. Present `workflow_checklist` first for path-driven parsing requests. Then summarize the current state in beginner-friendly language and explain the next choices.
6a. When `workflow_checklist_table`, `workflow_checklist_table_zh`, or `workflow_checklist[].missing_items_zh` is present, use it as the authoritative checklist rendering source and keep the `Still missing` / `仍缺少` column visible.
6b. For the `immune_deconvolution` checklist row, BayesPrism is an immune deconvolution method, but it is a standalone optional branch outside the default `tme_profile` / `runall` six-method bundle; when it is missing, keep it in the `仍缺少` column as `BayesPrism（独立可选，不由 tme_profile/runall 自动运行）`.
6c. Do not use `deconvolution_methods.default_bundle_missing_labels` to render the checklist `Still missing` / `仍缺少` column. Those labels only describe what `tme_profile` or `runall` can backfill automatically, and intentionally exclude standalone BayesPrism.
6aa. In English conversations, do not pull `missing_items_zh`, `missing_column_value_zh`, or other `_zh` checklist fields into English table cells. For the `immune_deconvolution` row in English, use `missing_column_value` or `missing_items` so BayesPrism keeps the English standalone-optional note.
7. Once `map_path` has succeeded, do not call `analyze_path_task` merely to render `workflow_checklist`, current stage, evidence paths, missing items, or next-step summaries. Answer directly from the map payload.
8. Keep the default `response_detail=summary` for normal user-facing MCP calls. Use `response_detail=structured` only when the user explicitly asks for technical/debug details or the map summary is insufficient for a real command-planning decision.
9. Do not say "the scan result is large" and then call `analyze_path_task` to extract fields. A large or folded MCP display is a rendering issue, not a reason to call another large helper.
10. Do not run or ask approval for helper shell commands whose only purpose is extracting fields from a previous `map_path` or MCP JSON result. In normal summary responses, render `default_user_facing_checklist` directly and use compact `missing_downstream_analysis_suggestions` from the tool result. Use structured fields such as `workflow_checklist` or `suggested_command_details` only when `response_detail=structured` explicitly returned them. In Chinese conversations, use the corresponding `_zh` fields when present.
11. If a real action still needs user confirmation, translate the intent into a concise natural-language question rather than exposing raw commands, temp-file paths, or parser snippets. For example, ask `Should I summarize the checklist and next-step suggestions from the scan result?` in English or `是否要我读取刚才扫描结果中的检查表和下一步建议？` in Chinese instead of showing a `python3 -c` command or a client `tool-results` path.
12. If an MCP helper response is folded or saved to a client `tool-results` file because it is too large, do not read that file just to recover checklist fields and do not run `python3 -c`, `jq`, `cat`, or similar JSON parser commands. Fall back to the most recent visible `map_path` payload fields, or rerun a compact `map_path` scan with the correct `language` value if genuinely necessary.
13. Do not limit yourself to three suggestions. For every recommendation you present, explain why it fits the current state, what input data it needs, and what result it will produce.
13a. When `missing_downstream_analysis_suggestions` is present, use it as the recommendation source after the checklist. It is expected to include every missing downstream analysis, including standalone BayesPrism; do not collapse it to only the default `tme_profile` bundle.
13b. When `suggested_command_details` contains parameter hints or placeholders, show concise natural-language prompts in the user's language, such as `path_salmon_index="替换你的 Salmon index 的路径"` for Chinese or `path_salmon_index="Replace with your Salmon index path."` for English. These hints are display guidance layered on top of `src/iobrpy/RAG_MCP/iobrpy_required_params.json`; do not invent flags outside that JSON.
14. For clustering requests, prefer `tme_cluster` before `nmf` unless the user explicitly asks for NMF or latent-factor discovery.
15. When clustering inputs already exist, prefer a detected CIBERSORT matrix first if one is present. Do not silently switch to MCPcounter or a very wide signature-score matrix without explaining the choice.
16. Do not auto-run clustering on a signature-score matrix unless the user clearly chose the signature subset or column family. Those files can be extremely wide.
17. For `tme_cluster` and `nmf`, remember that `--features` counts feature columns after excluding the first sample-ID column. Tutorial-style ranges therefore often start at `1`, such as `1:22`.
18. Before executing `tme_cluster` or `nmf`, make sure the output parent directory exists.
19. Mention matrix-preparation commands such as `anno_eset`, `mouse2human_eset`, and `log2_eset` only when the task really needs those preprocessing steps.
20. If a downstream matrix command fails and the payload includes `failure_guidance.upstream_recovery_paths`, use the preferred rollback path from existing Salmon or STAR outputs before giving up. Regenerated TPM should go through `log2_eset` before retrying downstream functions.
21. For BayesPrism specifically, do not jump straight to requiring a custom single-cell reference if a standardized TPM regeneration path is still available from detected upstream outputs.
22. Once the path state and parameters are clear, prefer top-level native commands such as `runall`, `tme_profile`, `hla_typing`, and `trust4`.
23. Before recommending, writing, or executing a native command invocation, validate the command and flags against `src/iobrpy/RAG_MCP/iobrpy_required_params.json` by using `iobrpy-cli commands --json` or `iobrpy-cli <native_command> --help`. Do not invent flags from general workflow concepts. In particular, `tme_profile` has no `--mode`; its JSON-defined parameters are `--input`, `--output`, and optional `--threads`. Reserve `--mode` for commands that declare it, such as `runall`.
24. When `concise_confirmation_plan` or `concise_confirmation_prompts` are present, use them instead of dumping raw parameter names. Ask only the most important 1-2 confirmation items first, and format them as short numbered natural-language choices when possible.
25. When a command has resource-sensitive parameters such as `threads`, `num_threads`, `parallel_size`, `num_processes`, `batch_size`, `-t`, or `-j`, do not blindly keep the native default or hard-code `8`. If the user did not specify the value, first check CPU cores, current load, and available memory on the same execution host, then ask a concise natural-language question such as `I checked the current CPU and memory load; can I continue with N threads?` in English or `我看了一下当前机器资源，建议用 N 个线程继续，是否可以？` in Chinese.

## Checklist Cell Rendering

- If `workflow_checklist[].detected_column_value` or `workflow_checklist[].missing_column_value` is present, use those pre-rendered cell values directly when rebuilding the checklist table.
- If a client strips HTML line breaks inside table cells, preserve a visible separator such as `•` between adjacent evidence paths or missing items instead of concatenating them into one string.
- For the `immune_deconvolution` row, prefer `detected_column_value` or `workflow_checklist[].evidence_display_paths` over method-name summaries.

- For the `tpm_matrix_ready` row, prefer `detected_column_value` or `workflow_checklist[].evidence_display_paths` directly and do not add downstream signature-scoring files such as `calculate_sig_score.csv`.

## Guardrails

- Keep the work inside `iobrpy-cli` or the registered iobrpy MCP tools unless the user explicitly asks for source inspection or code changes.
- Do not substitute the R `IOBR` package when Python `iobrpy` is available.
