# Source Classification For Mixed Directories

Classify every candidate file before visualization or interpretation. A shared directory,
similar filename, or related biological topic does not establish IOBRpy provenance.

## Mandatory Source Classes

Use exactly these reasoning classes:

| Source class | Admission rule | Allowed use |
|---|---|---|
| `iobrpy_confirmed` | Native method-specific content is verified, normally represented by `confirmed_by_content` and `result_source_id=iobrpy_confirmed_results` | Interpret and plot as an IOBRpy result |
| `iobrpy_likely` | Name/layout suggests IOBRpy but method-specific content or companions are incomplete; normally `likely_iobrpy_result` | Inspect further; label as likely, never confirmed |
| `compatible_reusable` | Generic content could be reused by an IOBRpy function but native execution is unproven; normally `reusable_result` | Treat as compatible existing data, not IOBRpy-produced output |
| `external_result` | Content, tool markers, or investigation evidence indicates a non-IOBRpy tool; normally `external_analysis_hints` or `external_tool_results` | Interpret separately when relevant; never claim native equivalence |
| `metadata_or_auxiliary` | Phenotype, clinical annotation, sample sheet, script, log, figure, checksum, or documentation | Use only for its declared supporting role |
| `unknown_or_ambiguous` | Evidence is insufficient or conflicting | Exclude from IOBRpy-specific claims and plots until resolved |

## Evidence Priority

Resolve conflicts from strongest to weakest:

1. Method-specific headers, row/column structure, value constraints, embedded tool/version
   markers, and companion-file consistency.
2. Native execution evidence such as expected child outputs, logs, and completion markers.
3. Exact IOBRpy output basename combined with compatible content.
4. Standard IOBRpy directory layout combined with compatible content.
5. Filename keywords, directory names, or broad biological vocabulary.

Levels 4 and 5 alone never establish confirmed native provenance. If a filename claims one
method but the content signature indicates another tool or an incompatible schema, classify
the file as `external_result` or `unknown_or_ambiguous`, and record the conflict.

## Required Map Fields

For a directory request, first consume:

- `function_detections[].status`;
- `function_detections[].matched_by_content`;
- `function_detections[].matched_by_name`;
- `function_detections[].result_source_id`;
- `content_verified_functions`;
- `likely_iobrpy_functions`;
- `reusable_result_functions`;
- `external_analysis_hints`;
- `agent_fallback.recommended`;
- `existing_result_investigation`, when present.

Use `matched_by_content` as the primary native evidence path list. Do not replace it with
all files from the same directory. A sibling file is not automatically part of the same
method.

When the initial map reports mixed, likely, reusable, external, or fallback evidence, run:

```text
iobrpy-cli map --path <path> --json --investigate-existing
```

Use `existing_result_investigation.target_summaries`, `bucket_counts`, `sample_paths`, and
external evidence profiles to refine the ledger. Investigation can improve classification,
but it must not upgrade an external result to native merely because an equivalent IOBRpy
function exists.

## Provenance Ledger

Create a working table before analysis:

```text
path | role | candidate function | source class | decisive evidence | confidence | use
```

Rules:

- One path can have only one primary source class.
- Preserve full paths; do not collapse evidence to a parent directory.
- Record companion relationships separately rather than merging files into one row.
- Add a separate completion field such as `complete`, `partial`, `failed`, or `unknown`.
  Source class answers who produced a file; completion answers whether the expected result
  set is present. Never infer one from the other.
- Mark inferred roles and unresolved conflicts explicitly.
- Keep files not used in the analysis in the ledger with `use=excluded`.

## Integration Rules

- Native and external results may appear in one figure only when they answer the same user
  question and their provenance remains explicit in panel titles, legends, captions, and
  exported source data.
- Add a `source` or `provenance` column before combining rows from different classes.
- Never concatenate same-named score columns from native and external tools without method
  and source suffixes.
- Metadata can be joined to results by verified sample IDs, but it remains metadata rather
  than an IOBRpy output.
- Scripts, notebooks, images, and prose reports are not quantitative result tables merely
  because their names contain an IOBRpy method.
- Unknown files must not influence feature selection, statistics, or biological conclusions.

## Common Collision Cases

- A generic `signature_scores.tsv` beside `calculate_sig_score.csv`: only the latter is
  native when its content signature is verified; the generic file remains external or
  ambiguous.
- A CellChat/CellPhoneDB-style communication table beside `lr_cal.csv`: it is an external
  communication result, not an `LR_cal` output.
- MiXCR/AIRR clonotype files beside TRUST4 summaries: keep external repertoire evidence
  separate from native TRUST4 evidence.
- A generic HLA genotype table beside `hla_result_merged.txt`: classify by header, version,
  tool markers, and companions; do not infer native SpecHLA from HLA vocabulary alone.
- A counts or TPM table with a user-defined filename: it may be compatible reusable data,
  but native `count2tpm`, `prepare_salmon`, or merge provenance remains unproven.
