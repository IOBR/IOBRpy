---
name: iobrpy-result
description: Adaptively visualize, audit, and interpret outputs from every registered IOBRpy result-producing command except the non-result ai orchestrator. Builds claim-led, reviewer-ready scientific figures from content-verified IOBRpy evidence while preserving method semantics, provenance, uncertainty, source data, and editable exports. Covers FASTQ QC, Salmon/STAR quantification, matrix annotation and transformation, signature scoring, TME methods, ligand-receptor scores, clustering, TRUST4, HLA processing and typing, runall, and tme_profile. Use when the user asks an agent to explore, plot, compare, summarize, explain, biologically interpret, or prepare publication-ready figures from IOBRpy result directories or tables.
---

# IOBRpy Result

Design a result-specific analysis rather than applying a fixed plotting script. Treat the
skill as a scientific reasoning framework: IOBRpy semantics and QA rules are guardrails;
figure structure, feature selection, visual encoding, and narrative should respond to the
actual question and data.

## Startup Banner

At the start of each new `iobrpy-result` request, render this banner once before language
selection questions, backend questions, scanning, planning, interpretation, or plotting:

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

Do not repeat the banner before every tool call, figure, or result section. Do not use the
banner text as language-selection evidence.

## Language And Localization

Support both English and Simplified Chinese interactions.

- Select the runtime language from the user's current request only. For a slash-command
  invocation, treat the text after the command name as the sole language-selection input.
- Do not inherit the language from earlier conversation turns, project instructions,
  memories, startup banners, paths, filenames, tool output, or assistant messages.
- An explicit output-language request in the current input takes precedence. Otherwise,
  ignore paths, code, command names, and identifiers, then choose Simplified Chinese when
  the remaining natural-language request is predominantly CJK Han text; choose English
  when it is predominantly Latin-script English text.
- Store the decision as one request-scoped language code: `zh` or `en`. Use that same code
  for every question, progress update, explanation, figure label, caption, interpretation,
  and final response in the task.
- When calling an IOBRpy MCP helper such as `map_path`, pass `language: "zh"` for a Chinese
  current request and `language: "en"` for an English current request. Do not mix `_zh`
  fields into English output or English fields into Chinese output.
- Localize explanations, figure labels, captions, and interpretation prose dynamically.
- Keep commands, code identifiers, file paths, column names, method names, and raw error
  messages in their original form.
- If the current request is genuinely balanced or contains no natural-language text, ask
  one concise language question before any backend question.
- Keep packaged instructions in English. Runtime questions, explanations, labels, captions,
  and interpretation remain dynamically localized from the current request.

## Request Routing Before Backend Selection

First decide whether the user explicitly requested creation or modification of a
visualization.

Visualization intent is explicit when the user asks to plot, draw, visualize, make a
figure, create a chart, revise an existing figure, or prepare publication graphics.

The following requests are interpretation-only unless they also contain explicit
visualization intent:

- interpret, explain, audit, review, summarize, compare, or inspect results;
- identify IOBRpy files or classify mixed-directory provenance;
- report completed functions, quality issues, biological meaning, or limitations;
- inspect a table, directory, existing image, or report without asking to modify it.

For interpretation-only requests:

- do not ask the user to choose Python or R;
- inspect and classify the evidence immediately;
- provide the requested interpretation without creating plotting code or figure artifacts;
- do not infer visualization intent merely because the skill can create figures.

If intent is ambiguous, default to interpretation-only behavior. Ask about Python or R only
after the user explicitly requests a visualization.

## Plotting Backend Gate

Backend selection is a blocking gate for every request that creates or modifies a figure.
If the user has not explicitly selected Python or R in the current request, ask one concise
question in the request-scoped language. The question should naturally ask which of Python
or R the user prefers. An English example is:

```text
Would you prefer Python or R?
```

This is an example, not required fixed wording. Emit only one natural question in the
request-scoped language, then stop and wait. Do not add an explanation, translation, or
bilingual pair unless the user explicitly requests bilingual output. For an English current
request, ask in English even if earlier turns were Chinese. For a Chinese current request,
ask naturally in Simplified Chinese even if earlier turns were English. Do not inspect result
files, run `iobrpy-cli map`, draft a figure contract, write plotting code, or create previews
before the user answers.

Concrete gate example: for a current request such as
`please analyze /path/to/results and visualize`, select `en` and ask naturally in English,
for example, `Would you prefer Python or R?`. Do not prepend a Chinese or bilingual
explanation before the question.

For any current request classified as `zh`, ask naturally in Simplified Chinese without
emitting the English question or a bilingual pair.

A clearly Python- or R-specific script supplied by the user counts as an explicit backend
selection.

After selection, use only that backend for figure generation, preview rendering, export,
and visual QA:

- Python selection: use Python plotting libraries only.
- R selection: use R plotting libraries only.
- Do not use the other language to create a temporary preview or fallback figure.
- If the selected runtime or required packages are missing, report the exact blocker and
  stop before rendering instead of silently switching backend.

Write task-specific plotting code in the user's output directory and keep the original
result files unchanged.

## Figure Quality Contract

For every explicit visualization request, complete a claim-led figure contract after the
backend gate and evidence inventory but before plotting code. The contract must define:

- one-sentence core conclusion with an active verb;
- result-specific figure archetype and target journal/output context;
- final physical dimensions and required editable/raster formats;
- panel map with one hero panel and non-redundant context, validation, or robustness roles;
- statistics, source-data, uncertainty/fit, and image-integrity requirements;
- the main reviewer risk and the visual choice that addresses it.

Do not start from a favorite template. If hiding a panel would not weaken the conclusion,
remove or merge it. Set the final-size and export contract before styling so typography,
annotation density, and panel balance are judged at the delivered dimensions.

## Core Workflow

1. Establish the evidence inventory.
   - For an IOBRpy output directory, run `iobrpy-cli map --path <path> --json` once
     and use its concrete evidence paths.
   - Read [references/source-classification.md](references/source-classification.md) and
     classify every candidate file before reading values into one analysis.
   - If the directory is mixed, nonstandard, or contains likely/reusable/external evidence,
     rerun `iobrpy-cli map --path <path> --json --investigate-existing` and use the enriched
     provenance fields.
   - For explicit tables, inspect headers, shape, IDs, ranges, missingness, and companion
     files directly, then apply the same source-classification rules.
   - Read [references/result-layout.md](references/result-layout.md) when the result type,
     stage, or orientation is unclear.
   - Build a provenance ledger with one row per candidate path: detected role, candidate
     function, source class, decisive evidence, confidence, and whether it is eligible for
     IOBRpy-specific interpretation or plotting.
2. Build a result profile before interpretation or visual design.
   - Open the actual result table and inspect its delimiter, headers, shape, orientation,
     sample/feature identifiers, numeric columns, representative values, ranges, and
     missingness. Do not infer method semantics from filenames alone.
   - Record the method, unit, orientation, sample key, feature family, uncertainty or fit
     fields, metadata coverage, exact IOBRpy version, implementation variant, and effective
     transformation chain.
   - Distinguish the published method family from the current IOBRpy implementation. Do not
     assume numerical parity with an official R package, web service, or reference
     implementation when IOBRpy uses different preprocessing, solvers, defaults, or output
     fields.
   - Read [references/function-coverage.md](references/function-coverage.md) to identify
     every detected IOBRpy function and its interpretation reference.
   - Read [references/interpretation.md](references/interpretation.md), then only the
     relevant domain guide named by the coverage table.
   - Interpret every detected function that produced evidence. Do not skip upstream QC,
     quantification, annotation, normalization, or workflow-wrapper results merely because
     downstream biological tables are also present.
   - Treat only content-confirmed native evidence as confirmed IOBRpy output. Keep likely,
     reusable, external, metadata, and unknown files visibly separate.
   - Separate provenance confirmation from execution completeness. One content-confirmed
     native file can confirm its own role without proving that all expected companion files,
     samples, loci, chains, diagnostics, or child stages completed.
3. Branch by request type.
   - Interpretation-only: skip the plotting backend question, figure contract, plotting code,
     preview generation, and visual QA. Continue directly to evidence-layer interpretation.
   - Visualization: require the plotting backend gate, establish the evidence inventory,
     then complete the claim-led figure contract with
     [references/figure-framework.md](references/figure-framework.md). State the core
     conclusion, archetype, target output, final dimensions, panel map, comparison unit,
     evidence hierarchy, selected backend, statistics/source-data needs, reviewer risk, and
     at least two materially different visual designs. Choose the design that best exposes
     the evidence, not the most familiar chart. Read
     [references/design-examples.md](references/design-examples.md) when the first design
     still looks generic or when several result families could be integrated.
4. For visualization requests, implement a bespoke figure.
   - Use [references/visualization.md](references/visualization.md) as a visual grammar,
     not a template catalog.
   - Adapt layout, ordering, transformations, annotations, and panel balance to the result.
   - Give the hero evidence the clearest axis or largest panel. Keep context, validation,
     uncertainty, and robustness panels visually quieter.
   - Reuse one biological color vocabulary, typography hierarchy, panel-label system, and
     legend strategy across the figure.
   - Create one coherent visual argument; do not automatically emit one plot per table or
     divide the canvas into equal panels when the evidence importance is unequal.
   - Export from the selected backend at final dimensions, open the rendered preview, inspect
     hierarchy, clipping, overlap, legibility, color scales, and empty space, then revise and
     re-render until the applicable QA checks pass.
   - Track every generated image file and resolve it to a normalized absolute path.
5. Interpret in evidence layers.
   - Separate direct observations, method-aware meaning, biological hypotheses, and limits.
   - Explain why the selected visual encoding supports the conclusion.
6. Run the applicable sections of [references/qa.md](references/qa.md). For visualization
   requests, complete its render-inspect-revise loop at final size and verify editable text,
   source-data traceability, and the export bundle before delivery.

## Coverage Boundary

The supported result-producing command set is defined by
`src/iobrpy/RAG_MCP/iobrpy_required_params.json` and mirrored in
[references/function-coverage.md](references/function-coverage.md).

- Cover all 27 registered commands in that table.
- Exclude `ai`: it orchestrates tools but does not define a standalone biological result
  format to interpret.
- A workflow wrapper such as `runall` or `tme_profile` requires both wrapper-level
  completeness/provenance interpretation and interpretation of each concrete child result.

## Freedom And Guardrails

The agent may freely choose:

- single-panel or multi-panel composition;
- heatmap, distribution, rank, composition, scatter, ordination, alluvial, tile, forest,
  or another defensible visual form;
- feature selection and ordering rules;
- whether several result tables should be integrated;
- labels, palette details, annotation density, and figure dimensions;
- exploratory statistics appropriate to the supplied design.

The agent must not freely change:

- the biological meaning or scale of an IOBRpy result;
- sample IDs, group labels, or missing values;
- fractions into percentages without disclosure;
- method-specific scores into cell counts;
- relative abundance into absolute abundance;
- LR co-expression scores into proven interactions;
- exploratory associations into causal or clinical claims.
- external or ambiguous files into native IOBRpy results;
- files from different provenance classes into one unlabeled matrix.

Do not add a generic `iobrpy-result plot` command or depend on a rigid packaged plotting
script. Reusable snippets are acceptable only as low-level styling or export helpers.

## Minimum Data Audit

Before analysis, verify:

- delimiter, table orientation, and sample-ID column;
- duplicate IDs and metadata join losses;
- intended numeric columns versus labels or fit statistics;
- NaN, Inf, all-zero, constant, and out-of-range values;
- whether a composition sums to one and whether a residual category is explicit;
- whether uncertainty, fit, depth, pairing, batch, or covariate information is available.

Do not silently normalize, impute, filter, or drop samples. Make every such choice explicit
in code and interpretation.

## Interpretation Structure

Use this order:

1. **Identified result**: method, exact file, table shape, sample orientation, and completion.
2. **Field meanings**: explain the important columns, units, score direction, and quality or
   uncertainty fields actually present.
3. **Key observations**: report concrete patterns, ranges, missing values, outliers, or group
   differences visible in the supplied values.
4. **Method-aware meaning**: explain what those observations can mean biologically or
   technically on the method's actual scale.
5. **Necessary limitations**: state only the main constraint needed to avoid overclaiming.
6. **Optional visualization handoff**: when useful, name the comparison, fields, quality
   overlay, and visual form that would best show the finding. Do not create a figure or ask
   about Python/R unless the user explicitly requests visualization.

Prefer cross-method concordance in direction or rank. Do not compare raw magnitudes across
different deconvolution methods.

## Deliverables

Choose deliverables that fit the request. For a publication-oriented task, normally include:

- editable SVG;
- PNG preview;
- PDF when the target workflow accepts vector PDF;
- TIFF when required by the target journal or image-heavy submission workflow;
- task-specific plotting script;
- exact plotted source data;
- concise figure contract and QA notes;
- concise interpretation and caveats;
- statistics table when comparisons were performed.

For every successful visualization request, the final response must include the complete
absolute path of every generated image file:

- list each SVG, PNG, PDF, TIFF, JPEG, or other image/figure export separately;
- verify each listed file exists before reporting completion;
- use the final normalized absolute path, not a basename, relative path, shell shorthand,
  or parent output directory;
- do not make the user infer filenames from an output directory;
- distinguish the editable figure from preview or alternate-format exports;
- if no image was successfully generated, state that clearly and do not report a planned
  or nonexistent path.

The absolute image-path list is mandatory even when the image is also embedded or linked in
the response.

Do not create empty boilerplate artifacts merely to satisfy a fixed folder layout.
