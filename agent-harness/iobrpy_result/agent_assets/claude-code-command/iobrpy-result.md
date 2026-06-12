---
description: Visualize and interpret IOBRpy result directories and tables.
argument-hint: [result path, files, and question]
---

Use the installed `~/.claude/iobrpy-result/SKILL.md` as the governing workflow.
Read only the referenced files needed for the supplied result type. Inspect, visualize,
and interpret the IOBRpy results described by:

$ARGUMENTS

Before any question, scan, planning text, interpretation, or plotting work, render this
startup banner once:

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

Do not repeat the banner, and do not use it as language-selection evidence.

Determine the response language only from `$ARGUMENTS`. Do not inherit language from
earlier conversation turns, project instructions, memories, startup banners, paths,
filenames, or tool output. An explicit output-language request in `$ARGUMENTS` wins;
otherwise ignore paths, code, command names, and identifiers, then select `zh` when the
remaining natural-language text is predominantly CJK Han text and `en` when it is
predominantly Latin-script English text. Use that request-scoped language for every
question, update, figure label, caption, interpretation, and final response. When calling
`map_path`, pass the same value as its `language` argument.

Ask about Python or R only when `$ARGUMENTS` explicitly requests creation or modification
of a plot, chart, figure, or other visualization. If visualization is explicitly requested
and no backend is selected, ask the localized equivalent of `Python or R?` in the user's
request-scoped language and stop before inspecting data or planning the figure. In
particular, English `$ARGUMENTS` must produce an English backend question even when prior
conversation turns were Chinese.

Concrete example: if `$ARGUMENTS` is
`please analyze /path/to/results and visualize`, select `en` and ask exactly
`Python or R?`. Do not prepend a Chinese or bilingual explanation.

For interpretation, explanation, audit, summary, comparison, provenance classification, or
result inspection without explicit visualization intent, do not ask about Python or R.
Inspect the evidence and answer directly. If intent is ambiguous, default to
interpretation-only behavior.

Keep the original result files unchanged. Design the analysis around the actual question
and data instead of applying a fixed plot template. Produce editable SVG figures, PNG
previews, source-data tables, reproducible task-specific plotting code, and method-aware
interpretation when the request calls for artifacts.

After successfully generating any figure, verify every image file exists and list the
complete normalized absolute path of each image in the final response. Do not report only a
basename, relative path, or parent output directory, even when the image is embedded.

Interpret every detected registered IOBRpy result-producing function, including upstream
QC, quantification, annotation, transformation, and workflow-wrapper evidence. Exclude
`deside` and exclude the non-result `ai` orchestration interface.

For mixed directories, classify each candidate path as confirmed IOBRpy, likely IOBRpy,
compatible reusable data, external result, metadata/auxiliary, or unknown before analysis.
Use content evidence and map provenance fields; never infer native provenance from a sibling
file, directory name, or filename keyword alone.

Support English and Simplified Chinese interactions while keeping commands, identifiers,
paths, and raw errors in their original form. All packaged source files remain English ASCII;
localization happens at runtime.
