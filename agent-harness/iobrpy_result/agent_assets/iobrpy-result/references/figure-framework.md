# Figure Contract And Design Framework

Complete this contract in working notes before plotting. Keep it concise and revise it after
the first data inspection. The figure is a visual argument, not a collection of available
tables.

Use this framework only after the user explicitly requests creation or modification of a
visualization. Do not use it for interpretation-only requests.

```text
User question:
Core conclusion:
Figure archetype:
Target journal/output:
Final physical size:
Result files:
Provenance ledger:
  confirmed native:
  likely or reusable:
  external:
  metadata:
  excluded or unresolved:
Selected backend: Python or R
Result profile:
  method and scale:
  rows / columns:
  sample key:
  metadata coverage:
  fit / uncertainty fields:
Comparison unit:
Panel map:
  a:
  b:
  c:
Evidence hierarchy:
  hero evidence:
  validation evidence:
  context / controls:
  robustness:
Main reviewer risk:
Candidate design A:
Candidate design B:
Chosen design and why:
Feature-selection rule:
Ordering rule:
Statistics:
Source data needed:
Fit / uncertainty treatment:
Image-integrity notes:
Final formats:
```

Do not complete this contract for a figure request until the user has explicitly selected the
backend. The selected backend is exclusive for plotting, previewing, exporting, and visual
QA.

## Core Conclusion Rules

- Write one sentence with an active verb, not a topic label. Prefer "Responders show a
  coordinated cytotoxic program across signatures and deconvolution methods" over "TME
  results".
- Give the primary evidence the hero panel or clearest quantitative axis.
- Make every panel answer one unique question. If covering a panel would not weaken the
  argument, remove or merge it.
- If the user supplied data without a claim, infer a provisional conclusion from the
  request and observed evidence, then keep it explicitly provisional until the analysis
  supports it.

## Result-Specific Archetypes

| Archetype | Use when | Hero evidence | Supporting evidence |
|---|---|---|---|
| `quantitative comparison grid` | The claim is mainly a set of aligned group or method effects | Raw distributions plus effect/interval summary | Fit, uncertainty, adjusted statistics |
| `landscape + focused inference` | A cohort-wide composition or state landscape must precede a specific comparison | Composition, heatmap, ordination, or ranked landscape | Focused group effects, metadata, fit strip |
| `workflow evidence composite` | `runall`, QC, quantification, and downstream completion form one reproducibility claim | Workflow/stage evidence or decisive QC/quantification panel | Compact downstream validation panels |
| `asymmetric mixed-result figure` | Several result families test one biological explanation | One dominant TME/signature/concordance panel spanning multiple cells | Smaller context, uncertainty, or robustness panels |

Choose the archetype before sizing subplots. Do not use equal panel areas when the evidence
is not equally important.

## Journal And Export Contract

- Use about `89 mm` for a single-column figure and about `183 mm` for a double-column figure
  unless the target journal or user specifies another size.
- Design typography, line weights, annotation density, and legend placement at the final
  physical size rather than on an enlarged development canvas.
- Treat editable SVG as the primary general-purpose figure output. Add PDF and TIFF only
  when required by the target workflow; always create a PNG preview for visual QA.
- Record the plotted source-data table and statistical output for every quantitative panel.
- Keep provenance visible when native IOBRpy, compatible reusable, and external evidence
  appear in one figure.

## Design Selection

Compare at least two materially different designs. A stacked bar with a second stacked bar
is not a second design. Consider whether the question is best answered by:

- sample composition versus group effect;
- global structure versus selected biology;
- one method in depth versus cross-method concordance;
- point estimates alone versus estimates paired with fit or uncertainty;
- one integrated figure versus several independent figures.

Choose the design that minimizes likely misinterpretation. A more complex layout is useful
only when each panel contributes non-redundant evidence.

## Panel Logic

Use this order unless the result story requires another:

1. Establish cohort, workflow stage, or global result landscape.
2. Show the main effect or primary comparison.
3. Add method-compatible biological support from a distinct result family.
4. Surface fit, uncertainty, depth, missingness, or model diagnostics.
5. Add robustness, subgroup, sensitivity, or alternative-method evidence.

Reuse the same group, method-family, and biological-compartment visual vocabulary across
all panels. Prefer one shared legend or direct labels over repeated legends.

## Result Profiles

Use a short result profile to stop visual design from outrunning semantics:

```text
result class:
measurement unit:
bounded or unbounded:
compositional:
sample orientation:
quality fields:
known transformation:
valid comparison:
invalid comparison:
```

Examples:

- CIBERSORT relative output: bounded, compositional, sample-level fit available, valid
  within-cell-type cohort comparison, invalid cross-method magnitude comparison.
- MCPcounter: unbounded method-specific abundance score, non-compositional, valid relative
  within-feature comparison, invalid percentage interpretation.
- LR score: minimum of paired bulk-expression signals after transformation, valid ranking
  of co-expression-compatible pairs, invalid proof of cell-cell interaction.

## Evidence Integration

When several outputs are available, integrate them only when they answer one shared
question. Useful integration patterns include:

- immune fraction plus immune signature plus ESTIMATE context;
- BayesPrism theta plus theta CV;
- cluster membership plus defining features plus phenotype metadata;
- TRUST4 diversity plus read depth plus top-clone concentration;
- HLA calls plus cohort allele frequency, without inferring antigen presentation.

Do not create a dashboard of every available result. Breadth without a claim is not a
scientific figure.

When a directory mixes IOBRpy and non-IOBRpy files, preserve provenance in every integrated
artifact. Panel titles, legends, captions, and exported source data must distinguish native
IOBRpy, compatible reusable, and external results. Unknown files stay excluded.

## Interpretation Prompts

After rendering, ask:

1. What exact visual pattern is robust to reasonable ordering or transformation choices?
2. Which parts come directly from data and which come from the method model?
3. What alternative technical or biological explanation remains?
4. Does fit, uncertainty, depth, or missingness weaken any highlighted sample?
5. What new evidence would most efficiently distinguish the leading explanations?
6. Does the conclusion remain obvious when the figure is viewed at final physical size?
