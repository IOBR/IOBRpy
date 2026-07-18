# Adaptive Visualization Grammar

Use this file to reason about visual encodings. It is not a lookup table that assigns one
mandatory chart to each IOBRpy method.

## Backend Exclusivity

Use the backend selected by the user for all visual work. Python figures must be generated,
previewed, exported, and revised in Python. R figures must be generated, previewed,
exported, and revised in R. Data inspection may use ordinary shell tools, but do not invoke
the non-selected language's graphics stack.

## Final Size And Typography

Choose final physical dimensions before styling:

- single-column target: about `89 mm` wide;
- double-column target: about `183 mm` wide;
- use another size only when the target journal, report, or user specifies it.

Judge typography at final size. For dense manuscript figures, body, tick, and legend text
should usually remain readable in the `5-7 pt` range, with bold lowercase panel labels near
`8 pt`. Use one sans-serif stack such as Arial, Helvetica, DejaVu Sans, then sans-serif.
Keep ordinary analytical panels on white; use dark backgrounds only when the input modality
requires them. Preview at the delivered dimensions rather than relying on an enlarged
development canvas.

## Start From The Analytical Question

Choose a visual grammar from the question:

| Question | Candidate visual grammar |
|---|---|
| What dominates each sample? | composition, rank, tile, small multiples |
| Which samples share a profile? | heatmap, ordination, distance map, cluster ribbon |
| What differs between groups? | raw points plus interval, paired slope, effect-size forest |
| Which features drive a state? | ranked effects, loading plot, annotated heatmap |
| Do methods agree? | within-method ranks, concordance scatter, sign matrix |
| How uncertain or reliable is the estimate? | estimate-uncertainty pairing, fit overlay, opacity/annotation |
| How does repertoire structure differ? | rank-frequency, diversity-concentration plane, clone tile |
| How are categorical calls distributed? | allele tile, frequency lollipop, mosaic |

Several grammars may be valid for the same table. Select the one that makes the intended
comparison easiest to perceive and hardest to misread.

## Compose A Visual Argument

Assign each panel one evidence role:

- **Hero**: the main pattern or comparison.
- **Context**: cohort structure, method fit, uncertainty, sequencing depth, or covariates.
- **Mechanism-compatible support**: a related signature or result family that tests the
  same biological explanation from another angle.
- **Robustness**: alternative method, sensitivity analysis, paired view, or residual check.

Use asymmetric layouts when evidence importance is asymmetric. Remove a panel when hiding
it would not weaken the conclusion.

## Layout Archetypes

Choose one archetype from the figure contract before allocating axes:

| Archetype | Practical layout |
|---|---|
| `quantitative comparison grid` | Aligned raw distributions and effect/interval panels with shared scales |
| `landscape + focused inference` | Wide composition/heatmap/ordination hero plus a smaller focused comparison row |
| `workflow evidence composite` | Stage/QC overview leading compact quantification and downstream validation panels |
| `asymmetric mixed-result figure` | One panel spans rows or columns while smaller context and robustness panels remain subordinate |

Do not force equal subplot sizes. A hero panel may occupy `40-60%` of the usable area when
it carries the primary evidence. Use bold lowercase panel labels near the top-left, real but
compact gutters, and alignment/whitespace rather than decorative panel boxes.

Prefer direct labels for stable identities. When legends are necessary, use one shared
legend strip or a legend-only axis instead of repeating the same legend in every panel.

## Method-Aware Design Constraints

### Fraction and proportion outputs

CIBERSORT relative mode, quanTIseq, EPIC, and BayesPrism theta can support compositional
encodings. Before stacking:

- verify the expected total and residual category;
- preserve observed totals when they do not equal one;
- avoid silently renormalizing;
- consider a log-ratio or selected-component view for group inference;
- surface fit or uncertainty beside the composition when available.

Stacked bars are useful for sample-level composition but weak for precise between-group
comparisons. Pair them with a focused effect or distribution panel when the question is
inferential.

### Method-specific score outputs

MCPcounter, ESTIMATE, IPS, signature scores, and LR scores are not parts of a whole. Do not
stack them. Use raw or within-feature standardized values depending on the question, and
label every transformation.

For cross-method views, compare direction, rank, standardized effect, or concordance.
Never imply that equal numeric values from different methods have equal biological meaning.

### Uncertainty and fit

- CIBERSORT: integrate `P-value`, `Correlation`, and `RMSE` into the figure or interpretation.
- BayesPrism: pair theta patterns with `theta_cv`; high-CV estimates should not look equally
  certain.
- TRUST4: show repertoire depth when interpreting richness or diversity.
- Clustering: show cluster size and defining features; use stability evidence if available.
- HLA: retain missing/no-call states rather than coloring them as an allele.

## Feature Selection

Choose features using a rule tied to the question, for example:

- prespecified biology;
- largest absolute group effect;
- highest cohort variance for exploratory structure;
- most abundant components for composition overview;
- strongest uncertainty-adjusted signal;
- cluster-defining loadings;
- concordant evidence across methods.

State the rule and tested feature family. Do not select by p-value and then present the same
p-value as confirmatory evidence.

## Ordering

Order samples and features by an interpretable variable:

- supplied group and paired subject;
- cluster;
- a declared score;
- hierarchical similarity;
- dominant component;
- genomic locus or cell lineage.

Avoid arbitrary sorting that creates a false trajectory. If clustering is used only for
display, say so.

## Palette

Use one restrained palette system per figure:

- one neutral family for context or baselines;
- one signal family for the main biological or clinical comparison;
- one accent family for the decisive feature, warning, or directional annotation.

Use the larger biological-role palette only when the encoding genuinely requires several
cell compartments, such as a composition landscape. Keep group colors distinct from
cell-type colors, and never remap the same group, method, or compartment between panels.
Do not use red and green as the only distinction.

Use stable biological roles rather than assigning colors from file order. The same hex
mapping can be represented as a Python dictionary or an R named vector:

```python
IOBRPY_COLORS = {
    "b_cell": "#4C78A8",
    "t_cell": "#59A14F",
    "nk_cell": "#76B7B2",
    "myeloid": "#E15759",
    "stromal": "#B07AA1",
    "tumor": "#F28E2B",
    "other": "#BAB0AC",
    "score": "#2F4B7C",
    "highlight": "#D45087",
}
```

Adapt shades and accents to the figure. Keep the same biological role consistent across
panels. Do not use red/green as the only distinction.

For an R implementation, preserve the same values explicitly:

```r
IOBRPY_COLORS <- c(
  b_cell = "#4C78A8", t_cell = "#59A14F", nk_cell = "#76B7B2",
  myeloid = "#E15759", stromal = "#B07AA1", tumor = "#F28E2B",
  other = "#BAB0AC", score = "#2F4B7C", highlight = "#D45087"
)
```

## Statistics

- Show individual samples where feasible.
- Respect pairing and repeated measures.
- Report group `n`, effect size, interval, test, and multiple-testing correction.
- Use permutation or rank-based approaches when assumptions are weak and sample size allows.
- Treat exploratory screens as exploratory.
- Do not run group tests without supplied or verifiable metadata.

## Publication Output

- Keep SVG text editable; for Python use `svg.fonttype = "none"`, and for R use `svglite`.
- Use final-size typography and inspect the rendered figure at target dimensions.
- Export editable SVG first, a PNG preview for QA, PDF when useful, and TIFF only when the
  target submission requires it.
- Export exact plotted source data and statistical results for quantitative panels.
- Use deterministic seeds for stochastic layout, clustering, jitter, or resampling.
- Keep plotting code specific to the actual analysis so another agent can audit every choice.
- Resolve every generated image export to an absolute path, verify it exists, and include
  each complete path in the final user response.

Selected-backend export minimum:

```python
# Python
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
fig.savefig("figure.svg", bbox_inches="tight")
fig.savefig("figure.pdf", bbox_inches="tight")
fig.savefig("figure.png", dpi=300, bbox_inches="tight")
```

```r
# R
svglite::svglite("figure.svg", width = width_mm / 25.4, height = height_mm / 25.4)
print(plot)
dev.off()
grDevices::cairo_pdf("figure.pdf", width = width_mm / 25.4, height = height_mm / 25.4)
print(plot)
dev.off()
ragg::agg_png("figure.png", width = width_mm, height = height_mm, units = "mm", res = 300)
print(plot)
dev.off()
```

Do not cross-render with the non-selected backend. If an export package is unavailable,
report or install that selected-backend dependency rather than substituting the other
language.
