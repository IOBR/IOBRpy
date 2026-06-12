# Adaptive Visualization Grammar

Use this file to reason about visual encodings. It is not a lookup table that assigns one
mandatory chart to each IOBRpy method.

## Backend Exclusivity

Use the backend selected by the user for all visual work. Python figures must be generated,
previewed, exported, and revised in Python. R figures must be generated, previewed,
exported, and revised in R. Data inspection may use ordinary shell tools, but do not invoke
the non-selected language's graphics stack.

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

Use stable biological roles rather than assigning colors from file order:

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

## Statistics

- Show individual samples where feasible.
- Respect pairing and repeated measures.
- Report group `n`, effect size, interval, test, and multiple-testing correction.
- Use permutation or rank-based approaches when assumptions are weak and sample size allows.
- Treat exploratory screens as exploratory.
- Do not run group tests without supplied or verifiable metadata.

## Publication Output

- Keep SVG text editable with `svg.fonttype = "none"`.
- Use final-size typography and inspect the rendered figure at target dimensions.
- Export a PNG preview and exact source data for quantitative panels.
- Use deterministic seeds for stochastic layout, clustering, jitter, or resampling.
- Keep plotting code specific to the actual analysis so another agent can audit every choice.
- Resolve every generated image export to an absolute path, verify it exists, and include
  each complete path in the final user response.
