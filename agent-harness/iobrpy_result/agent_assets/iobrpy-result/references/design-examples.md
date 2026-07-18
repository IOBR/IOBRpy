# Adaptive Design Examples

These examples illustrate decision-making, not prescribed layouts.

## One CIBERSORT Table, Different Questions

**Question: What dominates each tumor sample?**

Use a sample-level composition view as the hero. Add a compact fit-quality strip so poor-fit
samples are visibly qualified.

**Question: Which immune populations differ between responders and non-responders?**

Use raw group distributions and an effect-size panel as the hero. A composition overview is
supporting context, not the inferential centerpiece.

**Question: Are the apparent groups driven by global composition?**

Use an ordination or distance view with group overlays, then show the loadings or components
that explain separation. Do not begin with one box plot per cell type.

## Several TME Methods

**Weak design:** concatenate raw outputs into one heatmap and compare colors.

**Stronger design:** define a shared biological axis, standardize within each
method-feature, and show direction/rank concordance. Preserve method families and expose
disagreement rather than averaging it away.

**Publication architecture:** use the concordance view as an asymmetric hero panel. Place
one compact effect-size panel beside it and one fit/uncertainty strip below it. Do not give
three raw method heatmaps equal visual weight.

## TME Landscape Plus Focused Inference

**Question: What is the cohort landscape, and which populations differ by outcome?**

Use a wide sample composition, heatmap, or ordination panel to establish the landscape.
Follow with raw group distributions and effect intervals for a prespecified or reproducibly
selected subset. Add fit, uncertainty, or metadata coverage as a quiet context strip. This
is stronger than using the same relative-abundance matrix in both a stacked bar and an
absolute-value heatmap.

## Runall Or TME Profile Evidence

**Weak design:** create one equally sized panel for every output directory.

**Stronger design:** choose a `workflow evidence composite`. Let the decisive QC,
quantification, or downstream TME conclusion lead. Use a compact stage/completeness ribbon
for workflow context, then include only child-result panels that contribute unique evidence.
Workflow completion is not itself biological validation.

## BayesPrism

**Weak design:** plot theta alone.

**Stronger design:** let theta carry the biological pattern while theta CV controls visual
confidence through a companion panel, annotation, or selective emphasis. Discuss inferred
tumor expression separately because it is a different evidence type.

## TRUST4

**Question: Is repertoire diversity different between groups?**

Pair diversity metrics with read depth and top-clone concentration. If diversity changes
only with depth, make that limitation visually central rather than a footnote.

## HLA

**Question: What is the cohort call landscape?**

An allele-by-sample tile may expose individual genotypes; a frequency view may better serve
a cohort summary. Use both only when individual and cohort scales are both necessary.

## No Metadata

Do not manufacture group comparisons. Explore cohort structure, variance, correlations,
fit, uncertainty, and outliers, then state what metadata would unlock a stronger analysis.

## Final-Size Decision

For every example above, decide whether the figure is a single-column `89 mm` result or a
double-column `183 mm` composite before styling. At final size, remove redundant legends,
use direct labels when stable, and enlarge the hero evidence rather than shrinking every
panel equally.
