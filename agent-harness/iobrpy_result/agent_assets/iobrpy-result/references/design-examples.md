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
