# TME, Signature, LR, And Clustering Interpretation

These outputs are estimates derived from bulk expression. They are not direct single-cell
measurements.

The current modules are Python implementations in the IOBRpy codebase. Use the published
method name for scientific context, but do not imply bit-for-bit parity with an official R
package, web service, or reference implementation.

## `calculate_sig_score`

- Compare samples within the same signature and scoring method.
- Preprocessing conditionally applies a distribution-based `log2(x + 1)` heuristic only
  when the matrix has fewer than 5000 sample columns. `--adjust_eset` then removes
  nonnumeric columns and genes with nonfinite values or zero sample SD; it is filtering, not
  a general normalization step.
- PCA standardizes each matched gene and aligns PC1 so that it is non-negatively correlated
  with mean signature expression when that correlation is defined. Higher values therefore
  follow the mean-expression direction in the analyzed cohort, but the scale remains
  cohort-dependent and is not an absolute activation unit.
- Despite its name, the current `zscore` method does not calculate z-scores. It returns the
  arithmetic mean expression of matched signature genes after preprocessing. Describe it
  as a mean-expression signature score.
- ssGSEA is rank-based relative enrichment. Standalone ssGSEA requires at least five
  matched genes, while PCA and `zscore` require at least two.
- `integration` emits suffixed PCA, mean-expression `zscore`, and ssGSEA views when
  available. If `gseapy` is unavailable or no signature survives filtering, the current
  function can return only PCA and `zscore`; infer executed methods from actual suffixes.
- `TMEscore_CIR` and `TMEscore_plus` are A-minus-B contrasts, not independent signatures.
- Signature-name collisions with different gene sets are renamed with a `__<group>` suffix.
- Report signature collection, matched-gene counts, minimum-gene filtering, adjustment, and
  missing signatures.

## `cibersort`

- This is a CIBERSORT-like NuSVR implementation. Relative-mode cell columns are normalized
  to sum to one; absolute modes are method-derived scales, not cell counts.
- If the input maximum is below 50, the implementation applies `2**x` without subtracting
  one. This usually reverses log2-like input, but can misclassify a low-valued linear matrix.
  Report the input range and whether the heuristic fired.
- Quantile normalization defaults to enabled even though it is often disabled for RNA-seq.
  Report `QN` explicitly.
- Duplicate gene labels are made unique with suffixes rather than aggregated; only exact
  labels shared with LM22 enter the fit.
- Inspect `P-value_CIBERSORT`, `Correlation_CIBERSORT`, and `RMSE_CIBERSORT` when present.
- P-values use `count(null correlation >= observed) / permutations`, so zero is possible and
  resolution is limited by the permutation count. With `perm <= 0`, the output sentinel is
  `9999`, not a valid p-value.
- Poor fit weakens sample-level interpretation. Do not filter by a threshold unless it is
  stated and the number of excluded samples is reported.

## `quantiseq`

- Gene aliases are mapped to HGNC symbols and duplicate mapped genes are collapsed by
  median. If the input maximum is below 50, the implementation applies `2**x`; it then
  optionally quantile-normalizes array-mode input and always rescales each sample to one
  million.
- Final rows are normalized to sum to one. Values are compositional estimated fractions, not
  absolute cell counts.
- `Other` exists only for the LSEI path. It is calculated as a nonnegative residual before
  final row normalization, so interpret it as an implementation-specific unassigned
  compartment rather than direct tumor purity.
- Interpretation depends on tumor-gene removal, array handling, mRNA scaling, removed-gene
  rules, solver, signature coverage, and the low-Treg correction.
- Audit the input-scale heuristic, matched signature genes, row sums, and residual size
  before percentage displays.

## `epic`

- Values are cell fractions under the selected reference (`TRef`, `BRef`, or both). The
  default path uses weighted NNLS and CPM-like scaling; duplicate genes are collapsed by
  median.
- There is no automatic log-scale detection. Logged input must be explicitly back-transformed
  with `--unlog`; otherwise logged values are CPM-scaled as if they were linear expression.
- `otherCells` begins as a residual mRNA proportion, is adjusted by per-cell mRNA content,
  and is then renormalized with the other cell fractions. It is not automatically tumor
  purity.
- The API computes fit diagnostics internally, but the standard CLI writes only
  `cellFractions`. Do not claim that EPIC fit quality was audited unless a separate API-level
  result or diagnostic artifact is available.
- Report reference, solver path, scaling, unlog setting, signature coverage, fraction sums,
  and whether `otherCells` was enabled.

## `mcpcounter`

- The current implementation computes the arithmetic mean expression of matched marker
  genes for each population. Treat it as an IOBRpy marker-mean abundance score and do not
  assume numerical parity with another MCP-counter implementation.
- Gene IDs are version-stripped, uppercased, and duplicate rows are collapsed by per-sample
  maximum before marker matching.
- They are not constrained fractions and must not be stacked or read as percentages.
- Scores are strongly input-scale dependent. Compare within a population feature and
  consistently processed cohort; raw scales across populations, transformations, or studies
  are not directly interchangeable.
- Report feature-ID mode and matched marker count per population when logs or reference data
  make that information available.

## `estimate`

- `ImmuneSignature_estimate` and `StromalSignature_estimate` are expression-derived scores;
  `ESTIMATEScore_estimate` combines immune and stromal context.
- The implementation rank-normalizes genes within each sample before enrichment. Scores are
  relative to the retained common-gene universe; report common-gene overlap and do not read
  score differences as linear expression differences.
- The current implementation emits `TumorPurity_estimate` only for `affymetrix`; agilent and
  illumina runs return score rows without purity.
- Purity is modeled from ESTIMATEScore with a cosine transform. Negative transformed values
  become missing; treat missing/out-of-range values as quality issues.

## `IPS`

- `MHC_IPS`, `EC_IPS`, `SC_IPS`, and `CP_IPS` are component scores; `AZ_IPS` is their
  aggregate before transformation and `IPS_IPS` is the final score.
- Input is log2-transformed only when the matrix maximum exceeds 100. Record the input range;
  a low-valued linear matrix can otherwise be treated as if it were already logged.
- Each matched IPS gene is standardized within its sample against the mean and SD of all
  genes in that sample. Missing IPS genes are dropped, so report gene and group coverage.
- Components are assembled from positional slices of the ordered signature groups. If an
  entire expected IPS group is absent, later groups can shift between component slices and
  the result should be treated as structurally unreliable.
- The final value is `round(AZ * 10 / 3)` converted to an integer. The implementation does
  not clip to 0-10, so do not impose a canonical range that is absent from the code.
- Higher IPS can indicate a more immunogenic expression context. It is not a measured
  treatment response or validated prediction by itself.

## `bayesprism`

- `theta.csv` contains estimated proportions and `theta_cv.csv` contains uncertainty.
- The current entrypoint casts both the single-cell reference and transposed bulk matrix to
  `int32`, then runs BayesPrism with `input_type="count.matrix"`. Decimal TPM-like values are
  truncated toward zero before modeling. Report the original scale, noninteger fraction,
  and the amount of zero inflation introduced by this cast when the input is available.
- Do not describe the current run as preserving continuous TPM values. Logged input is also
  inappropriate for this count-matrix path because it is integer-truncated rather than
  back-transformed.
- High coefficient of variation must visibly temper sample- or cell-type conclusions.
- `Z_tumor.csv` is model-inferred tumor expression, not directly isolated tumor-cell data.
- Treat the expected companion set as `theta.csv`, `theta_cv.csv`, and `Z_tumor.csv`.
  Detection of one file confirms only that artifact, not complete BayesPrism output.
- Report reference source, tumor key, matched genes, proportion sums, uncertainty, and
  sensitivity to rare or closely related reference states.

## `LR_cal`

- For linear TPM input, each available pair is scored as the minimum of ligand and receptor
  `log2(TPM + 1)` expression. A high score means both genes are sufficiently expressed in
  the bulk sample.
- Use "co-expression-compatible LR signal" or "LR score", not proven interaction.
- The score does not identify expressing cell types, signaling direction, binding, downstream
  activation, or causality.
- Before scoring, the implementation removes genes with any missing/nonfinite value,
  nonnumeric values, and, for multi-sample input, zero variance. A pair is then dropped if
  either gene is unavailable in any sample. Report gene and pair-set loss before comparing
  cohorts or datasets.

## `tme_cluster`

- The function applies k-means-style clustering to selected, optionally standardized
  features and chooses `k` by the KL index over neighboring within-cluster sums.
- Scaling is enabled by default. With scaling, features containing any missing/nonfinite
  value or zero SD are removed; samples are not imputed. With `--no-scale`, this filtering
  path is bypassed and invalid values can make the clustering unusable.
- The output table contains the feature values actually clustered, which are z-scored when
  scaling was enabled rather than the original input values.
- Labels such as `TME1` are arbitrary identifiers.
- Report selected features, scaling, candidate range, chosen `k`, cluster sizes, missing
  samples, and cluster-defining feature patterns.
- KL scores, within-cluster sums, and the chosen-k message are not saved in the result table.
  Require captured logs or a rerun for those diagnostics instead of reconstructing them
  from cluster labels.
- KL selection does not establish stability. Use resampling, sensitivity checks, or external
  validation before naming a subtype.

## `nmf`

- NMF requires nonnegative input. Record log1p, row normalization, shifts, and feature
  selection because they change the factorization's meaning.
- For each successful `k`, labels are assigned from the largest W loading and silhouette is
  computed in W space; the maximum silhouette is selected, with smaller `k` winning ties.
- The default behavior can include `k=2`; exclude it only when `--skip_k_2` was used.
- `clusters.csv` contains labels and processed features,
  `top_features_per_cluster.csv` reports high H loadings, and `pca_plot.png` is a descriptive
  PCA view of W.
- `top_features_per_cluster.csv` stores feature names but not their H loading values.
  Silhouette and reconstruction-error values are printed during execution but are not saved
  in the standard result files. Do not claim they were audited without logs.
- `pca_plot.png` is a PCA projection of W and does not report explained-variance percentages.
  Ellipses are descriptive covariance summaries, not inferential confidence regions.
- Missing/nonfinite values are not imputed. A requested shift changes all values by
  `abs(minimum) + shift`, so record the resulting minimum and not only the user-supplied
  shift argument.
- Cluster labels and top features are exploratory. Report reconstruction error, silhouette,
  cluster sizes, seeds/stability, and phenotype associations separately.

## Cross-Method TME Interpretation

- Prefer concordance in direction, rank, or group separation.
- Do not average raw CIBERSORT, quanTIseq, EPIC, MCPcounter, ESTIMATE, and IPS values into a
  single "immune score" without a defined model.
- Differences can arise from signatures, references, composition constraints, scaling, and
  tumor-content assumptions rather than biology alone.
