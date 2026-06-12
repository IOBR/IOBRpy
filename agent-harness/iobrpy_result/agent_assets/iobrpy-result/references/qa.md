# Delivery QA

## Data

- Confirm every candidate file has a provenance-ledger row.
- Confirm only content-verified paths are labeled as native IOBRpy results.
- Confirm likely, reusable, external, metadata, and unknown files are not silently promoted
  to native IOBRpy evidence.
- Confirm sibling files were not included solely because they share a directory.
- Confirm provenance status and method-level completion status were assessed separately.
- For multi-file functions, compare the expected companion set with the observed files and
  report partial output explicitly.
- Record every input path, delimiter, table shape, and sample-ID column.
- Verify sample joins and report unmatched IDs.
- Confirm whether values are fractions, abundance scores, transformed expression, or
  method-specific scores.
- Preserve raw files and write transformed plotting data separately.
- Make all feature filtering and ordering reproducible in code.

## Figure

- Apply this section only when the user explicitly requested creation or modification of a
  visualization.
- Confirm interpretation-only requests did not trigger a Python/R question, plotting code,
  preview, or figure artifact.
- Confirm the user explicitly selected Python or R before any requested figure work began.
- Confirm the selected backend produced the figure, preview, exports, and visual revisions;
  no graphics fallback from the other language was used.
- Ensure every panel has one distinct evidence role.
- Confirm the figure responds to the supplied question and observed data. It should not be
  interchangeable with a figure produced from an unrelated IOBRpy result directory.
- Confirm at least one plausible alternative design was considered and record why the
  chosen encoding is less misleading or more informative.
- Check labels, units, sample counts, legends, and color consistency.
- Confirm small text remains readable at final dimensions.
- Open the PNG preview and inspect clipping, overlap, empty panels, and colorbar ranges.
- Confirm SVG text remains editable and the output contains no accidental rasterization.
- If sources were integrated, confirm panel labels, legends, captions, and source-data
  tables retain explicit provenance.
- Confirm every generated image exists at the reported path.
- Confirm the final response lists each generated image with its normalized absolute path.
- Reject basename-only, relative-path-only, output-directory-only, or nonexistent image
  path reporting.

## Interpretation

- Name the IOBRpy method and the exact result file.
- State the source class for every non-native or uncertain result that is discussed.
- Tie every substantive sentence to a visible pattern, a reported statistic, or a documented
  method property.
- Distinguish observed values from method-aware interpretation and biological hypotheses.
- Report fit or uncertainty columns when available.
- Do not claim clustering selection diagnostics, BayesPrism uncertainty, or HLA schema
  consistency from a final table when the required logs or companion artifacts are absent.
- State missing metadata, untested confounders, compositional constraints, and scale limits.
- Avoid cross-method raw-value comparisons and causal or clinical claims unsupported by the
  supplied data.

## Reproducibility

- Deliver the plotting script and source-data CSV with the figure.
- Use fixed random seeds for clustering, jitter, label placement, or resampling.
- Include package versions when statistical results or layout depend on them.
- Do not overwrite original IOBRpy results.
