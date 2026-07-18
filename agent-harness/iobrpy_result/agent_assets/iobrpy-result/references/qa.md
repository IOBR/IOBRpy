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
- Confirm the figure contract states one core conclusion with an active verb, a declared
  archetype, target output, final physical dimensions, panel map, and reviewer risk.
- Ensure every panel has one distinct evidence role.
- Confirm the hero evidence has the clearest axis or dominant panel and that supporting
  panels are visually subordinate rather than uniformly sized.
- Confirm the figure responds to the supplied question and observed data. It should not be
  interchangeable with a figure produced from an unrelated IOBRpy result directory.
- Confirm at least one plausible alternative design was considered and record why the
  chosen encoding is less misleading or more informative.
- Check labels, units, sample counts, legends, and color consistency.
- Confirm group, method-family, and biological-compartment colors remain stable across
  panels; red/green is not the only distinction.
- Confirm lowercase panel labels, sans-serif typography, line weights, and legend placement
  are consistent at the declared final size.
- Prefer direct labels or one shared legend over repeated legends when identities are stable.
- Confirm small text remains readable at final dimensions.
- Open the PNG preview and inspect clipping, overlap, empty panels, and colorbar ranges.
- Confirm SVG text remains editable and the output contains no accidental rasterization.
- Confirm quantitative panels have exact plotted source data and statistical output.
- Confirm `n`, center, spread/interval, test, correction, and exact comparison are documented
  whenever statistical claims appear.
- If sources were integrated, confirm panel labels, legends, captions, and source-data
  tables retain explicit provenance.
- Confirm every generated image exists at the reported path.
- Confirm the final response lists each generated image with its normalized absolute path.
- Reject basename-only, relative-path-only, output-directory-only, or nonexistent image
  path reporting.

## Render-Inspect-Revise Loop

Apply this loop only for an explicit visualization request:

1. Export SVG and PNG from the selected backend at the declared final dimensions.
2. Open the PNG preview and inspect the whole page, then zoom into the densest panel.
3. Check conclusion hierarchy, panel balance, clipping, overlap, small text, legend travel,
   colorbar limits, empty space, and whether uncertainty/fit warnings are visually visible.
4. Revise the selected-backend code when any check fails, re-export, and reopen the preview.
5. Inspect SVG/PDF structure separately for editable text and accidental rasterization.
6. Stop only when the figure reads at final size and every applicable Figure check passes.

Do not treat file existence as visual QA. Do not use the non-selected backend to create a
preview or repair an export.

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
- Preserve a concise figure-contract record and QA note with the delivered artifacts for
  publication-oriented requests.
