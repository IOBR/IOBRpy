# IOBRpy Interpretation Router

Interpretation is required for every detected IOBRpy result-producing function, not only
for final TME tables. Use `function-coverage.md` to route each result to one domain guide:

- `interpretation-preprocessing.md`: QC, Salmon/STAR, annotation, ortholog conversion,
  count-to-TPM, and log transformation.
- `interpretation-tme.md`: signatures, deconvolution, score-based TME methods,
  ligand-receptor analysis, and clustering.
- `interpretation-immune.md`: TRUST4, HLA-read extraction, SpecHLA, and batch HLA typing.
- `interpretation-workflows.md`: `runall` and `tme_profile` wrappers.

Before routing, apply `source-classification.md`. Method-aware interpretation does not
override provenance: an external or ambiguous table remains external or ambiguous even when
its columns resemble an IOBRpy-supported analysis family.

## Required Structure

For each function, report:

1. **Detected evidence**: exact files, sample count, shape, completion markers, and relevant
   parameters or modes.
2. **Observed result**: concrete patterns, missing outputs, range violations, or quality
   signals visible in the supplied evidence.
3. **Method-aware meaning**: what the function computed and the scale of its output.
4. **Bounded interpretation**: biological or technical meaning supported by the result.
5. **Limitations and next evidence**: assumptions, confounding, uncertainty, and the most
   useful validation when relevant.

Do not manufacture a biological narrative for a technical preprocessing function. Its
interpretation should instead explain data quality, transformation provenance, information
loss, and suitability for the next stage.

## Cross-Function Rules

- Preserve the distinction between raw counts, estimated counts, TPM, fractions, relative
  abundance scores, transformed expression, enrichment scores, and genotype calls.
- Reconstruct the effective transformation chain from the concrete entrypoint and flags.
  A filename containing `TPM` is not proof that the values remain linear TPM.
- Treat method names as method-family labels until the current IOBRpy implementation has
  been checked. State implementation-specific preprocessing, solvers, defaults, and missing
  output fields when they materially affect interpretation.
- Do not import canonical score ranges, fit statistics, or validation thresholds from an
  external implementation unless the IOBRpy output and code support them.
- State all transformations and duplicate-resolution rules that materially change meaning.
- Distinguish a completed wrapper from concrete evidence that every child stage succeeded.
- Distinguish native provenance from result-set completeness. For multi-file functions,
  enumerate expected and observed companions before calling the function complete.
- Treat `.done`, `task.complete`, and similar markers as execution evidence, not biological
  evidence.
- Do not compare raw magnitudes across different deconvolution or scoring methods.
- Prefer cross-method agreement in direction, rank, or cohort separation.
- Do not silently normalize, impute, filter, collapse alleles, or discard missing samples.

## Evidence Language

Prefer:

- "was higher in the observed cohort";
- "was associated with";
- "the method estimated";
- "is consistent with";
- "supports the hypothesis that";
- "the transformation retained/lost";
- "the available files indicate that the stage completed".

Avoid without additional evidence:

- "caused";
- "activated" from a score alone;
- "increased immune-cell number" from a relative fraction;
- "demonstrated cell-cell communication" from ligand-receptor scores;
- "predicts response" without a validated outcome model;
- "passed QC" when only output files or completion markers are present.
