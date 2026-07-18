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

## File Understanding Gate

Before applying a method guide, open the actual result file and inspect:

- delimiter, headers, table shape, and row/column orientation;
- sample and feature identifiers, duplicates, and companion-file consistency;
- numeric versus label columns, representative values, ranges, missingness, and nonfinite
  values;
- method-specific quality, fit, uncertainty, depth, residual, or completion fields that are
  actually present.

A filename or directory name identifies only a candidate result. Do not explain fields,
claim completion, or apply biological semantics until the content is compatible with the
method. If the values cannot be read, state that the file was identified but not
value-interpreted.

## Default User-Facing Interpretation

For each function, report:

1. **Identified result**: exact file, IOBRpy method, sample count, shape, orientation, and
   completion evidence.
2. **Field meanings**: explain important columns, units, score direction, and quality or
   uncertainty fields present in the file.
3. **Key observations**: report concrete patterns, ranges, missing outputs, outliers, range
   violations, or group differences visible in the supplied values.
4. **Method-aware meaning**: explain the supported biological or technical meaning on the
   method's actual scale.
5. **Necessary limitations**: state the main assumption, confounder, uncertainty, or scale
   constraint needed to avoid overclaiming.
6. **Optional visualization handoff**: when useful, name the finding, source fields,
   comparison unit, quality overlay, and suitable visual grammar. Do not create a figure or
   trigger the Python/R backend gate unless visualization was explicitly requested.

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
