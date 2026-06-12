# Preprocessing And Matrix Interpretation

These functions mainly establish technical quality and data provenance. Interpret whether
the output is complete and suitable for downstream use; do not force a biological claim.

## `fastq_qc`

- Evidence includes cleaned FASTQ files, per-sample `*_fastp.json`/HTML, task markers, and
  optional `multiqc_report/multiqc_fastp_report.html`.
- Compare before/after read counts, base quality, adapter content, N bases, duplication,
  insert-size patterns, and sample outliers. A completion marker alone is not a QC pass.
- The current command disables fastp length and quality filtering. It performs adapter/N-base
  processing but the output must not be described generically as "quality-filtered reads".
- N-base filtering remains active with `--n_base_limit 6`; distinguish reads removed for
  excessive N bases from disabled quality- and length-filtering.
- The current `length_required` argument is not passed to fastp, so do not attribute an
  observed length cutoff to that parameter.
- MultiQC summarizes technical metrics across samples; it does not establish biological
  comparability or absence of batch effects.

## `batch_salmon`

- Interpret each `quant.sf` as transcript-level abundance estimation under Salmon's model.
- `TPM` is relative abundance; `NumReads` is an estimated count and may be fractional. It is
  not a direct raw-alignment count.
- Record sample completeness, index and Salmon version compatibility, mapping/log warnings,
  and whether all expected samples have non-empty `quant.sf` plus completion evidence.
- The current workflow uses paired-end `ISF`, `--gcBias`, and `--validateMappings`; library
  orientation mismatch can bias estimates.

## `merge_salmon`

- The function recursively finds `quant.sf`, derives sample names from parent directories,
  and outer-merges `Name`, `TPM`, and `NumReads`.
- Expected outputs are `<project>_salmon_tpm.tsv.gz` and
  `<project>_salmon_count.tsv.gz`.
- Report sample-name collisions, missing transcripts, NaN introduced by the outer merge,
  and differences in detected transcript sets.
- Do not call the NumReads matrix integer counts or the TPM matrix gene-level until its
  identifiers and aggregation state are verified.

## `prepare_salmon`

- The function parses composite `Name` annotation fields and retains ENST, ENSG, or symbol.
- Duplicate retained identifiers are averaged across numeric values, not summed.
- Because duplicate aggregation uses a mean, postprocessed sample totals may no longer sum
  to one million. Describe the result as a cleaned/aggregated TPM-like matrix rather than
  untouched Salmon TPM.
- Report malformed names, missing annotations, duplicate-collapse rate, selected identifier,
  version stripping, and row loss.

## `batch_star_count`

- Evidence includes coordinate-sorted BAM files, STAR logs, and
  `*_ReadsPerGene.out.tab`.
- The current workflow uses two-pass STAR, sorted BAM output, and GeneCounts. Interpret
  mapping rate, uniquely/multi-mapped reads, splice metrics, sample outliers, and file
  completeness before gene counts.
- Confirm that the STAR index, genome build, annotation, paired-read suffixes, and actual
  compression match the inputs.
- A BAM file proves alignment output exists; it does not by itself prove acceptable mapping
  quality.

## `merge_star_count`

- The function always reads column 1 from `ReadsPerGene.out.tab`, which is the unstranded
  count column.
- The merged file is `<project>.STAR.count.tsv.gz`.
- The first STAR summary rows such as `N_unmapped` are retained by the current implementation.
  Identify and exclude them before treating rows as genes.
- For stranded libraries, verify whether column 1 is appropriate; the current merge does not
  select the forward/reverse stranded columns.
- Report missing samples, duplicated sample names, count depth, and large library-size
  differences. These are counts, not TPM.

## `count2tpm`

- The function converts counts to RPK using effective length, then normalizes to TPM.
- Report effective-length source, organism, ID type, version stripping, unmatched genes,
  invalid lengths, and row loss.
- Duplicate symbols are resolved by retaining the row with the highest selected row score;
  duplicate rows are not summed.
- TPM totals are one million before representative-symbol deduplication. Totals may be lower
  afterward, so audit totals at both stages when possible.
- The output is normalized relative abundance and is not appropriate for raw-count
  differential-expression models.

## `anno_eset`

- Report annotation match rate, unmatched IDs, rows removed as invalid, and the number of
  symbols affected by duplicates.
- `mean`, `sd`, and `sum` choose a row-level score used to retain one representative row per
  symbol; they do not aggregate all duplicated probes into one expression profile.
- Interpretation must preserve the selected annotation resource, ID type, version stripping,
  and representative-selection method because each can alter downstream features.

## `mouse2human_eset`

- Treat the output as an ortholog-mapped expression matrix, not as direct human measurement.
- Report mapping coverage, one-to-many/many-to-one collapse, duplicate handling, unmatched
  genes, and row loss.
- Human ortholog mapping supports cross-species feature alignment but does not establish
  identical regulation, cell composition, or biological function.

## `log2_eset`

- The output scale is `log2(x + 1)` and must not be described as counts or linear TPM.
- Zero remains zero; the pseudocount has the greatest relative effect on low values.
- Nonnumeric values are coerced to missing. Values below `-1` are rejected; exactly `-1`
  produces negative infinity, while values between `-1` and zero remain finite but require
  explicit justification.
- Report input/output range, missingness introduced, and whether downstream methods expect
  linear or logged expression.
