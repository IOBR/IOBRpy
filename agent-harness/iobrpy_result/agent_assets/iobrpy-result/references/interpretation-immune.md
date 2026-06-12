# Repertoire And HLA Interpretation

## `trust4`

- Interpret raw `*_report.tsv` evidence together with `trust4_immdata.csv` and
  `trust4_immune_indices.csv`.
- `Nreads` is repertoire evidence depth; `Nclones` is the number of detected non-zero clones.
- In the current summary code, `Nclones` counts positive report rows; it does not collapse
  identical CDR3 sequences into unique clonotypes.
- `Length_CDR3` is read-weighted. Shannon and Gini-Simpson increase with diversity;
  Evenness describes balance; Gini and top-clone frequency increase with concentration.
- Evenness divides Shannon by `log(number of positive clones)` and is therefore undefined
  when only one positive clone is present. Rare-clone fields can be zero because zero-count
  report rows are not removed before sorting.
- Metrics are grouped only by sample. Different receptor chains, loci, and productive states
  present in a report are pooled unless the input reports were filtered beforehand. Do not
  interpret the summary as chain-specific diversity.
- Missing `CDR3_dna` values are converted to the string `nan` before length calculation,
  which contributes a length of three when its read count is positive. Audit missing CDR3
  rows before using `Length_CDR3`.
- Compare richness/diversity only with sequencing depth, sample type, chain, and productive
  status in view. Low depth can mimic low diversity or high clonality.
- Report missing chains/samples, clone filtering, dominant clonotypes, and whether metrics
  are read- or clone-weighted. Repertoire concentration is not equivalent to antigen
  specificity.

## `extract_hla_read`

- This is a preprocessing result, not an HLA genotype.
- Interpret paired HLA-region FASTQ existence, read yield, pair integrity, extraction
  completeness, BAM/CRAM indexing, and reference-build match (`hg19` versus `hg38`).
- Low yield or missing mates weakens downstream typing confidence. Do not infer alleles from
  extracted-read files alone.

## `spechla`

- Preserve per-locus two-allele calls, allele resolution, no-calls, ambiguity, and the
  reported SpecHLA version.
- Major merged fields can include HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, and DRB1 allele pairs.
- Record whether exon/RNA mode (`use_exon=1`) or full-length/WGS mode (`use_exon=0`) was used.
- The standalone wrapper merger treats the first physical nonempty line of each
  `hla.result.txt` as the header. If the file starts with a version comment, the true header
  can be copied as a data row. Inspect per-sample files before trusting a standalone merged
  table.
- An allele call is not evidence of expression, antigen presentation, peptide binding,
  immune escape, transplant compatibility, or clinical response.

## `hla_typing`

- Interpret the wrapper at three levels: extraction success, per-sample SpecHLA completion,
  and integrity of `hla_result_merged.txt`.
- Report expected versus typed samples, missing loci, no-calls, duplicate samples, mixed
  headers or versions, reference build, and exon/full-length mode.
- The batch merger writes the first available version line and header once, then appends
  later data rows without validating later headers or versions. A well-formed merged file
  therefore does not prove cross-sample schema or database-version consistency.
- A SpecHLA done marker requires a row with the expected sample ID, but does not require
  complete or high-confidence allele calls at every locus.
- Completion markers support resume and execution auditing but do not replace inspection of
  allele calls and missingness.
- Cohort allele-frequency summaries require appropriate ancestry/population context and
  should retain ambiguity rather than silently collapsing calls.
