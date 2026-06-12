# Workflow Wrapper Interpretation

Wrapper interpretation combines stage completeness, data provenance, and the semantics of
every concrete child output. Never replace child-result interpretation with a single
pipeline-level summary.

## `runall`

The expected branch depends on mode:

- Salmon: `fastq_qc -> batch_salmon -> merge_salmon -> prepare_salmon -> log2_eset`.
- STAR: `fastq_qc -> batch_star_count -> merge_star_count -> count2tpm -> log2_eset`.
- Both then continue through signature scoring, the six default TME methods, merged TME
  output, ligand-receptor scoring, and TRUST4.

Interpretation checklist:

- Identify the selected branch, index/reference, expected samples, and exact completed,
  partial, failed, skipped, or resumed stages.
- Require concrete files for each stage; do not infer full completion from the root directory
  or a final-stage file.
- The final `03-tpm/tpm_matrix.csv` is currently produced after `log2(x + 1)`. It is not a
  linear TPM matrix even though the directory name contains TPM.
- Downstream children do not all receive the same effective scale:
  - CIBERSORT and quanTIseq usually reverse the logged matrix with `2**x`, yielding
    approximately `TPM + 1`; both also use quantile normalization under current wrapper
    defaults.
  - IPS and signature scoring consume the logged matrix directly under their usual
    heuristics; MCPcounter computes marker means on that logged scale.
  - ESTIMATE is rank-based and is less sensitive to a monotone log transform.
  - EPIC is called without `--unlog`, so it CPM-scales logged values as if they were linear.
  - `LR_cal` applies another `log2(x + 1)`, so the wrapper LR output is based on
    `log2(log2(TPM + 1) + 1)`, not the standalone linear-TPM formula.
- Treat EPIC and LR results from the current wrapper as implementation-specific transformed
  outputs and make this provenance visible in biological interpretation and cross-method
  comparisons.
- The wrapper always adds quanTIseq `--arrays`, `--tumor`, and `--scale_mrna` defaults unless
  the implementation is changed; record these settings rather than assuming standalone
  defaults.
- Trace row/sample loss and transformations through QC, quantification, merge, identifier
  handling, log transformation, and downstream analyses.
- Interpret each detected signature, TME, LR, and repertoire result with its domain guide.
- State whether missing outputs reflect an intentional branch, optional standalone method,
  failed child stage, or incomplete evidence. BayesPrism and clustering are not implied by
  the default wrapper.

## `tme_profile`

Expected components are:

- `01-signatures/`: `calculate_sig_score`;
- `02-tme/`: CIBERSORT, IPS, ESTIMATE, MCPcounter, quanTIseq, EPIC, and merged output;
- `03-LR_cal/`: ligand-receptor scores.

Interpretation checklist:

- Confirm that the input was a gene-by-sample TPM-like matrix on the scale required by the
  child methods. The intended input is linear TPM; already logged input creates inconsistent
  child transformations, especially for EPIC and `LR_cal`.
- The wrapper enables quanTIseq array quantile normalization, tumor-gene removal, and mRNA
  scaling by default. CIBERSORT also defaults to quantile normalization.
- Audit expected versus observed samples and features for each child table.
- Interpret every produced method separately, then synthesize concordance, disagreement,
  missingness, residual compartments, and fit/uncertainty fields.
- `deconvo_merged.csv` is an outer join, not an independent algorithm or a validated composite
  immune score.
- The wrapper does not include BayesPrism, `tme_cluster`, or `nmf` by default. Do not imply
  those analyses ran without their own output evidence.
