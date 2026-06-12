# IOBRpy Result Layout

Use concrete file evidence. The same method may also be run as a standalone command,
so accept explicit files outside these standard directories.

File patterns in this guide identify candidates, not provenance by themselves. Apply
`source-classification.md` and verify content before calling a candidate an IOBRpy result.

## Standard workflows

| Workflow | Stage | Typical result |
|---|---|---|
| `runall` | `03-tpm/` | `tpm_matrix.csv` |
| `runall` | `04-signatures/` | `calculate_sig_score.csv` |
| `runall` | `05-tme/` | six method tables and `deconvo_merged.csv` |
| `runall` | `06-LR_cal/` | `lr_cal.csv` |
| `runall` | `07-TCRBCR/` | TRUST4 reports and immune summaries |
| `tme_profile` | `01-signatures/` | `calculate_sig_score.csv` |
| `tme_profile` | `02-tme/` | six method tables and `deconvo_merged.csv` |
| `tme_profile` | `03-LR_cal/` | `lr_cal.csv` |

The `runall` TPM file is currently produced after `log2(x + 1)` transformation. Confirm
the actual input scale before interpreting or reusing it.

## Result signatures

| File pattern | Result class |
|---|---|
| `*_fastp.json`, `*_fastp.html` | per-sample fastp QC metrics and report |
| `multiqc_fastp_report.html` | cross-sample fastp/MultiQC summary |
| cleaned `*.fastq.gz` plus task marker | `fastq_qc` cleaned reads and execution evidence |
| `quant.sf` | per-sample Salmon transcript estimates |
| `*_salmon_tpm.tsv.gz` | merged Salmon transcript TPM matrix |
| `*_salmon_count.tsv.gz` | merged Salmon estimated NumReads matrix |
| user-named `prepare_salmon` output | cleaned/aggregated TPM-like matrix |
| `*_Aligned.sortedByCoord.out.bam` | per-sample STAR aligned reads |
| `*_ReadsPerGene.out.tab` | STAR per-sample unstranded/stranded count columns |
| `*.STAR.count.tsv.gz` | merged STAR column-1 unstranded count matrix |
| user-named `count2tpm` output | gene-level TPM matrix |
| user-named `anno_eset` output | annotated representative-deduplicated matrix |
| user-named `mouse2human_eset` output | human-ortholog expression matrix |
| user-named `log2_eset` output | `log2(x + 1)` expression matrix |
| `calculate_sig_score.csv` | pathway and biological signature scores |
| `cibersort_results.csv` | CIBERSORT fractions plus fit statistics |
| `quantiseq_results.csv` | quanTIseq fractions, often including `Other` |
| `epic_results.csv` | EPIC fractions, optionally including `otherCells` |
| `estimate_results.csv` | stromal, immune, ESTIMATE, and tumor-purity scores |
| `mcpcounter_results.csv` | population abundance scores, not fractions |
| `IPS_results.csv` | immunophenoscore components and total IPS |
| `deconvo_merged.csv` | outer join of produced TME result tables by `ID` |
| `theta.csv` | BayesPrism cell-type proportions |
| `theta_cv.csv` | BayesPrism coefficient of variation for proportions |
| `Z_tumor.csv` | BayesPrism tumor expression estimates |
| `lr_cal.csv` | ligand-receptor pair scores by sample |
| `clusters.csv` or `tme_cluster.csv` | sample cluster labels and used features |
| `top_features_per_cluster.csv` | NMF component-leading features |
| `trust4_immdata.csv` | clone-level repertoire table |
| `trust4_immune_indices.csv` | per-sample repertoire summary metrics |
| extracted paired HLA FASTQ | HLA-region reads for typing, not allele calls |
| per-sample `hla.result.txt` | SpecHLA allele calls |
| `hla_result_merged.txt` | merged SpecHLA allele calls |

## Table handling

- In mixed directories, use only `function_detections[].matched_by_content` as the default
  confirmed-native path set. Do not include every sibling file from the stage directory.
- Treat `ID`, `Sample`, and `sample` as candidate sample-ID columns.
- Result tables are usually samples by features; expression matrices are usually genes by
  samples. Confirm orientation from headers and values instead of guessing.
- `deconvo_merged.csv` is an outer join. Missing cells can mean a method lacked a sample,
  not that the biological value is zero.
- For BayesPrism, inventory `theta.csv`, `theta_cv.csv`, and `Z_tumor.csv` separately.
  A single detected companion does not establish complete execution.
- STAR `ReadsPerGene.out.tab` contains summary pseudo-rows before gene rows, and
  `merge_star_count` currently preserves them.
- `prepare_salmon` averages duplicate retained identifiers; `anno_eset` retains one
  representative duplicate rather than aggregating all rows.
- Preserve method suffixes such as `_CIBERSORT`, `_quantiseq`, `_EPIC`, `_estimate`,
  `_MCPcounter`, and `_IPS` in source data. Strip them only for display labels.
- Keep a provenance column when native, compatible, and external tables are intentionally
  integrated.
- Record `provenance_status` and `completion_status` separately. A file may be confirmed
  native while its method-level companion set remains partial.
