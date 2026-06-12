# IOBRpy Result Interpretation Coverage

This table is the complete coverage contract for `iobrpy-result`. It mirrors the 27 keys in
`src/iobrpy/RAG_MCP/iobrpy_required_params.json`.

| Function | Primary result or evidence | Interpretation guide |
|---|---|---|
| `runall` | end-to-end staged output directory | `interpretation-workflows.md` |
| `prepare_salmon` | cleaned/aggregated Salmon TPM-like matrix | `interpretation-preprocessing.md` |
| `count2tpm` | gene-level TPM matrix | `interpretation-preprocessing.md` |
| `anno_eset` | annotated and representative-deduplicated expression matrix | `interpretation-preprocessing.md` |
| `calculate_sig_score` | signature scores by sample | `interpretation-tme.md` |
| `cibersort` | cell estimates and fit statistics | `interpretation-tme.md` |
| `IPS` | IPS component and aggregate scores | `interpretation-tme.md` |
| `estimate` | immune, stromal, ESTIMATE, and optional purity scores | `interpretation-tme.md` |
| `mcpcounter` | immune/stromal population abundance scores | `interpretation-tme.md` |
| `quantiseq` | immune cell fractions and residual compartment | `interpretation-tme.md` |
| `epic` | immune/stromal fractions and residual compartment | `interpretation-tme.md` |
| `tme_cluster` | TME cluster labels and selected features | `interpretation-tme.md` |
| `LR_cal` | bulk ligand-receptor co-expression scores | `interpretation-tme.md` |
| `nmf` | NMF clusters, component-leading features, and PCA view | `interpretation-tme.md` |
| `mouse2human_eset` | human-ortholog expression matrix | `interpretation-preprocessing.md` |
| `batch_salmon` | per-sample `quant.sf` and logs | `interpretation-preprocessing.md` |
| `merge_salmon` | merged TPM and estimated-count matrices | `interpretation-preprocessing.md` |
| `merge_star_count` | merged STAR unstranded count matrix | `interpretation-preprocessing.md` |
| `batch_star_count` | BAM and `ReadsPerGene.out.tab` per sample | `interpretation-preprocessing.md` |
| `fastq_qc` | cleaned FASTQ, fastp JSON/HTML, and optional MultiQC report | `interpretation-preprocessing.md` |
| `log2_eset` | `log2(x + 1)` expression matrix | `interpretation-preprocessing.md` |
| `trust4` | repertoire reports, clone table, and immune indices | `interpretation-immune.md` |
| `spechla` | per-sample HLA allele calls | `interpretation-immune.md` |
| `extract_hla_read` | paired HLA-region FASTQ files | `interpretation-immune.md` |
| `hla_typing` | per-sample SpecHLA outputs and merged allele table | `interpretation-immune.md` |
| `tme_profile` | signature, six-method TME, merged, and LR outputs | `interpretation-workflows.md` |
| `bayesprism` | proportions, uncertainty, and inferred tumor expression | `interpretation-tme.md` |

## Explicit Exclusions

- `deside` is intentionally not covered.
- `ai` is intentionally not covered because it is an orchestration interface rather than a
  result-producing analysis function with its own biological output semantics.

When a result directory contains files from more than one row, interpret every evidenced
row. Do not infer a function solely from a directory name when a concrete output signature
is absent.
