# IOBRpy CLI Index

This page lists all available CLI commands with one‑line summaries and quick links.
For full parameter details, see **[Commands & common options](./commands-and-options.md)**.
For an example end‑to‑end run and expected folder layout, see **[Typical end‑to‑end workflow](./workflow.md)**.

---

## Core pipeline
- **runall** — Orchestrate QC → quantification (Salmon/STAR) → TPM → signatures → deconvolution → clustering → LR scoring. → [details](./commands-and-options.md#runall--orchestrate-the-full-pipeline)

## Preprocessing
- **fastq_qc** — FASTQ QC/trimming via fastp; optional MultiQC report. → [details](./commands-and-options.md#fastq_qc)

## Quantification — Salmon line
- **batch_salmon** — Batch `salmon quant` on paired‑end FASTQs. → [details](./commands-and-options.md#batch_salmon)
- **merge_salmon** — Merge per‑sample `quant.sf` to project‑level TPM/NumReads. → [details](./commands-and-options.md#merge_salmon)
- **prepare_salmon** — Clean combined Salmon TPMs (ID mapping, version stripping). → [details](./commands-and-options.md#prepare_salmon)

## Quantification — STAR line
- **batch_star_count** — STAR alignment with `--quantMode GeneCounts`. → [details](./commands-and-options.md#batch_star_count)
- **merge_star_count** — Merge multiple `*_ReadsPerGene.out.tab` into one count matrix. → [details](./commands-and-options.md#merge_star_count)
- **count2tpm** — Convert counts to TPM (supports Ensembl/Entrez/Symbol/MGI). → [details](./commands-and-options.md#count2tpm)

## Expression utilities
- **log2_eset** — Apply `log2(x+1)` to an expression matrix. → [details](./commands-and-options.md#log2_eset)
- **anno_eset** — Annotate / de‑duplicate an expression matrix. → [details](./commands-and-options.md#anno_eset)
- **mouse2human_eset** — Map mouse symbols to human symbols (matrix/table modes). → [details](./commands-and-options.md#mouse2human_eset)

## Signature scoring
- **calculate_sig_score** — Compute pathway/signature scores (`pca`/`zscore`/`ssgsea`/`integration`). → [details](./commands-and-options.md#calculate_sig_score)
- **IPS** — Immunophenoscore (subscores + total). → [details](./commands-and-options.md#ips)

## Immune deconvolution
- **cibersort** — CIBERSORT wrapper (permutations, QN, absolute mode). → [details](./commands-and-options.md#cibersort)
- **quantiseq** — quanTIseq fractions with robust options. → [details](./commands-and-options.md#quantiseq)
- **epic** — EPIC fractions (TRef/BRef). → [details](./commands-and-options.md#epic)
- **estimate** — ESTIMATE stromal/immune/purity scores. → [details](./commands-and-options.md#estimate)
- **mcpcounter** — MCPcounter infiltration scores. → [details](./commands-and-options.md#mcpcounter)

## Clustering / decomposition
- **tme_cluster** — k‑means with KL index auto‑k; feature selection. → [details](./commands-and-options.md#tme_cluster)
- **nmf** — Non‑negative Matrix Factorization (auto‑k, optional skip k=2). → [details](./commands-and-options.md#nmf)
- **deside** — Deep learning deconvolution; requires pre‑downloaded model; runs on Python 3.9. → [details](./commands-and-options.md#deside)

## Ligand–receptor
- **LR_cal** — Ligand–receptor interaction scoring. → [details](./commands-and-options.md#lr_cal)

---

### Quick links
- **End‑to‑end workflow & file layout:** [workflow.md](./workflow.md)
- **All commands & options:** [commands-and-options.md](./commands-and-options.md)