---
title: Commands & common options
layout: default
nav_order: 3
---

# Commands & common options

## `runall` — From FASTQ to TME
- **runall**
  - `--mode {salmon|star}` (required)
  - `--outdir <DIR>` (required): root output directory
  - `--fastq <DIR>` (required): forwarded to `fastq_qc --path1_fastq`
  - `--threads <INT>` (per-block): CPU/concurrency control set via block-level flags (e.g., fastq_qc --num_threads, batch_salmon --num_threads, batch_star_count --num_threads, merge_salmon --num_processes, cibersort --threads, calculate_sig_score --parallel_size).
  - `--batch_size <INT>` (per-block): batching size set via block-level flags (e.g., fastq_qc --batch_size, batch_salmon --batch_size, batch_star_count --batch_size).
  - `--resume`: skip steps if outputs already exist
  - `--dry_run`: print planned commands without executing

## From FASTQ through FASTQ Quality Control and Salmon/STAR to TPM
- **fastq_qc**
  - `--path1_fastq <DIR>` (required): raw FASTQ directory
  - `--path2_fastp <DIR>` (required): output directory for fastp results (`01-qc/`)
  - `--num_threads <int>` (default: `8`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`): forward read suffix
  - `--batch_size <int>` (default: `5`)
  - `--se`: single-end mode
  - `--length_required <int>` (default: `50`)
  - Notes: Writes per-sample `*_fastp.html/json`; if **multiqc** is present, also writes `01-qc/multiqc_report/multiqc_fastp_report.html`.  
  *(Implementation: automatic MultiQC invocation and output path)*

### Salmon mode
- **batch_salmon**
  - `--index <DIR>` (required): salmon index
  - `--path_fq <DIR>` (required): directory of FASTQs (after `fastq_qc`)
  - `--path_out <DIR>` (required): output root (e.g., `02-salmon/`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`)
  - `--batch_size <int>` (default: `1`): concurrent samples (processes)
  - `--num_threads <int>` (default: `8`): threads per salmon
  - `--gtf <FILE>`: optional GTF for `-g` gene-level quant
  - Behavior: safe R1 to R2 inference; per-sample `task.complete`; progress; preflight prints salmon version & index meta keys.

- **merge_salmon**
  - `--path_salmon <DIR>` (required): root containing per-sample salmon outputs (searched recursively)
  - `--project <STR>` (required): prefix for outputs
  - `--num_processes <int>`: I/O threads (default: CPU count)
  - Output: `<project>_salmon_tpm.tsv.gz`, `<project>_salmon_count.tsv.gz` under `--path_salmon` with progress and head preview.

- **prepare_salmon**
  - `-i/--input <TSV|TSV.GZ>` (required): Salmon-combined gene TPM table
  - `-o/--output <CSV/TSV>` (required): cleaned TPM matrix (genes × samples)
  - `-r/--return_feature {ENST|ENSG|symbol}` (default: `symbol`): which identifier to keep
  - `--remove_version`: strip version suffix from gene IDs (e.g., `ENSG000001.12 to ENSG000001`)

### STAR mode
- **batch_star_count**
  - `--index <DIR>` (required): STAR genomeDir
  - `--path_fq <DIR>` (required): directory of FASTQs (after `fastq_qc`)
  - `--path_out <DIR>` (required): outputs (e.g., `02-star/`)
  - `--suffix1 <str>` (default: `_1.fastq.gz`)
  - `--batch_size <int>` (default: `1`)
  - `--num_threads <int>` (default: `8`)
  - Notes: generates sorted BAM and `_ReadsPerGene.out.tab` per sample and a summary of paths.

- **merge_star_count**
  - `--path <DIR>` (required): directory containing multiple `*_ReadsPerGene.out.tab`
  - `--project <STR>` (required): output prefix
  - Output: `<project>.STAR.count.tsv.gz` (gzipped TSV with gene IDs as rows and samples as columns)

- **count2tpm**
  - `-i/--input <CSV/TSV[.gz]>` (required): raw count matrix (genes × samples)
  - `-o/--output <CSV/TSV>` (required): output TPM matrix
  - `--effLength_csv <CSV>`: optional effective-length file with columns `id`, `eff_length`, `symbol`
  - `--idtype {ensembl|entrez|symbol|mgi}` (default: `ensembl`)
  - `--org {hsa|mmus}` (default: `hsa`)
  - `--id <str>` (default: `id`): ID column name in `--effLength_csv`
  - `--length <str>` (default: `eff_length`): length column
  - `--gene_symbol <str>` (default: `symbol`): gene symbol column
  - `--check_data`: check & drop missing/invalid entries before conversion
  - `--remove_version`: strip version suffix from gene IDs

## (Optional) Mouse to Human symbol mapping
- **mouse2human_eset**
  - `-i/--input <CSV|TSV|TXT[.gz]>` (required): input expression **matrix** or **table**
  - `-o/--output <CSV|TSV|TXT[.gz]>` (required): converted matrix indexed by **human** symbols (genes × samples)
  - `--is_matrix`: treat input as a matrix (rows = **mouse** gene symbols, columns = samples); if omitted, runs in **table mode**
  - `--column_of_symbol <str>` (required in table mode): column name that contains **mouse** gene symbols
  - `--sep <,|\t>`: override input separator; if omitted, inferred by extension.
  - `--out_sep <,|\t>`: override output separator; if omitted, inferred by **output** path extension
  - `--verbose`: print shapes and basic run info

## (Optional) Annotate / de‑duplicate
- **anno_eset**
  - `-i/--input <CSV/TSV/TXT>` (required)
  - `-o/--output <CSV/TSV/TXT>` (required)
  - `--annotation {anno_hug133plus2|anno_rnaseq|anno_illumina|anno_grch38}` (required unless using external file)
  - `--annotation-file <pkl/csv/tsv/xlsx>`: external annotation (overrides built-in)
  - `--annotation-key <str>`: key to pick a table if external `.pkl` stores a dict of DataFrames
  - `--symbol <str>` (default: `symbol`): column used as gene symbol
  - `--probe  <str>` (default: `id`): column used as probe/feature ID
  - `--method {mean|sd|sum}` (default: `mean`): duplicate-ID aggregation
  - `--remove_version`: strip version suffix from gene IDs

## (Optional) Log2 transform
- **log2_eset**
  - `-i/--input <CSV/TSV/TXT>` (required)
  - `-o/--output <CSV/TSV/TXT>` (required)

## Signature scoring
- **calculate_sig_score**
  - `-i/--input <CSV/TSV/TXT>` (required), `-o/--output <CSV/TSV/TXT>` (required)
  - `--signature <one or more groups>` (required; space- or comma-separated; `all` uses every group)  
    Groups: `go_bp`, `go_cc`, `go_mf`, `signature_collection`, `signature_tme`, `signature_sc`, `signature_tumor`, `signature_metabolism`, `kegg`, `hallmark`, `reactome`
  - `--method {pca|zscore|ssgsea|integration}` (default: `pca`)
  - `--mini_gene_count <int>` (default: `3`)
  - `--adjust_eset`: apply extra filtering after log2 transform
  - `--parallel_size <int>` (default: `1`; threads for scoring (`PCA`/`zscore`/`ssGSEA`))

## Deconvolution / scoring
- **cibersort**
  - `-i/--input <CSV/TSV>` (required), `-o/--output <CSV/TSV>` (required)
  - `--perm <int>` (default: `100`)
  - `--QN <True|False>` (default: `True`): quantile normalization
  - `--absolute <True|False>` (default: `False`): absolute mode
  - `--abs_method {sig.score|no.sumto1}` (default: `sig.score`)
  - `--threads <int>` (default: `1`)  
  *Output: columns are suffixed with `_CIBERSORT`, index name is `ID`, separator inferred from output extension.*

- **quantiseq**
  - `-i/--input <CSV/TSV>` (required; genes × samples), `-o/--output <TSV>` (required)
  - `--arrays`: perform quantile normalization for arrays
  - `--signame <str>` (default: `TIL10`)
  - `--tumor`: remove genes highly expressed in tumors
  - `--scale_mrna`: enable mRNA scaling (otherwise raw signature proportions)
  - `--method {lsei|hampel|huber|bisquare}` (default: `lsei`)
  - `--rmgenes <str>` (default: `unassigned`; allowed: `default`, `none`, or comma-separated list)

- **epic**
  - `-i/--input <CSV/TSV>` (required; genes × samples)
  - `-o/--output <CSV/TSV>` (required)
  - `--reference {TRef|BRef|both}` (default: `TRef`)

- **estimate**
  - `-i/--input <CSV/TSV/TXT>` (required; genes × samples)
  - `-p/--platform {affymetrix|agilent|illumina}` (default: `affymetrix`)
  - `-o/--output <CSV/TSV/TXT>` (required)  
  *Output is transposed; columns are suffixed with `_estimate`; index label is `ID`; separator inferred from extension.*

- **mcpcounter**
  - `-i/--input <TSV>` (required; genes × samples)
  - `-f/--features {affy133P2_probesets|HUGO_symbols|ENTREZ_ID|ENSEMBL_ID}` (required)
  - `-o/--output <CSV/TSV>` (required)  
  *Output: suffixed with `_MCPcounter`; index label `ID`; separator inferred from extension.*

- **IPS**
  - `-i/--input <matrix>` (required), `-o/--output <file>` (required)  
  *No extra flags (the expression matrix yields IPS sub-scores and a total score).*

- **deside** (deep learning–based deconvolution)
  - `-m/--model_dir <dir>` (required): path to the pre-downloaded DeSide model directory
  - `-i/--input <CSV/TSV>` (required): rows = genes, columns = samples
  - `-o/--output <CSV>` (required)
  - `--exp_type {TPM|log_space|linear}` (default: `TPM`)  
    - `TPM`: already log2 processed  
    - `log_space`: `log2(TPM+1)`  
    - `linear`: linear space (TPM/counts)
  - `--gmt <file1.gmt file2.gmt ...>`: optional one or more GMT files for pathway masking
  - `--method_adding_pathway {add_to_end|convert}` (default: `add_to_end`)
  - `--scaling_by_constant`, `--scaling_by_sample`, `--one_minus_alpha`: optional scaling/transforms
  - `--print_info`: verbose logs
  - `--add_cell_type`: append predicted cell-type labels
  - `--transpose`: use if your file is *samples × genes*
  - `-r/--result_dir <dir>`: optional directory to save result plots/logs


## Clustering / decomposition
- **tme_cluster**
  - `-i/--input <CSV/TSV/TXT>` (required): input table for clustering.
    - Expected shape: first column = sample ID (use `--id` if not first), remaining columns = features.
  - `-o/--output <CSV/TSV/TXT>` (required): output file for clustering results.
  - `--features <spec>`: select feature columns by 1-based inclusive range, e.g. `1:22` (intended for CIBERSORT outputs; **exclude** the sample ID column when counting).
  - `--pattern <regex>`: alternatively select features by a regex on column names (e.g. `^CD8|^NK`).  
    *Tip: use one of `--features` or `--pattern`.*
  - `--id <str>` (default: first column): column name containing sample IDs.
  - `--scale` / `--no-scale`: toggle z-score scaling of features (help text: default = **True**).
  - `--min_nc <int>` (default: `2`): minimum number of clusters to try.
  - `--max_nc <int>` (default: `6`): maximum number of clusters to try.
  - `--max_iter <int>` (default: `10`): maximum iterations for k-means.
  - `--tol <float>` (default: `1e-4`): convergence tolerance for centroid updates.
  - `--print_result`: print intermediate KL scores and cluster counts.
  - `--input_sep <str>` (default: auto): input delimiter (e.g. `,` or `\t`); auto-detected if unset.
  - `--output_sep <str>` (default: auto): output delimiter; inferred from filename if unset.

- **nmf**
  - `-i/--input <CSV/TSV>` (required): matrix to factorize; first column should be sample names (index).
  - `-o/--output <DIR>` (required): directory to save results.
  - `--kmin <int>` (default: `2`): minimum `k` (inclusive).
  - `--kmax <int>` (default: `8`): maximum `k` (inclusive).
  - `--features <spec>`: 1-based inclusive selection of feature columns (e.g. `2-10` or `1:5`), typically cell-type columns.
  - `--log1p`: apply `log1p` to the input (useful for counts).
  - `--normalize`: L1 row normalization (each sample sums to 1).
  - `--shift <float>` (default: `None`): if data contain negatives, add a constant to make all values non-negative.
  - `--random-state <int>` (default: `42`): random seed for NMF.
  - `--max-iter <int>` (default: `1000`): NMF max iterations.
  - `--skip_k_2`: skip evaluating `k = 2` when searching for the best `k`.

## Ligand–receptor
- **LR_cal**
  - `-i/--input <CSV/TSV>` (required): expression matrix (genes × samples).
  - `-o/--output <CSV/TSV>` (required): file to save LR scores.
  - `--data_type {count|tpm}` (default: `tpm`): type of the input matrix.
  - `--id_type <str>` (default: `ensembl`): gene ID type expected by the LR backend.Choices: `ensembl`, `entrez`, `symbol`, `mgi`.
  - `--cancer_type <str>` (default: `pancan`): cancer-type network to use.
  - `--verbose`: verbose logging.