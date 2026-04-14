"""
Detailed semantic profiles for native iobrpy commands exposed by the harness.
"""

from __future__ import annotations

from typing import Any, Dict, List


COMMAND_PROFILES: Dict[str, Dict[str, Any]] = {
    "prepare_salmon": {
        "summary": "Turn a merged Salmon expression table into the cleaned TPM matrix used by downstream IOBRpy analyses.",
        "detailed_description": "Use this after Salmon quantification has already been merged. It rewrites the merged Salmon matrix into a TPM-style expression matrix and can keep ENST, ENSG, or gene-symbol identifiers.",
        "input_expectations": [
            "Input is a merged Salmon TSV/TSV.GZ table, not raw FASTQ and not a count matrix.",
            "Rows should represent transcript or gene identifiers and columns should represent samples.",
        ],
        "outputs": [
            "A cleaned TPM matrix ready for commands such as calculate_sig_score, tme_profile, and LR_cal.",
        ],
        "when_to_use": [
            "You already ran batch_salmon and merge_salmon and now need the final TPM matrix.",
            "You want to convert Salmon outputs into symbol-, ENSG-, or ENST-based expression space.",
        ],
        "use_another_command_when": [
            "Use runall when you want the full FASTQ-to-TME workflow instead of a single post-Salmon step.",
            "Use count2tpm when your starting point is a gene-level count matrix rather than Salmon output.",
        ],
        "project_context_notes": [
            "In a standard Salmon-mode runall tree, this sits between `02-salmon/` and `03-tpm/`.",
            "Typical upstream files are `<project>_salmon_tpm.tsv.gz` under the Salmon output directory.",
            "Typical downstream outputs are a cleaned TPM matrix such as `03-tpm/prepare_salmon.csv`, which can then feed `log2_eset`, `calculate_sig_score`, `tme_profile`, or `LR_cal`.",
            "The CLI keeps `symbol` as the normal downstream identifier target unless the user explicitly asks to retain ENST or ENSG identifiers.",
            "In `runall`, this stage is normally followed immediately by `log2_eset`, which creates the standard `03-tpm/tpm_matrix.csv` handoff matrix.",
        ],
    },
    "count2tpm": {
        "summary": "Convert a gene-level count matrix into TPM using built-in or supplied gene-length annotation.",
        "detailed_description": "This is the STAR/count-matrix entrypoint into TPM space. It takes counts, resolves gene lengths, optionally strips version suffixes, and writes the TPM matrix used by downstream TME workflows.",
        "input_expectations": [
            "Input is a counts matrix with genes as rows and samples as columns.",
            "Gene IDs should match the selected idtype or the supplied effective-length CSV.",
        ],
        "outputs": [
            "A TPM matrix suitable for tme_profile or any of the downstream scoring and deconvolution commands.",
        ],
        "when_to_use": [
            "You already have counts from STAR or another quantifier and need TPM before TME analysis.",
            "You want to control gene-length annotation explicitly with --effLength_csv.",
        ],
        "use_another_command_when": [
            "Use prepare_salmon when the input comes from Salmon quantification instead of counts.",
            "Use runall when you want FASTQ QC, quantification, TPM generation, and downstream TME steps in one command.",
        ],
        "project_context_notes": [
            "In a standard STAR-mode runall tree, this sits after `merge_star_count` and before downstream TME analysis.",
            "Typical upstream input is `<project>.STAR.count.tsv.gz` from `merge_star_count`.",
            "Typical downstream output is a TPM matrix such as `03-tpm/count2tpm.csv`, which is often followed by `log2_eset` to create `03-tpm/tpm_matrix.csv`.",
            "The implementation can use built-in human or mouse annotation resources, so the `idtype` and organism assumptions need to match the count matrix rather than being guessed loosely.",
            "When users provide an `effLength_csv`, it overrides the built-in reference lengths and becomes the key determinant of the TPM denominator.",
        ],
    },
    "anno_eset": {
        "summary": "Annotate an expression matrix and collapse duplicate probes or identifiers into one expression set.",
        "detailed_description": "This command harmonizes identifiers using a built-in or external annotation table and then resolves duplicates with a configurable aggregation rule such as mean, sd, or sum.",
        "input_expectations": [
            "Input is an expression matrix or expression set that needs identifier annotation.",
            "If you use an external annotation file, it must contain the symbol and probe columns you specify.",
        ],
        "outputs": [
            "An annotated expression matrix with duplicate handling already applied.",
        ],
        "when_to_use": [
            "You need to map probe IDs or raw identifiers onto a stable gene-symbol representation.",
            "You need to harmonize a matrix before signature scoring or deconvolution.",
        ],
        "use_another_command_when": [
            "Use mouse2human_eset when the main task is species conversion from mouse to human gene symbols.",
            "Use prepare_salmon or count2tpm when the main task is TPM generation rather than annotation.",
        ],
        "project_context_notes": [
            "Built-in annotation keys such as `anno_hug133plus2`, `anno_rnaseq`, `anno_illumina`, and `anno_grch38` come from `iobrpy.resources/anno_eset.pkl`.",
            "An external annotation file overrides the bundled annotations and must provide the symbol and probe columns named by `--symbol` and `--probe`.",
            "The command keeps one row per final symbol by scoring duplicated rows with `mean`, `sd`, or `sum`, then retaining the highest-scoring representative.",
            "It can optionally strip Ensembl version suffixes before annotation, which is important when the matrix index looks like `ENSG000001.12` rather than bare Ensembl IDs.",
        ],
    },
    "calculate_sig_score": {
        "summary": "Compute pathway or signature scores from an expression matrix using PCA, z-score, ssGSEA, or integration.",
        "detailed_description": "This is the signature-scoring entrypoint used directly or inside tme_profile and runall. It reads a genes-by-samples expression matrix and produces per-sample signature scores for the requested signature groups.",
        "input_expectations": [
            "Input should be a TPM-like or otherwise normalized expression matrix with genes as rows and samples as columns.",
            "You must choose one or more signature groups, or use all to score the full built-in collection.",
        ],
        "outputs": [
            "A per-sample signature score table for the requested pathway or gene-signature groups.",
        ],
        "when_to_use": [
            "You want pathway or immune-signature scores without running the full downstream TME workflow.",
            "You want to compare PCA, z-score, ssGSEA, or the integrated signature-scoring strategy.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you also want immune deconvolution and LR_cal from the same TPM matrix.",
            "Use runall when you are starting from FASTQ rather than an existing expression matrix.",
        ],
        "project_context_notes": [
            "In standard `tme_profile` or `runall` outputs, this usually writes `01-signatures/calculate_sig_score.csv` or `04-signatures/calculate_sig_score.csv` depending on the workflow root.",
            "The `integration` method combines PCA, z-score, and ssGSEA-style information into one merged signature matrix.",
            "This command is often the best direct input source for downstream `tme_cluster` or `nmf` when the user wants subtype-style analysis from feature scores.",
        ],
    },
    "cibersort": {
        "summary": "Run CIBERSORT cell-fraction deconvolution on a bulk expression matrix.",
        "detailed_description": "This command estimates immune cell fractions from a genes-by-samples mixture matrix. It is one of the default downstream TME methods used inside tme_profile and runall.",
        "input_expectations": [
            "Input is an expression mixture matrix, typically TPM-like, with genes as rows and samples as columns.",
            "Use --QN, --absolute, and --abs_method only when you intentionally need those CIBERSORT modes.",
        ],
        "outputs": [
            "A cell-fraction result table with CIBERSORT-derived estimates for each sample.",
        ],
        "when_to_use": [
            "You want only CIBERSORT output instead of the full tme_profile bundle.",
            "You need to tune permutations or absolute scoring behavior explicitly.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want CIBERSORT plus the other default TME methods and LR_cal.",
            "Use runall when you are starting from FASTQ rather than an existing expression matrix.",
        ],
        "project_context_notes": [
            "The implementation uses the bundled LM22 signature matrix from `iobrpy.resources/lm22.txt`.",
            "In standard `tme_profile` or `runall` outputs, this usually writes `02-tme/cibersort_results.csv` or `05-tme/cibersort_results.csv`.",
            "The output matrix is sample-by-celltype with an `ID` index label, and cell-type columns are suffixed with `_CIBERSORT`.",
            "The native default for `QN` is still important to surface because array-style normalization assumptions do not always match TPM-like RNA-seq input.",
        ],
    },
    "IPS": {
        "summary": "Calculate Immunophenoscore-related outputs from an expression matrix.",
        "detailed_description": "IPS is a score-based TME command rather than a cell-fraction deconvolution method. It summarizes immunophenoscore-related immune activity from an existing expression matrix.",
        "input_expectations": [
            "Input is an expression matrix file rather than raw sequencing data.",
        ],
        "outputs": [
            "An IPS result table for the input samples.",
        ],
        "when_to_use": [
            "You specifically need IPS output as a standalone score table.",
            "You want to inspect the score-based component separately from deconvolution outputs.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want IPS to run alongside the rest of the standard downstream TME methods.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "IPS is a score-based branch, not a cell-fraction deconvolution method, so its outputs should not be described as fractions or cell abundances.",
            "In the bundled downstream workflows, this usually writes `02-tme/IPS_results.csv` or `05-tme/IPS_results.csv`.",
            "The standalone CLI keeps the interface intentionally simple: it mainly needs the expression matrix and output path, without the larger option surface used by several other deconvolution tools.",
            "This makes IPS a good low-friction add-on when the user wants immunophenoscore-style context rather than a new fraction matrix.",
        ],
    },
    "estimate": {
        "summary": "Compute ESTIMATE stromal and immune scores, plus purity-related context, from an expression matrix.",
        "detailed_description": "ESTIMATE is a score-based tumor-microenvironment command. It produces stromal and immune context scores from a genes-by-samples expression matrix and is also included inside tme_profile and runall.",
        "input_expectations": [
            "Input is an expression matrix with genes as rows and samples as columns.",
            "Choose the platform that matches the matrix when ESTIMATE platform handling matters.",
        ],
        "outputs": [
            "A score table with ESTIMATE-derived tumor microenvironment context per sample.",
        ],
        "when_to_use": [
            "You need only ESTIMATE output or want to tune the platform explicitly.",
            "You want a score-based TME view rather than cell-fraction deconvolution alone.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want ESTIMATE together with the rest of the default TME bundle.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "In the bundled downstream workflows, this usually writes `02-tme/estimate_results.csv` or `05-tme/estimate_results.csv`.",
            "The orchestrators default ESTIMATE to `affymetrix` unless the user explicitly overrides the platform.",
            "The standalone output is transposed into a sample-by-score table, with score columns suffixed by `_estimate` and an `ID` index label.",
            "ESTIMATE is part of the default downstream TME bundle, but it remains a score-based context command rather than a cell-fraction estimator.",
        ],
    },
    "mcpcounter": {
        "summary": "Compute MCPcounter immune and stromal abundance scores from an expression matrix.",
        "detailed_description": "MCPcounter estimates abundance-style immune and stromal signals rather than normalized cell fractions. It requires you to tell IOBRpy which identifier system the matrix uses.",
        "input_expectations": [
            "Input is an expression matrix with genes as rows and samples as columns.",
            "You must provide the correct --features identifier type for the input matrix.",
        ],
        "outputs": [
            "A sample-by-feature MCPcounter score table.",
        ],
        "when_to_use": [
            "You want marker-based abundance scores rather than fraction-based deconvolution.",
            "You need MCPcounter as a standalone result outside the full tme_profile bundle.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want MCPcounter together with the other default downstream TME methods.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "The critical command-specific choice is `--features`, which must reflect the row-identifier namespace such as `HUGO_symbols`, `ENSEMBL_ID`, or `ENTREZ_ID`.",
            "The downstream orchestrators default this to `HUGO_symbols` unless the user explicitly supplies another namespace.",
            "In standard bundled outputs, this usually writes `02-tme/mcpcounter_results.csv` or `05-tme/mcpcounter_results.csv`.",
            "The output is a sample-by-feature abundance matrix whose columns are suffixed with `_MCPcounter` and whose index label is `ID`.",
        ],
    },
    "quantiseq": {
        "summary": "Run quanTIseq deconvolution on a bulk expression matrix to estimate immune-cell composition.",
        "detailed_description": "quanTIseq is one of the default downstream TME methods in IOBRpy. It works on an existing expression matrix and can optionally enable array normalization, tumor-gene filtering, and mRNA scaling.",
        "input_expectations": [
            "Input is a genes-by-samples expression matrix, typically TPM-like.",
            "Optional flags such as --arrays, --tumor, and --scale_mrna change preprocessing and interpretation.",
        ],
        "outputs": [
            "A quanTIseq deconvolution result table for each sample.",
        ],
        "when_to_use": [
            "You want quanTIseq output only, or you want to tune its method and preprocessing options directly.",
            "You want to compare quanTIseq against CIBERSORT, EPIC, or MCPcounter.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want quanTIseq together with the rest of the standard TME bundle.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "The bundled downstream workflows usually write this to `02-tme/quantiseq_results.csv` or `05-tme/quantiseq_results.csv`.",
            "Inside `tme_profile`, the orchestrator turns on `--arrays`, `--tumor`, and `--scale_mrna` unless the user explicitly overrides those behaviors.",
            "The implementation supports solver choices such as `lsei`, `hampel`, `huber`, and `bisquare`, so `--method` changes more than just a cosmetic setting.",
            "This command returns a deconvolution matrix, making it a much cleaner clustering candidate than a huge all-signature score table when the user wants standard TME subtyping.",
        ],
    },
    "epic": {
        "summary": "Run EPIC deconvolution on a bulk expression matrix using TRef, BRef, or both reference modes.",
        "detailed_description": "EPIC estimates immune, stromal, and tumor-context signals from an existing bulk expression matrix. It is one of the standard deconvolution methods included in tme_profile and runall.",
        "input_expectations": [
            "Input is a genes-by-samples expression matrix, typically TPM-like.",
            "Choose the EPIC reference mode that matches the analysis context.",
        ],
        "outputs": [
            "An EPIC result table with per-sample cell-fraction estimates.",
        ],
        "when_to_use": [
            "You want EPIC output only or want to control the EPIC reference explicitly.",
            "You want to compare EPIC with other deconvolution strategies on the same matrix.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want EPIC together with the rest of the default TME bundle.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "The standard bundled workflows usually write this to `02-tme/epic_results.csv` or `05-tme/epic_results.csv`.",
            "The orchestrators default EPIC to the `TRef` reference unless the user explicitly asks for `BRef` or `both`.",
            "Although the internal implementation still retains a slower range-based optimization path, the exposed CLI surface in `iobrpy-cli` centers on the reference choice rather than those solver internals.",
            "EPIC is part of the default downstream TME bundle, so it is usually already present when a directory looks like a completed `tme_profile` or `runall` output.",
        ],
    },
    "tme_cluster": {
        "summary": "Cluster samples using existing TME feature matrices such as signature-score tables or deconvolution outputs.",
        "detailed_description": "This is the primary clustering route for downstream TME subtype analysis. It runs clustering on an existing sample-feature matrix after you already computed immune-deconvolution outputs or signature scores. In the tutorials, the default clustering examples use a CIBERSORT result matrix first; other deconvolution outputs or a curated feature table can also work, while very wide signature-score matrices should usually be user-selected rather than auto-chosen by the agent.",
        "input_expectations": [
            "Input should be a sample-by-feature matrix whose first column is the sample ID and whose remaining columns are TME-related features.",
            "The tutorial-style default input is a CIBERSORT result table such as `cibersort.csv` or `cibersort_results.csv`.",
            "Feature ranges in `--features` are counted after excluding the sample ID column. For example, `--features 1:22` usually means the second through twenty-third actual CSV columns.",
            "Signature-score matrices can also be used, but they may be very wide; if the file contains many pathways or signatures, the feature subset should usually be chosen by the user rather than guessed automatically.",
        ],
        "outputs": [
            "A clustering result table or output directory with sample cluster assignments and subtype-oriented summaries.",
        ],
        "when_to_use": [
            "You already have deconvolution or selected signature-feature outputs and want to stratify samples.",
            "You want a clustering-focused downstream analysis rather than TPM generation or deconvolution.",
        ],
        "use_another_command_when": [
            "Use nmf when you want an NMF-based clustering/decomposition workflow instead of k-means-style TME clustering.",
            "Use tme_profile or runall first if you still need to generate the TME feature matrix itself.",
        ],
        "project_context_notes": [
            "When the user asks for clustering in general, prefer `tme_cluster` before `nmf` unless they explicitly want NMF or factorization-style latent programs.",
            "The tutorial examples use CIBERSORT outputs as the default clustering input, so a detected `cibersort_results.csv` is usually the safest first choice.",
            "For `tme_cluster`, the sample ID column is excluded before counting `--features`, so left bounds commonly start at `1` rather than `2`.",
            "If only signature-score outputs are available, or if the signature matrix is very wide, it is usually better to ask the user which feature subset to cluster on instead of auto-running on all columns.",
            "The parent directory of the `--output` file should exist before execution, because the workflow does not create missing output parents automatically.",
        ],
    },
    "LR_cal": {
        "summary": "Compute ligand-receptor interaction scores from an expression matrix.",
        "detailed_description": "LR_cal is the ligand-receptor analysis step used in downstream TME workflows. It expects an existing expression matrix and produces interaction scores using a cancer-type-specific ligand-receptor network.",
        "input_expectations": [
            "Input is an expression matrix with genes as rows and samples as columns.",
            "You should provide the correct data_type and id_type when the defaults do not match the matrix.",
        ],
        "outputs": [
            "A ligand-receptor interaction score table for the input samples.",
        ],
        "when_to_use": [
            "You already have TPM-like expression data and want ligand-receptor scoring only.",
            "You want to run LR_cal separately from tme_profile or runall.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you want LR_cal to run together with the standard downstream TME bundle.",
            "Use runall when you are starting from FASTQ.",
        ],
        "project_context_notes": [
            "In standard downstream layouts, this usually writes `03-LR_cal/lr_cal.csv` under `tme_profile` or `06-LR_cal/lr_cal.csv` under `runall`.",
            "The default assumptions are `data_type=tpm`, `id_type=ensembl`, and `cancer_type=pancan`, so those defaults should be called out when the matrix or disease context differs.",
            "This step belongs after TPM-ready expression analysis; it is not designed to consume deconvolution matrices or clustering tables.",
            "When a path already has signatures and TME outputs, LR_cal is often the last standard bundled branch before optional clustering or repertoire-specific follow-up work.",
        ],
    },
    "nmf": {
        "summary": "Run NMF-based clustering and decomposition on an existing sample-feature matrix.",
        "detailed_description": "This is the decomposition-oriented clustering alternative. It performs non-negative matrix factorization on a downstream sample-feature matrix such as CIBERSORT outputs, other TME matrices, or a curated signature matrix. In normal clustering requests, `tme_cluster` should usually be suggested first and `nmf` should stay as the factorization-focused alternative.",
        "input_expectations": [
            "Input is a sample-by-feature matrix whose first column is usually the sample name or sample ID and whose remaining columns are clustering features.",
            "The tutorial examples use a CIBERSORT result matrix with `--features 1:22`, where feature counting excludes the first sample-ID column.",
            "Use normalization, shifting, and feature selection options only after confirming they fit the matrix scale and the intended factorization workflow.",
            "Very wide signature-score matrices may need user-selected feature subsets instead of a blind all-column run.",
        ],
        "outputs": [
            "An output directory containing NMF clustering results and related artifacts.",
        ],
        "when_to_use": [
            "You want latent-feature discovery or NMF-based cluster selection on TME-related matrices.",
            "You need a downstream clustering stage after deconvolution or signature scoring.",
        ],
        "use_another_command_when": [
            "Use tme_cluster when you want the dedicated TME clustering command instead of NMF decomposition.",
            "Use tme_profile or runall first if you still need to create the TME feature matrix itself.",
        ],
        "project_context_notes": [
            "When the user asks for ordinary TME clustering or subtype assignment, keep `nmf` behind `tme_cluster` unless they explicitly ask for NMF or latent-program discovery.",
            "The tutorial examples also use CIBERSORT outputs as the default NMF input, not a raw TPM matrix.",
            "As with `tme_cluster`, `--features` counts only feature columns and excludes the first sample-ID column when interpreting ranges such as `1:22`.",
            "The parent directory of the `--output` directory path should exist before execution planning, because the workflow expects a real writable location.",
        ],
    },
    "mouse2human_eset": {
        "summary": "Convert mouse gene symbols to human gene symbols in matrix or table form.",
        "detailed_description": "This command rewrites mouse identifiers into human symbols so downstream human-focused scoring or deconvolution workflows can consume the matrix more easily.",
        "input_expectations": [
            "Input can be a genes-by-samples matrix or a tabular file with a gene-symbol column.",
            "Use --is_matrix for matrix mode or --column_of_symbol for table mode.",
        ],
        "outputs": [
            "A converted file with mouse gene symbols mapped onto human symbols.",
        ],
        "when_to_use": [
            "You have mouse-symbol input but need human-symbol output for downstream workflows.",
            "You need identifier harmonization before annotation, signature scoring, or deconvolution.",
        ],
        "use_another_command_when": [
            "Use anno_eset when the task is annotation and duplicate handling rather than species conversion.",
            "Use tme_profile or runall only after the expression matrix is already in the correct species space.",
        ],
        "project_context_notes": [
            "The conversion uses the bundled `mus_human.pkl` mapping resource from `iobrpy.resources`.",
            "Matrix mode expects mouse symbols on the row index, while table mode expects a named symbol column and then collapses duplicates before conversion.",
            "After mapping, the workflow reuses `anno_eset(..., method='mean')` so duplicate handling stays consistent with the annotation pipeline.",
            "This is a preprocessing bridge into human-focused downstream commands such as `calculate_sig_score`, `cibersort`, `tme_profile`, and `LR_cal`.",
        ],
    },
    "batch_salmon": {
        "summary": "Run Salmon quantification across a directory of paired-end FASTQ files.",
        "detailed_description": "This is the Salmon quantification step only. It produces per-sample Salmon outputs and does not itself merge quantifications, create the final TPM matrix, or run TME analyses.",
        "input_expectations": [
            "Input is a FASTQ directory plus a Salmon transcriptome index.",
            "FASTQ naming must allow paired-end inference from the supplied suffix1 pattern.",
        ],
        "outputs": [
            "Per-sample Salmon quantification directories containing quant.sf and related outputs.",
        ],
        "when_to_use": [
            "You want only the Salmon quantification stage and will handle merging or downstream analysis separately.",
            "You are building a custom Salmon-based pipeline step by step.",
        ],
        "use_another_command_when": [
            "Use runall when you want the full FASTQ-to-TME workflow rather than quantification alone.",
            "Use merge_salmon next when you need a combined matrix from the per-sample Salmon outputs.",
        ],
        "project_context_notes": [
            "Typical output is one Salmon folder per sample containing `quant.sf`.",
            "In a standard runall tree this corresponds to the `02-salmon/` stage after fastp/QC has already populated `01-qc/`.",
            "This stage does not create the merged project matrix on its own; `merge_salmon` is the next normal handoff point.",
            "The harness should describe this as a quantification stage, not as a finished TPM-preparation or downstream-TME result.",
        ],
    },
    "merge_salmon": {
        "summary": "Merge multiple Salmon quant.sf outputs into project-level TPM and NumReads matrices.",
        "detailed_description": "This command collects per-sample Salmon quantification outputs and builds the project-wide matrices that feed into prepare_salmon or custom downstream analysis.",
        "input_expectations": [
            "Input is the root directory containing Salmon quantification subdirectories.",
            "The path should already contain per-sample quant.sf files produced by batch_salmon or another Salmon workflow.",
        ],
        "outputs": [
            "Merged project-level TPM and NumReads matrices.",
        ],
        "when_to_use": [
            "You already finished Salmon quantification and now need a combined matrix across samples.",
            "You want to prepare a project-wide matrix before prepare_salmon or other downstream work.",
        ],
        "use_another_command_when": [
            "Use batch_salmon first when you still need to quantify raw FASTQ files.",
            "Use runall when you want Salmon quantification, merge, TPM preparation, and downstream TME analysis in one workflow.",
        ],
        "project_context_notes": [
            "Typical outputs are `<project>_salmon_tpm.tsv.gz` and `<project>_salmon_count.tsv.gz` under the Salmon result root.",
            "Those merged matrices are the normal handoff point into `prepare_salmon`.",
            "The implementation searches recursively for `quant.sf` files, so the input path is the Salmon result root rather than a single sample directory.",
            "This merge step is still upstream of the cleaned TPM matrix that downstream TME commands normally expect.",
        ],
    },
    "merge_star_count": {
        "summary": "Merge multiple STAR GeneCounts outputs into a single count matrix.",
        "detailed_description": "This is the STAR merge step only. It collects `_ReadsPerGene.out.tab` files from STAR outputs and builds a project-level count matrix for later TPM conversion.",
        "input_expectations": [
            "Input is the directory that contains STAR GeneCounts outputs.",
            "The folder should already contain per-sample STAR results from batch_star_count or a compatible workflow.",
        ],
        "outputs": [
            "A merged count matrix across samples.",
        ],
        "when_to_use": [
            "You already completed STAR counting and now need a project-level matrix.",
            "You want to convert STAR counts to TPM with count2tpm.",
        ],
        "use_another_command_when": [
            "Use batch_star_count first when you still need to align raw FASTQ files.",
            "Use runall when you want the full STAR-mode FASTQ-to-TME workflow.",
        ],
        "project_context_notes": [
            "Typical upstream evidence is a STAR result directory containing many `*_ReadsPerGene.out.tab` files.",
            "Typical output is `<project>.STAR.count.tsv.gz`, which is the normal input to `count2tpm`.",
            "This command only merges count-like STAR outputs; it does not itself apply gene-length normalization or create the final TPM matrix.",
            "Because the merged output is gzipped TSV with genes as rows and samples as columns, it is the matrix-level handoff point before `count2tpm` and `log2_eset`.",
        ],
    },
    "batch_star_count": {
        "summary": "Run STAR alignment plus GeneCounts over a directory of paired-end FASTQ files.",
        "detailed_description": "This is the STAR quantification/counting step only. It produces per-sample STAR outputs with gene counts, but it does not itself merge counts, convert them to TPM, or run downstream TME analysis.",
        "input_expectations": [
            "Input is a FASTQ directory plus a STAR genome index.",
            "FASTQ naming must allow paired-end inference from the supplied suffix1 pattern.",
        ],
        "outputs": [
            "Per-sample STAR outputs, including GeneCounts files and aligned BAM outputs.",
        ],
        "when_to_use": [
            "You want only the STAR alignment/counting stage and will manage merging or downstream analysis separately.",
            "You are building a custom STAR-based pipeline step by step.",
        ],
        "use_another_command_when": [
            "Use runall with --mode star when you want the full FASTQ-to-TME workflow rather than just STAR counting.",
            "Use merge_star_count next when you need a combined count matrix across samples.",
        ],
        "project_context_notes": [
            "Typical outputs include aligned BAMs, STAR logs, and per-sample `*_ReadsPerGene.out.tab` files.",
            "In a standard runall tree this corresponds to the `02-star/` stage before `merge_star_count` and `count2tpm`.",
            "These STAR BAM outputs can later feed immune-focused branches such as TRUST4 or HLA workflows, but this command itself is still only the quantification/counting stage.",
            "The harness should therefore treat it as upstream evidence, not as proof that TPM generation or downstream TME analysis is already complete.",
        ],
    },
    "fastq_qc": {
        "summary": "Run fastp-based QC and trimming over a FASTQ directory and organize cleaned reads plus reports.",
        "detailed_description": "This is the FASTQ preprocessing stage used inside runall. It prepares cleaned reads, per-sample QC artifacts, and a standardized QC output directory before quantification.",
        "input_expectations": [
            "Input is a directory of raw FASTQ files.",
            "Use the paired-end suffix and single-end options to match the sequencing layout.",
        ],
        "outputs": [
            "A cleaned FASTQ output directory plus fastp QC reports.",
        ],
        "when_to_use": [
            "You want only FASTQ QC and trimming before running the rest of the pipeline yourself.",
            "You need cleaned reads for a custom downstream workflow.",
        ],
        "use_another_command_when": [
            "Use runall when you want FASTQ QC followed immediately by quantification, TPM generation, and TME analysis.",
            "Use batch_salmon or batch_star_count only after QC if your workflow requires cleaned reads first.",
        ],
        "project_context_notes": [
            "Per-sample outputs include cleaned FASTQ files plus `<sample>_fastp.html` and `<sample>_fastp.json` reports.",
            "The implementation writes per-sample task markers and then attempts to build `multiqc_report/multiqc_fastp_report.html` from the fastp JSON outputs.",
            "In a standard runall tree this is the `01-qc/` stage and becomes the direct input root for `batch_salmon` or `batch_star_count`.",
            "Single-end mode is explicitly supported through `--se`, so the agent should not assume paired-end reads if the user says the data are single-end.",
        ],
    },
    "log2_eset": {
        "summary": "Apply log2(x+1) transformation to an existing expression matrix.",
        "detailed_description": "This is a simple matrix transformation step for workflows that need log-transformed expression values before downstream analysis or visualization.",
        "input_expectations": [
            "Input is an existing expression matrix with genes as rows and samples as columns.",
        ],
        "outputs": [
            "A transformed matrix in the same broad shape as the input, but with log2(x+1) applied.",
        ],
        "when_to_use": [
            "You specifically need a log-transformed matrix for a later analysis step.",
            "You want to normalize presentation or scoring input outside the end-to-end workflows.",
        ],
        "use_another_command_when": [
            "Use prepare_salmon or count2tpm when the main task is TPM generation rather than log transformation.",
            "Use runall or tme_profile when you want the complete workflow rather than a single matrix transform.",
        ],
        "project_context_notes": [
            "The implementation auto-detects common delimiters and chooses the output delimiter from the target file extension, mirroring the input when the extension is ambiguous.",
            "It errors on values below -1 because `log2(x+1)` is undefined there, and it warns when negative-but-valid values are present.",
            "In the standard runall workflow, this converts `03-tpm/prepare_salmon.csv` or `03-tpm/count2tpm.csv` into the canonical `03-tpm/tpm_matrix.csv` handoff matrix.",
            "This is a matrix-preparation step, not a downstream TME analysis result by itself.",
        ],
    },
    "trust4": {
        "summary": "Run TRUST4 for TCR/BCR repertoire reconstruction from BAM or FASTQ inputs.",
        "detailed_description": "This command wraps TRUST4 and focuses on immune-repertoire reconstruction rather than tumor-microenvironment deconvolution. It supports BAM input, paired FASTQ input, single-end FASTQ input, or FASTQ-directory batch mode.",
        "input_expectations": [
            "Input can be a BAM file, BAM directory, paired FASTQ files, single-end FASTQ, or a FASTQ directory depending on mode.",
            "Reference arguments and barcode/UMI options matter when you are using advanced TRUST4 modes.",
        ],
        "outputs": [
            "TRUST4 repertoire outputs for TCR/BCR reconstruction.",
        ],
        "when_to_use": [
            "You want immune repertoire analysis rather than expression-based TME deconvolution.",
            "You want TRUST4 as a standalone step or as the repertoire branch of a larger workflow.",
        ],
        "use_another_command_when": [
            "Use runall when you want FASTQ-to-TME plus the standard downstream branch, including TRUST4 outputs.",
            "Use hla_typing or extract_hla_read/spechla for HLA typing rather than TCR/BCR reconstruction.",
        ],
        "project_context_notes": [
            "Besides raw TRUST4 assembly outputs, iobrpy also derives merged repertoire summaries such as `trust4_immdata.csv` and `trust4_immune_indices.csv` after the run finishes.",
            "Batch FASTQ mode expects paired files under `--fqdir` with suffixes such as `_1.fastq.gz` and `_2.fastq.gz`.",
            "This is an immune-repertoire branch, not part of the six default TME deconvolution methods.",
        ],
    },
    "spechla": {
        "summary": "Run SpecHLA on a single sample from paired FASTQ inputs for exon/RNA or WGS-style HLA typing.",
        "detailed_description": "This is the per-sample SpecHLA command. It expects paired FASTQ reads for one sample and produces HLA typing output for that sample only.",
        "input_expectations": [
            "Input is a single sample represented by read1 and read2 FASTQ files.",
            "Choose --use-exon 1 for exon/RNA style runs and 0 for WGS mode.",
        ],
        "outputs": [
            "Per-sample SpecHLA typing results in the requested output directory.",
        ],
        "when_to_use": [
            "You already have the FASTQ inputs for one sample and want to run SpecHLA directly.",
            "You are debugging or customizing a single-sample HLA typing run.",
        ],
        "use_another_command_when": [
            "Use extract_hla_read first when your starting point is BAM/CRAM rather than FASTQ.",
            "Use hla_typing when you want the batch BAM-directory workflow rather than one sample at a time.",
        ],
        "project_context_notes": [
            "This wrapper auto-detects the bundled SpecHLA root, checks Python and external dependencies, and can trigger helper build steps such as SpecHap installation when needed.",
            "It is the FASTQ-level HLA typing stage that normally follows `extract_hla_read` when the starting point was BAM or CRAM.",
            "Typical per-sample outputs include `hla.result.txt` or `hla.results.txt` inside the requested sample output directory.",
            "The main exposed choice is whether to run exon/RNA mode or WGS mode through `--use_exon`, not to reinterpret the command as a matrix-based analysis step.",
        ],
    },
    "extract_hla_read": {
        "summary": "Extract HLA-related reads from a BAM/CRAM file and convert them into FASTQ for downstream SpecHLA typing.",
        "detailed_description": "This is the preparation step for HLA typing when you start from aligned reads. It wraps the SpecHLA extraction helper and produces the HLA-focused FASTQ files that SpecHLA consumes.",
        "input_expectations": [
            "Input is a sorted and indexed BAM or CRAM file plus a reference selection.",
            "This command is about HLA-read extraction, not general RNA-seq quantification or TME analysis.",
        ],
        "outputs": [
            "An output directory containing extracted HLA-focused FASTQ reads for the sample.",
        ],
        "when_to_use": [
            "You start from BAM/CRAM and need FASTQ-like HLA reads before running SpecHLA.",
            "You want to validate or customize the extraction step separately from full batch HLA typing.",
        ],
        "use_another_command_when": [
            "Use spechla after extraction when you are ready to type a single sample from FASTQ.",
            "Use hla_typing when you want the full batch BAM-directory workflow that chains extraction and SpecHLA for you.",
        ],
        "project_context_notes": [
            "Typical inputs are sorted and indexed BAM or CRAM files plus the matching hg19/hg38 reference choice.",
            "Typical outputs are paired HLA-read FASTQ files for one sample, which then become the direct input to `spechla`.",
            "IOBRpy can lazily install missing helper dependencies for this workflow unless the user disables auto-install behavior.",
        ],
    },
    "hla_typing": {
        "summary": "Run the batch HLA typing workflow over a BAM directory by chaining ExtractHLAread and SpecHLA.",
        "detailed_description": "This is the directory-scale HLA typing entrypoint. It scans a BAM directory, extracts HLA reads for each sample, runs SpecHLA, and then merges the per-sample outputs.",
        "input_expectations": [
            "Input is a directory containing BAM files, not a single FASTQ pair or a TPM matrix.",
            "You must provide the reference genome and output directory for the batch workflow.",
        ],
        "outputs": [
            "A batch HLA typing output tree containing extraction outputs, SpecHLA results, and merged summaries.",
        ],
        "when_to_use": [
            "You want end-to-end HLA typing from a BAM directory.",
            "You want the batch wrapper instead of manually running extract_hla_read and spechla for each sample.",
        ],
        "use_another_command_when": [
            "Use extract_hla_read and spechla separately when you need finer control over each stage.",
            "Use trust4 when the task is TCR/BCR repertoire reconstruction rather than HLA typing.",
        ],
        "project_context_notes": [
            "This workflow is resume-friendly and writes per-sample completion markers for extraction and SpecHLA stages.",
            "A successful batch run normally produces a merged `hla_result_merged.txt` alongside the per-sample result folders.",
            "This workflow is for BAM-directory HLA genotyping, not for TPM matrices or general TME profiling.",
        ],
    },
    "tme_profile": {
        "summary": "Run the downstream TME workflow from an existing TPM-like matrix: signature scoring, six default TME methods, merged deconvolution, and LR_cal.",
        "detailed_description": "This is the TPM-to-TME entrypoint. It does not quantify FASTQ files and it does not accept a runall output directory as its main input; it expects an existing genes-by-samples TPM-like expression matrix.",
        "input_expectations": [
            "Input must be an expression matrix file, typically TPM-like, with genes as rows and samples as columns.",
            "This command is for downstream TME analysis only and is not the entrypoint for raw FASTQ data.",
        ],
        "outputs": [
            "A standardized output directory containing signature scores, deconvolution results, merged TME outputs, and LR_cal output.",
        ],
        "when_to_use": [
            "You already have a TPM matrix and want the standard downstream TME workflow.",
            "You want signature scoring, CIBERSORT, IPS, ESTIMATE, MCPcounter, quanTIseq, EPIC, merged deconvolution, and LR_cal from one matrix.",
        ],
        "use_another_command_when": [
            "Use runall when you are starting from raw FASTQ data.",
            "Use standalone commands such as calculate_sig_score or cibersort when you want only one downstream step instead of the whole bundle.",
        ],
        "project_context_notes": [
            "A standard `tme_profile` output tree contains `01-signatures/`, `02-tme/`, and `03-LR_cal/`.",
            "The default bundled TME branch includes `cibersort`, `IPS`, `estimate`, `mcpcounter`, `quantiseq`, and `epic`, plus a merged deconvolution table such as `deconvo_merged.csv`.",
            "This bundled workflow does not include BayesPrism, `tme_cluster`, or `nmf`; those are optional follow-up branches.",
            "The orchestrator also applies method-specific defaults unless the user overrides them: `estimate --platform affymetrix`, `mcpcounter --features HUGO_symbols`, `quantiseq --arrays --tumor --scale_mrna`, and `epic --reference TRef`.",
            "Because the bundled workflow already writes concrete per-method files such as `cibersort_results.csv`, `quantiseq_results.csv`, and `deconvo_merged.csv`, those should be named explicitly when the agent recommends downstream reuse or clustering input candidates.",
        ],
    },
    "runall": {
        "summary": "Run the full end-to-end FASTQ workflow from raw reads through QC, quantification, TPM generation, downstream TME analysis, LR_cal, and TRUST4 outputs.",
        "detailed_description": "runall is not just a quantification wrapper. In salmon mode it performs fastq_qc, batch_salmon, merge_salmon, prepare_salmon, downstream signature scoring, six default TME methods, merged deconvolution, LR_cal, and TRUST4 outputs; in star mode it swaps in the STAR/count-to-TPM branch before the same downstream analysis.",
        "input_expectations": [
            "Input is a raw FASTQ directory plus a Salmon or STAR index, depending on mode.",
            "This command is the FASTQ entrypoint and not the right command when you already have a TPM matrix.",
        ],
        "outputs": [
            "A standardized analysis tree with QC, quantification, TPM, signature scores, TME outputs, LR_cal results, and TRUST4 outputs.",
        ],
        "when_to_use": [
            "You want the main FASTQ-to-TME workflow with as little manual orchestration as possible.",
            "You want IOBRpy to take you from raw reads to TPM, downstream TME analysis, and repertoire-related outputs in one pipeline.",
        ],
        "use_another_command_when": [
            "Use tme_profile when you already have a TPM-like matrix and only need the downstream TME stage.",
            "Use batch_salmon or batch_star_count when you only want the quantification/counting stage rather than the full pipeline.",
        ],
        "project_context_notes": [
            "Salmon-mode runall typically builds `01-qc/`, `02-salmon/`, `03-tpm/`, `04-signatures/`, `05-tme/`, and `06-LR_cal/`.",
            "STAR-mode runall uses the same downstream layout but replaces the quantification branch with `02-star/` and `count2tpm` before the log2 TPM matrix is created.",
            "The final TPM matrix in a standard runall tree is typically `03-tpm/tpm_matrix.csv`, which is the log2(x+1) matrix produced after TPM generation.",
            "The downstream TME bundle still excludes BayesPrism and clustering, so those remain optional follow-up branches even after runall completes.",
            "Top-level `--threads` and `--batch_size` are propagated into multiple substeps, so the harness should treat them as workflow-wide concurrency controls rather than one-step-only flags.",
            "The workflow is resume-aware and uses per-stage completion markers, including TRUST4 completion under `07-TCRBCR/`, so the agent should prefer resume/continue reasoning when outputs already exist.",
            "The final TRUST4 branch is mode-aware: STAR-mode runall feeds TRUST4 from the STAR result root, while Salmon-mode runall feeds TRUST4 from the cleaned FASTQ/QC directory.",
        ],
    },
    "bayesprism": {
        "summary": "Run BayesPrism Bayesian deconvolution on a bulk expression matrix, using the built-in reference by default or a custom single-cell reference when needed.",
        "detailed_description": "BayesPrism is a standalone Bayesian deconvolution workflow. Unlike the default six-method bundle in tme_profile and runall, it writes its own dedicated output directory. IOBRpy can use its bundled BayesPrism reference setup by default; a custom single-cell reference is only needed when you intentionally want to override that built-in reference.",
        "input_expectations": [
            "Input is a bulk expression matrix with genes as rows and samples as columns.",
            "The bulk matrix is usually a TPM-like matrix that is already ready for downstream analysis.",
            "When you intentionally use a custom single-cell reference, the related cell-state, cell-type, and key arguments become important.",
        ],
        "outputs": [
            "A BayesPrism result directory containing theta-style outputs and related deconvolution artifacts.",
        ],
        "when_to_use": [
            "You specifically want Bayesian deconvolution rather than the default six-method downstream bundle.",
            "You want to add an optional BayesPrism branch on top of an existing TPM-ready analysis state.",
            "You want to use the bundled reference by default or override it with a cohort-specific single-cell reference.",
        ],
        "use_another_command_when": [
            "Use tme_profile or runall when you want the standard six-method downstream TME bundle instead of BayesPrism alone.",
            "Use quantiseq, EPIC, or CIBERSORT directly when you want one of the lighter-weight default deconvolution methods.",
        ],
        "project_context_notes": [
            "BayesPrism is an optional standalone branch rather than one of the default bundled TME methods in `tme_profile` or `runall`.",
            "IOBRpy already has a bundled BayesPrism reference setup, so a custom `sc_dat` reference is optional rather than mandatory.",
            "Typical outputs include `theta.csv`, `theta_cv.csv`, and `Z_tumor.csv` under the BayesPrism output directory.",
            "This makes BayesPrism a good optional backfill step when the main directory already has the six default TME methods but still lacks the Bayesian branch.",
        ],
    },
}


def _important(
    parameter: str,
    why_it_matters: str,
    native_default_behavior: str = "",
    ask_when: str = "",
) -> Dict[str, Any]:
    return {
        "parameter": parameter,
        "why_it_matters": why_it_matters,
        "native_default_behavior": native_default_behavior,
        "ask_when": ask_when,
    }


IMPORTANT_OPTIONAL_PARAMETERS: Dict[str, List[Dict[str, Any]]] = {
    "prepare_salmon": [
        _important(
            "return_feature",
            "This decides which identifier space the downstream TPM matrix uses, such as ENST, ENSG, or gene symbols.",
            "If omitted, the native command keeps its built-in default feature choice.",
            "Confirm when downstream analysis expects a specific identifier type.",
        ),
        _important(
            "remove_version",
            "Version suffix handling changes whether identifiers like ENSG000001.1 are normalized before later steps.",
            "If omitted, version stripping stays at the native default.",
            "Confirm when the input matrix mixes versioned and non-versioned identifiers.",
        ),
    ],
    "count2tpm": [
        _important(
            "effLength_csv",
            "A custom effective-length table changes the TPM denominator and can be necessary for non-default references.",
            "If omitted, the command falls back to built-in annotation sources.",
            "Confirm when the counts come from a custom annotation or non-standard reference build.",
        ),
        _important(
            "idtype",
            "The gene identifier type controls how lengths and annotations are resolved.",
            "If omitted, the command keeps its native default identifier assumption.",
            "Confirm when the count matrix uses Ensembl IDs, symbols, or another identifier system.",
        ),
        _important(
            "remove_version",
            "Removing version suffixes can be required for Ensembl-style identifiers to match annotation resources.",
            "If omitted, version handling stays at the native default.",
            "Confirm when the input has identifiers like ENSG000001.1.",
        ),
    ],
    "anno_eset": [
        _important(
            "annotation",
            "This chooses the built-in annotation source and changes how identifiers are mapped.",
            "If omitted, the command uses its native default annotation mode.",
            "Confirm when the matrix comes from a specific platform or probe set.",
        ),
        _important(
            "annotation-file",
            "A custom annotation file overrides the bundled annotation mapping and can materially change the output matrix.",
            "If omitted, only built-in annotation sources are used.",
            "Confirm when the data uses custom probes or project-specific annotations.",
        ),
        _important(
            "method",
            "Duplicate-resolution strategy affects how multiple probes or identifiers collapse into one gene value.",
            "If omitted, the command keeps its native duplicate-handling default.",
            "Confirm when mean, sum, or variability-based collapsing changes interpretation.",
        ),
    ],
    "calculate_sig_score": [
        _important(
            "signature",
            "This controls which signature families are scored and can dramatically change runtime and output breadth.",
            "If omitted, the command uses only its native default signature selection.",
            "Confirm when the user wants a subset instead of the full built-in collection.",
        ),
        _important(
            "method",
            "PCA, z-score, ssGSEA, and integration produce meaningfully different signature summaries.",
            "If omitted, the command keeps its native scoring default.",
            "Confirm when the user cares about the scoring formulation or comparability with earlier analyses.",
        ),
        _important(
            "mini_gene_count",
            "Low gene-count thresholds can retain sparse signatures, while stricter thresholds drop incomplete signatures.",
            "If omitted, the command keeps its native minimum-gene threshold.",
            "Confirm when signatures are expected to be partially represented in the matrix.",
        ),
    ],
    "cibersort": [
        _important(
            "perm",
            "Permutation count changes runtime and the robustness of significance estimates.",
            "If omitted, the command uses its native permutation default.",
            "Confirm when the user cares about statistical stability versus turnaround time.",
        ),
        _important(
            "QN",
            "Quantile normalization can be appropriate for some array-like inputs and inappropriate for RNA-seq TPM-like inputs.",
            "If omitted, the command uses its native QN default.",
            "Confirm when the matrix origin is not obvious or when matching a published workflow matters.",
        ),
        _important(
            "absolute",
            "Absolute mode changes the interpretation of the output from relative fractions toward absolute scores.",
            "If omitted, the command stays in its native default mode.",
            "Confirm when the user needs absolute rather than relative abundance output.",
        ),
    ],
    "IPS": [],
    "estimate": [
        _important(
            "platform",
            "ESTIMATE platform handling affects how the matrix is normalized for scoring.",
            "If omitted, the command uses its native default platform assumption.",
            "Confirm when the data did not come from the default expression platform.",
        ),
    ],
    "mcpcounter": [
        _important(
            "features",
            "MCPcounter needs the correct feature namespace so it can interpret row identifiers correctly.",
            "If omitted, the command uses its native default feature namespace.",
            "Confirm when the matrix uses Ensembl IDs, gene symbols, or another non-default identifier set.",
        ),
    ],
    "quantiseq": [
        _important(
            "arrays",
            "Array-style preprocessing changes how the expression matrix is normalized before deconvolution.",
            "If omitted, the command keeps its native preprocessing default.",
            "Confirm when the input is not a standard RNA-seq TPM-like matrix.",
        ),
        _important(
            "tumor",
            "Tumor-specific handling changes how quanTIseq interprets and rescales the mixture.",
            "If omitted, the command uses its native tumor-mode default.",
            "Confirm when the samples are non-tumor tissue or mixed study types.",
        ),
        _important(
            "scale_mrna",
            "mRNA scaling affects the final reported composition values.",
            "If omitted, the command uses its native scaling default.",
            "Confirm when the user needs outputs aligned with a particular publication or prior pipeline.",
        ),
    ],
    "epic": [
        _important(
            "reference",
            "The EPIC reference selection changes which cell reference profiles are used during deconvolution.",
            "If omitted, the command uses its native default reference.",
            "Confirm when the cohort matches a specific tumor or blood-reference setting.",
        ),
    ],
    "tme_cluster": [
        _important(
            "features",
            "Feature selection changes which TME signals drive the clustering result.",
            "If omitted, the command uses all eligible features. Feature counting excludes the sample ID column, so a common tutorial-style range such as `1:22` refers to the second through twenty-third actual CSV columns.",
            "Confirm when the user wants to cluster on a targeted subset of signatures or cell scores, or when the matrix has many hundreds or thousands of feature columns.",
        ),
        _important(
            "pattern",
            "Pattern-based feature selection can include or exclude large parts of the matrix automatically.",
            "If omitted, no pattern filter is applied.",
            "Confirm when only a prefix or feature family should define the clustering space.",
        ),
        _important(
            "min_nc",
            "The lower bound for cluster-number search changes the solution space.",
            "If omitted, the command keeps its native clustering-range default.",
            "Confirm when the user has prior expectations for the number of TME subtypes.",
        ),
        _important(
            "max_nc",
            "The upper bound for cluster-number search changes runtime and the possible subtype granularity.",
            "If omitted, the command keeps its native clustering-range default.",
            "Confirm when broad versus fine-grained stratification matters.",
        ),
    ],
    "LR_cal": [
        _important(
            "data_type",
            "Ligand-receptor scoring depends on whether the matrix should be interpreted as TPM or another expression scale.",
            "If omitted, the command uses its native default data-type assumption.",
            "Confirm when the input is not standard TPM-like expression.",
        ),
        _important(
            "id_type",
            "Identifier type determines how ligand and receptor genes are matched.",
            "If omitted, the command uses its native default identifier type.",
            "Confirm when the matrix rows are Ensembl IDs instead of symbols.",
        ),
        _important(
            "cancer_type",
            "Cancer type changes the ligand-receptor context and can alter the interaction network used for scoring.",
            "If omitted, the command uses its native default cancer-type context.",
            "Confirm when the user wants cohort-specific rather than pan-cancer behavior.",
        ),
    ],
    "nmf": [
        _important(
            "kmin",
            "The lower rank bound changes the factorization search space and the smallest number of latent programs considered.",
            "If omitted, the command keeps its native lower-rank default.",
            "Confirm when the user has a hypothesis about the minimum number of clusters or latent programs.",
        ),
        _important(
            "kmax",
            "The upper rank bound changes runtime and the largest factorization model explored.",
            "If omitted, the command keeps its native upper-rank default.",
            "Confirm when a narrower or wider model-selection sweep is desired.",
        ),
        _important(
            "features",
            "Feature selection decides which variables drive the NMF decomposition.",
            "If omitted, the command uses its native feature-selection behavior. Feature counting excludes the first sample-ID column, so tutorial-style ranges such as `1:22` refer to the second through twenty-third actual CSV columns.",
            "Confirm when only selected signatures or TME scores should be included, or when the matrix is so wide that all-column NMF would be hard to justify automatically.",
        ),
    ],
    "mouse2human_eset": [
        _important(
            "is_matrix",
            "Matrix mode changes how the command reads the file and which columns are interpreted as identifiers.",
            "If omitted, the command keeps its native input-mode default.",
            "Confirm when the input is a general table rather than a strict genes-by-samples matrix.",
        ),
        _important(
            "column_of_symbol",
            "For table mode, this decides which column contains the source gene symbols to convert.",
            "If omitted, the command uses its native default symbol column handling.",
            "Confirm when the symbol column is not in the usual position or name.",
        ),
    ],
    "batch_salmon": [
        _important(
            "suffix1",
            "Read-pair suffix handling determines how FASTQ files are matched into samples.",
            "If omitted, the command uses its native suffix assumption.",
            "Confirm when the FASTQ naming scheme is not the default one expected by IOBRpy.",
        ),
        _important(
            "num_threads",
            "Thread count directly affects runtime and cluster-resource usage.",
            "If omitted, the command uses its native threading default.",
            "Confirm when running on shared compute or under scheduler limits.",
        ),
        _important(
            "gtf",
            "A custom GTF can change how quantification outputs are interpreted or annotated.",
            "If omitted, the command runs without an explicit custom GTF override.",
            "Confirm when the reference annotation version must match the index exactly.",
        ),
    ],
    "merge_salmon": [
        _important(
            "project",
            "Project naming affects output filenames and how merged matrices are labeled.",
            "If omitted, the command uses its native project-name default.",
            "Confirm when the user wants stable project-specific naming.",
        ),
        _important(
            "num_processes",
            "Process count changes merge throughput and compute usage.",
            "If omitted, the command uses its native multiprocessing default.",
            "Confirm when the directory is large or the environment has process-count limits.",
        ),
    ],
    "merge_star_count": [
        _important(
            "project",
            "Project naming affects output filenames and merged count-matrix labels.",
            "If omitted, the command uses its native project-name default.",
            "Confirm when the user wants stable project-specific naming.",
        ),
    ],
    "batch_star_count": [
        _important(
            "suffix1",
            "Read-pair suffix handling determines how FASTQ files are matched into samples.",
            "If omitted, the command uses its native suffix assumption.",
            "Confirm when the FASTQ naming scheme is not the default expected by the STAR wrapper.",
        ),
        _important(
            "num_threads",
            "Thread count directly affects runtime and cluster-resource usage.",
            "If omitted, the command uses its native threading default.",
            "Confirm when running on shared compute or under scheduler limits.",
        ),
        _important(
            "batch_size",
            "Batch size changes how many samples are processed together and can affect memory usage.",
            "If omitted, the command uses its native batching default.",
            "Confirm when the FASTQ directory is large or per-node memory is limited.",
        ),
    ],
    "fastq_qc": [
        _important(
            "suffix1",
            "Read naming affects how paired FASTQ files are discovered and matched.",
            "If omitted, the command uses its native suffix assumption.",
            "Confirm when the FASTQ naming scheme is not the default expected by fastq_qc.",
        ),
        _important(
            "se",
            "Single-end mode changes how FASTQ files are interpreted and prevents paired-end matching.",
            "If omitted, the command stays in paired-end mode.",
            "Confirm when the input is single-end rather than paired-end sequencing.",
        ),
        _important(
            "length_required",
            "Minimum post-trim read length affects how aggressively reads are filtered.",
            "If omitted, the command uses its native quality-filtering default.",
            "Confirm when short-read retention versus strict QC matters.",
        ),
    ],
    "log2_eset": [],
    "trust4": [
        _important(
            "ref",
            "TRUST4 reference selection changes repertoire reconstruction behavior and reference matching.",
            "If omitted, the command uses its native reference default.",
            "Confirm when the cohort requires a specific TRUST4 reference bundle.",
        ),
        _important(
            "fqdir",
            "FASTQ-directory mode changes the input discovery strategy for batch repertoire analysis.",
            "If omitted, the command relies on whichever direct BAM or FASTQ inputs were supplied.",
            "Confirm when the user wants directory-scale processing instead of one sample at a time.",
        ),
        _important(
            "t",
            "Thread count changes runtime and scheduler footprint.",
            "If omitted, the command uses its native threading default.",
            "Confirm when running on shared compute or when the user has a specific CPU allocation.",
        ),
        _important(
            "barcode",
            "Barcode handling is critical for single-cell or barcode-aware repertoire workflows.",
            "If omitted, the command runs without explicit barcode configuration.",
            "Confirm when the input is barcode-tagged or UMI-aware data.",
        ),
        _important(
            "UMI",
            "UMI handling changes molecule collapsing and repertoire quantification behavior.",
            "If omitted, the command uses its native non-UMI default.",
            "Confirm when the data comes from UMI-tagged protocols.",
        ),
    ],
    "spechla": [
        _important(
            "threads",
            "Thread count changes runtime and compute usage for per-sample HLA typing.",
            "If omitted, the command uses its native threading default.",
            "Confirm when the user has a specific CPU allocation.",
        ),
        _important(
            "use-exon",
            "SpecHLA mode changes whether the workflow behaves like exon/RNA typing or WGS-style typing.",
            "If omitted, the command uses its native mode default.",
            "Confirm when the input data type is RNA/exome rather than WGS.",
        ),
    ],
    "extract_hla_read": [
        _important(
            "ref",
            "Reference selection changes how HLA reads are extracted from BAM or CRAM input.",
            "If omitted, the command uses its native reference default.",
            "Confirm when the alignment reference differs from the usual setup.",
        ),
        _important(
            "no-auto-install",
            "Disabling auto-install changes whether missing helper assets are provisioned automatically.",
            "If omitted, the command follows its native installation behavior.",
            "Confirm when running on locked-down systems without install permissions.",
        ),
    ],
    "hla_typing": [
        _important(
            "threads",
            "Thread count affects throughput across batch HLA typing runs.",
            "If omitted, the command uses its native threading default.",
            "Confirm when the batch runs on shared compute or scheduler-managed nodes.",
        ),
        _important(
            "use-exon",
            "This switches the SpecHLA mode used for the batch workflow and can change biological interpretation.",
            "If omitted, the command uses its native typing-mode default.",
            "Confirm when the input data is exon/RNA-like rather than WGS-like.",
        ),
    ],
    "tme_profile": [
        _important(
            "threads",
            "Thread count affects downstream scoring and deconvolution runtime.",
            "If omitted, the command keeps its native threading default.",
            "Confirm when running on shared compute or when the user already knows the CPU allocation.",
        ),
        _important(
            "signature",
            "This controls which signature collections are scored before the deconvolution bundle runs.",
            "If omitted, tme_profile forwards its native default signature selection.",
            "Confirm when the user wants only selected signatures instead of the full built-in collection.",
        ),
        _important(
            "method",
            "Signature scoring method affects the meaning of the signature output even when the rest of the TME bundle stays the same.",
            "If omitted, tme_profile keeps its native default scoring method.",
            "Confirm when the user needs PCA, z-score, ssGSEA, or integration specifically.",
        ),
        _important(
            "platform",
            "ESTIMATE platform handling matters when the matrix does not match the native default platform assumption.",
            "If omitted, tme_profile uses its native ESTIMATE platform default.",
            "Confirm when the input expression platform is known and not the usual default.",
        ),
        _important(
            "features",
            "MCPcounter feature namespace must match the matrix row identifiers for correct abundance scoring.",
            "If omitted, tme_profile uses its native feature-namespace default.",
            "Confirm when the matrix uses Ensembl IDs rather than gene symbols.",
        ),
        _important(
            "reference",
            "EPIC reference choice changes one of the default downstream deconvolution branches.",
            "If omitted, tme_profile uses its native EPIC reference default.",
            "Confirm when the cohort fits a specific EPIC reference setting.",
        ),
        _important(
            "cancer_type",
            "Cancer type changes LR_cal context and can alter ligand-receptor scoring behavior.",
            "If omitted, tme_profile uses its native pan-cancer context.",
            "Confirm when the user wants tumor-type-specific LR_cal behavior.",
        ),
    ],
    "runall": [
        _important(
            "threads",
            "Thread count affects QC, quantification, and downstream analysis runtime across the whole pipeline.",
            "If omitted, runall uses its native unified threading default.",
            "Confirm when running on shared compute or when the user already has a CPU allocation in mind.",
        ),
        _important(
            "batch_size",
            "Batch size affects how many samples move through QC and quantification together, which changes memory pressure and scheduling behavior.",
            "If omitted, runall keeps its native batch-processing default.",
            "Confirm when the FASTQ directory is large or the environment has memory constraints.",
        ),
        _important(
            "suffix1",
            "Read-pair suffix handling determines how FASTQ files are matched into samples before QC and quantification.",
            "If omitted, runall assumes its native FASTQ naming pattern.",
            "Confirm when filenames use `_1.fastq.gz`, `_R1_001.fastq.gz`, or another non-default pattern.",
        ),
        _important(
            "se",
            "Single-end mode changes the entire pipeline's FASTQ interpretation.",
            "If omitted, runall stays in paired-end mode.",
            "Confirm when the sequencing data is single-end rather than paired-end.",
        ),
        _important(
            "signature",
            "This controls which signature families are scored in the downstream TME stage.",
            "If omitted, runall forwards its native default signature selection.",
            "Confirm when the user wants only selected signatures instead of the full collection.",
        ),
        _important(
            "platform",
            "ESTIMATE platform handling matters for one branch of the downstream TME analysis bundle.",
            "If omitted, runall uses its native ESTIMATE platform default.",
            "Confirm when the expression platform is known and differs from the usual default assumption.",
        ),
        _important(
            "features",
            "MCPcounter feature namespace must match the matrix row identifiers produced by the TPM step.",
            "If omitted, runall keeps its native MCPcounter feature default.",
            "Confirm when downstream matrices use Ensembl IDs rather than symbols.",
        ),
        _important(
            "reference",
            "EPIC reference choice changes one of the default downstream deconvolution branches.",
            "If omitted, runall uses its native EPIC reference default.",
            "Confirm when the cohort fits a specific EPIC reference setting.",
        ),
        _important(
            "cancer_type",
            "Cancer type changes LR_cal context and can materially change ligand-receptor scoring outputs.",
            "If omitted, runall uses its native pan-cancer context.",
            "Confirm when the user wants disease-specific LR_cal behavior rather than pan-cancer defaults.",
        ),
    ],
    "bayesprism": [
        _important(
            "threads",
            "Thread count affects runtime for the Bayesian deconvolution workflow.",
            "If omitted, the command uses its native threading default.",
            "Confirm when the user has explicit CPU or scheduler limits.",
        ),
        _important(
            "sc_dat",
            "A custom single-cell reference fundamentally changes the BayesPrism model input.",
            "If omitted, the command uses its native built-in reference behavior.",
            "Confirm when the user has a cohort-specific single-cell reference.",
        ),
        _important(
            "cell_state_labels",
            "Cell-state labels decide how the single-cell reference is interpreted and aggregated.",
            "If omitted, the command uses its native label handling.",
            "Confirm when the reference object contains multiple state-resolution label columns.",
        ),
        _important(
            "cell_type_labels",
            "Cell-type labels determine the biological grouping used by the deconvolution model.",
            "If omitted, the command uses its native label handling.",
            "Confirm when the reference object has several candidate cell-type label fields.",
        ),
    ],
}


for _command_name, _important_parameters in IMPORTANT_OPTIONAL_PARAMETERS.items():
    COMMAND_PROFILES.setdefault(_command_name, {})["important_optional_parameters"] = _important_parameters
