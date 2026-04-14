"""
Pipeline stage mapping helpers for agent-facing directory introspection.
"""

from __future__ import annotations

import fnmatch
import gzip
import json
import os
import pickle
import re
from collections import deque
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Set, Tuple

from .pipeline_map_render import (
    _agent_rendering_hints_payload,
    _enrich_workflow_checklist_payload,
    _filter_item_specific_evidence,
    _helper_or_nonresult_file,
    _looks_like_raw_fastq_path,
    _recommended_deep_focus_roots,
    _scan_checklist_display_hints,
    _scan_limit_messages,
    _terminal_stage_map,
    format_pipeline_map_report,
)
from .pipeline_map_scenario import (
    _classify_scenario,
    _next_nodes,
    _recommended_action,
    _roadmap_progress_text,
    _scenario_payload,
)
from .pipeline_map_scan import (
    _DEFAULT_MAP_MAX_DEPTH,
    _DEFAULT_MAP_MAX_ENTRIES,
    _QUICK_MIN_MAX_DEPTH,
    _QUICK_MIN_MAX_ENTRIES,
    _build_entry_matcher,
    _collect_stage_matches,
    _is_hla_intermediate_read_path,
    _matches_any_pattern,
    _normalize_pattern,
    _normalize_relpath,
    _normalize_token,
    _path_has_exact_part,
    _path_has_keyword,
    _path_has_part_token,
    _path_has_part_tokens,
    _path_parts_lower,
    _part_tokens,
    _scan_path_entries,
    _suggest_initial_scan_limits,
)


_STAGE_DISPLAY_ORDER = [
    "raw_fastq",
    "fastq_qc",
    "salmon_quant",
    "salmon_merge",
    "prepare_salmon",
    "star_quant",
    "star_merge",
    "count2tpm",
    "tpm_matrix",
    "signature_scoring",
    "deconvolution",
    "ligand_receptor",
    "clustering",
    "trust4",
    "spechla",
    "hla_typing",
]

_RESUME_CHOICES = [
    "continue_downstream",
    "rerun_current_stage",
    "rerun_full_pipeline",
]

_STAGE_TO_NODE = {
    "raw_fastq": "Input Data & QC",
    "fastq_qc": "Input Data & QC",
    "salmon_quant": "Salmon Mode",
    "salmon_merge": "Salmon Mode",
    "prepare_salmon": "Salmon Mode",
    "star_quant": "STAR Mode",
    "star_merge": "STAR Mode",
    "count2tpm": "STAR Mode",
    "tpm_matrix": "TPM Matrix",
    "signature_scoring": "Signature Scoring",
    "deconvolution": "Deconvolution",
    "ligand_receptor": "Ligand-Receptor",
    "clustering": "Clustering",
    "trust4": "TRUST4",
    "spechla": "SpecHLA",
    "hla_typing": "HLA Genotypes",
}

_NON_IMPLIED_STAGE_IDS = {"raw_fastq", "fastq_qc", "salmon_quant"}
_WEAK_STAGE_REQUIRES_CONTENT = {"trust4"}

_NEXT_COMMANDS = {
    "fastq_qc": "iobrpy-cli fastq_qc --path1_fastq <path_fastq_dir> --path2_fastp <path_qc_outdir> --num_threads <num_threads> --batch_size <batch_size>",
    "salmon_quant": "iobrpy-cli batch_salmon --index <path_salmon_index> --path_fq <path_fastq_dir> --path_out <path_salmon_outdir> --batch_size <batch_size> --num_threads <num_threads>",
    "salmon_merge": "iobrpy-cli merge_salmon --path_salmon <path_salmon_outdir> --project <project_name> --num_processes <num_processes>",
    "prepare_salmon": "iobrpy-cli prepare_salmon --input <path_merged_salmon_tpm> --output <path_tpm_matrix>",
    "star_quant": "iobrpy-cli batch_star_count --index <path_star_index> --path_fq <path_fastq_dir> --path_out <path_star_outdir> --batch_size <batch_size> --num_threads <num_threads>",
    "star_merge": "iobrpy-cli merge_star_count --path <path_star_outdir> --project <project_name>",
    "count2tpm": "iobrpy-cli count2tpm --input <path_count_matrix> --output <path_tpm_matrix> --idtype <ensembl|entrez|symbol|mgi> --org <hsa|mmus>",
    "tme_profile": "iobrpy-cli tme_profile --input <path_tpm_matrix> --output <path_tme_profile_outdir> --threads <threads>",
    "signature_scoring": "iobrpy-cli calculate_sig_score --input <path_tpm_matrix> --output <path_signature_scores> --signature <signature_group>",
    "deconvolution": "iobrpy-cli cibersort --input <path_tpm_matrix> --output <path_cibersort_results>",
    "ligand_receptor": "iobrpy-cli LR_cal --input <path_tpm_matrix> --output <path_lr_scores> --cancer_type <cancer_type>",
    "clustering": "iobrpy-cli tme_cluster --input <path_tme_matrix> --output <path_tme_cluster_outdir> --features <feature_columns>",
    "trust4": "iobrpy-cli trust4 --fqdir <path_fastq_dir> -o <path_trust4_outdir> -t <threads>",
    "hla_typing": "iobrpy-cli hla_typing -b <path_bam_dir> -r <hg19|hg38> -o <path_hla_outdir> -j <threads>",
}

_SCENARIO_LABELS = {
    "raw_fastq_only": "Raw FASTQ only",
    "upstream_partial": "Upstream processing already started",
    "mixed_raw_and_processed": "Raw input plus existing processed outputs",
    "tpm_ready": "TPM-ready for downstream TME analysis",
    "downstream_partial": "Partial downstream TME outputs already exist",
    "downstream_complete": "Most downstream TME outputs already exist",
    "immune_only": "Immune-specific outputs detected",
    "unknown": "Unclassified directory",
}

_CHOICE_LABELS = {
    "continue_downstream": "Continue downstream analysis",
    "rerun_current_stage": "Rerun the current stage",
    "rerun_full_pipeline": "Rerun the full pipeline in a fresh output directory",
    "start_pipeline": "Start an IOBRpy workflow from the beginning",
    "inspect_commands": "Inspect commands and choose a workflow manually",
}

_QUICK_MAX_PREVIEW_FILES = 192
_FOCUSED_DEEP_MAX_PREVIEW_FILES = 1024
_FULL_MAX_PREVIEW_FILES = 2048

_EXTERNAL_ANALYSIS_HINTS = [
    {
        "id": "external_tcr_bcr",
        "label": "External TCR/BCR result hint",
        "patterns": [
            "*mixcr*",
            "*clonotype*",
            "*clones*.txt",
            "*clones*.tsv",
            "*vdj*",
            "*tcr*",
            "*bcr*",
            "*airr*",
            "*cdr3*",
            "*contig*annotation*",
            "*consensus_annotation*",
            "*rearrangement*",
            "*immunarch*",
            "*immunearch*",
        ],
        "description": "Detected external TCR/BCR result-like files or directories that should not be assumed to be iobrpy TRUST4 outputs.",
    },
    {
        "id": "external_hla",
        "label": "External HLA analysis hint",
        "patterns": [
            "*-hla",
            "*-hla-*",
            "*_hla",
            "*_hla_*",
            "*optitype*",
            "*arcas*",
            "*hlahd*",
            "*xhla*",
            "*hlab*",
            "*genotype*h*la*",
        ],
        "description": "Detected HLA-related files or directories that should not be assumed to be iobrpy HLA outputs.",
    },
]

_DECONVOLUTION_METHOD_SPECS: List[Dict[str, Any]] = [
    {
        "id": "cibersort",
        "label": "CIBERSORT",
        "patterns": ["*cibersort*.csv", "*cibersort*.txt"],
        "marker_groups": [
            ["p-value", "rmse"],
            ["b cells naive", "t cells cd8"],
        ],
    },
    {
        "id": "epic",
        "label": "EPIC",
        "patterns": ["*epic*.csv", "*epic*.txt"],
        "marker_groups": [
            ["bcells", "cd4_tcells", "cd8_tcells"],
            ["endothelial", "macrophages", "tumor"],
        ],
    },
    {
        "id": "quantiseq",
        "label": "quanTIseq",
        "patterns": ["*quantiseq*.csv", "*quantiseq*.txt"],
        "marker_groups": [
            ["m1 macrophages", "m2 macrophages", "tregs"],
            ["b cells", "neutrophils", "nk cells"],
        ],
    },
    {
        "id": "mcpcounter",
        "label": "MCPcounter",
        "patterns": ["*mcpcounter*.csv", "*mcpcounter*.txt"],
        "marker_groups": [
            ["cytotoxic lymphocytes", "b lineage", "fibroblasts"],
            ["cd8 t cells", "nk cells", "endothelial cells"],
        ],
    },
    {
        "id": "estimate",
        "label": "ESTIMATE",
        "patterns": ["*estimate*.csv", "*estimate*.txt"],
        "marker_groups": [
            ["stromalscore", "immunescore", "estimatescore"],
            ["tumorpurity"],
        ],
    },
    {
        "id": "ips",
        "label": "IPS",
        "patterns": ["*ips*.csv", "*ips*.txt"],
        "marker_groups": [
            ["az", "cp", "ec", "mhc", "sc"],
            ["immunophenoscore"],
        ],
    },
    {
        "id": "bayesprism",
        "label": "BayesPrism",
        "patterns": ["*bayesprism*.csv", "*bayesprism*.txt", "*bayesprism*"],
        "marker_groups": [
            ["bayesprism"],
            ["cell_type_fraction", "theta"],
        ],
    },
]

_WRAPPER_DEFAULT_DECONVOLUTION_IDS: Tuple[str, ...] = (
    "cibersort",
    "ips",
    "estimate",
    "mcpcounter",
    "quantiseq",
    "epic",
)
_WRAPPER_DEFAULT_DECONVOLUTION_ID_SET = set(_WRAPPER_DEFAULT_DECONVOLUTION_IDS)
_IPS_SUFFIX_COLUMNS = {"az_ips", "cp_ips", "ec_ips", "mhc_ips", "sc_ips", "ips_ips"}
_DECONVOLUTION_FUNCTION_ID_BY_CANONICAL = {
    "cibersort": "cibersort",
    "ips": "IPS",
    "estimate": "estimate",
    "mcpcounter": "mcpcounter",
    "quantiseq": "quantiseq",
    "epic": "epic",
    "bayesprism": "bayesprism",
}

_FUNCTION_STATUS_ORDER = {
    "not_detected": 0,
    "reusable_result": 1,
    "likely_iobrpy_result": 2,
    "confirmed_by_content": 3,
}

_STRICT_IOBRPY_STAGE_IDS = {
    "signature_scoring",
    "deconvolution",
    "ligand_receptor",
    "clustering",
    "trust4",
    "spechla",
    "hla_typing",
}

_STRICT_IOBRPY_INVESTIGATION_STAGE_IDS = {
    "signature_scoring",
    "deconvolution",
    "ligand_receptor",
    "clustering",
}

_STRICT_IOBRPY_INVESTIGATION_ITEM_IDS = {
    "feature_scoring",
    "immune_deconvolution",
    "ligand_receptor_analysis",
    "clustering",
}
_EXTERNAL_STATUS_CHECKLIST_ITEM_IDS = {"hla_typing_summary", "tcr_bcr_summary"}

_FUNCTION_DETECTION_SPECS: List[Dict[str, Any]] = [
    {
        "id": "runall",
        "label": "runall",
        "policy": "iobrpy_wrapper",
        "stage_ids": [],
    },
    {
        "id": "fastq_qc",
        "label": "fastq_qc",
        "policy": "iobrpy",
        "stage_ids": ["fastq_qc"],
        "patterns": ["*_fastp.html", "*_fastp.json", "01-qc/*fastp*"],
        "content_checker": "fastp_report",
    },
    {
        "id": "batch_salmon",
        "label": "batch_salmon",
        "policy": "reusable",
        "stage_ids": ["salmon_quant"],
        "patterns": ["*/quant.sf", "02-salmon/*/quant.sf"],
        "content_checker": "salmon_quant",
    },
    {
        "id": "merge_salmon",
        "label": "merge_salmon",
        "policy": "reusable",
        "stage_ids": ["salmon_merge"],
        "patterns": [
            "*_salmon_tpm.tsv.gz",
            "*_salmon_tpm.tsv",
            "*_salmon_count.tsv.gz",
            "*_salmon_count.tsv",
            "*merged*salmon*tpm*.tsv*",
            "*merged*salmon*count*.tsv*",
            "*tximport*tpm*.rdata",
            "*tximport*count*.rdata",
            "*salmon*tpm*symbol*.rdata",
            "*salmon*count*symbol*.rdata",
        ],
    },
    {
        "id": "prepare_salmon",
        "label": "prepare_salmon",
        "policy": "iobrpy",
        "stage_ids": ["prepare_salmon", "tpm_matrix"],
        "patterns": [
            "prepare_salmon.csv",
            "03-tpm/prepare_salmon.csv",
            "*prepare_salmon*.csv",
            "*prepare_salmon*.tsv",
        ],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["prepare_salmon", "prepare-salmon"],
    },
    {
        "id": "batch_star_count",
        "label": "batch_star_count",
        "policy": "reusable",
        "stage_ids": ["star_quant"],
        "patterns": ["*_ReadsPerGene.out.tab", "*_Aligned.sortedByCoord.out.bam", "02-star/*_ReadsPerGene.out.tab"],
        "content_checker": "star_reads_per_gene",
    },
    {
        "id": "merge_star_count",
        "label": "merge_star_count",
        "policy": "reusable",
        "stage_ids": ["star_merge"],
        "patterns": [
            "*.STAR.count.tsv.gz",
            "*.STAR.count.tsv",
            "*_star_ReadsPerGene.tsv.gz",
            "*_star_ReadsPerGene.tsv",
            "featureCounts.txt",
            "*featureCounts*.txt",
            "*gene_count*.tsv",
            "*gene_count*.csv",
            "*count_matrix*.tsv",
            "*count_matrix*.csv",
            "*counts_matrix*.tsv",
            "*counts_matrix*.csv",
            "*raw_count*.tsv",
            "*raw_count*.csv",
            "*read_count*.tsv",
            "*read_count*.csv",
        ],
        "content_checker": "count_matrix",
        "content_path_keywords": ["readspergene", "featurecount", "featurecounts", "star", "star_count", "star.count"],
    },
    {
        "id": "count2tpm",
        "label": "count2tpm",
        "policy": "iobrpy",
        "stage_ids": ["count2tpm", "tpm_matrix"],
        "patterns": ["count2tpm.csv", "03-tpm/count2tpm.csv", "*count2tpm*.csv", "*count2tpm*.tsv"],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["count2tpm"],
    },
    {
        "id": "anno_eset",
        "label": "anno_eset",
        "policy": "iobrpy",
        "stage_ids": [],
        "patterns": [
            "*anno_eset*.csv",
            "*anno_eset*.tsv",
            "*anno_eset*.txt",
            "*annotated*.csv",
            "*annotated*.tsv",
            "*annotated_eset*.csv",
            "*annotated_eset*.tsv",
        ],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["anno_eset", "annotated"],
    },
    {
        "id": "calculate_sig_score",
        "label": "calculate_sig_score",
        "policy": "iobrpy",
        "stage_ids": ["signature_scoring"],
        "patterns": [
            "calculate_sig_score.csv",
            "*calculate_sig_score*.csv",
            "*calculate_sig_score*.tsv",
            "*signature_score*.csv",
            "*signature_score*.tsv",
            "*sig_score*.csv",
            "*sig_score*.tsv",
            "*ssgsea*.csv",
            "*ssgsea*.tsv",
        ],
        "native_filenames": ["calculate_sig_score.csv", "calculate_sig_score.tsv"],
        "require_content_match": True,
        "name_context": "tme_profile",
        "content_checker": "signature_scoring",
        "content_context": "tme_profile",
    },
    {
        "id": "cibersort",
        "label": "cibersort",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["cibersort_results.csv", "*cibersort*.csv", "*cibersort*.txt"],
        "native_filenames": [
            "cibersort_results.csv",
            "cibersort_results.tsv",
            "cibersort_results.txt",
            "cibersort_results.csv.gz",
            "cibersort_results.tsv.gz",
        ],
        "name_context": "tme_profile",
        "content_checker": "cibersort_output",
        "content_context": "tme_profile",
    },
    {
        "id": "IPS",
        "label": "IPS",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["IPS_results.csv", "*IPS_results*.csv", "*ips_results*.csv", "*immunophenoscore*.csv", "*IPS*.txt"],
        "native_filenames": [
            "IPS_results.csv",
            "ips_results.csv",
            "IPS_results.tsv",
            "ips_results.tsv",
            "IPS_results.txt",
            "IPS_results.csv.gz",
            "IPS_results.tsv.gz",
        ],
        "content_checker": "ips_output",
        "name_context": "tme_profile",
        "content_context": "tme_profile",
    },
    {
        "id": "estimate",
        "label": "estimate",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["estimate_results.csv", "*estimate_results*.csv", "*estimate*.txt"],
        "native_filenames": [
            "estimate_results.csv",
            "estimate_results.tsv",
            "estimate_results.txt",
            "estimate_results.csv.gz",
            "estimate_results.tsv.gz",
        ],
        "require_content_match": True,
        "name_context": "tme_profile",
        "content_checker": "estimate_output",
        "content_context": "tme_profile",
    },
    {
        "id": "mcpcounter",
        "label": "mcpcounter",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["mcpcounter_results.csv", "*mcpcounter_results*.csv", "*mcpcounter*.csv", "*mcpcounter*.txt"],
        "native_filenames": [
            "mcpcounter_results.csv",
            "mcpcounter_results.tsv",
            "mcpcounter_results.txt",
            "mcpcounter_results.csv.gz",
            "mcpcounter_results.tsv.gz",
        ],
        "require_content_match": True,
        "name_context": "tme_profile",
        "content_checker": "mcpcounter_output",
        "content_context": "tme_profile",
    },
    {
        "id": "quantiseq",
        "label": "quantiseq",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["quantiseq_results.csv", "*quantiseq_results*.csv", "*quantiseq*.csv", "*quantiseq*.txt"],
        "native_filenames": [
            "quantiseq_results.csv",
            "quantiseq_results.tsv",
            "quantiseq_results.txt",
            "quantiseq_results.csv.gz",
            "quantiseq_results.tsv.gz",
        ],
        "name_context": "tme_profile",
        "content_checker": "quantiseq_output",
        "content_context": "tme_profile",
    },
    {
        "id": "epic",
        "label": "epic",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["epic_results.csv", "*epic_results*.csv", "*epic*.csv", "*epic*.txt"],
        "native_filenames": [
            "epic_results.csv",
            "epic_results.tsv",
            "epic_results.txt",
            "epic_results.csv.gz",
            "epic_results.tsv.gz",
        ],
        "name_context": "tme_profile",
        "content_checker": "epic_output",
        "content_context": "tme_profile",
    },
    {
        "id": "bayesprism",
        "label": "bayesprism",
        "policy": "iobrpy",
        "stage_ids": ["deconvolution"],
        "patterns": ["*bayesprism*.csv", "*bayesprism*.txt", "*bayesprism*"],
        "native_filenames": ["bayesprism_results.csv", "bayesprism_results.txt"],
        "name_context": "tme_profile",
        "content_checker": "bayesprism_output",
        "content_context": "tme_profile",
    },
    {
        "id": "LR_cal",
        "label": "LR_cal",
        "policy": "iobrpy",
        "stage_ids": ["ligand_receptor"],
        "patterns": [
            "lr_cal.csv",
            "lr_scores.csv",
            "*lr_cal*.csv",
            "*LR_cal*.csv",
            "*lr_scores*.csv",
            "*LR_scores*.csv",
            "*ligand*receptor*.csv",
            "*ligand*receptor*.tsv",
        ],
        "native_filenames": ["lr_cal.csv", "lr_cal.tsv", "LR_cal.csv", "LR_cal.tsv", "lr_scores.csv", "LR_scores.csv"],
        "require_content_match": True,
        "name_context": "tme_profile",
        "content_checker": "ligand_receptor_output",
        "content_context": "tme_profile",
    },
    {
        "id": "tme_cluster",
        "label": "tme_cluster",
        "policy": "iobrpy",
        "stage_ids": ["clustering"],
        "patterns": ["tme_cluster.csv", "*tme_cluster*.csv", "*tme_cluster*.tsv", "*cluster_assignments*.csv"],
        "native_filenames": ["tme_cluster.csv", "tme_cluster.tsv"],
        "name_context": "clustering",
        "content_checker": "clustering_output",
        "content_context": "clustering",
        "content_path_keywords": ["tme_cluster", "cluster_assignments"],
    },
    {
        "id": "nmf",
        "label": "nmf",
        "policy": "iobrpy",
        "stage_ids": ["clustering"],
        "patterns": [
            "*nmf*/clusters.csv",
            "*nmf*clusters*.csv",
            "*top_features_per_cluster*.csv",
            "*pca_plot.png",
            "*nmf*/top_features_per_cluster.csv",
        ],
        "native_filenames": ["clusters.csv", "top_features_per_cluster.csv"],
        "name_context": "clustering",
        "content_checker": "nmf_clusters",
        "content_context": "clustering",
        "content_path_keywords": ["nmf", "top_features_per_cluster"],
    },
    {
        "id": "mouse2human",
        "label": "mouse2human",
        "policy": "iobrpy_wrapper",
        "stage_ids": [],
        "patterns": ["*mouse2human*.csv", "*mouse2human*.tsv", "*mouse2human*.txt"],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["mouse2human"],
    },
    {
        "id": "mouse2human_eset",
        "label": "mouse2human_eset",
        "policy": "iobrpy",
        "stage_ids": [],
        "patterns": ["*mouse2human*.csv", "*mouse2human*.tsv", "*mouse2human*.txt"],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["mouse2human"],
    },
    {
        "id": "log2_eset",
        "label": "log2_eset",
        "policy": "iobrpy",
        "stage_ids": [],
        "patterns": ["*log2_eset*.csv", "*log2_eset*.tsv", "*log2_eset*.txt", "*log2*.csv", "*log2*.tsv"],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["log2_eset", "log2-eset"],
    },
    {
        "id": "remove_version",
        "label": "remove_version",
        "policy": "iobrpy",
        "stage_ids": [],
        "patterns": [
            "*remove_version*.csv",
            "*remove_version*.tsv",
            "*remove-version*.csv",
            "*remove-version*.tsv",
            "*version_stripped*.csv",
            "*version_stripped*.tsv",
            "*stripped_version*.csv",
            "*stripped_version*.tsv",
            "*noversion*.csv",
            "*noversion*.tsv",
        ],
        "content_checker": "expression_matrix",
        "content_path_keywords": ["remove_version", "remove-version", "version_stripped", "stripped_version", "noversion"],
    },
    {
        "id": "extract_hla_read",
        "label": "extract_hla_read",
        "policy": "iobrpy",
        "stage_ids": [],
        "patterns": ["ExtractHLAread/*/*_1.fq.gz", "ExtractHLAread/*/*_2.fq.gz", "*.ExtractHLAread.done"],
    },
    {
        "id": "spechla",
        "label": "spechla",
        "policy": "iobrpy",
        "stage_ids": ["spechla"],
        "patterns": ["*HLAfinal.type.txt", "SpecHLA/*/hla.result.txt", "SpecHLA/*/hla.results.txt", "*.SpecHLA.done"],
        "content_checker": "spechla_output",
        "content_context": "hla",
        "content_path_exact_parts": ["spechla"],
        "content_path_keywords": ["hlafinal.type", "hla.result", "hla.results"],
    },
    {
        "id": "hla_typing",
        "label": "hla_typing",
        "policy": "iobrpy",
        "stage_ids": ["hla_typing"],
        "patterns": [
            "hla_result_merged.txt",
            "*hla_result*.txt",
            "*hla_result*.tsv",
            "*hla_result*.csv",
            "*hla*merged*.txt",
            "*hla*merged*.tsv",
            "*hla*merged*.csv",
            "*cohort*genotype*.tsv",
            "*cohort*genotype*.csv",
        ],
        "native_filenames": ["hla_result_merged.txt", "hla_result_merged.tsv", "hla_result_merged.csv"],
        "name_context": "hla",
        "content_checker": "hla_merged",
        "content_context": "hla",
        "content_path_keywords": ["hla_result", "merged", "cohort"],
    },
    {
        "id": "trust4",
        "label": "trust4",
        "policy": "iobrpy",
        "stage_ids": ["trust4"],
        "patterns": [
            "*TRUST*_report.tsv",
            "*TRUST*_airr.tsv",
            "*TRUST*_airr_align.tsv",
            "*TRUST*_cdr3.out",
            "*.TRUST4.done",
            "trust4_immdata.csv",
            "trust4_immune_indices.csv",
            "07-TCRBCR",
            "*-TCRBCR",
            "*_TCRBCR",
        ],
        "content_checker": "trust4_output",
        "content_context": "trust4",
    },
    {
        "id": "tme_profile",
        "label": "tme_profile",
        "policy": "iobrpy_wrapper",
        "stage_ids": ["signature_scoring", "deconvolution", "ligand_receptor"],
        "patterns": [
            "01-signatures/*",
            "02-tme/*",
            "03-LR_cal/*",
            "*deconvo_merged*.csv",
            "*deconvo_merged*.tsv",
        ],
    },
]

_STAGE_USER_COPY: Dict[str, Dict[str, str]] = {
    "raw_fastq": {
        "label": "Original sequencing data",
        "label_zh": "原始测序数据",
        "done_text": "Original FASTQ sequencing files are available",
        "done_text_zh": "已检测到原始 FASTQ 测序文件",
        "pending_text": "Original FASTQ sequencing files have not been confirmed",
        "pending_text_zh": "尚未确认原始 FASTQ 测序文件",
        "description": "Raw sequencing reads are present and can be used as the starting point for a full analysis run.",
        "description_zh": "目录中存在原始测序 reads，可作为完整分析流程的起点。",
        "next_text": "Start from the original sequencing files",
        "next_text_zh": "从原始测序文件开始",
    },
    "fastq_qc": {
        "label": "FASTQ quality control and cleaning with fastp",
        "label_zh": "使用 fastp 进行 FASTQ 质量控制与清洗",
        "done_text": "fastp-based FASTQ quality control and cleaning has been completed",
        "done_text_zh": "已完成使用 fastp 的 FASTQ 质量控制与清洗",
        "pending_text": "fastp FASTQ quality-control outputs have not been confirmed yet",
        "pending_text_zh": "尚未确认 fastp 的 FASTQ 质量控制结果",
        "description": "Reads have been checked and cleaned with fastp so they are ready for expression quantification.",
        "description_zh": "reads 已使用 fastp 完成质控与清洗，可继续进行表达定量。",
        "next_text": "Perform fastp-based FASTQ quality control",
        "next_text_zh": "使用 fastp 进行 FASTQ 质量控制",
    },
    "salmon_quant": {
        "label": "Expression quantification with Salmon",
        "label_zh": "Salmon 表达定量",
        "done_text": "Expression quantification with Salmon has been completed",
        "done_text_zh": "Salmon 表达定量已完成",
        "pending_text": "Expression quantification with Salmon has not been confirmed yet",
        "pending_text_zh": "尚未确认 Salmon 表达定量结果",
        "description": "Sequencing reads have been quantified to estimate gene or transcript abundance.",
        "description_zh": "测序 reads 已完成定量，可估计基因或转录本丰度。",
        "next_text": "Quantify expression with Salmon",
        "next_text_zh": "使用 Salmon 进行表达定量",
    },
    "salmon_merge": {
        "label": "Merged Salmon quantification results",
        "label_zh": "合并 Salmon 定量结果",
        "done_text": "Salmon quantification results have been merged across samples",
        "done_text_zh": "Salmon 定量结果已在样本间完成合并",
        "pending_text": "Merged Salmon quantification results have not been confirmed yet",
        "pending_text_zh": "尚未确认合并后的 Salmon 定量结果",
        "description": "Per-sample Salmon results have been combined into matrix-style outputs.",
        "description_zh": "各样本的 Salmon 结果已合并为矩阵形式输出。",
        "next_text": "Merge the per-sample Salmon quantification results",
        "next_text_zh": "合并各样本的 Salmon 定量结果",
    },
    "prepare_salmon": {
        "label": "TPM matrix prepared from Salmon results",
        "label_zh": "由 Salmon 结果生成 TPM 矩阵",
        "done_text": "A TPM expression matrix has been prepared from Salmon quantification",
        "done_text_zh": "已根据 Salmon 定量结果生成 TPM 表达矩阵",
        "pending_text": "A TPM expression matrix has not been prepared from Salmon results yet",
        "pending_text_zh": "尚未根据 Salmon 结果生成 TPM 表达矩阵",
        "description": "The Salmon quantification results have been converted into a TPM matrix suitable for downstream analysis.",
        "description_zh": "Salmon 定量结果已转换为适合下游分析的 TPM 矩阵。",
        "next_text": "Prepare a TPM matrix from Salmon results",
        "next_text_zh": "根据 Salmon 结果生成 TPM 矩阵",
    },
    "star_quant": {
        "label": "Expression quantification with STAR",
        "label_zh": "STAR 表达定量",
        "done_text": "Expression quantification with STAR has been completed",
        "done_text_zh": "STAR 表达定量已完成",
        "pending_text": "Expression quantification with STAR has not been confirmed yet",
        "pending_text_zh": "尚未确认 STAR 表达定量结果",
        "description": "Sequencing reads have been aligned and counted with STAR.",
        "description_zh": "测序 reads 已通过 STAR 完成比对和计数。",
        "next_text": "Quantify expression with STAR",
        "next_text_zh": "使用 STAR 进行表达定量",
    },
    "star_merge": {
        "label": "Merged gene count results",
        "label_zh": "合并基因计数矩阵",
        "done_text": "Gene count results have been merged across samples",
        "done_text_zh": "基因计数结果已在样本间完成合并",
        "pending_text": "Merged gene count results have not been confirmed yet",
        "pending_text_zh": "尚未确认合并后的基因计数矩阵",
        "description": "Per-sample count outputs have been combined into one matrix for downstream conversion.",
        "description_zh": "各样本的 count 输出已合并为一个矩阵，可继续下游转换。",
        "next_text": "Merge the per-sample gene count results",
        "next_text_zh": "合并各样本的基因计数结果",
    },
    "count2tpm": {
        "label": "TPM matrix prepared from gene counts",
        "label_zh": "由基因计数转换得到 TPM 矩阵",
        "done_text": "A TPM expression matrix has been prepared from gene counts",
        "done_text_zh": "已根据 gene count 生成 TPM 表达矩阵",
        "pending_text": "A TPM expression matrix has not been prepared from gene counts yet",
        "pending_text_zh": "尚未根据 gene count 生成 TPM 表达矩阵",
        "description": "Gene counts have been converted into a TPM matrix suitable for downstream analysis.",
        "description_zh": "gene count 已转换为适合下游分析的 TPM 矩阵。",
        "next_text": "Convert gene counts into a TPM matrix",
        "next_text_zh": "将 gene count 转换为 TPM 矩阵",
    },
    "tpm_matrix": {
        "label": "TPM expression matrix ready for downstream analysis",
        "label_zh": "可用于下游分析的 TPM 表达矩阵",
        "done_text": "A TPM expression matrix is ready for downstream analysis",
        "done_text_zh": "TPM 表达矩阵已可用于下游分析",
        "pending_text": "A TPM expression matrix has not been confirmed yet",
        "pending_text_zh": "尚未确认可用于下游分析的 TPM 表达矩阵",
        "description": "The data is at the stage where tumor microenvironment analysis can begin.",
        "description_zh": "数据已达到可以开始肿瘤微环境分析的阶段。",
        "next_text": "Use the TPM matrix for downstream tumor microenvironment analysis",
        "next_text_zh": "使用 TPM 矩阵开展下游肿瘤微环境分析",
    },
    "signature_scoring": {
        "label": "Pathway and signature scoring",
        "label_zh": "通路与特征集评分",
        "done_text": "Pathway and signature scoring has been completed",
        "done_text_zh": "通路与特征集评分已完成",
        "pending_text": "Pathway and signature scoring has not been completed yet",
        "pending_text_zh": "尚未确认通路与特征集评分结果",
        "description": "Pathway activity or biological signature scores have been calculated from the expression matrix.",
        "description_zh": "已根据表达矩阵计算通路活性或生物学特征集评分。",
        "next_text": "Calculate pathway and signature scores",
        "next_text_zh": "计算通路与特征集评分",
    },
    "deconvolution": {
        "label": "Immune and stromal composition analysis",
        "label_zh": "免疫与基质组成分析",
        "done_text": "Immune and stromal composition analysis has been completed",
        "done_text_zh": "免疫与基质组成分析已完成",
        "pending_text": "Immune and stromal composition analysis has not been completed yet",
        "pending_text_zh": "尚未确认免疫与基质组成分析结果",
        "description": "Cell composition in the tumor microenvironment has been estimated from bulk RNA-seq data.",
        "description_zh": "已根据 bulk RNA-seq 数据估计肿瘤微环境中的细胞组成。",
        "next_text": "Estimate immune and stromal composition",
        "next_text_zh": "估计免疫与基质组成",
    },
    "ligand_receptor": {
        "label": "Cell-cell communication analysis",
        "label_zh": "细胞间通讯分析",
        "done_text": "Cell-cell communication analysis has been completed",
        "done_text_zh": "细胞间通讯分析已完成",
        "pending_text": "Cell-cell communication analysis has not been completed yet",
        "pending_text_zh": "尚未确认细胞间通讯分析结果",
        "description": "Ligand-receptor activity has been evaluated to study interactions in the tumor microenvironment.",
        "description_zh": "已评估配体-受体活性，用于研究肿瘤微环境中的相互作用。",
        "next_text": "Run cell-cell communication analysis",
        "next_text_zh": "运行细胞间通讯分析",
    },
    "clustering": {
        "label": "TME subtype / clustering analysis",
        "label_zh": "TME 亚型 / 聚类分析",
        "done_text": "TME subtype or clustering analysis has been completed",
        "done_text_zh": "TME 亚型或聚类分析已完成",
        "pending_text": "TME subtype or clustering analysis has not been completed yet",
        "pending_text_zh": "尚未确认 TME 亚型或聚类分析结果",
        "description": "Samples have been grouped into tumor microenvironment subtypes or clusters.",
        "description_zh": "样本已根据肿瘤微环境特征分为不同亚型或聚类。",
        "next_text": "Run TME subtype or clustering analysis",
        "next_text_zh": "运行 TME 亚型或聚类分析",
    },
    "trust4": {
        "label": "TCR/BCR repertoire analysis",
        "label_zh": "TCR/BCR 库分析",
        "done_text": "TCR/BCR repertoire analysis has been completed",
        "done_text_zh": "TCR/BCR 库分析已完成",
        "pending_text": "TCR/BCR repertoire analysis has not been completed yet",
        "pending_text_zh": "尚未确认 TCR/BCR 库分析结果",
        "description": "Immune receptor repertoire information has been derived from the sequencing data.",
        "description_zh": "已根据测序数据解析免疫受体库信息。",
        "next_text": "Run TCR/BCR repertoire analysis",
        "next_text_zh": "运行 TCR/BCR 库分析",
    },
    "spechla": {
        "label": "Per-sample HLA typing",
        "label_zh": "单样本 HLA 分型",
        "done_text": "Per-sample HLA typing has been completed",
        "done_text_zh": "单样本 HLA 分型已完成",
        "pending_text": "Per-sample HLA typing has not been completed yet",
        "pending_text_zh": "尚未确认单样本 HLA 分型结果",
        "description": "Sample-level HLA genotypes have been inferred.",
        "description_zh": "已推断单样本层面的 HLA 基因型。",
        "next_text": "Run per-sample HLA typing",
        "next_text_zh": "运行单样本 HLA 分型",
    },
    "hla_typing": {
        "label": "Merged HLA typing results",
        "label_zh": "合并 HLA 分型结果",
        "done_text": "Merged HLA typing results are available",
        "done_text_zh": "已检测到合并后的 HLA 分型结果",
        "pending_text": "Merged HLA typing results have not been confirmed yet",
        "pending_text_zh": "尚未确认合并后的 HLA 分型结果",
        "description": "HLA typing results from multiple samples have been merged into a summary result.",
        "description_zh": "多样本的 HLA 分型结果已合并为汇总输出。",
        "next_text": "Review the merged HLA typing results",
        "next_text_zh": "查看合并后的 HLA 分型结果",
    },
}

_WORKFLOW_CHECKLIST_GROUPS: List[Dict[str, Any]] = [
    {
        "id": "raw_data",
        "label": "Original sequencing data",
        "label_zh": "原始测序数据",
        "details": "Raw FASTQ/FQ files that can be used as the starting point for a full IOBRpy analysis run.",
        "details_zh": "可作为完整 IOBRpy 分析起点的原始 FASTQ/FQ 文件。",
        "stage_ids": ["raw_fastq"],
        "completed_text": "Raw FASTQ/FQ sequencing files were detected, so the dataset can be started from the original reads if needed.",
        "completed_text_zh": "已检测到原始 FASTQ/FQ 测序文件，如有需要可从原始 reads 开始分析。",
        "pending_text": "No raw FASTQ/FQ sequencing files were detected in the scanned path.",
        "pending_text_zh": "在扫描路径下未检测到原始 FASTQ/FQ 测序文件。",
    },
    {
        "id": "quality_control",
        "label": "Read quality control and cleaning with fastp",
        "label_zh": "使用 fastp 进行 Read 质量控制与清洗",
        "details": "Use fastp to perform quality checking, trimming, and cleaning of FASTQ reads before expression quantification.",
        "details_zh": "在表达定量前使用 fastp 对 FASTQ reads 进行质量检查、剪切和清洗。",
        "stage_ids": ["fastq_qc"],
        "completed_text": "fastp outputs were detected, which suggests the FASTQ files have already been checked and cleaned.",
        "completed_text_zh": "已检测到 fastp 输出，说明 FASTQ 文件已完成检查和清洗。",
        "pending_text": "No fastp quality-control outputs were detected yet.",
        "pending_text_zh": "尚未检测到 fastp 质量控制输出。",
    },
    {
        "id": "salmon_quant_merge",
        "label": "Salmon quantification and result merging",
        "label_zh": "Salmon 定量与结果合并",
        "details": "Sample-level Salmon quantification plus cross-sample merged abundance/count tables for the whole cohort.",
        "details_zh": "样本级 Salmon 定量，以及队列层面的跨样本 abundance/count 合并表。",
        "stage_ids": ["salmon_quant", "salmon_merge"],
        "completed_text": "Salmon quantification and/or merged Salmon result tables were detected, so the Salmon branch has already been started.",
        "completed_text_zh": "已检测到 Salmon 定量和/或合并后的 Salmon 结果表，说明 Salmon 分支已经开始。",
        "partial_text": "Some Salmon-related outputs were detected, but not every Salmon step is clearly complete yet.",
        "partial_text_zh": "已检测到部分 Salmon 相关输出，但并非每个 Salmon 步骤都已明确完成。",
        "pending_text": "No Salmon quantification or merged Salmon result tables were detected.",
        "pending_text_zh": "未检测到 Salmon 定量或合并后的 Salmon 结果表。",
    },
    {
        "id": "star_quant_merge",
        "label": "STAR quantification and merged count matrix",
        "label_zh": "STAR 表达定量和合并计数矩阵",
        "details": "STAR-based alignment/counting plus a merged cross-sample gene-count matrix for the cohort.",
        "details_zh": "基于 STAR 的比对/计数，以及队列层面的跨样本基因 count 合并矩阵。",
        "stage_ids": ["star_quant", "star_merge"],
        "completed_text": "STAR quantification and/or merged gene-count outputs were detected, so the STAR branch has already been started.",
        "completed_text_zh": "已检测到 STAR 定量和/或合并后的基因 count 输出，说明 STAR 分支已经开始。",
        "partial_text": "Some STAR-related outputs were detected, but not every STAR step is clearly complete yet.",
        "partial_text_zh": "已检测到部分 STAR 相关输出，但并非每个 STAR 步骤都已明确完成。",
        "pending_text": "No STAR quantification outputs were detected. If the directory follows the Salmon route, this STAR branch is usually not needed.",
        "pending_text_zh": "未检测到 STAR 定量输出。如果该目录走的是 Salmon 路线，通常不需要这一 STAR 分支。",
    },
    {
        "id": "tpm_matrix_ready",
        "label": "TPM expression matrix for downstream analysis",
        "label_zh": "TPM 表达矩阵",
        "details": "A cross-sample TPM matrix that can be used directly for downstream tumor microenvironment analysis.",
        "details_zh": "可直接用于下游肿瘤微环境分析的跨样本 TPM 矩阵。",
        "stage_ids": ["prepare_salmon", "count2tpm", "tpm_matrix"],
        "completed_text": "A TPM-like expression matrix was detected and appears ready, or nearly ready, for downstream TME analysis.",
        "completed_text_zh": "已检测到 TPM 类表达矩阵，并且看起来已经准备好或接近可用于下游 TME 分析。",
        "pending_text": "No TPM matrix that is clearly ready for downstream analysis was confirmed yet.",
        "pending_text_zh": "尚未确认明确可用于下游分析的 TPM 矩阵。",
    },
    {
        "id": "feature_scoring",
        "label": "Feature-set scoring",
        "label_zh": "特征集评分",
        "details": "Pathway or gene-signature scoring from the TPM matrix, such as curated feature sets, ssGSEA-style scores, or related signature summaries.",
        "details_zh": "基于 TPM 矩阵进行通路或基因特征集评分，例如 curated feature sets、ssGSEA 风格评分或相关 signature 汇总。",
        "stage_ids": ["signature_scoring"],
        "completed_text": "Feature-set or pathway-signature scoring outputs were detected.",
        "completed_text_zh": "已检测到特征集或通路 signature 评分输出。",
        "pending_text": "No feature-set or pathway-signature scoring outputs were confirmed yet.",
        "pending_text_zh": "尚未确认特征集或通路 signature 评分输出。",
    },
    {
        "id": "immune_deconvolution",
        "label": "Immune deconvolution",
        "label_zh": "免疫反卷积",
        "details": "Estimate immune or stromal cell proportions and related abundance or score outputs from bulk RNA-seq, such as CIBERSORT, EPIC, quanTIseq, MCPcounter, ESTIMATE, or IPS.",
        "details_zh": "根据 bulk RNA-seq 估计免疫或基质细胞比例及相关 abundance/score 输出，例如 CIBERSORT、EPIC、quanTIseq、MCPcounter、ESTIMATE 或 IPS。",
        "stage_ids": ["deconvolution"],
        "completed_text": "Immune-deconvolution outputs were detected, including cell-fraction or immune/stromal score results.",
        "completed_text_zh": "已检测到免疫反卷积输出，包括细胞比例或免疫/基质评分结果。",
        "pending_text": "No immune-deconvolution outputs were confirmed yet.",
        "pending_text_zh": "尚未确认免疫反卷积输出。",
    },
    {
        "id": "ligand_receptor_analysis",
        "label": "Ligand-receptor analysis",
        "label_zh": "配体-受体分析",
        "details": "Infer potential cell-cell communication signals from ligand-receptor pairs based on the expression matrix.",
        "details_zh": "基于表达矩阵和配体-受体对推断潜在的细胞间通讯信号。",
        "stage_ids": ["ligand_receptor"],
        "completed_text": "Ligand-receptor analysis outputs were detected.",
        "completed_text_zh": "已检测到配体-受体分析输出。",
        "pending_text": "No ligand-receptor analysis outputs were confirmed yet.",
        "pending_text_zh": "尚未确认配体-受体分析输出。",
    },
    {
        "id": "clustering",
        "label": "TME clustering or subtype analysis",
        "label_zh": "TME 聚类/亚型分析",
        "details": "Clustering or subtype discovery based on tumor-microenvironment features derived from the cohort.",
        "details_zh": "基于队列肿瘤微环境特征进行聚类或亚型发现。",
        "stage_ids": ["clustering"],
        "completed_text": "TME clustering or subtype-analysis outputs were detected.",
        "completed_text_zh": "已检测到 TME 聚类或亚型分析输出。",
        "pending_text": "No TME clustering or subtype-analysis outputs were detected yet.",
        "pending_text_zh": "尚未检测到 TME 聚类或亚型分析输出。",
    },
    {
        "id": "hla_typing_summary",
        "label": "HLA typing",
        "label_zh": "HLA 分型",
        "details": "Sample-level or cohort-level HLA genotype inference and merged HLA summary results.",
        "details_zh": "样本级或队列级 HLA 基因型推断，以及合并后的 HLA 汇总结果。",
        "stage_ids": ["spechla", "hla_typing"],
        "external_hint_ids": ["external_hla"],
        "completed_text": "HLA typing outputs were detected.",
        "completed_text_zh": "已检测到 HLA 分型输出。",
        "external_text": "HLA-related result files or directories were detected, but they may come from another tool rather than an IOBRpy HLA workflow. IOBRpy can also perform HLA typing if you want to rerun it inside the same pipeline.",
        "external_text_zh": "已检测到 HLA 相关结果文件或目录，但它们可能来自其他工具而不是 IOBRpy 的 HLA 工作流。如果你希望在同一条流程中重跑，IOBRpy 也可以执行 HLA 分型。",
        "partial_text": "Some HLA-related outputs were detected, but not every HLA step is clearly complete yet.",
        "partial_text_zh": "已检测到部分 HLA 相关输出，但并非每个 HLA 步骤都已明确完成。",
        "pending_text": "No HLA typing outputs were confirmed in the scanned path.",
        "pending_text_zh": "在扫描路径中尚未确认 HLA 分型输出。",
    },
    {
        "id": "tcr_bcr_summary",
        "label": "TCR/BCR repertoire analysis",
        "label_zh": "TCR/BCR 库分析",
        "details": "Immune-receptor repertoire analysis for clonotypes, clonal expansion, and related TCR/BCR summary outputs.",
        "details_zh": "用于克隆型、克隆扩增及相关 TCR/BCR 汇总输出的免疫受体库分析。",
        "stage_ids": ["trust4"],
        "external_hint_ids": ["external_tcr_bcr"],
        "completed_text": "TCR/BCR repertoire-analysis outputs were detected.",
        "completed_text_zh": "已检测到 TCR/BCR 库分析输出。",
        "external_text": "TCR/BCR-related result files or directories were detected, but they appear to come from an external tool rather than IOBRpy TRUST4. IOBRpy can also perform TCR/BCR repertoire analysis if you want to rerun it inside the same pipeline.",
        "external_text_zh": "已检测到 TCR/BCR 相关结果文件或目录，但它们看起来来自外部工具而不是 IOBRpy TRUST4。如果你希望在同一条流程中重跑，IOBRpy 也可以执行 TCR/BCR 库分析。",
        "pending_text": "No TCR/BCR repertoire-analysis outputs were confirmed in the scanned path.",
        "pending_text_zh": "在扫描路径中尚未确认 TCR/BCR 库分析输出。",
    },
]

_AGENT_FALLBACK_RESULT_BUCKETS: List[Dict[str, str]] = [
    {
        "id": "iobrpy_confirmed_results",
        "label": "IOBRpy-confirmed results",
        "label_zh": "IOBRpy 已确认结果",
        "definition": "Use only map-confirmed native results here, such as content_verified_functions, content_verified_stages, and checklist sections that have direct IOBRpy evidence.",
        "definition_zh": "这一类只放 map 已确认的 IOBRpy 原生结果，例如 content_verified_functions、content_verified_stages，以及具有直接 IOBRpy 证据的 checklist 已完成部分。",
    },
    {
        "id": "agent_inferred_existing_results",
        "label": "Agent-inferred existing results",
        "label_zh": "Agent 推断的已有结果",
        "definition": "Use this bucket for results the agent finds by targeted search or file preview when map did not confirm them as native IOBRpy outputs.",
        "definition_zh": "这一类用于 agent 通过定向搜索或文件预览发现、但 map 尚未确认其为 IOBRpy 原生输出的结果。",
    },
    {
        "id": "external_tool_results",
        "label": "External-tool results",
        "label_zh": "外部工具结果",
        "definition": "Use this bucket for outputs that appear to come from non-IOBRpy tools or generic external analyses.",
        "definition_zh": "这一类用于看起来来自非 IOBRpy 工具或通用外部分析的软件结果。",
    },
]

_AGENT_FALLBACK_TERM_OVERRIDES: Dict[str, List[str]] = {
    "raw_data": ["fastq"],
    "quality_control": ["fastp"],
    "salmon_quant_merge": ["salmon", "quant.sf", "merge_salmon", "abundance", "numreads"],
    "star_quant_merge": ["star", "readspergene", "merge_star_count", "count", "featurecounts"],
    "tpm_matrix_ready": ["tpm", "tpm_matrix", "prepare_salmon", "count2tpm", "expression"],
    "feature_scoring": ["signature", "signatures", "score", "scores", "ssgsea", "zscore", "pca", "calculate_sig_score"],
    "immune_deconvolution": ["deconvolution", "cibersort", "epic", "quantiseq", "mcpcounter", "estimate", "ips", "bayesprism"],
    "ligand_receptor_analysis": ["lr", "ligand", "receptor", "lr_cal"],
    "clustering": ["cluster", "clustering", "subtype", "tme_cluster", "nmf"],
    "hla_typing_summary": ["hla", "spechla", "optitype", "arcas", "hla_typing"],
    "tcr_bcr_summary": ["tcr", "bcr", "trust4", "mixcr", "airr", "clonotype", "vdj"],
}

_GENERIC_INVESTIGATION_TERMS = {
    "analysis",
    "cleaning",
    "data",
    "expression",
    "matrix",
    "merged",
    "original",
    "quality",
    "quantification",
    "read",
    "reads",
    "result",
    "results",
    "or",
    "set",
    "subtype",
    "summary",
    "typing",
}


@lru_cache(maxsize=1)
def load_pipeline_map_rules() -> Dict[str, Any]:
    raw = (files("iobrpy.RAG_MCP") / "iobrpy_pipeline_map.json").read_text(encoding="utf-8")
    return json.loads(raw)


@lru_cache(maxsize=1)
def load_iobrpy_required_params() -> Dict[str, Any]:
    raw = (files("iobrpy.RAG_MCP") / "iobrpy_required_params.json").read_text(encoding="utf-8")
    return json.loads(raw)


@lru_cache(maxsize=1)
def load_agent_parameter_hints() -> Dict[str, Any]:
    try:
        raw = (files("iobrpy.RAG_MCP") / "iobrpy_agent_parameter_hints.json").read_text(encoding="utf-8")
    except FileNotFoundError:
        return {"version": 1, "commands": {}}
    return json.loads(raw)


def _closure_with_implies(stage_ids: Iterable[str], stage_index: Dict[str, Dict[str, Any]]) -> Set[str]:
    completed = set(stage_ids)
    changed = True
    while changed:
        changed = False
        for stage_id in list(completed):
            for implied in stage_index.get(stage_id, {}).get("implies", []):
                if implied in _NON_IMPLIED_STAGE_IDS:
                    continue
                if implied not in completed:
                    completed.add(implied)
                    changed = True
    return completed


def _ordered_stage_ids(stage_ids: Iterable[str]) -> List[str]:
    stage_set = set(stage_ids)
    return [stage_id for stage_id in _STAGE_DISPLAY_ORDER if stage_id in stage_set]


def _stage_copy(stage_id: str, stage_index: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    copy = dict(_STAGE_USER_COPY.get(stage_id, {}))
    stage = stage_index.get(stage_id, {})
    if "label" not in copy:
        copy["label"] = stage.get("label", stage_id.replace("_", " ").title())
    if "done_text" not in copy:
        copy["done_text"] = f"{copy['label']} has been completed"
    if "pending_text" not in copy:
        copy["pending_text"] = f"{copy['label']} has not been completed yet"
    if "description" not in copy:
        copy["description"] = copy["label"]
    if "next_text" not in copy:
        copy["next_text"] = copy["label"]
    if "label_zh" not in copy:
        copy["label_zh"] = copy["label"]
    if "done_text_zh" not in copy:
        copy["done_text_zh"] = copy["done_text"]
    if "pending_text_zh" not in copy:
        copy["pending_text_zh"] = copy["pending_text"]
    if "description_zh" not in copy:
        copy["description_zh"] = copy["description"]
    if "next_text_zh" not in copy:
        copy["next_text_zh"] = copy["next_text"]
    return copy


def _localized_copy_value(copy: Dict[str, str], key: str, language: str = "en") -> str:
    if language == "zh":
        return copy.get(f"{key}_zh", copy.get(key, ""))
    return copy.get(key, "")


def _is_text_like_result(relpath: str) -> bool:
    suffixes = [suffix.lower() for suffix in Path(relpath).suffixes]
    if not suffixes:
        return False
    if suffixes[-1] in {".csv", ".tsv", ".txt", ".json", ".html", ".out", ".sf"}:
        return True
    if suffixes[-1] == ".gz" and len(suffixes) >= 2 and suffixes[-2] in {".csv", ".tsv", ".txt", ".out"}:
        return True
    return False


def _read_text_preview(path: Path, *, max_chars: int = 64000) -> str:
    try:
        if path.suffix.lower() == ".gz":
            with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
                return handle.read(max_chars)
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            return handle.read(max_chars)
    except OSError:
        return ""


def _matches_marker_groups(text: str, marker_groups: Iterable[Iterable[str]]) -> bool:
    normalized = text.lower()
    return any(all(marker in normalized for marker in group) for group in marker_groups)


def _first_nonempty_lines(text: str, *, limit: int = 6) -> List[str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines[:limit]


def _split_table_line(line: str) -> List[str]:
    if "\t" in line:
        return [item.strip() for item in line.split("\t")]
    return [item.strip() for item in line.split(",")]


def _preview_table_rows(preview: str, *, limit: int = 3) -> List[List[str]]:
    return [_split_table_line(line) for line in _first_nonempty_lines(preview, limit=limit)]


def _preview_header(preview: str) -> List[str]:
    rows = _preview_table_rows(preview, limit=1)
    return rows[0] if rows else []


def _preview_feature_columns(preview: str) -> List[str]:
    header = _preview_header(preview)
    return header[1:] if len(header) >= 2 else []


def _lowered_columns(columns: Iterable[str]) -> List[str]:
    return [str(item).strip().lower() for item in columns if str(item).strip()]


def _preview_deconvolution_feature_columns(preview: str) -> List[str]:
    return _lowered_columns(_preview_feature_columns(preview))


def _looks_like_delimited_result_name(filename: str) -> bool:
    lowered = filename.lower()
    return lowered.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz"))


def _canonical_deconvolution_method_id(method_id: str) -> str:
    lowered = str(method_id).strip().lower()
    if lowered == "ips":
        return "ips"
    return lowered


def _infer_deconvolution_method_ids_from_columns(columns: Iterable[str]) -> List[str]:
    feature_columns = _lowered_columns(columns)
    feature_set = set(feature_columns)
    detected: List[str] = []
    if sum(column.endswith("_cibersort") for column in feature_columns) >= 2:
        detected.append("cibersort")
    if len(feature_set.intersection(_IPS_SUFFIX_COLUMNS)) >= 2 or sum(column.endswith("_ips") for column in feature_columns) >= 5:
        detected.append("ips")
    estimate_suffix_count = sum(column.endswith("_estimate") for column in feature_columns)
    if (
        {"stromalscore_estimate", "immunescore_estimate", "estimatescore_estimate"}.issubset(feature_set)
        or {"stromalsignature_estimate", "immunesignature_estimate", "estimatescore_estimate"}.issubset(feature_set)
        or estimate_suffix_count >= 3
    ):
        detected.append("estimate")
    if sum(column.endswith("_mcpcounter") for column in feature_columns) >= 2:
        detected.append("mcpcounter")
    if sum(column.endswith("_quantiseq") for column in feature_columns) >= 2:
        detected.append("quantiseq")
    if sum(column.endswith("_epic") for column in feature_columns) >= 2:
        detected.append("epic")
    if any(column.endswith("_bayesprism") for column in feature_columns):
        detected.append("bayesprism")
    return detected


def _collect_merged_deconvolution_method_matches(preview_texts: Dict[str, str]) -> Dict[str, List[str]]:
    matches: Dict[str, List[str]] = {}
    for relpath, preview in preview_texts.items():
        filename = Path(relpath).name.lower()
        if "deconvo_merged" not in filename or not _looks_like_delimited_result_name(filename):
            continue
        for method_id in _infer_deconvolution_method_ids_from_columns(_preview_deconvolution_feature_columns(preview)):
            matches.setdefault(method_id, [])
            if relpath not in matches[method_id]:
                matches[method_id].append(relpath)
    return matches


def _looks_like_star_count_matrix_path(relpath: str) -> bool:
    lower = relpath.lower()
    filename = Path(relpath).name.lower()
    if not _looks_like_delimited_result_name(filename):
        return False
    if filename in {"biosample_count.txt", "biosample_counts.txt"}:
        return False
    if any(token in lower for token in ("star.count", "star_count", "readspergene", "featurecounts", "featurecount")):
        return True
    if not (_path_has_part_token(relpath, "star") or _path_has_exact_part(relpath, "02-star", "03-star", "04-star", "star")):
        return False
    return any(token in filename for token in ("count_matrix", "counts_matrix", "gene_count", "read_count"))


def _looks_like_hla_per_sample_result_path(relpath: str) -> bool:
    filename = Path(relpath).name.lower()
    if not _looks_like_delimited_result_name(filename):
        return False
    if not filename.endswith(("_result.tsv", "_result.csv", "_result.txt")):
        return False
    return _path_has_part_token(relpath, "hla", "spechla", "hlatyping", "hla_typing")


def _should_force_preview_result(relpath: str) -> bool:
    lower = relpath.lower()
    parts = set(_path_parts_lower(relpath))
    filename = Path(relpath).name.lower()
    native_filenames = {
        str(name).strip().lower()
        for spec in _FUNCTION_DETECTION_SPECS
        for name in spec.get("native_filenames", [])
        if str(name).strip()
    }
    if filename in native_filenames:
        return True
    if _looks_like_star_count_matrix_path(relpath):
        return True
    if _looks_like_hla_per_sample_result_path(relpath):
        return True
    if _looks_like_delimited_result_name(filename) and any(
        token in lower
        for token in (
            "cibersort",
            "estimate",
            "mcpcounter",
            "quantiseq",
            "epic",
            "ips",
            "bayesprism",
            "deconvo",
            "tme_profile",
            "tme-profile",
            "lr_cal",
            "lr_scores",
            "hla_result",
            "hla.results",
            "hla.result",
            "trust4",
            "trust_",
            "tcrbcr",
            "airr",
            "cdr3",
            "immune_indices",
            "immdata",
        )
    ):
        return True
    if "tme_profile" in lower and _looks_like_delimited_result_name(filename):
        return True
    in_deconvolution_dir = bool(parts.intersection({"02-tme", "05-tme"}))
    if in_deconvolution_dir:
        if _looks_like_delimited_result_name(filename):
            return True
    if "deconvo_merged" in lower and _looks_like_delimited_result_name(Path(relpath).name):
        return True
    return False


def _looks_like_confident_tpm_path_match(relpath: str) -> bool:
    filename = Path(relpath).name.lower()
    return (
        "tpm" in filename
        or filename in {"prepare_salmon.csv", "prepare_salmon.tsv", "count2tpm.csv", "count2tpm.tsv"}
    )


@lru_cache(maxsize=1)
def _load_pickled_resource(filename: str) -> Any:
    try:
        with files("iobrpy.resources").joinpath(filename).open("rb") as handle:
            return pickle.load(handle)
    except Exception:
        repo_resource_path = Path(__file__).resolve().parents[3] / "src" / "iobrpy" / "resources" / filename
        try:
            with repo_resource_path.open("rb") as handle:
                return pickle.load(handle)
        except Exception:
            return None


def _looks_like_gene_list(value: Any) -> bool:
    if isinstance(value, (list, tuple, set, frozenset)):
        return True
    if isinstance(value, (str, bytes, dict)):
        return False
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        try:
            return isinstance(tolist(), list)
        except Exception:
            return False
    return False


def _collect_signature_name_variants(obj: Any, *, group_name: str = "", names: Optional[Set[str]] = None) -> Set[str]:
    collected = names if names is not None else set()
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_str = str(key).strip()
            if isinstance(value, dict):
                _collect_signature_name_variants(value, group_name=key_str or group_name, names=collected)
                continue
            if not key_str or not _looks_like_gene_list(value):
                continue
            collected.add(key_str.lower())
            if group_name:
                collected.add(f"{key_str}__{group_name}".lower())
    return collected


@lru_cache(maxsize=1)
def _known_signature_score_columns() -> Set[str]:
    raw = _load_pickled_resource("calculate_data.pkl")
    names = _collect_signature_name_variants(raw)
    if "tmescorea_cir" in names and "tmescoreb_cir" in names:
        names.add("tmescore_cir")
    if "tmescorea_plus" in names and "tmescoreb_plus" in names:
        names.add("tmescore_plus")
    derived = set(names)
    for suffix in ("_pca", "_zscore", "_ssgsea"):
        derived.update(f"{name}{suffix}" for name in names)
    return derived


def _coerce_string_list(value: Any) -> List[str]:
    if isinstance(value, str):
        return [value] if value else []
    if isinstance(value, (list, tuple, set, frozenset)):
        return [str(item) for item in value if str(item)]
    return []


def _lr_output_columns_for_network(network: Dict[str, Any], group_lrpairs: Iterable[Dict[str, Any]]) -> List[str]:
    ligands = _coerce_string_list(network.get("ligands", []))
    receptors = _coerce_string_list(network.get("receptors", []))
    raw_pairs: List[str] = []
    seen_pairs: Set[str] = set()
    for ligand, receptor in zip(ligands, receptors):
        pair = f"{ligand}_{receptor}"
        if pair not in seen_pairs:
            seen_pairs.add(pair)
            raw_pairs.append(pair)

    cols = list(raw_pairs)
    for group in group_lrpairs:
        if not isinstance(group, dict):
            continue
        main = str(group.get("main", "")).strip()
        if not main:
            continue
        combo = str(group.get("combo_name", "")).strip() or main
        cols = [combo if column == main else column for column in cols]
        remove_items = set(_coerce_string_list(group.get("involved_pairs", [])))
        if remove_items:
            cols = [column for column in cols if column not in remove_items]

    deduped: List[str] = []
    seen: Set[str] = set()
    for column in cols:
        if column not in seen:
            seen.add(column)
            deduped.append(column)
    return deduped


@lru_cache(maxsize=1)
def _known_ligand_receptor_columns() -> Set[str]:
    raw = _load_pickled_resource("lr_data.pkl")
    if not isinstance(raw, dict):
        return set()
    networks = raw.get("intercell_networks", {})
    group_lrpairs = raw.get("group_lrpairs", [])
    if not isinstance(networks, dict):
        return set()

    columns: Set[str] = set()
    for network in networks.values():
        if not isinstance(network, dict):
            continue
        columns.update(column.lower() for column in _lr_output_columns_for_network(network, group_lrpairs))
    return columns


_ENSEMBL_GENE_ID_PATTERN = re.compile(r"^ENS(?:[A-Z]{0,6})?(?:G|T|P|R)?\d+(?:\.\d+)?$", flags=re.IGNORECASE)
_ENTREZ_GENE_ID_PATTERN = re.compile(r"^\d{2,9}$")
_GENE_SYMBOL_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{1,31}$")
_SAMPLE_LIKE_IDENTIFIER_PATTERN = re.compile(
    r"^(?:sample(?:[-_].*)?|patient(?:[-_].*)?|tcga[-_].*|gsm\d+$|srr\d+$|err\d+$|drr\d+$|srx\d+$|erx\d+$|drx\d+$|[sp]\d{1,4}$)",
    flags=re.IGNORECASE,
)


def _looks_like_sample_identifier(value: str) -> bool:
    token = str(value).strip()
    return bool(token and _SAMPLE_LIKE_IDENTIFIER_PATTERN.match(token))


def _looks_like_gene_identifier(value: str) -> bool:
    token = str(value).strip()
    if not token:
        return False
    if _looks_like_sample_identifier(token):
        return False
    if token.count("-") >= 3 or token.count("_") >= 3:
        return False
    if _ENSEMBL_GENE_ID_PATTERN.match(token) or _ENTREZ_GENE_ID_PATTERN.match(token):
        return True
    if not _GENE_SYMBOL_PATTERN.match(token):
        return False
    letters = sum(char.isalpha() for char in token)
    return letters >= 2


def _looks_like_fastp_report(preview: str) -> bool:
    lower = preview.lower()
    return (
        ("before_filtering" in lower and "after_filtering" in lower)
        or ("fastp report" in lower)
        or ("duplication" in lower and "filtering_result" in lower)
    )


def _looks_like_salmon_quant(preview: str) -> bool:
    lines = _first_nonempty_lines(preview, limit=2)
    if not lines:
        return False
    header = [item.lower() for item in _split_table_line(lines[0])]
    return header[:5] == ["name", "length", "effectivelength", "tpm", "numreads"]


def _looks_like_star_reads_per_gene(preview: str) -> bool:
    lines = _first_nonempty_lines(preview, limit=6)
    if len(lines) < 4:
        return False
    return lines[0].startswith("N_unmapped") and any(line.startswith("N_multimapping") for line in lines[:4])


def _looks_like_featurecounts_matrix(preview: str) -> bool:
    lines = _first_nonempty_lines(preview, limit=3)
    if not lines:
        return False
    lower = lines[0].lower()
    return lower.startswith("geneid") and "strand" in lower and "length" in lower


def _looks_like_expression_matrix(preview: str) -> bool:
    rows = _preview_table_rows(preview, limit=3)
    if len(rows) < 2:
        return False
    header = rows[0]
    row = rows[1]
    if len(header) < 2 or len(row) < 2:
        return False
    first = header[0].lower().replace(" ", "").replace("_", "")
    if first not in {"gene", "genes", "geneid", "genesymbol", "symbol", "name", "feature", "ensembl", "id", "genename"}:
        return False
    numeric_like = 0
    for item in row[1: min(len(row), 5)]:
        try:
            float(item)
            numeric_like += 1
        except ValueError:
            pass
    required_numeric = 2 if len(row) >= 3 else 1
    return numeric_like >= required_numeric


def _looks_like_tpm_like_expression_matrix(preview: str) -> bool:
    rows = _preview_table_rows(preview, limit=4)
    if len(rows) < 2 or not _looks_like_expression_matrix(preview):
        return False
    data_rows = [row for row in rows[1:] if len(row) >= 2]
    if not data_rows:
        return False
    if sum(_looks_like_sample_identifier(row[0]) for row in data_rows[:3]) >= 2:
        return False
    gene_like_rows = sum(_looks_like_gene_identifier(row[0]) for row in data_rows[:3])
    required_gene_like_rows = 2 if len(data_rows) >= 2 else 1
    return gene_like_rows >= required_gene_like_rows


def _looks_like_count_matrix(preview: str) -> bool:
    rows = _preview_table_rows(preview, limit=3)
    if len(rows) < 2:
        return False
    header = rows[0]
    row = rows[1]
    if len(header) < 2 or len(row) < 2:
        return False
    first = _normalize_token(header[0])
    if first not in {"gene", "genes", "geneid", "genesymbol", "symbol", "name", "feature", "ensembl", "id", "genename"}:
        return False
    numeric_values: List[float] = []
    for item in row[1: min(len(row), 6)]:
        try:
            numeric_values.append(float(item))
        except ValueError:
            continue
    required_numeric = 2 if len(row) >= 3 else 1
    if len(numeric_values) < required_numeric:
        return False
    integer_like = sum(abs(value - round(value)) < 1e-9 for value in numeric_values)
    return integer_like >= required_numeric


def _looks_like_signature_scoring(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    known_columns = _known_signature_score_columns()
    recognized = [column for column in feature_columns if column in known_columns]
    return bool(recognized and (len(recognized) >= 2 or len(feature_columns) == 1))


def _looks_like_cibersort_output(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    if sum(column.endswith("_cibersort") for column in feature_columns) >= 2:
        return True
    lower = preview.lower()
    return (
        _matches_marker_groups(preview, [["b cells naive", "t cells cd8"], ["p-value", "rmse"]])
        or "_cibersort" in lower
    )


def _looks_like_ips_output(preview: str) -> bool:
    feature_columns = set(_lowered_columns(_preview_feature_columns(preview)))
    if {"az_ips", "cp_ips", "ec_ips", "mhc_ips", "sc_ips"}.issubset(feature_columns):
        return True
    if sum(column.endswith("_ips") for column in feature_columns) >= 5:
        return True
    lower = preview.lower()
    lines = _first_nonempty_lines(preview, limit=2)
    if not lines:
        return False
    header_tokens = {item.strip().lower() for item in _split_table_line(lines[0])}
    return (
        {"az", "cp", "ec", "mhc", "sc"}.issubset(header_tokens)
        or any("immunophenoscore" in item for item in header_tokens)
        or "_ips" in lower
    )


def _looks_like_estimate_output(preview: str) -> bool:
    feature_columns = _preview_deconvolution_feature_columns(preview)
    return "estimate" in _infer_deconvolution_method_ids_from_columns(feature_columns)


def _looks_like_mcpcounter_output(preview: str) -> bool:
    feature_columns = _preview_deconvolution_feature_columns(preview)
    return "mcpcounter" in _infer_deconvolution_method_ids_from_columns(feature_columns)


def _looks_like_quantiseq_output(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    if sum(column.endswith("_quantiseq") for column in feature_columns) >= 2:
        return True
    lower = preview.lower()
    return (
        _matches_marker_groups(preview, [["m1 macrophages", "m2 macrophages", "tregs"], ["b cells", "nk cells"]])
        or "_quantiseq" in lower
    )


def _looks_like_epic_output(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    if sum(column.endswith("_epic") for column in feature_columns) >= 2:
        return True
    lower = preview.lower()
    return (
        _matches_marker_groups(preview, [["bcells", "cd4_tcells", "cd8_tcells"], ["endothelial", "macrophages", "tumor"]])
        or "_epic" in lower
    )


def _looks_like_bayesprism_output(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    if any(column.endswith("_bayesprism") for column in feature_columns):
        return True
    lower = preview.lower()
    return (
        _matches_marker_groups(preview, [["bayesprism"], ["cell_type_fraction", "theta"]])
        or "_bayesprism" in lower
    )


def _looks_like_clustering_output(preview: str) -> bool:
    rows = _preview_table_rows(preview, limit=3)
    if not rows:
        return False
    header_tokens = set(_lowered_columns(rows[0]))
    if "cluster" in header_tokens and header_tokens.intersection({"id", "sample", "samples", "name"}):
        if len(rows) >= 2 and len(rows[1]) >= 2 and re.search(r"\btme\d+\b", rows[1][1], flags=re.IGNORECASE):
            return True
        return True
    lower = ",".join(rows[0]).lower()
    return ("cluster" in lower or "subtype" in lower or "nmf" in lower) and "sample" in lower


def _looks_like_ligand_receptor_output(preview: str) -> bool:
    feature_columns = _lowered_columns(_preview_feature_columns(preview))
    known_columns = _known_ligand_receptor_columns()
    recognized = [column for column in feature_columns if column in known_columns]
    return bool(recognized and (len(recognized) >= 2 or len(feature_columns) == 1))


_HLA_ALLELE_PATTERN = re.compile(
    r"\b(?:HLA[-_ ]?)?(?:A|B|C|E|F|G|DPA1|DPB1|DQA1|DQB1|DRA|DRB1|DRB3|DRB4|DRB5)\*?\d{1,3}(?::\d{1,3}){1,3}[A-Z]?\b",
    flags=re.IGNORECASE,
)
_HLA_SAMPLE_HEADER_TOKENS = {
    "sample",
    "samples",
    "sampleid",
    "sample_id",
    "id",
    "patient",
    "patientid",
    "patient_id",
}
_HLA_CLASS_I_HEADER_GROUPS = [
    {"a1", "hlaa1"},
    {"a2", "hlaa2"},
    {"b1", "hlab1"},
    {"b2", "hlab2"},
]
_HLA_LOCUS_HEADER_TOKENS = {
    "a",
    "b",
    "c",
    "a1",
    "a2",
    "b1",
    "b2",
    "c1",
    "c2",
    "hlaa",
    "hlab",
    "hlac",
    "hlaa1",
    "hlaa2",
    "hlab1",
    "hlab2",
    "hlac1",
    "hlac2",
    "drb1",
    "drb3",
    "drb4",
    "drb5",
    "dqa1",
    "dqb1",
    "dpa1",
    "dpb1",
    "classi",
    "classii",
}
_TCR_BCR_HEADER_TOKEN_GROUPS = {
    "airr_calls": [{"vcall", "vgene"}, {"jcall", "jgene"}, {"junctionaa", "cdr3aa", "cdr3"}],
    "clonotype_table": [{"rawclonotypeid", "clonotypeid", "cloneid", "clonotype"}, {"cdr3", "cdr3aa"}, {"productive", "count", "frequency", "fraction"}],
    "contig_annotation": [{"barcode", "cellid", "cellbarcode"}, {"chain", "locus"}, {"cdr3", "cdr3aa", "junctionaa"}],
}
_TCR_BCR_CHAIN_TOKENS = {
    "tra",
    "trb",
    "trg",
    "trd",
    "igh",
    "igk",
    "igl",
    "tcra",
    "tcrb",
    "tcrg",
    "tcrd",
    "bcr",
}
_EXTERNAL_HINT_CONFIDENCE_RANK = {"low": 1, "medium": 2, "high": 3}


def _has_tcr_bcr_repertoire_context(relpath: str) -> bool:
    lower = relpath.lower()
    if _path_has_part_token(
        relpath,
        "mixcr",
        "trust4",
        "tcr",
        "bcr",
        "airr",
        "vdj",
        "clonotype",
        "clonotypes",
        "cdr3",
        "immunarch",
        "immunearch",
    ):
        return True
    return any(
        token in lower
        for token in (
            ".clones_",
            ".clonotypes",
            "contig_annotation",
            "consensus_annotation",
            "rearrangement",
            "immune_indices",
            "immdata",
        )
    )


def _has_explicit_hla_result_context(relpath: str) -> bool:
    lower = relpath.lower()
    return (
        _path_has_part_token(relpath, "hla", "spechla", "extracthlaread", "hla_typing")
        or "hla_result" in lower
        or "hlafinal.type" in lower
        or "hla.result" in lower
        or "hla.results" in lower
    )


def _has_non_hla_downstream_context(relpath: str) -> bool:
    lower = relpath.lower()
    return any(
        token in lower
        for token in (
            "01-signatures",
            "04-signatures",
            "signature",
            "calculate_sig_score",
            "02-tme",
            "05-tme",
            "deconvo",
            "cibersort",
            "estimate",
            "mcpcounter",
            "quantiseq",
            "epic",
            "ips",
            "bayesprism",
            "03-lr_cal",
            "06-lr_cal",
            "lr_cal",
            "lr_scores",
            "tme_cluster",
            "cluster",
            "03-tpm",
            "tpm_matrix",
            "prepare_salmon",
            "count2tpm",
        )
    )


def _preview_header_tokens(preview: str) -> List[str]:
    lines = _first_nonempty_lines(preview, limit=1)
    if not lines:
        return []
    return [_normalize_token(item) for item in _split_table_line(lines[0]) if item.strip()]


def _contains_token_group(header_tokens: Set[str], groups: Iterable[Set[str]]) -> bool:
    for options in groups:
        if not header_tokens.intersection(options):
            return False
    return True


def _external_hla_hint_profile(
    relpath: str,
    preview: str,
    *,
    include_path_signals: bool = True,
) -> Optional[Dict[str, Any]]:
    if _has_tcr_bcr_repertoire_context(relpath):
        return None
    if _has_non_hla_downstream_context(relpath) and not _has_explicit_hla_result_context(relpath):
        return None

    header_tokens = set(_preview_header_tokens(preview))
    lower_preview = preview.lower()
    lower_path = relpath.lower()
    matched_signal_ids: List[str] = []
    matched_signal_labels: List[str] = []
    inference_basis: List[str] = []
    likely_tool_families: List[str] = []

    def add_signal(signal_id: str, label: str, basis: str) -> None:
        matched_signal_ids.append(signal_id)
        matched_signal_labels.append(label)
        inference_basis.append(basis)

    has_sample_header = bool(header_tokens.intersection(_HLA_SAMPLE_HEADER_TOKENS))
    has_class_i_pairs = _contains_token_group(header_tokens, _HLA_CLASS_I_HEADER_GROUPS)
    if has_sample_header and has_class_i_pairs:
        add_signal(
            "merged_hla_header_signature",
            "Header looks like a merged per-sample HLA typing table.",
            "header_signature",
        )

    locus_header_tokens = sorted(header_tokens.intersection(_HLA_LOCUS_HEADER_TOKENS))
    if len(locus_header_tokens) >= 3:
        add_signal(
            "hla_locus_header_signature",
            "Header contains multiple HLA locus columns.",
            "header_signature",
        )

    allele_matches = _dedupe_terms(_HLA_ALLELE_PATTERN.findall(preview), limit=6)
    if len(allele_matches) >= 2:
        add_signal(
            "hla_allele_pattern",
            "Content contains HLA allele-like strings such as A*02:01.",
            "allele_pattern",
        )

    if "class ii" in lower_preview or "classi" in lower_preview or "classii" in lower_preview:
        add_signal(
            "hla_class_annotation",
            "Content mentions HLA class I/II style annotations.",
            "content_signature",
        )

    if include_path_signals:
        if "optitype" in lower_path:
            likely_tool_families.append("OptiType-like HLA typing output")
            inference_basis.append("filename_pattern")
        elif "arcas" in lower_path:
            likely_tool_families.append("arcasHLA-like HLA typing output")
            inference_basis.append("filename_pattern")
        elif "xhla" in lower_path:
            likely_tool_families.append("xHLA-like HLA typing output")
            inference_basis.append("filename_pattern")
        elif "hlahd" in lower_path:
            likely_tool_families.append("HLA-HD-like HLA typing output")
            inference_basis.append("filename_pattern")
        elif any(token in lower_path for token in ("hla", "genotype")):
            likely_tool_families.append("generic HLA typing output")
            inference_basis.append("filename_pattern")

    if not matched_signal_ids and not likely_tool_families:
        return None

    if has_sample_header and has_class_i_pairs and allele_matches:
        confidence = "high"
        result_kind_id = "merged_sample_level_hla_table"
        result_kind_label = "merged sample-level HLA typing table"
    elif allele_matches or len(locus_header_tokens) >= 3:
        confidence = "medium"
        result_kind_id = "hla_typing_result_table"
        result_kind_label = "HLA typing result table"
    else:
        confidence = "low"
        result_kind_id = "hla_related_output"
        result_kind_label = "HLA-related output"

    if not likely_tool_families:
        if result_kind_id == "merged_sample_level_hla_table":
            likely_tool_families.append("merged HLA typing summary from an external workflow")
        else:
            likely_tool_families.append("external HLA typing family output")

    return {
        "analysis_family_id": "external_hla_family",
        "analysis_family_label": "External HLA typing family",
        "confidence": confidence,
        "inference_basis": _dedupe_terms(inference_basis, limit=6),
        "matched_signal_ids": _dedupe_terms(matched_signal_ids, limit=8),
        "matched_signal_labels": _dedupe_terms(matched_signal_labels, limit=8),
        "result_kind_id": result_kind_id,
        "result_kind_label": result_kind_label,
        "likely_tool_families": _dedupe_terms(likely_tool_families, limit=6),
        "header_tokens": locus_header_tokens[:8],
        "allele_examples": allele_matches[:4],
        "family_summary": (
            "This looks like an external HLA typing result family. It may be reusable for interpretation, but it does not prove that native iobrpy hla_typing or spechla was run."
        ),
        "preview_excerpt": _preview_excerpt(preview) if preview else "",
        "path": relpath,
    }


def _external_tcr_bcr_hint_profile(
    relpath: str,
    preview: str,
    *,
    include_path_signals: bool = True,
) -> Optional[Dict[str, Any]]:
    header_tokens = set(_preview_header_tokens(preview))
    lower_preview = preview.lower()
    lower_path = relpath.lower()
    matched_signal_ids: List[str] = []
    matched_signal_labels: List[str] = []
    inference_basis: List[str] = []
    likely_tool_families: List[str] = []

    def add_signal(signal_id: str, label: str, basis: str) -> None:
        matched_signal_ids.append(signal_id)
        matched_signal_labels.append(label)
        inference_basis.append(basis)

    if _contains_token_group(header_tokens, _TCR_BCR_HEADER_TOKEN_GROUPS["airr_calls"]):
        add_signal(
            "airr_style_rearrangement_header",
            "Header looks like an AIRR-style rearrangement table.",
            "header_signature",
        )

    if _contains_token_group(header_tokens, _TCR_BCR_HEADER_TOKEN_GROUPS["clonotype_table"]):
        add_signal(
            "clonotype_table_header",
            "Header looks like a clonotype or clone abundance table.",
            "header_signature",
        )

    if _contains_token_group(header_tokens, _TCR_BCR_HEADER_TOKEN_GROUPS["contig_annotation"]):
        add_signal(
            "contig_annotation_header",
            "Header looks like a contig or cell-level VDJ annotation table.",
            "header_signature",
        )

    if any(metric in lower_preview for metric in ("shannon", "simpson", "clonality", "evenness")):
        add_signal(
            "repertoire_diversity_metrics",
            "Content contains repertoire diversity metrics such as Shannon or Simpson.",
            "content_signature",
        )

    if any(token in lower_preview for token in ("cdr3", "cdr3aa", "junction_aa")) and any(
        token in lower_preview for token in ("v_gene", "j_gene", "v_call", "j_call", "productive")
    ):
        add_signal(
            "cdr3_vj_signature",
            "Content contains CDR3 plus V/J gene style repertoire annotations.",
            "content_signature",
        )

    chain_examples = _dedupe_terms(
        [token for token in _TCR_BCR_CHAIN_TOKENS if token in lower_preview or token in header_tokens],
        limit=6,
    )
    if chain_examples:
        add_signal(
            "immune_receptor_chain_tokens",
            "Content contains TCR/BCR chain or locus tokens.",
            "content_signature",
        )

    if include_path_signals:
        if "mixcr" in lower_path:
            likely_tool_families.append("MiXCR-like clonotype output")
            inference_basis.append("filename_pattern")
        elif "immunarch" in lower_path or "immunarch" in lower_preview:
            likely_tool_families.append("immunarch-style repertoire summary")
            inference_basis.append("filename_pattern")
        elif "airr" in lower_path:
            likely_tool_families.append("AIRR-format rearrangement output")
            inference_basis.append("filename_pattern")
        elif any(token in lower_path for token in ("vdj", "contig", "clonotype", "cdr3")):
            likely_tool_families.append("generic repertoire-analysis output")
            inference_basis.append("filename_pattern")

    if {"count", "frequency"}.issubset(header_tokens) and header_tokens.intersection({"cdr3aa", "cdr3nt", "cdr3"}) and header_tokens.intersection({"v", "d", "j", "c"}):
        add_signal(
            "trust4_report_header",
            "Header looks like a TRUST4 repertoire report table.",
            "header_signature",
        )

    if not matched_signal_ids and not likely_tool_families:
        return None

    if "airr_style_rearrangement_header" in matched_signal_ids or "clonotype_table_header" in matched_signal_ids or "trust4_report_header" in matched_signal_ids:
        confidence = "high"
    elif "cdr3_vj_signature" in matched_signal_ids or "repertoire_diversity_metrics" in matched_signal_ids:
        confidence = "medium"
    else:
        confidence = "low"

    if "airr_style_rearrangement_header" in matched_signal_ids:
        result_kind_id = "airr_rearrangement_table"
        result_kind_label = "AIRR-style rearrangement table"
    elif "trust4_report_header" in matched_signal_ids:
        result_kind_id = "trust4_repertoire_report"
        result_kind_label = "TRUST4 repertoire report table"
    elif "clonotype_table_header" in matched_signal_ids:
        result_kind_id = "clonotype_abundance_table"
        result_kind_label = "clonotype abundance table"
    elif "contig_annotation_header" in matched_signal_ids:
        result_kind_id = "vdj_contig_annotation_table"
        result_kind_label = "VDJ contig annotation table"
    elif "repertoire_diversity_metrics" in matched_signal_ids:
        result_kind_id = "repertoire_diversity_summary"
        result_kind_label = "repertoire diversity summary"
    else:
        result_kind_id = "tcr_bcr_related_output"
        result_kind_label = "TCR/BCR-related output"

    if not likely_tool_families:
        likely_tool_families.append("external TCR/BCR repertoire family output")

    return {
        "analysis_family_id": "external_tcr_bcr_family",
        "analysis_family_label": "External TCR/BCR repertoire family",
        "confidence": confidence,
        "inference_basis": _dedupe_terms(inference_basis, limit=6),
        "matched_signal_ids": _dedupe_terms(matched_signal_ids, limit=8),
        "matched_signal_labels": _dedupe_terms(matched_signal_labels, limit=8),
        "result_kind_id": result_kind_id,
        "result_kind_label": result_kind_label,
        "likely_tool_families": _dedupe_terms(likely_tool_families, limit=6),
        "header_tokens": sorted(header_tokens.intersection({"barcode", "chain", "cdr3", "cdr3aa", "junctionaa", "vcall", "vgene", "jcall", "jgene", "productive", "rawclonotypeid", "clonotypeid"}))[:8],
        "chain_examples": chain_examples[:4],
        "family_summary": (
            "This looks like an external TCR/BCR repertoire-analysis result family. It may be reusable for interpretation, but it does not prove that native iobrpy trust4 was run."
        ),
        "preview_excerpt": _preview_excerpt(preview) if preview else "",
        "path": relpath,
    }


def _external_hint_profile(
    relpath: str,
    preview: str,
    hint_id: str,
    *,
    include_path_signals: bool = True,
    respect_native_context: bool = True,
) -> Optional[Dict[str, Any]]:
    if hint_id == "external_hla":
        if respect_native_context and _has_iobrpy_context(relpath, "hla"):
            return None
        return _external_hla_hint_profile(relpath, preview, include_path_signals=include_path_signals)
    if hint_id == "external_tcr_bcr":
        if respect_native_context and _has_iobrpy_context(relpath, "trust4"):
            return None
        return _external_tcr_bcr_hint_profile(relpath, preview, include_path_signals=include_path_signals)
    return None


def _looks_like_hla_merged(preview: str) -> bool:
    header_tokens = set(_preview_header_tokens(preview))
    has_sample_header = bool(header_tokens.intersection(_HLA_SAMPLE_HEADER_TOKENS))
    has_locus_header = bool(header_tokens.intersection(_HLA_LOCUS_HEADER_TOKENS))
    allele_matches = _dedupe_terms(_HLA_ALLELE_PATTERN.findall(preview), limit=6)
    if has_sample_header and has_locus_header and allele_matches:
        return True
    if has_sample_header and allele_matches and any(token.startswith("hla") for token in header_tokens):
        return True
    profile = _external_hla_hint_profile("external_hla_preview.tsv", preview, include_path_signals=False)
    return bool(
        profile
        and profile.get("result_kind_id") in {"merged_sample_level_hla_table", "hla_typing_result_table"}
        and profile.get("confidence") in {"medium", "high"}
    )


def _looks_like_spechla_result(preview: str) -> bool:
    lower = preview.lower()
    return "hla-a" in lower or "hlafinal" in lower or ("a*" in lower and "b*" in lower and "c*" in lower)


def _looks_like_trust4_output(preview: str) -> bool:
    profile = _external_tcr_bcr_hint_profile("external_tcr_bcr_preview.tsv", preview, include_path_signals=False)
    return bool(profile and profile.get("confidence") in {"medium", "high"})


_FUNCTION_CONTENT_CHECKERS: Dict[str, Callable[[str], bool]] = {
    "fastp_report": _looks_like_fastp_report,
    "salmon_quant": _looks_like_salmon_quant,
    "star_reads_per_gene": _looks_like_star_reads_per_gene,
    "featurecounts_matrix": _looks_like_featurecounts_matrix,
    "count_matrix": _looks_like_count_matrix,
    "expression_matrix": _looks_like_expression_matrix,
    "signature_scoring": _looks_like_signature_scoring,
    "cibersort_output": _looks_like_cibersort_output,
    "ips_output": _looks_like_ips_output,
    "estimate_output": _looks_like_estimate_output,
    "mcpcounter_output": _looks_like_mcpcounter_output,
    "quantiseq_output": _looks_like_quantiseq_output,
    "epic_output": _looks_like_epic_output,
    "bayesprism_output": _looks_like_bayesprism_output,
    "ligand_receptor_output": _looks_like_ligand_receptor_output,
    "clustering_output": _looks_like_clustering_output,
    "nmf_clusters": _looks_like_clustering_output,
    "hla_merged": _looks_like_hla_merged,
    "spechla_output": _looks_like_spechla_result,
    "trust4_output": _looks_like_trust4_output,
}


def _has_iobrpy_context(relpath: str, family: str) -> bool:
    lower = relpath.lower()
    parts = set(_path_parts_lower(relpath))
    if family == "tme_profile":
        return (
            bool(parts.intersection({"01-signatures", "02-tme", "03-lr_cal", "04-signatures", "05-tme", "06-lr_cal"}))
            or _path_has_part_tokens(relpath, "tme", "profile")
            or "tme_profile_output" in lower
            or _path_has_part_token(relpath, "runall")
        )
    if family == "hla":
        return (
            _path_has_part_token(relpath, "hla", "spechla", "extracthlaread", "hla_typing")
            or "iobrpy" in lower
            or "hla_result_merged" in lower
            or "hlafinal.type" in lower
        )
    if family == "trust4":
        if any(token in lower for token in ("mixcr", "immunarch", "immunearch")):
            return False
        return any(
            token in lower
            for token in (
                "07-tcrbcr",
                "tcrbcr",
                "trust4",
                "trust_",
                "cdr3.out",
            )
        ) or bool(parts.intersection({"07-tcrbcr", "tcrbcr", "trust4"}))
    if family == "clustering":
        return any(
            token in lower
            for token in (
                "tme_cluster",
                "cluster_assignments",
                "top_features_per_cluster",
                "nmf",
                "iobrpy",
            )
        ) or bool(parts.intersection({"tme_cluster", "nmf", "iobrpy"}))
    return False


def _dedupe_matches(matches: Iterable[str], *, limit: int = 5) -> List[str]:
    return sorted(dict.fromkeys(matches))[:limit]


def _function_status_label(status: str) -> str:
    return {
        "not_detected": "Not detected",
        "reusable_result": "Reusable result detected",
        "likely_iobrpy_result": "Likely iobrpy result",
        "confirmed_by_content": "Confirmed by content",
    }.get(status, status.replace("_", " ").title())


def _function_status_label_zh(status: str) -> str:
    return {
        "not_detected": "未检测到",
        "reusable_result": "检测到可复用结果",
        "likely_iobrpy_result": "疑似已有 IOBRpy 结果",
        "confirmed_by_content": "已通过内容确认",
    }.get(status, status.replace("_", " ").title())


def _merge_function_status(*statuses: str) -> str:
    ranked = [status for status in statuses if status]
    if not ranked:
        return "not_detected"
    return max(ranked, key=lambda item: _FUNCTION_STATUS_ORDER.get(item, -1))


def _default_function_note(spec: Dict[str, Any]) -> str:
    if spec.get("policy") == "reusable":
        return "This command is treated as reusable when compatible upstream quantification or TPM results are present, even if they were produced outside iobrpy."
    return "This command is only counted as an iobrpy result when the filename, directory context, or file content looks consistent with iobrpy output."


def _function_stage_matches(detections: Iterable[Dict[str, Any]]) -> Dict[str, List[str]]:
    stage_matches: Dict[str, List[str]] = {}
    for detection in detections:
        if detection.get("status") == "not_detected":
            continue
        evidence = detection.get("evidence", [])
        for stage_id in detection.get("stage_ids", []):
            stage_matches.setdefault(stage_id, [])
            for match in evidence:
                if match not in stage_matches[stage_id]:
                    stage_matches[stage_id].append(match)
    for stage_id, matches in list(stage_matches.items()):
        stage_matches[stage_id] = _dedupe_matches(matches)
    return stage_matches


def _function_detection_payload(
    spec: Dict[str, Any],
    *,
    status: str,
    evidence: List[str],
    matched_by_name: List[str],
    matched_by_content: List[str],
) -> Dict[str, Any]:
    return {
        "id": spec["id"],
        "label": spec["label"],
        "status": status,
        "status_label": _function_status_label(status),
        "status_label_zh": _function_status_label_zh(status),
        "policy": spec.get("policy", "iobrpy"),
        "stage_ids": spec.get("stage_ids", []),
        "matched_by_name": matched_by_name,
        "matched_by_content": matched_by_content,
        "evidence": evidence,
        "note": spec.get("note", _default_function_note(spec)),
    }


def _collect_function_name_matches(
    entries: Iterable[str],
    spec: Dict[str, Any],
    *,
    ignore_context: bool = False,
    entry_matcher: Optional[Any] = None,
) -> List[str]:
    matches = _collect_stage_matches(entries, spec.get("patterns", []), matcher=entry_matcher)
    name_context = None if ignore_context else spec.get("name_context")
    native_filenames = {name.lower() for name in spec.get("native_filenames", [])}
    if name_context:
        matches = [
            match
            for match in matches
            if _has_iobrpy_context(match, name_context) or Path(match).name.lower() in native_filenames
        ]
    return _dedupe_matches(matches)


def _collect_function_content_matches(
    preview_texts: Dict[str, str],
    spec: Dict[str, Any],
    *,
    ignore_context: bool = False,
) -> List[str]:
    checker_name = spec.get("content_checker")
    if not checker_name:
        return []
    checker = _FUNCTION_CONTENT_CHECKERS[checker_name]
    matches: List[str] = []
    for relpath, preview in preview_texts.items():
        content_context = None if ignore_context else spec.get("content_context")
        allow_content_without_context = bool(spec.get("allow_content_without_context"))
        native_filenames = {name.lower() for name in spec.get("native_filenames", [])}
        if content_context and not (
            allow_content_without_context
            or _has_iobrpy_context(relpath, content_context)
            or Path(relpath).name.lower() in native_filenames
        ):
            continue
        content_path_exact_parts = spec.get("content_path_exact_parts", [])
        if content_path_exact_parts and not _path_has_exact_part(relpath, *content_path_exact_parts):
            continue
        content_path_keywords = spec.get("content_path_keywords", [])
        if content_path_keywords and not _path_has_keyword(relpath, content_path_keywords):
            continue
        if not checker(preview):
            continue
        matches.append(relpath)
    return _dedupe_matches(matches)


def _detect_generic_stage_content_matches(preview_texts: Dict[str, str]) -> Dict[str, List[str]]:
    matches: Dict[str, List[str]] = {}
    excluded_tpm_keywords = (
        "signature",
        "ssgsea",
        "score",
        "cibersort",
        "estimate",
        "mcpcounter",
        "quantiseq",
        "epic",
        "ips",
        "bayesprism",
        "ligand",
        "receptor",
        "hla",
        "trust4",
        "tcr",
        "bcr",
        "cluster",
        "nmf",
    )

    for relpath, preview in preview_texts.items():
        lower = relpath.lower()
        if _looks_like_count_matrix(preview) and _path_has_keyword(
            relpath,
            ("readspergene", "featurecount", "featurecounts", "star", "star_count", "star.count"),
        ):
            matches.setdefault("star_merge", []).append(relpath)

        if not _looks_like_tpm_like_expression_matrix(preview):
            continue
        if any(keyword in lower for keyword in excluded_tpm_keywords):
            continue
        if "count" in lower and "tpm" not in lower:
            continue
        if (
            "03-tpm" in lower
            or "03_tpm" in lower
            or "tpm" in lower
            or ("expression" in lower and "matrix" in lower)
            or ("expr" in lower and "matrix" in lower)
        ):
            matches.setdefault("tpm_matrix", []).append(relpath)

    return {stage_id: _dedupe_matches(stage_matches) for stage_id, stage_matches in matches.items()}


def _external_hint_content_matches(
    preview_texts: Dict[str, str],
    hint_id: str,
) -> List[str]:
    matches: List[str] = []
    for relpath, preview in preview_texts.items():
        profile = _external_hint_profile(
            relpath,
            preview,
            hint_id,
            include_path_signals=False,
            respect_native_context=False,
        )
        if profile and profile.get("confidence") == "low":
            continue
        if profile:
            matches.append(relpath)
    return _dedupe_matches(matches)


def _looks_like_external_tcr_bcr_result_path(relpath: str) -> bool:
    lower = relpath.lower()
    name = Path(relpath.rstrip("/")).name.lower()
    if relpath.endswith("/") or _has_iobrpy_context(relpath, "trust4"):
        return False
    return any(
        token in lower
        for token in (
            "mixcr",
            "immunarch",
            "immunearch",
            "clonotype",
            "clonotypes",
            ".clones_",
            "airr",
            "cdr3",
            "vdj",
            "rearrangement",
            "align.report.json",
        )
    ) or name.endswith((".align.report.json", "_airr.tsv", "_report.tsv"))


def _external_hint_match_rank(
    relpath: str,
    hint_id: str,
    preview_texts: Optional[Dict[str, str]] = None,
) -> int:
    preview_texts = preview_texts or {}
    preview = preview_texts.get(relpath, "")
    profile = _external_hint_profile(
        relpath,
        preview,
        hint_id,
        respect_native_context=True,
    )
    rank = 0
    if profile:
        rank = _EXTERNAL_HINT_CONFIDENCE_RANK.get(str(profile.get("confidence", "low")), 0)

    if hint_id == "external_hla":
        if _looks_like_hla_per_sample_result_path(relpath) or Path(relpath.rstrip("/")).name.lower().endswith(".task.complete"):
            rank = max(rank, 2)
        if relpath.endswith("/") and not preview.strip():
            return 0
        return rank

    if hint_id == "external_tcr_bcr":
        if _looks_like_external_tcr_bcr_result_path(relpath):
            rank = max(rank, 2)
        if relpath.endswith("/") and not preview.strip():
            return 0
        return rank

    return rank


def _filter_external_hint_matches(
    matches: List[str],
    hint_id: str,
    preview_texts: Optional[Dict[str, str]] = None,
) -> List[str]:
    preview_texts = preview_texts or {}
    deduped = _dedupe_matches(matches, limit=32)
    ranked = [
        (match, _external_hint_match_rank(match, hint_id, preview_texts))
        for match in deduped
    ]
    kept = [match for match, rank in ranked if rank > 0]
    if not kept:
        return []

    pruned: List[str] = []
    normalized_kept = [item.rstrip("/\\") for item in kept]
    for match in kept:
        normalized = match.rstrip("/\\")
        if any(other != normalized and other.startswith(f"{normalized}/") for other in normalized_kept):
            continue
        pruned.append(match)
    return _dedupe_matches(pruned, limit=16)


def _has_native_hla_result_evidence(
    entries: Iterable[str],
    preview_texts: Optional[Dict[str, str]] = None,
) -> bool:
    preview_texts = preview_texts or {}
    merged_filenames = {"hla_result_merged.txt", "hla_result_merged.tsv", "hla_result_merged.csv"}
    for relpath in entries:
        name = Path(relpath.rstrip("/")).name.lower()
        if name in merged_filenames:
            return True
        if name == "hlafinal.type.txt" and _path_has_part_token(relpath, "spechla"):
            return True
        if name.endswith(".spechla.done"):
            return True
        preview = preview_texts.get(relpath, "")
        if preview and _looks_like_hla_merged(preview) and ("merged" in name or "cohort" in name):
            return True
    return False


def _resolve_function_status(
    spec: Dict[str, Any],
    *,
    matched_by_name: List[str],
    matched_by_content: List[str],
) -> str:
    if matched_by_content:
        return "confirmed_by_content"
    if matched_by_name:
        if spec.get("require_content_match"):
            return "not_detected"
        if spec.get("policy") == "reusable":
            return "reusable_result"
        return "likely_iobrpy_result"
    return "not_detected"


def _aggregate_wrapper_detections(
    entries: Iterable[str],
    detections: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    entry_set = set(entries)
    aggregated: List[Dict[str, Any]] = []

    def _merge_status(*statuses: str) -> str:
        ranked = [status for status in statuses if status]
        if not ranked:
            return "not_detected"
        return max(ranked, key=lambda item: _FUNCTION_STATUS_ORDER.get(item, -1))

    def _is_merged_deconvolution_name(entry: str) -> bool:
        filename = Path(entry).name.lower()
        return "deconvo_merged" in filename and _looks_like_delimited_result_name(filename)

    runall_spec = next(spec for spec in _FUNCTION_DETECTION_SPECS if spec["id"] == "runall")
    runall_matches: List[str] = []
    for token in ("01-qc", "02-salmon", "02-star", "03-tpm", "04-signatures", "05-tme", "06-LR_cal", "07-TCRBCR", "runall"):
        for entry in entries:
            if token.lower() in entry.lower():
                runall_matches.append(entry)
    for function_id in (
        "fastq_qc",
        "batch_salmon",
        "merge_salmon",
        "prepare_salmon",
        "batch_star_count",
        "merge_star_count",
        "count2tpm",
        "calculate_sig_score",
        "cibersort",
        "IPS",
        "estimate",
        "mcpcounter",
        "quantiseq",
        "epic",
        "bayesprism",
        "LR_cal",
        "trust4",
    ):
        runall_matches.extend(detections.get(function_id, {}).get("evidence", []))
    runall_matches = _dedupe_matches(runall_matches, limit=20)

    has_qc_layout = any("01-qc" in item.lower() for item in runall_matches)
    has_tpm_layout = any("03-tpm" in item.lower() for item in runall_matches)
    has_standard_layout = has_qc_layout and has_tpm_layout
    has_downstream_layout = any(
        token in item.lower()
        for item in runall_matches
        for token in ("04-signatures", "05-tme", "06-lr_cal")
    )
    has_quant_layout = any(
        token in item.lower()
        for item in runall_matches
        for token in ("02-salmon", "02-star")
    )
    if has_standard_layout and has_quant_layout and has_downstream_layout:
        runall_status = "confirmed_by_content"
    elif len(runall_matches) >= 3 and has_quant_layout:
        runall_status = "likely_iobrpy_result"
    else:
        runall_status = "not_detected"
    aggregated.append(
        _function_detection_payload(
            runall_spec,
            status=runall_status,
            evidence=_dedupe_matches(runall_matches),
            matched_by_name=_dedupe_matches(runall_matches),
            matched_by_content=_dedupe_matches(runall_matches if runall_status == "confirmed_by_content" else []),
        )
    )

    tme_profile_spec = next(spec for spec in _FUNCTION_DETECTION_SPECS if spec["id"] == "tme_profile")
    tme_matches: List[str] = []
    for token in ("01-signatures", "02-tme", "03-LR_cal", "tme_profile", "tme-profile"):
        for entry in entries:
            if token.lower() in entry.lower():
                tme_matches.append(entry)
    legacy_sig_matches = [
        entry
        for entry in entries
        if entry.lower().endswith(".rdata")
        and "signature-score" in entry.lower()
    ]
    legacy_tme_bundle_matches = [
        entry
        for entry in entries
        if (
            entry.lower().endswith("pdata.rdata")
            or "count2tpm-symbol.rdata" in entry.lower()
            or _is_merged_deconvolution_name(entry)
        )
    ]
    legacy_lr_matches = [
        entry
        for entry in entries
        if (
            _path_has_exact_part(entry, "03-lr_cal", "06-lr_cal")
            or _path_has_part_tokens(entry, "lr", "cal")
            or _path_has_part_token(entry, "ligand", "receptor")
        )
        and any(entry.lower().endswith(suffix) for suffix in (".csv", ".tsv", ".txt", ".rdata", ".rds"))
    ]
    tme_matches.extend(legacy_sig_matches)
    tme_matches.extend(legacy_tme_bundle_matches)
    tme_matches.extend(legacy_lr_matches)
    deconvo_merged_matches = [
        entry
        for entry in entries
        if _is_merged_deconvolution_name(entry)
    ]
    tme_matches.extend(deconvo_merged_matches)
    for function_id in ("calculate_sig_score", "cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic", "LR_cal"):
        tme_matches.extend(detections.get(function_id, {}).get("evidence", []))
    tme_matches = _dedupe_matches(tme_matches, limit=20)
    has_tme_dirs = any("01-signatures" in item.lower() for item in tme_matches) and any("02-tme" in item.lower() for item in tme_matches)
    has_lr = bool(detections.get("LR_cal", {}).get("evidence"))
    has_sig = bool(detections.get("calculate_sig_score", {}).get("evidence"))
    has_deconv = bool(deconvo_merged_matches) or any(
        detections.get(item, {}).get("evidence")
        for item in _WRAPPER_DEFAULT_DECONVOLUTION_IDS
    )
    has_legacy_bundle = bool(legacy_sig_matches) and bool(legacy_tme_bundle_matches)
    has_legacy_lr = bool(legacy_lr_matches)
    downstream_family_count = sum([has_sig, has_deconv, has_lr])
    if has_tme_dirs and has_sig and has_lr and has_deconv:
        tme_status = "confirmed_by_content"
    elif downstream_family_count >= 2 or has_legacy_bundle or (bool(legacy_sig_matches) and has_legacy_lr):
        tme_status = "likely_iobrpy_result"
    else:
        tme_status = "not_detected"
    aggregated.append(
        _function_detection_payload(
            tme_profile_spec,
            status=tme_status,
            evidence=tme_matches,
            matched_by_name=tme_matches if tme_status != "not_detected" else [],
            matched_by_content=tme_matches if tme_status == "confirmed_by_content" else [],
        )
    )

    log2_spec = next(spec for spec in _FUNCTION_DETECTION_SPECS if spec["id"] == "log2_eset")
    log2_existing = detections.get("log2_eset", {})
    log2_matches: List[str] = []
    tpm_outputs = [
        entry
        for entry in entries
        if Path(entry).name.lower() == "tpm_matrix.csv"
        and ("03-tpm/" in entry.lower() or entry.lower().startswith("03-tpm"))
    ]
    has_upstream_tpm_prep = bool(detections.get("prepare_salmon", {}).get("evidence")) or bool(
        detections.get("count2tpm", {}).get("evidence")
    )
    if tpm_outputs and has_upstream_tpm_prep:
        log2_matches.extend(tpm_outputs)
    log2_name_matches = _dedupe_matches(log2_existing.get("matched_by_name", []))
    log2_content_matches = _dedupe_matches(log2_existing.get("matched_by_content", []) + log2_matches)
    log2_status = _merge_status(
        log2_existing.get("status", "not_detected"),
        "confirmed_by_content" if log2_matches else "not_detected",
    )
    aggregated.append(
        _function_detection_payload(
            log2_spec,
            status=log2_status,
            evidence=_dedupe_matches(log2_name_matches + log2_content_matches),
            matched_by_name=log2_name_matches,
            matched_by_content=log2_content_matches,
        )
    )

    mouse2human_spec = next(spec for spec in _FUNCTION_DETECTION_SPECS if spec["id"] == "mouse2human")
    mouse2human_existing = detections.get("mouse2human", {})
    mouse2human_eset_detection = detections.get("mouse2human_eset", {})
    mouse2human_name_matches = _dedupe_matches(
        mouse2human_existing.get("matched_by_name", []) + mouse2human_eset_detection.get("matched_by_name", [])
    )
    mouse2human_content_matches = _dedupe_matches(
        mouse2human_existing.get("matched_by_content", []) + mouse2human_eset_detection.get("matched_by_content", [])
    )
    mouse2human_status = _merge_status(
        mouse2human_existing.get("status", "not_detected"),
        mouse2human_eset_detection.get("status", "not_detected"),
    )
    aggregated.append(
        _function_detection_payload(
            mouse2human_spec,
            status=mouse2human_status,
            evidence=_dedupe_matches(mouse2human_name_matches + mouse2human_content_matches),
            matched_by_name=mouse2human_name_matches,
            matched_by_content=mouse2human_content_matches,
        )
    )
    return aggregated


def _detect_content_signatures(
    root: Path,
    entries: Iterable[str],
    *,
    priority_patterns: Optional[Iterable[str]] = None,
    max_preview_files: Optional[int] = None,
    entry_matcher: Optional[Any] = None,
) -> Dict[str, Any]:
    stage_matches: Dict[str, List[str]] = {}
    method_matches: Dict[str, List[str]] = {}
    method_labels = {spec["id"]: spec["label"] for spec in _DECONVOLUTION_METHOD_SPECS}
    preview_texts: Dict[str, str] = {}
    preview_files = 0
    normalized_priority_patterns = list(priority_patterns or [])

    for relpath in entries:
        if not _is_text_like_result(relpath):
            continue
        path_matches = bool(normalized_priority_patterns and _matches_any_pattern(relpath, normalized_priority_patterns))
        force_preview = _should_force_preview_result(relpath)
        if max_preview_files is not None and preview_files >= max_preview_files and not path_matches:
            if not force_preview:
                continue
        absolute = root / relpath
        if not absolute.is_file():
            continue
        preview = _read_text_preview(absolute)
        if preview:
            preview_texts[relpath] = preview
            if not path_matches and not force_preview:
                preview_files += 1

    for stage_id, matches in _detect_generic_stage_content_matches(preview_texts).items():
        for match in matches:
            stage_matches.setdefault(stage_id, [])
            if match not in stage_matches[stage_id]:
                stage_matches[stage_id].append(match)

    function_detections: List[Dict[str, Any]] = []
    detections_by_id: Dict[str, Dict[str, Any]] = {}
    for spec in _FUNCTION_DETECTION_SPECS:
        if spec["id"] in {"runall", "tme_profile"}:
            continue
        matched_by_name = _collect_function_name_matches(entries, spec, entry_matcher=entry_matcher)
        matched_by_content = _collect_function_content_matches(preview_texts, spec)
        status = _resolve_function_status(
            spec,
            matched_by_name=matched_by_name,
            matched_by_content=matched_by_content,
        )
        evidence = _dedupe_matches(matched_by_name + matched_by_content)
        if status == "not_detected" and spec.get("require_content_match"):
            evidence = []
        detection = _function_detection_payload(
            spec,
            status=status,
            evidence=evidence,
            matched_by_name=matched_by_name,
            matched_by_content=matched_by_content,
        )
        function_detections.append(detection)
        detections_by_id[spec["id"]] = detection

    function_detections.extend(_aggregate_wrapper_detections(entries, detections_by_id))
    detections_by_id = {item["id"]: item for item in function_detections}
    spec_index = {spec["id"]: spec for spec in _FUNCTION_DETECTION_SPECS}
    for method_id, merged_matches in _collect_merged_deconvolution_method_matches(preview_texts).items():
        function_id = _DECONVOLUTION_FUNCTION_ID_BY_CANONICAL.get(method_id, method_id)
        spec = spec_index.get(function_id)
        if spec is None:
            continue
        existing = detections_by_id.get(function_id, _empty_function_detection(spec))
        merged_detection = _function_detection_payload(
            spec,
            status="confirmed_by_content",
            evidence=merged_matches,
            matched_by_name=[],
            matched_by_content=merged_matches,
        )
        detections_by_id[function_id] = _build_function_detection_from_sources(
            spec,
            existing,
            merged_detection,
        )
    function_detections = [
        detections_by_id[spec["id"]]
        for spec in _FUNCTION_DETECTION_SPECS
        if spec["id"] in detections_by_id
    ]
    detections_by_id = {item["id"]: item for item in function_detections}

    for detection in function_detections:
        if detection["status"] == "not_detected":
            continue
        if detection.get("id") in {"runall", "tme_profile"}:
            continue
        has_confirmed_content = bool(detection.get("matched_by_content"))
        for stage_id in detection.get("stage_ids", []):
            if stage_id in _STRICT_IOBRPY_STAGE_IDS and not has_confirmed_content:
                continue
            stage_matches.setdefault(stage_id, [])
            evidence_paths = (
                detection.get("matched_by_content", [])
                if stage_id in _STRICT_IOBRPY_STAGE_IDS
                else detection.get("evidence", [])
            )
            for evidence in evidence_paths:
                if evidence not in stage_matches[stage_id]:
                    stage_matches[stage_id].append(evidence)

    for method_spec in _DECONVOLUTION_METHOD_SPECS:
        detection = detections_by_id.get(_DECONVOLUTION_FUNCTION_ID_BY_CANONICAL.get(method_spec["id"], method_spec["id"]), {})
        status = detection.get("status")
        has_confirmed_content = bool(detection.get("matched_by_content"))
        if status and status != "not_detected" and has_confirmed_content:
            method_matches[method_spec["id"]] = detection.get("evidence", [])
            if method_spec["id"] not in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET:
                continue
            for evidence in detection.get("evidence", []):
                stage_matches.setdefault("deconvolution", [])
                if evidence not in stage_matches["deconvolution"]:
                    stage_matches["deconvolution"].append(evidence)

    for stage_id, matches in list(stage_matches.items()):
        stage_matches[stage_id] = _dedupe_matches(matches)
    for method_id, matches in list(method_matches.items()):
        method_matches[method_id] = _dedupe_matches(matches)

    default_method_specs = [
        spec for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET
    ]
    optional_method_specs = [
        spec for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] not in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET
    ]
    detected_method_ids = [spec["id"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] in method_matches]
    missing_method_ids = [spec["id"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] not in method_matches]
    default_detected_method_ids = [spec["id"] for spec in default_method_specs if spec["id"] in method_matches]
    default_missing_method_ids = [spec["id"] for spec in default_method_specs if spec["id"] not in method_matches]
    optional_detected_method_ids = [spec["id"] for spec in optional_method_specs if spec["id"] in method_matches]
    optional_missing_method_ids = [spec["id"] for spec in optional_method_specs if spec["id"] not in method_matches]
    reusable_function_ids = [
        item["id"]
        for item in function_detections
        if item["status"] == "reusable_result"
    ]
    likely_iobrpy_function_ids = [
        item["id"]
        for item in function_detections
        if item["status"] == "likely_iobrpy_result"
    ]
    content_verified_function_ids = [
        item["id"]
        for item in function_detections
        if item["status"] == "confirmed_by_content"
    ]

    return {
        "stage_matches": stage_matches,
        "content_verified_stage_ids": sorted(
            {
                stage_id
                for item in function_detections
                if item["status"] == "confirmed_by_content"
                for stage_id in item.get("stage_ids", [])
            }
        ),
        "function_detections": function_detections,
        "detected_function_ids": [item["id"] for item in function_detections if item["status"] != "not_detected"],
        "reusable_function_ids": reusable_function_ids,
        "likely_iobrpy_function_ids": likely_iobrpy_function_ids,
        "content_verified_function_ids": content_verified_function_ids,
        "deconvolution_methods": {
            "total_supported": len(_DECONVOLUTION_METHOD_SPECS),
            "total_available": len(_DECONVOLUTION_METHOD_SPECS),
            "default_bundle_total": len(default_method_specs),
            "detected_ids": detected_method_ids,
            "detected_labels": [method_labels[item] for item in detected_method_ids],
            "detected_count": len(detected_method_ids),
            "missing_ids": missing_method_ids,
            "missing_labels": [method_labels[item] for item in missing_method_ids],
            "default_bundle_ids": list(_WRAPPER_DEFAULT_DECONVOLUTION_IDS),
            "default_bundle_detected_ids": default_detected_method_ids,
            "default_bundle_detected_labels": [method_labels[item] for item in default_detected_method_ids],
            "default_bundle_detected_count": len(default_detected_method_ids),
            "default_bundle_missing_ids": default_missing_method_ids,
            "default_bundle_missing_labels": [method_labels[item] for item in default_missing_method_ids],
            "optional_ids": [spec["id"] for spec in optional_method_specs],
            "optional_labels": [spec["label"] for spec in optional_method_specs],
            "optional_detected_ids": optional_detected_method_ids,
            "optional_detected_labels": [method_labels[item] for item in optional_detected_method_ids],
            "optional_missing_ids": optional_missing_method_ids,
            "optional_missing_labels": [method_labels[item] for item in optional_missing_method_ids],
            "evidence": method_matches,
        },
        "preview_texts": preview_texts,
        "preview_file_count": len(preview_texts),
        "preview_file_limit": max_preview_files,
    }


def _find_latest_stage(stage_ids: Set[str]) -> Optional[str]:
    for stage_id in reversed(_STAGE_DISPLAY_ORDER):
        if stage_id in stage_ids:
            return stage_id
    return None


def _content_signature_preview_limit(
    *,
    quick: bool,
    focus: Optional[Iterable[str]] = None,
) -> Optional[int]:
    if quick:
        return _QUICK_MAX_PREVIEW_FILES
    if list(focus or []):
        return _FOCUSED_DEEP_MAX_PREVIEW_FILES
    return _FULL_MAX_PREVIEW_FILES


def _find_current_stage(stage_ids: Set[str], stage_index: Dict[str, Dict[str, Any]]) -> Optional[str]:
    primary_candidates = [stage_id for stage_id in _STAGE_DISPLAY_ORDER if stage_id in stage_ids]
    if not primary_candidates:
        return None

    mainline_groups = {"input", "salmon", "star", "core"}
    nontrivial_mainline = [
        stage_id
        for stage_id in primary_candidates
        if stage_index.get(stage_id, {}).get("group") in mainline_groups
        and stage_id not in {"raw_fastq", "fastq_qc"}
    ]
    if nontrivial_mainline:
        return nontrivial_mainline[-1]

    return primary_candidates[-1]


def _supplementary_stage_ids(
    completed: Set[str],
    stage_index: Dict[str, Dict[str, Any]],
    current_stage: Optional[str],
) -> List[str]:
    if current_stage and stage_index.get(current_stage, {}).get("group") != "immune":
        return [
            stage_id
            for stage_id in _ordered_stage_ids(completed)
            if stage_index.get(stage_id, {}).get("group") == "immune"
        ]
    return []


def _stage_summary_list(
    stage_ids: Iterable[str],
    stage_index: Dict[str, Dict[str, Any]],
    *,
    completed: bool,
) -> List[Dict[str, str]]:
    summaries: List[Dict[str, str]] = []
    for stage_id in _ordered_stage_ids(stage_ids):
        copy = _stage_copy(stage_id, stage_index)
        summaries.append(
            {
                "id": stage_id,
                "label": copy["label"],
                "text": copy["done_text"] if completed else copy["next_text"],
                "description": copy["description"],
                "label_zh": copy["label_zh"],
                "text_zh": copy["done_text_zh"] if completed else copy["next_text_zh"],
                "description_zh": copy["description_zh"],
            }
        )
    return summaries


def _has_salmon_tpm_like_evidence(strong_evidence: Dict[str, List[str]]) -> bool:
    for evidence in strong_evidence.get("salmon_merge", []):
        if "tpm" in evidence.lower():
            return True
    return False


def _localized_group_value(group: Dict[str, Any], key: str, language: str = "en") -> str:
    if language == "zh":
        return group.get(f"{key}_zh", group.get(key, ""))
    return group.get(key, "")


def _agent_term_tokens(text: str) -> List[str]:
    tokens: List[str] = []
    current: List[str] = []
    for char in text:
        if char.isalnum():
            current.append(char.lower())
        else:
            if current:
                tokens.append("".join(current))
                current = []
    if current:
        tokens.append("".join(current))
    return tokens


def _dedupe_terms(terms: Iterable[str], *, limit: int = 12) -> List[str]:
    deduped: List[str] = []
    seen: Set[str] = set()
    for term in terms:
        cleaned = str(term).strip()
        if not cleaned:
            continue
        lowered = cleaned.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        deduped.append(cleaned)
        if len(deduped) >= limit:
            break
    return deduped


def _checklist_marker(checked: bool, status: str) -> str:
    if status == "external":
        return "⚠"
    if status in {"partial", "inferred_existing"}:
        return "△"
    return "✓" if checked else "[ ]"


def _workflow_checklist(
    completed: Set[str],
    strong_evidence: Dict[str, List[str]],
    external_analysis_hints: List[Dict[str, Any]],
    deconvolution_methods: Dict[str, Any],
    stage_index: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    checklist: List[Dict[str, Any]] = []
    external_ids = {item.get("id") for item in external_analysis_hints}

    for group in _WORKFLOW_CHECKLIST_GROUPS:
        stage_ids = group.get("stage_ids", [])
        matched_stage_ids = [stage_id for stage_id in stage_ids if stage_id in completed]
        pending_stage_ids = [stage_id for stage_id in stage_ids if stage_id not in completed]
        external_hint_ids = group.get("external_hint_ids", [])
        matched_external_hint_ids = [hint_id for hint_id in external_hint_ids if hint_id in external_ids]
        matched_stage_labels = _ordered_labels(matched_stage_ids, stage_index)
        matched_stage_labels_zh = _ordered_labels(matched_stage_ids, stage_index, language="zh")
        pending_stage_labels = _ordered_labels(pending_stage_ids, stage_index)
        pending_stage_labels_zh = _ordered_labels(pending_stage_ids, stage_index, language="zh")

        checked = bool(matched_stage_ids or matched_external_hint_ids)
        if group["id"] == "tpm_matrix_ready" and _has_salmon_tpm_like_evidence(strong_evidence):
            checked = True

        status = "completed" if checked else "pending"
        if group["id"] == "immune_deconvolution":
            detected_labels = deconvolution_methods.get("detected_labels", [])
            missing_labels = deconvolution_methods.get("missing_labels", [])
            detected_count = deconvolution_methods.get("detected_count", 0)
            total_supported = deconvolution_methods.get("total_supported", 0)
            default_missing_labels = deconvolution_methods.get("default_bundle_missing_labels", missing_labels)
            optional_detected_labels = deconvolution_methods.get("optional_detected_labels", [])
            optional_missing_labels = deconvolution_methods.get("optional_missing_labels", [])
            optional_detected_text = ", ".join(optional_detected_labels)
            optional_missing_text = ", ".join(optional_missing_labels)
            optional_note = ""
            optional_note_zh = ""
            if optional_detected_labels:
                optional_note = f" Optional standalone deconvolution already detected: {optional_detected_text}."
                optional_note_zh = f" 已检测到独立可选反卷积方法：{'、'.join(optional_detected_labels)}。"
            elif optional_missing_labels:
                optional_note = (
                    f" Optional standalone deconvolution such as {optional_missing_text} "
                    "must be run separately; it is not run by tme_profile/runall."
                )
                optional_note_zh = (
                    f" 独立可选反卷积方法（如 {'、'.join(optional_missing_labels)}）需要单独运行；"
                    "不由 tme_profile/runall 自动运行。"
                )
            if detected_count:
                checked = True
                status = "completed" if not missing_labels else "partial"
                detected_text = ", ".join(detected_labels)
                detected_text_zh = "、".join(detected_labels)
                missing_text_zh = "、".join(missing_labels)
                default_missing_text = ", ".join(default_missing_labels)
                default_missing_text_zh = "、".join(default_missing_labels)
                default_gap_note = (
                    f" Default bundled methods still missing: {default_missing_text}."
                    if default_missing_labels
                    else ""
                )
                default_gap_note_zh = (
                    f" 默认 bundle 方法仍缺少：{default_missing_text_zh}。"
                    if default_missing_labels
                    else ""
                )
                if missing_labels:
                    text = (
                        f"Detected {detected_count} of {total_supported} immune-deconvolution methods: "
                        f"{detected_text}. Still available to run: {', '.join(missing_labels)}."
                        f"{default_gap_note}"
                        f"{optional_note}"
                    )
                    text_zh = (
                        f"已检测到 {detected_count}/{total_supported} 个免疫反卷积方法："
                        f"{detected_text_zh}。仍可继续运行：{missing_text_zh}。"
                        f"{default_gap_note_zh}"
                        f"{optional_note_zh}"
                    )
                else:
                    text = (
                        f"All {total_supported} immune-deconvolution methods were detected: "
                        f"{detected_text}.{optional_note}"
                    )
                    text_zh = (
                        f"已检测到全部 {total_supported} 个免疫反卷积方法："
                        f"{detected_text_zh}。{optional_note_zh}"
                    )
            else:
                if optional_detected_labels:
                    checked = True
                    status = "partial"
                    text = (
                        "Only optional standalone deconvolution output was detected: "
                        f"{optional_detected_text}. Default bundled methods still available to run: "
                        f"{', '.join(missing_labels)}."
                    )
                    text_zh = (
                        f"仅检测到独立可选反卷积输出：{'、'.join(optional_detected_labels)}。"
                        f"默认 bundle 方法仍可运行：{'、'.join(missing_labels)}。"
                    )
                else:
                    checked = False
                    status = "pending"
                    text = (
                        "No immune-deconvolution outputs were confirmed yet. "
                        f"Supported methods include: {', '.join(missing_labels)}."
                        f" Default bundled methods include: {', '.join(default_missing_labels)}."
                        f"{optional_note}"
                    )
                    text_zh = (
                        "尚未确认免疫反卷积输出。"
                        f"支持的方法包括：{'、'.join(missing_labels)}。"
                        f"默认 bundle 方法包括：{'、'.join(default_missing_labels)}。"
                        f"{optional_note_zh}"
                    )
        elif matched_external_hint_ids and not matched_stage_ids and group.get("external_text"):
            status = "external"
            text = group["external_text"]
            text_zh = _localized_group_value(group, "external_text", "zh")
        elif checked and group.get("partial_text") and 0 < len(matched_stage_ids) < len(stage_ids):
            status = "partial"
            text = group["partial_text"]
            text_zh = _localized_group_value(group, "partial_text", "zh")
        elif checked:
            status = "completed"
            text = group["completed_text"]
            text_zh = _localized_group_value(group, "completed_text", "zh")
        else:
            status = "pending"
            text = group["pending_text"]
            text_zh = _localized_group_value(group, "pending_text", "zh")

        if status == "partial":
            detail_bits: List[str] = []
            detail_bits_zh: List[str] = []
            if matched_stage_labels:
                detail_bits.append(f"Completed in this section: {', '.join(matched_stage_labels)}.")
            if matched_stage_labels_zh:
                detail_bits_zh.append(f"这一部分已完成：{'、'.join(matched_stage_labels_zh)}。")
            if pending_stage_labels:
                detail_bits.append(f"Still not confirmed: {', '.join(pending_stage_labels)}.")
            if pending_stage_labels_zh:
                detail_bits_zh.append(f"这一部分尚未确认：{'、'.join(pending_stage_labels_zh)}。")
            if detail_bits:
                text = f"{text} {' '.join(detail_bits)}"
            if detail_bits_zh:
                text_zh = f"{text_zh} {' '.join(detail_bits_zh)}"

        checklist.append(
            {
                "id": group["id"],
                "checked": checked,
                "marker": "✓" if checked else "[ ]",
                "label": group["label"],
                "label_zh": _localized_group_value(group, "label", "zh"),
                "status": status,
                "text": text,
                "text_zh": text_zh,
                "details": group["details"],
                "details_zh": _localized_group_value(group, "details", "zh"),
                "matched_stage_ids": matched_stage_ids,
                "pending_stage_ids": pending_stage_ids,
                "matched_stage_labels": matched_stage_labels,
                "matched_stage_labels_zh": matched_stage_labels_zh,
                "pending_stage_labels": pending_stage_labels,
                "pending_stage_labels_zh": pending_stage_labels_zh,
                "matched_external_hint_ids": matched_external_hint_ids,
            }
        )
    return checklist

def _resolve_supported_weak_stages(
    strong_detected: Set[str],
    weak_detected: Set[str],
    stage_index: Dict[str, Dict[str, Any]],
) -> Set[str]:
    supported: Set[str] = set()
    changed = True
    while changed:
        changed = False
        support_base = _closure_with_implies(strong_detected | supported, stage_index)
        for stage_id in weak_detected:
            if stage_id in strong_detected or stage_id in supported:
                continue
            if stage_id in _WEAK_STAGE_REQUIRES_CONTENT:
                continue
            implied = set(stage_index.get(stage_id, {}).get("implies", []))
            if implied and implied.intersection(support_base):
                supported.add(stage_id)
                changed = True
    return supported


def _stage_confidence(
    stage_id: str,
    strong_detected: Set[str],
    supported_weak: Set[str],
    completed: Set[str],
) -> str:
    if stage_id in strong_detected:
        return "high"
    if stage_id in supported_weak:
        return "medium"
    if stage_id in completed:
        return "inferred"
    return "none"


def _stage_completion_rows(
    rules: Dict[str, Any],
    stage_index: Dict[str, Dict[str, Any]],
    strong_detected: Set[str],
    supported_weak: Set[str],
    completed: Set[str],
    strong_evidence: Dict[str, List[str]],
    weak_evidence: Dict[str, List[str]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for stage in rules.get("stages", []):
        stage_id = stage["id"]
        confidence = _stage_confidence(stage_id, strong_detected, supported_weak, completed)
        copy = _stage_copy(stage_id, stage_index)
        rows.append(
            {
                "id": stage_id,
                "label": copy["label"],
                "label_zh": copy["label_zh"],
                "technical_label": stage.get("technical_label", stage["id"]),
                "group": stage.get("group"),
                "completed": stage_id in completed,
                "detected_directly": stage_id in strong_detected or stage_id in supported_weak,
                "confidence": confidence,
                "done_text": copy["done_text"],
                "done_text_zh": copy["done_text_zh"],
                "pending_text": copy["pending_text"],
                "pending_text_zh": copy["pending_text_zh"],
                "description": copy["description"],
                "description_zh": copy["description_zh"],
                "next_text": copy["next_text"],
                "next_text_zh": copy["next_text_zh"],
                "evidence": {
                    "strong": strong_evidence.get(stage_id, []),
                    "weak": weak_evidence.get(stage_id, []),
                },
                "roadmap_node": _STAGE_TO_NODE.get(stage_id, copy["label"]),
            }
        )
    return rows


def _infer_next_stage_ids(completed: Set[str]) -> List[str]:
    if "tpm_matrix" in completed:
        next_ids: List[str] = []
        if "signature_scoring" not in completed:
            next_ids.append("signature_scoring")
        if "deconvolution" not in completed:
            next_ids.append("deconvolution")
        if "ligand_receptor" not in completed:
            next_ids.append("ligand_receptor")
        if "deconvolution" in completed and "clustering" not in completed:
            next_ids.append("clustering")
        return next_ids
    if (
        "raw_fastq" in completed
        and "fastq_qc" not in completed
        and not completed.intersection(
            {
                "salmon_quant",
                "salmon_merge",
                "prepare_salmon",
                "star_quant",
                "star_merge",
                "count2tpm",
                "tpm_matrix",
                "signature_scoring",
                "deconvolution",
                "ligand_receptor",
                "clustering",
                "trust4",
                "spechla",
                "hla_typing",
            }
        )
    ):
        return ["fastq_qc", "trust4"]
    if "salmon_quant" in completed and "salmon_merge" not in completed:
        return ["salmon_merge"]
    if "salmon_merge" in completed and "prepare_salmon" not in completed:
        return ["prepare_salmon"]
    if "star_quant" in completed and "star_merge" not in completed:
        return ["star_merge"]
    if "star_merge" in completed and "count2tpm" not in completed:
        return ["count2tpm"]
    if "fastq_qc" in completed and "tpm_matrix" not in completed:
        if "salmon_quant" not in completed:
            return ["salmon_quant", "star_quant", "trust4"]
        if "star_quant" not in completed:
            return ["star_quant", "trust4"]
    return []


def _infer_suggested_commands(completed: Set[str], next_stage_ids: List[str]) -> List[str]:
    if "tpm_matrix" in completed and not any(
        stage in completed for stage in ("signature_scoring", "deconvolution", "ligand_receptor")
    ):
        return [_NEXT_COMMANDS["tme_profile"]]

    commands: List[str] = []
    for stage_id in next_stage_ids:
        if stage_id == "signature_scoring":
            commands.append(_NEXT_COMMANDS["tme_profile"])
            continue
        command = _NEXT_COMMANDS.get(stage_id)
        if command and command not in commands:
            commands.append(command)
    return commands


_RESOURCE_SENSITIVE_PARAMETER_NAMES = {
    "threads",
    "num_threads",
    "parallel_size",
    "num_processes",
    "batch_size",
    "t",
    "j",
}


def _command_id_from_invocation(invocation: str) -> Optional[str]:
    tokens = str(invocation).split()
    for token in tokens:
        if token in {"iobrpy-cli", "iobrpy"}:
            continue
        if token.startswith("-"):
            continue
        return token
    return None


def _agent_command_hint(command_id: str) -> Dict[str, Any]:
    return dict(load_agent_parameter_hints().get("commands", {}).get(command_id, {}))


def _agent_command_template(command_id: str) -> str:
    hint = _agent_command_hint(command_id)
    template = hint.get("command_template")
    if isinstance(template, str) and template.strip():
        return template
    for command in _NEXT_COMMANDS.values():
        if _command_id_from_invocation(command) == command_id:
            return command
    return f"iobrpy-cli {command_id} <parameters>"


def _fallback_parameter_placeholder(name: str, param_type: Optional[str]) -> str:
    if param_type == "path":
        return f"path_{name}"
    if param_type == "choice":
        return f"{name}_choice"
    return name


def _parameter_hints_for_command(command_id: str) -> List[Dict[str, Any]]:
    spec = load_iobrpy_required_params().get(command_id, {})
    if not isinstance(spec, dict):
        return []
    command_hint = _agent_command_hint(command_id)
    overrides = command_hint.get("parameter_overrides", {})
    if not isinstance(overrides, dict):
        overrides = {}

    required = list(spec.get("required", []) or [])
    optional = list(spec.get("optional", []) or [])
    known_params = required + [name for name in optional if name not in required]
    priority = command_hint.get("priority_parameters")
    if not isinstance(priority, list) or not priority:
        priority = required + [
            name
            for name in optional
            if name in _RESOURCE_SENSITIVE_PARAMETER_NAMES or name in spec.get("choices", {})
        ]

    cli_flags = spec.get("cli_flags", {}) if isinstance(spec.get("cli_flags"), dict) else {}
    param_types = spec.get("param_types", {}) if isinstance(spec.get("param_types"), dict) else {}
    param_hints = spec.get("param_hints", {}) if isinstance(spec.get("param_hints"), dict) else {}
    choices = spec.get("choices", {}) if isinstance(spec.get("choices"), dict) else {}
    defaults = spec.get("optional_defaults", {}) if isinstance(spec.get("optional_defaults"), dict) else {}

    hints: List[Dict[str, Any]] = []
    for name in priority:
        if name not in known_params:
            continue
        override = overrides.get(name, {})
        if not isinstance(override, dict):
            override = {}
        base_hint = param_hints.get(name, {}) if isinstance(param_hints.get(name), dict) else {}
        param_type = param_types.get(name)
        hint = {
            "name": name,
            "flag": cli_flags.get(name),
            "required": name in required,
            "type": param_type,
            "placeholder": override.get("placeholder") or _fallback_parameter_placeholder(name, param_type),
            "prompt_zh": override.get("prompt_zh") or base_hint.get("zh"),
            "prompt_en": override.get("prompt_en") or base_hint.get("en"),
            "choices": choices.get(name),
            "default": defaults.get(name),
            "resource_sensitive": name in _RESOURCE_SENSITIVE_PARAMETER_NAMES
            or bool(override.get("resource_sensitive")),
            "contextual_placeholders": override.get("contextual_placeholders"),
        }
        hints.append({key: value for key, value in hint.items() if value not in (None, [], {})})
    return hints


def _suggested_command_detail(command: str, command_id: Optional[str] = None) -> Dict[str, Any]:
    resolved_command_id = command_id or _command_id_from_invocation(command) or ""
    spec = load_iobrpy_required_params().get(resolved_command_id, {})
    command_hint = _agent_command_hint(resolved_command_id)
    detail = {
        "command_id": resolved_command_id,
        "command": command,
        "command_template": command_hint.get("command_template") or command,
        "summary": spec.get("function_summary", {}).get("en") if isinstance(spec, dict) else None,
        "summary_zh": spec.get("function_summary", {}).get("zh") if isinstance(spec, dict) else None,
        "parameter_hints": _parameter_hints_for_command(resolved_command_id),
        "parameter_groups": command_hint.get("parameter_groups"),
        "parameter_source": "src/iobrpy/RAG_MCP/iobrpy_required_params.json",
        "agent_hint_source": "src/iobrpy/RAG_MCP/iobrpy_agent_parameter_hints.json",
    }
    return {key: value for key, value in detail.items() if value not in (None, [], {})}


def _suggested_command_details(commands: List[str]) -> List[Dict[str, Any]]:
    details: List[Dict[str, Any]] = []
    seen: Set[str] = set()
    for command in commands:
        command_id = _command_id_from_invocation(command)
        if not command_id:
            continue
        key = f"{command_id}\0{command}"
        if key in seen:
            continue
        seen.add(key)
        details.append(_suggested_command_detail(command, command_id))
    return details


def _missing_downstream_suggestion(
    suggestion_id: str,
    command_id: str,
    *,
    label: str,
    label_zh: str,
    reason: str,
    reason_zh: str,
    missing_analysis_ids: List[str],
    missing_analysis_labels: List[str],
    missing_analysis_labels_zh: Optional[List[str]] = None,
    notes: Optional[List[str]] = None,
    notes_zh: Optional[List[str]] = None,
) -> Dict[str, Any]:
    command = _agent_command_template(command_id)
    return {
        "id": suggestion_id,
        "command_id": command_id,
        "label": label,
        "label_zh": label_zh,
        "reason": reason,
        "reason_zh": reason_zh,
        "missing_analysis_ids": missing_analysis_ids,
        "missing_analysis_labels": missing_analysis_labels,
        "missing_analysis_labels_zh": missing_analysis_labels_zh or missing_analysis_labels,
        "suggested_commands": [command],
        "suggested_command_details": [_suggested_command_detail(command, command_id)],
        "notes": notes or [],
        "notes_zh": notes_zh or [],
    }


def _missing_downstream_analysis_suggestions(
    completed: Set[str],
    deconvolution_methods: Dict[str, Any],
) -> List[Dict[str, Any]]:
    suggestions: List[Dict[str, Any]] = []
    if "tpm_matrix" in completed:
        default_missing_ids = list(deconvolution_methods.get("default_bundle_missing_ids", []) or [])
        default_missing_labels = list(deconvolution_methods.get("default_bundle_missing_labels", []) or [])
        optional_missing_ids = list(deconvolution_methods.get("optional_missing_ids", []) or [])
        optional_missing_labels = list(deconvolution_methods.get("optional_missing_labels", []) or [])

        tme_profile_missing_ids: List[str] = []
        tme_profile_missing_labels: List[str] = []
        tme_profile_missing_labels_zh: List[str] = []
        if "signature_scoring" not in completed:
            tme_profile_missing_ids.append("signature_scoring")
            tme_profile_missing_labels.append("Signature scoring")
            tme_profile_missing_labels_zh.append("特征集评分")
        if default_missing_ids:
            tme_profile_missing_ids.extend(default_missing_ids)
            tme_profile_missing_labels.extend(default_missing_labels)
            tme_profile_missing_labels_zh.extend(default_missing_labels)
        if "ligand_receptor" not in completed:
            tme_profile_missing_ids.append("ligand_receptor")
            tme_profile_missing_labels.append("Ligand-receptor analysis")
            tme_profile_missing_labels_zh.append("配体-受体分析")

        if tme_profile_missing_ids:
            suggestions.append(
                _missing_downstream_suggestion(
                    "fill_tme_profile_default_bundle",
                    "tme_profile",
                    label="Fill missing default TME profile analyses",
                    label_zh="补齐默认 TME profile 下游分析",
                    reason="A TPM matrix is available and one or more default tme_profile/runall downstream outputs are missing.",
                    reason_zh="已检测到 TPM 矩阵，且默认 tme_profile/runall 下游结果仍有缺失。",
                    missing_analysis_ids=tme_profile_missing_ids,
                    missing_analysis_labels=tme_profile_missing_labels,
                    missing_analysis_labels_zh=tme_profile_missing_labels_zh,
                    notes_zh=[
                        "tme_profile 默认补齐 signature scoring、CIBERSORT、EPIC、quanTIseq、MCPcounter、ESTIMATE、IPS 和 LR_cal 相关结果。",
                        "BayesPrism 是独立可选分支，不由 tme_profile/runall 自动运行；如缺失会单独列为建议。",
                    ],
                    notes=[
                        "tme_profile fills the default signature scoring, CIBERSORT, EPIC, quanTIseq, MCPcounter, ESTIMATE, IPS, and LR_cal outputs.",
                        "BayesPrism is a standalone optional branch and is not run automatically by tme_profile/runall; list it separately when missing.",
                    ],
                )
            )

        if "bayesprism" in optional_missing_ids:
            bayes_label = "BayesPrism"
            bayes_index = optional_missing_ids.index("bayesprism")
            if bayes_index < len(optional_missing_labels):
                bayes_label = optional_missing_labels[bayes_index]
            suggestions.append(
                _missing_downstream_suggestion(
                    "run_bayesprism",
                    "bayesprism",
                    label="Run standalone BayesPrism deconvolution",
                    label_zh="运行独立 BayesPrism 反卷积",
                    reason="BayesPrism is an immune deconvolution method, but it is a standalone optional branch outside the default tme_profile/runall bundle.",
                    reason_zh="BayesPrism 也是免疫反卷积方法，但属于独立可选分支，不在默认 tme_profile/runall 六方法包里。",
                    missing_analysis_ids=["bayesprism"],
                    missing_analysis_labels=[bayes_label],
                    missing_analysis_labels_zh=[bayes_label],
                    notes_zh=[
                        "如果不使用自定义单细胞参考，只需要确认 input、output 和 threads。",
                        "如果使用自定义单细胞参考，需要同时确认 sc_dat、cell_state_labels、cell_type_labels 和 key。",
                    ],
                    notes=[
                        "If you are not using a custom single-cell reference, only confirm input, output, and threads.",
                        "If you use a custom single-cell reference, confirm sc_dat, cell_state_labels, cell_type_labels, and key together.",
                    ],
                )
            )

        if "clustering" not in completed:
            suggestions.append(
                _missing_downstream_suggestion(
                    "run_tme_cluster",
                    "tme_cluster",
                    label="Run TME subtype / clustering analysis",
                    label_zh="运行 TME 亚型 / 聚类分析",
                    reason="A TPM-ready downstream context exists, but no TME clustering result was detected.",
                    reason_zh="已进入可做下游分析的状态，但未检测到 TME 聚类/亚型结果。",
                    missing_analysis_ids=["clustering"],
                    missing_analysis_labels=["TME subtype / clustering analysis"],
                    missing_analysis_labels_zh=["TME 亚型 / 聚类分析"],
                    notes_zh=["聚类输入通常应选择已有的 deconvo_merged.csv、CIBERSORT 结果或 signature score 矩阵。"],
                    notes=[
                        "The clustering input is usually an existing deconvo_merged.csv file, a CIBERSORT result, or a signature score matrix."
                    ],
                )
            )

    if "hla_typing" not in completed:
        suggestions.append(
            _missing_downstream_suggestion(
                "run_hla_typing",
                "hla_typing",
                label="Run HLA typing",
                label_zh="运行 HLA 分型",
                reason="No merged iobrpy HLA typing result was detected.",
                reason_zh="未检测到 iobrpy HLA 分型合并结果。",
                missing_analysis_ids=["hla_typing"],
                missing_analysis_labels=["HLA typing"],
                missing_analysis_labels_zh=["HLA 分型"],
                notes_zh=["需要确认 BAM 输入目录、参考基因组版本 hg19/hg38、输出目录和线程数。"],
                notes=["Confirm the BAM input directory, hg19/hg38 reference genome version, output directory, and thread count."],
            )
        )

    if "trust4" not in completed:
        suggestions.append(
            _missing_downstream_suggestion(
                "run_trust4",
                "trust4",
                label="Run iobrpy TRUST4 TCR/BCR repertoire analysis",
                label_zh="运行 iobrpy TRUST4 TCR/BCR 库分析",
                reason="No iobrpy TRUST4 repertoire result was detected.",
                reason_zh="未检测到 iobrpy TRUST4 TCR/BCR 库分析结果。",
                missing_analysis_ids=["trust4"],
                missing_analysis_labels=["TCR/BCR repertoire analysis"],
                missing_analysis_labels_zh=["TCR/BCR 库分析"],
                notes_zh=["TRUST4 输入模式四选一：BAM、R1/R2、单端 reads，或 FASTQ 目录。"],
                notes=["Choose exactly one TRUST4 input mode: BAM, R1/R2, single-end reads, or a FASTQ directory."],
            )
        )

    return suggestions


def _extend_suggested_commands(
    base_commands: List[str],
    missing_downstream_suggestions: List[Dict[str, Any]],
) -> List[str]:
    commands = list(base_commands)
    for suggestion in missing_downstream_suggestions:
        for command in suggestion.get("suggested_commands", []):
            if command not in commands:
                commands.append(command)
    return commands


def _enrich_scenario_with_command_details(scenario: Dict[str, Any]) -> Dict[str, Any]:
    enriched = dict(scenario)
    choice_details: List[Dict[str, Any]] = []
    for item in enriched.get("choice_details", []):
        if not isinstance(item, dict):
            continue
        updated = dict(item)
        suggested_commands = updated.get("suggested_commands", [])
        if suggested_commands:
            updated["suggested_command_details"] = _suggested_command_details(suggested_commands)
        choice_details.append(updated)
    enriched["choice_details"] = choice_details
    return enriched


def _infer_decision_prompt(completed: Set[str], next_stage_ids: List[str], current_label: str) -> str:
    if "tpm_matrix" in completed:
        if any(stage in completed for stage in ("signature_scoring", "deconvolution", "ligand_receptor", "clustering")):
            return (
                f"Detected existing downstream TME outputs around {current_label}. "
                "Continue from the current results by filling in missing analyses unless you explicitly want a rerun."
            )
        return (
            "Detected a TPM-ready state. Continue with downstream TME analysis rather than rerunning FASTQ quantification "
            "unless you need a different quantification mode or reference."
        )
    if next_stage_ids and any(stage in completed for stage in ("fastq_qc", "salmon_quant", "salmon_merge", "star_quant", "star_merge")):
        return (
            f"Detected intermediate analysis outputs up to {current_label}. "
            "Choose whether to continue from the current stage, rerun this stage, or rerun the full pipeline in a fresh output directory."
        )
    if "raw_fastq" in completed:
        return (
            "Detected raw FASTQ input. Start from runall, or choose a side workflow such as TRUST4 if you only want immune repertoire outputs."
        )
    return (
        "The directory could not be placed on the IOBRpy roadmap confidently. Use recommend plus commands/help to choose the next workflow."
    )


def _infer_resume_recommendation(completed: Set[str]) -> bool:
    return any(
        stage in completed
        for stage in (
            "fastq_qc",
            "salmon_quant",
            "salmon_merge",
            "prepare_salmon",
            "star_quant",
            "star_merge",
            "count2tpm",
            "tpm_matrix",
            "signature_scoring",
            "deconvolution",
            "ligand_receptor",
            "clustering",
        )
    )


def _roadmap_links(rules: Dict[str, Any]) -> Dict[str, Optional[str]]:
    roadmap = rules.get("roadmap", {})
    local_doc = roadmap.get("local_doc")
    local_path: Optional[str] = None
    if local_doc:
        candidate = Path(__file__).resolve().parents[3] / local_doc
        if candidate.exists():
            local_path = str(candidate)
    return {
        "title": roadmap.get("title", "IOBRpy workflow"),
        "url": roadmap.get("url"),
        "local_doc": local_path,
    }

def _agent_fallback_targets(
    workflow_checklist: List[Dict[str, Any]],
    function_detections: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    detection_by_id = {item["id"]: item for item in function_detections}
    targets: List[Dict[str, Any]] = []

    for item in workflow_checklist:
        if item.get("status") == "completed":
            continue

        related_specs = [
            spec
            for spec in _FUNCTION_DETECTION_SPECS
            if set(spec.get("stage_ids", [])).intersection(item.get("matched_stage_ids", []) + item.get("pending_stage_ids", []))
        ]
        related_function_ids = [spec["id"] for spec in related_specs]
        related_function_labels = [spec["label"] for spec in related_specs]
        confirmed_function_ids = [
            function_id
            for function_id in related_function_ids
            if detection_by_id.get(function_id, {}).get("status") == "confirmed_by_content"
        ]
        likely_function_ids = [
            function_id
            for function_id in related_function_ids
            if detection_by_id.get(function_id, {}).get("status") in {"likely_iobrpy_result", "reusable_result"}
        ]
        filename_patterns = _dedupe_terms(
            pattern
            for spec in related_specs
            for pattern in spec.get("patterns", [])
        )
        search_terms = _dedupe_terms(
            list(_AGENT_FALLBACK_TERM_OVERRIDES.get(item["id"], []))
            + item.get("matched_stage_ids", [])
            + item.get("pending_stage_ids", [])
            + related_function_ids
            + related_function_labels
            + _agent_term_tokens(item.get("label", ""))
        )
        bucket_hint = "external_tool_results" if item.get("status") == "external" else "agent_inferred_existing_results"
        targets.append(
            {
                "id": item["id"],
                "label": item.get("label"),
                "label_zh": item.get("label_zh", item.get("label")),
                "status": item.get("status"),
                "bucket_hint": bucket_hint,
                "search_terms": search_terms,
                "filename_patterns": filename_patterns,
                "related_stage_ids": _dedupe_matches(item.get("matched_stage_ids", []) + item.get("pending_stage_ids", [])),
                "related_function_ids": related_function_ids,
                "related_function_labels": related_function_labels,
                "confirmed_related_function_ids": confirmed_function_ids,
                "likely_related_function_ids": likely_function_ids,
            }
        )

    return targets


def _agent_fallback_payload(
    workflow_checklist: List[Dict[str, Any]],
    function_detections: List[Dict[str, Any]],
    external_analysis_hints: List[Dict[str, Any]],
    scan_meta: Dict[str, Any],
    likely_iobrpy_functions: List[str],
    reusable_result_functions: List[str],
) -> Dict[str, Any]:
    reason_ids: List[str] = []
    if scan_meta.get("scan_retry_recommended"):
        reason_ids.append("scan_limited")
    if external_analysis_hints:
        reason_ids.append("external_results_present")
    if likely_iobrpy_functions:
        reason_ids.append("likely_iobrpy_results_need_review")
    if reusable_result_functions:
        reason_ids.append("reusable_results_need_review")
    if any(item.get("status") in {"partial", "external"} for item in workflow_checklist):
        reason_ids.append("ambiguous_checklist_sections")

    recommended = bool(reason_ids)
    summary = (
        "Use agent fallback investigation when map still appears to miss existing outputs, especially after the automatic scan retry. "
        "Search targeted subdirectories, inspect a few candidate file headers, and separate findings into confirmed IOBRpy results, agent-inferred existing results, and external-tool results."
    )
    summary_zh = (
        "当 map 在自动重扫后仍然看起来漏掉已有结果时，使用 agent 的兜底调查模式。"
        "对目标子目录做定向搜索，抽样查看候选文件头部，并把发现分成 IOBRpy 已确认结果、Agent 推断的已有结果、外部工具结果 三类。"
    )

    return {
        "recommended": recommended,
        "reason_ids": reason_ids,
        "summary": summary,
        "summary_zh": summary_zh,
        "trigger_when": [
            "map output still conflicts with the user's statement about existing TPM, HLA, TCR/BCR, or tme_profile-like outputs",
            "scan_warning remains after the automatic retry",
            "the directory contains external or likely-iobrpy-like results that are not yet confirmed natively",
        ],
        "trigger_when_zh": [
            "map 结果仍然和用户声明的已有 TPM、HLA、TCR/BCR 或 tme_profile 类输出相冲突",
            "自动重扫后仍然存在 scan_warning",
            "目录中存在外部结果或疑似 IOBRpy 结果，但尚未被原生规则确认",
        ],
        "actions": [
            "Search the path for target-specific search_terms and filename_patterns instead of relying only on the checklist negative result",
            "Preview a few candidate files to decide whether they are native IOBRpy outputs, nonstandard existing results, or external-tool results",
            "Report findings under the three result buckets without mixing inferred or external results into the confirmed IOBRpy bucket",
        ],
        "actions_zh": [
            "不要只依赖 checklist 的未完成结论，而是按 investigation_targets 里的 search_terms 和 filename_patterns 去定向搜索",
            "抽样预览几个候选文件，判断它们属于 IOBRpy 原生输出、非标准已有结果，还是外部工具结果",
            "汇报时严格使用三类结果桶，不要把推断结果或外部结果混进 IOBRpy 已确认结果里",
        ],
        "result_buckets": list(_AGENT_FALLBACK_RESULT_BUCKETS),
        "investigation_targets": _agent_fallback_targets(workflow_checklist, function_detections),
    }


def _auto_investigation_reason_ids(payload: Dict[str, Any]) -> List[str]:
    fallback_reason_ids = list(payload.get("agent_fallback", {}).get("reason_ids", []))
    current_confidence = str(payload.get("current_stage_confidence", "none"))
    external_hints_present = bool(payload.get("external_analysis_hints"))
    auto_reason_ids: List[str] = []
    for reason_id in fallback_reason_ids:
        if reason_id in {"scan_limited", "likely_iobrpy_results_need_review", "reusable_results_need_review"}:
            auto_reason_ids.append(reason_id)
        elif (
            reason_id == "ambiguous_checklist_sections"
            and current_confidence in {"low", "none"}
            and not external_hints_present
        ):
            auto_reason_ids.append(reason_id)
    return _dedupe_terms(auto_reason_ids)


def _function_result_source_labels(source_id: str) -> Tuple[str, str]:
    labels = {
        "iobrpy_confirmed_results": ("IOBRpy-confirmed", "IOBRpy 已确认"),
        "agent_inferred_existing_results": ("Existing results inferred", "已推断存在结果"),
        "external_tool_results": ("External-tool results", "外部工具结果"),
        "pending": ("Not yet confirmed", "尚未确认"),
    }
    return labels.get(source_id, (source_id.replace("_", " ").title(), source_id))


def _empty_function_detection(spec: Dict[str, Any]) -> Dict[str, Any]:
    return _function_detection_payload(
        spec,
        status="not_detected",
        evidence=[],
        matched_by_name=[],
        matched_by_content=[],
    )


def _function_investigation_records(
    investigation: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    records: Dict[str, Dict[str, Any]] = {}
    known_ids = {spec["id"] for spec in _FUNCTION_DETECTION_SPECS}
    for finding in investigation.get("findings", []):
        relpath = finding.get("relpath", "")
        bucket_id = finding.get("bucket_id", "agent_inferred_existing_results")
        matched_ids = {
            function_id
            for function_id in finding.get("matched_function_ids", [])
            if function_id in known_ids
        }
        confirmed_ids = {
            function_id
            for function_id in finding.get("confirmed_related_function_ids", [])
            if function_id in known_ids
        }
        for function_id in matched_ids | confirmed_ids:
            record = records.setdefault(
                function_id,
                {
                    "finding_count": 0,
                    "bucket_counts": {},
                    "bucket_paths": {
                        "iobrpy_confirmed_results": [],
                        "agent_inferred_existing_results": [],
                        "external_tool_results": [],
                    },
                    "confirmed_paths": [],
                },
            )
            record["finding_count"] += 1
            record["bucket_counts"][bucket_id] = record["bucket_counts"].get(bucket_id, 0) + 1
            if relpath and relpath not in record["bucket_paths"].setdefault(bucket_id, []):
                record["bucket_paths"][bucket_id].append(relpath)
            if function_id in confirmed_ids and relpath and relpath not in record["confirmed_paths"]:
                record["confirmed_paths"].append(relpath)

    for record in records.values():
        for bucket_id, paths in list(record.get("bucket_paths", {}).items()):
            record["bucket_paths"][bucket_id] = _dedupe_matches(paths, limit=max(len(paths), 1))
        confirmed_paths = record.get("confirmed_paths", [])
        record["confirmed_paths"] = _dedupe_matches(confirmed_paths, limit=max(len(confirmed_paths), 1))
    return records


def _function_result_source_id(
    status: str,
    investigation_record: Optional[Dict[str, Any]] = None,
) -> str:
    bucket_counts = (investigation_record or {}).get("bucket_counts", {})
    if status == "confirmed_by_content" or bucket_counts.get("iobrpy_confirmed_results"):
        return "iobrpy_confirmed_results"
    if status in {"likely_iobrpy_result", "reusable_result"} or bucket_counts.get("agent_inferred_existing_results"):
        return "agent_inferred_existing_results"
    if bucket_counts.get("external_tool_results"):
        return "external_tool_results"
    return "pending"


def _build_function_detection_from_sources(
    spec: Dict[str, Any],
    *sources: Dict[str, Any],
    investigation_record: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    source_items = [item for item in sources if item]
    bucket_paths = (investigation_record or {}).get("bucket_paths", {})
    inferred_paths = list(bucket_paths.get("agent_inferred_existing_results", []))
    confirmed_paths = list((investigation_record or {}).get("confirmed_paths", []))
    external_paths = list(bucket_paths.get("external_tool_results", []))
    bucket_counts = dict((investigation_record or {}).get("bucket_counts", {}))
    finding_count = int((investigation_record or {}).get("finding_count", 0))
    strict_iobrpy_only = bool(set(spec.get("stage_ids", [])).intersection(_STRICT_IOBRPY_INVESTIGATION_STAGE_IDS))

    if strict_iobrpy_only:
        inferred_paths = []
        external_paths = []
        if not confirmed_paths:
            bucket_counts.pop("agent_inferred_existing_results", None)
            bucket_counts.pop("external_tool_results", None)

    inferred_status = "not_detected"
    if confirmed_paths:
        inferred_status = "confirmed_by_content"
    elif inferred_paths:
        inferred_status = "reusable_result" if spec.get("policy") == "reusable" else "likely_iobrpy_result"

    status = _merge_function_status(
        *(item.get("status", "not_detected") for item in source_items),
        inferred_status,
    )
    matched_by_name = _dedupe_matches(
        match
        for item in source_items
        for match in item.get("matched_by_name", [])
    )
    matched_by_content = _dedupe_matches(
        match
        for item in source_items
        for match in item.get("matched_by_content", [])
    )
    for relpath in inferred_paths:
        if relpath not in matched_by_name:
            matched_by_name.append(relpath)
    for relpath in confirmed_paths:
        if relpath not in matched_by_content:
            matched_by_content.append(relpath)

    evidence = _dedupe_matches(
        match
        for item in source_items
        for match in item.get("evidence", [])
    )
    for relpath in confirmed_paths + inferred_paths:
        if relpath not in evidence:
            evidence.append(relpath)
    evidence = _dedupe_matches(evidence, limit=max(len(evidence), 1))
    matched_by_name = _dedupe_matches(matched_by_name, limit=max(len(matched_by_name), 1))
    matched_by_content = _dedupe_matches(matched_by_content, limit=max(len(matched_by_content), 1))

    payload = _function_detection_payload(
        spec,
        status=status,
        evidence=evidence,
        matched_by_name=matched_by_name,
        matched_by_content=matched_by_content,
    )
    source_id = _function_result_source_id(status, investigation_record)
    source_label, source_label_zh = _function_result_source_labels(source_id)
    investigation_evidence = confirmed_paths + inferred_paths + external_paths
    payload.update(
        {
            "result_source_id": source_id,
            "result_source_label": source_label,
            "result_source_label_zh": source_label_zh,
            "investigation_applied": bool(investigation_record),
            "investigation_finding_count": finding_count,
            "investigation_bucket_counts": bucket_counts,
            "investigation_evidence": _dedupe_matches(investigation_evidence, limit=max(len(investigation_evidence), 1)),
            "investigation_confirmed_evidence": _dedupe_matches(confirmed_paths, limit=max(len(confirmed_paths), 1)),
            "investigation_inferred_evidence": _dedupe_matches(inferred_paths, limit=max(len(inferred_paths), 1)),
            "investigation_external_evidence": _dedupe_matches(external_paths, limit=max(len(external_paths), 1)),
        }
    )
    return payload


def _wrapper_entry_hints_from_detections(function_detections: List[Dict[str, Any]]) -> List[str]:
    entries: Set[str] = set()
    for detection in function_detections:
        for relpath in (
            detection.get("evidence", [])
            + detection.get("matched_by_name", [])
            + detection.get("matched_by_content", [])
            + detection.get("investigation_evidence", [])
        ):
            if not relpath:
                continue
            normalized = _normalize_relpath(Path(relpath))
            if normalized == ".":
                continue
            entries.add(normalized)
            parts = Path(normalized).parts
            for depth in range(1, len(parts)):
                entries.add(_normalize_relpath(Path(*parts[:depth])))
    return sorted(entries)


def _function_has_confirmed_native_evidence(detection: Dict[str, Any]) -> bool:
    return bool(
        detection.get("matched_by_content")
        or detection.get("investigation_confirmed_evidence")
    )


def _deconvolution_methods_from_function_detections(function_detections: List[Dict[str, Any]]) -> Dict[str, Any]:
    raw_detection_by_id = {item["id"]: item for item in function_detections}
    detection_by_id = {
        _canonical_deconvolution_method_id(item["id"]): item
        for item in function_detections
    }
    default_specs = [
        spec for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET
    ]
    optional_specs = [
        spec for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] not in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET
    ]
    detected_ids: List[str] = []
    detected_labels: List[str] = []
    optional_detected_ids: List[str] = []
    optional_detected_labels: List[str] = []
    evidence: Dict[str, List[str]] = {}
    for spec in _DECONVOLUTION_METHOD_SPECS:
        detection = detection_by_id.get(spec["id"], {})
        if detection.get("status") in {None, "", "not_detected"}:
            continue
        if not _function_has_confirmed_native_evidence(detection):
            continue
        if spec["id"] in _WRAPPER_DEFAULT_DECONVOLUTION_ID_SET:
            detected_ids.append(spec["id"])
            detected_labels.append(spec["label"])
        else:
            optional_detected_ids.append(spec["id"])
            optional_detected_labels.append(spec["label"])
        evidence[spec["id"]] = _dedupe_matches(
            detection.get("matched_by_content", []) + detection.get("investigation_confirmed_evidence", []),
            limit=max(
                len(detection.get("matched_by_content", [])) + len(detection.get("investigation_confirmed_evidence", [])),
                1,
            ),
        )

    default_detected_ids = [spec["id"] for spec in default_specs if spec["id"] in detected_ids]
    default_detected_labels = [spec["label"] for spec in default_specs if spec["id"] in detected_ids]
    default_missing_ids = [spec["id"] for spec in default_specs if spec["id"] not in detected_ids]
    default_missing_labels = [spec["label"] for spec in default_specs if spec["id"] not in detected_ids]
    all_detected_ids = [spec["id"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] in set(detected_ids + optional_detected_ids)]
    all_detected_labels = [spec["label"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] in set(all_detected_ids)]
    all_missing_ids = [spec["id"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] not in set(all_detected_ids)]
    all_missing_labels = [spec["label"] for spec in _DECONVOLUTION_METHOD_SPECS if spec["id"] not in set(all_detected_ids)]

    wrapper_evidence = _dedupe_matches(
        raw_detection_by_id.get("tme_profile", {}).get("evidence", [])
        + raw_detection_by_id.get("runall", {}).get("evidence", []),
        limit=12,
    )
    wrapper_statuses = {
        raw_detection_by_id.get("tme_profile", {}).get("status"),
        raw_detection_by_id.get("runall", {}).get("status"),
    }
    wrapper_bundle_confirmed = bool(wrapper_evidence and "confirmed_by_content" in wrapper_statuses)
    return {
        "total_supported": len(_DECONVOLUTION_METHOD_SPECS),
        "total_available": len(_DECONVOLUTION_METHOD_SPECS),
        "default_bundle_total": len(default_specs),
        "detected_ids": all_detected_ids,
        "detected_labels": all_detected_labels,
        "detected_count": len(all_detected_ids),
        "missing_ids": all_missing_ids,
        "missing_labels": all_missing_labels,
        "default_bundle_ids": list(_WRAPPER_DEFAULT_DECONVOLUTION_IDS),
        "default_bundle_detected_ids": default_detected_ids,
        "default_bundle_detected_labels": default_detected_labels,
        "default_bundle_detected_count": len(default_detected_ids),
        "default_bundle_missing_ids": default_missing_ids,
        "default_bundle_missing_labels": default_missing_labels,
        "optional_ids": [spec["id"] for spec in optional_specs],
        "optional_labels": [spec["label"] for spec in optional_specs],
        "optional_detected_ids": optional_detected_ids,
        "optional_detected_labels": optional_detected_labels,
        "optional_missing_ids": [spec["id"] for spec in optional_specs if spec["id"] not in optional_detected_ids],
        "optional_missing_labels": [
            spec["label"] for spec in optional_specs if spec["id"] not in optional_detected_ids
        ],
        "evidence": evidence,
        "wrapper_bundle_confirmed": wrapper_bundle_confirmed,
        "wrapper_bundle_evidence": wrapper_evidence,
        "inferred_from_wrapper_ids": [],
    }


def _refresh_payload_from_function_detections(payload: Dict[str, Any]) -> Dict[str, Any]:
    if not payload.get("success"):
        return payload

    rules = load_pipeline_map_rules()
    stage_index = {stage["id"]: stage for stage in rules.get("stages", [])}
    function_detections = list(payload.get("function_detections", []))
    strong_detected = set(payload.get("strong_detected_stages", []))
    supported_weak = set(payload.get("supported_weak_stages", []))
    weak_detected = set(payload.get("weak_hint_stages", [])) | supported_weak
    strong_evidence = {
        stage_id: list(values.get("strong", []))
        for stage_id, values in payload.get("evidence", {}).items()
    }
    weak_evidence = {
        stage_id: list(values.get("weak", []))
        for stage_id, values in payload.get("evidence", {}).items()
    }
    content_verified_stage_ids = set(payload.get("content_verified_stages", []))

    for detection in function_detections:
        if detection.get("status") == "confirmed_by_content":
            for stage_id in detection.get("stage_ids", []):
                strong_detected.add(stage_id)
                content_verified_stage_ids.add(stage_id)
                strong_evidence.setdefault(stage_id, [])
                for evidence in detection.get("evidence", []):
                    if evidence not in strong_evidence[stage_id]:
                        strong_evidence[stage_id].append(evidence)
        elif detection.get("status") in {"likely_iobrpy_result", "reusable_result"}:
            for stage_id in detection.get("stage_ids", []):
                weak_detected.add(stage_id)
                weak_evidence.setdefault(stage_id, [])
                for evidence in detection.get("evidence", []):
                    if evidence not in weak_evidence[stage_id]:
                        weak_evidence[stage_id].append(evidence)

    supported_weak = _resolve_supported_weak_stages(strong_detected, weak_detected, stage_index)
    completed = _closure_with_implies(strong_detected | supported_weak, stage_index)
    ordered_completed = _ordered_stage_ids(completed)
    current_stage = _find_current_stage(completed, stage_index)
    current_label = _localized_copy_value(_stage_copy(current_stage, stage_index), "label") if current_stage else "Unclassified"
    current_label_zh = _localized_copy_value(_stage_copy(current_stage, stage_index), "label", "zh") if current_stage else "未分类"
    next_stage_ids = _infer_next_stage_ids(completed)
    deconvolution_methods = _deconvolution_methods_from_function_detections(function_detections)
    missing_downstream_analysis_suggestions = _missing_downstream_analysis_suggestions(completed, deconvolution_methods)
    suggested_commands = _extend_suggested_commands(
        _infer_suggested_commands(completed, next_stage_ids),
        missing_downstream_analysis_suggestions,
    )
    suggested_command_details = _suggested_command_details(suggested_commands)
    current_confidence = _stage_confidence(current_stage, strong_detected, supported_weak, completed) if current_stage else "none"
    supplementary_stages = _supplementary_stage_ids(completed, stage_index, current_stage)
    choices = _RESUME_CHOICES if _infer_resume_recommendation(completed) else ["start_pipeline", "inspect_commands"]
    scenario_id = _classify_scenario(completed, strong_detected | supported_weak, current_stage)
    scenario = _scenario_payload(
        scenario_id,
        completed_labels=_ordered_labels(completed, stage_index),
        next_labels=_ordered_labels(next_stage_ids, stage_index),
        choices=choices,
        current_stage=current_stage,
        suggested_commands=suggested_commands,
        scenario_labels=_SCENARIO_LABELS,
        choice_labels=_CHOICE_LABELS,
        next_commands=_NEXT_COMMANDS,
    )
    scenario = _enrich_scenario_with_command_details(scenario)
    recommended_action = _recommended_action(
        scenario,
        scenario.get("choice_details", []),
        choice_labels=_CHOICE_LABELS,
    )
    workflow_checklist = _workflow_checklist(
        completed,
        strong_evidence,
        payload.get("external_analysis_hints", []),
        deconvolution_methods,
        stage_index,
    )
    detected_functions = _dedupe_matches(
        [item["id"] for item in function_detections if item.get("status") != "not_detected"],
        limit=len(_FUNCTION_DETECTION_SPECS),
    )
    reusable_result_functions = _dedupe_matches(
        [item["id"] for item in function_detections if item.get("status") == "reusable_result"],
        limit=len(_FUNCTION_DETECTION_SPECS),
    )
    likely_iobrpy_functions = _dedupe_matches(
        [item["id"] for item in function_detections if item.get("status") == "likely_iobrpy_result"],
        limit=len(_FUNCTION_DETECTION_SPECS),
    )
    content_verified_functions = _dedupe_matches(
        [item["id"] for item in function_detections if item.get("status") == "confirmed_by_content"],
        limit=len(_FUNCTION_DETECTION_SPECS),
    )
    roadmap = _roadmap_links(rules)
    roadmap_progress = _roadmap_progress_text(
        completed,
        current_stage,
        next_stage_ids,
        ordered_stage_ids=_ordered_stage_ids,
        stage_to_node=_STAGE_TO_NODE,
    )
    roadmap_payload = {
        **roadmap,
        "current_node": _STAGE_TO_NODE.get(current_stage, "Unclassified"),
        "next_nodes": _next_nodes(next_stage_ids, _STAGE_TO_NODE),
        "completed_nodes_text": roadmap_progress["completed_text"],
        "next_nodes_text": roadmap_progress["next_text"],
        "position_summary": roadmap_progress["sentence"],
    }
    stage_rows = _stage_completion_rows(
        rules,
        stage_index,
        strong_detected,
        supported_weak,
        completed,
        strong_evidence,
        weak_evidence,
    )
    evidence = {
        stage_id: {
            "strong": _dedupe_matches(strong_evidence.get(stage_id, []), limit=max(len(strong_evidence.get(stage_id, [])), 1)),
            "weak": _dedupe_matches(weak_evidence.get(stage_id, []), limit=max(len(weak_evidence.get(stage_id, [])), 1)),
        }
        for stage_id in _ordered_stage_ids(strong_detected | weak_detected)
    }
    scan_meta = {
        "scan_retry_recommended": payload.get("scan_retry_recommended"),
        "scan_warning": payload.get("scan_warning"),
        "scan_retry_command": payload.get("scan_retry_command"),
        "scan_retry_limits": payload.get("scan_retry_limits", {}),
    }
    agent_fallback = _agent_fallback_payload(
        workflow_checklist,
        function_detections,
        payload.get("external_analysis_hints", []),
        scan_meta,
        likely_iobrpy_functions,
        reusable_result_functions,
    )

    payload = dict(payload)
    payload.update(
        {
            "strong_detected_stages": _ordered_stage_ids(strong_detected),
            "supported_weak_stages": _ordered_stage_ids(supported_weak),
            "weak_hint_stages": _ordered_stage_ids(weak_detected - supported_weak - strong_detected),
            "detected_stages": _ordered_stage_ids(strong_detected | supported_weak),
            "completed_stages": ordered_completed,
            "completed_stage_labels": _ordered_labels(ordered_completed, stage_index),
            "completed_stage_labels_zh": _ordered_labels(ordered_completed, stage_index, language="zh"),
            "completed_stage_summaries": _stage_summary_list(ordered_completed, stage_index, completed=True),
            "current_stage": current_stage,
            "current_label": current_label,
            "current_label_zh": current_label_zh,
            "current_stage_confidence": current_confidence,
            "supplementary_stages": supplementary_stages,
            "supplementary_stage_labels": _ordered_labels(supplementary_stages, stage_index),
            "supplementary_stage_labels_zh": _ordered_labels(supplementary_stages, stage_index, language="zh"),
            "supplementary_stage_summaries": _stage_summary_list(supplementary_stages, stage_index, completed=True),
            "next_stages": next_stage_ids,
            "next_stage_labels": _ordered_labels(next_stage_ids, stage_index),
            "next_stage_labels_zh": _ordered_labels(next_stage_ids, stage_index, language="zh"),
            "next_step_summaries": _stage_summary_list(next_stage_ids, stage_index, completed=False),
            "resume_candidates": next_stage_ids,
            "rerun_candidates": [current_stage, "runall"] if current_stage else ["runall"],
            "resume_recommended": _infer_resume_recommendation(completed),
            "rerun_recommended": not _infer_resume_recommendation(completed),
            "decision_prompt": _infer_decision_prompt(completed, next_stage_ids, current_label),
            "choices": choices,
            "suggested_commands": suggested_commands,
            "suggested_command_details": suggested_command_details,
            "missing_downstream_analysis_suggestions": missing_downstream_analysis_suggestions,
            "scenario": scenario,
            "recommended_action": recommended_action,
            "deconvolution_methods": deconvolution_methods,
            "content_verified_stages": _ordered_stage_ids(content_verified_stage_ids),
            "detected_functions": detected_functions,
            "reusable_result_functions": reusable_result_functions,
            "likely_iobrpy_functions": likely_iobrpy_functions,
            "content_verified_functions": content_verified_functions,
            "agent_fallback": agent_fallback,
            "agent_rendering_hints": _agent_rendering_hints_payload(payload),
            "workflow_checklist": workflow_checklist,
            "evidence": evidence,
            "roadmap": roadmap_payload,
            "stage_completion": stage_rows,
            "terminal_stage_map": _terminal_stage_map(workflow_checklist),
            "terminal_stage_map_zh": _terminal_stage_map(workflow_checklist, language="zh"),
        }
    )
    return payload


def _enrich_function_detections_with_investigation(
    payload: Dict[str, Any],
    investigation: Dict[str, Any],
) -> Dict[str, Any]:
    if not payload.get("success"):
        return payload

    existing_by_id = {
        item["id"]: dict(item)
        for item in payload.get("function_detections", [])
    }
    investigation_records = _function_investigation_records(investigation)
    updated_by_id: Dict[str, Dict[str, Any]] = {}

    for spec in _FUNCTION_DETECTION_SPECS:
        existing = existing_by_id.get(spec["id"], _empty_function_detection(spec))
        updated_by_id[spec["id"]] = _build_function_detection_from_sources(
            spec,
            existing,
            investigation_record=investigation_records.get(spec["id"]),
        )

    wrapper_entries = _wrapper_entry_hints_from_detections(list(updated_by_id.values()))
    wrapper_detections = _aggregate_wrapper_detections(wrapper_entries, updated_by_id)
    spec_index = {spec["id"]: spec for spec in _FUNCTION_DETECTION_SPECS}
    for wrapper_detection in wrapper_detections:
        spec = spec_index[wrapper_detection["id"]]
        existing = updated_by_id.get(spec["id"], _empty_function_detection(spec))
        updated_by_id[spec["id"]] = _build_function_detection_from_sources(
            spec,
            existing,
            wrapper_detection,
            investigation_record=investigation_records.get(spec["id"]),
        )

    payload = dict(payload)
    payload["function_detections"] = [
        updated_by_id[spec["id"]]
        for spec in _FUNCTION_DETECTION_SPECS
        if spec["id"] in updated_by_id
    ]
    return _refresh_payload_from_function_detections(payload)


def _candidate_has_iobrpy_context(relpath: str, spec: Dict[str, Any]) -> bool:
    for key in ("content_context", "name_context"):
        family = spec.get(key)
        if family and _has_iobrpy_context(relpath, family):
            return True

    lower = relpath.lower()
    stage_ids = set(spec.get("stage_ids", []))
    if stage_ids.intersection({"signature_scoring", "deconvolution", "ligand_receptor", "clustering"}):
        if (
            any(
                token in lower
                for token in (
                    "04-signatures",
                    "05-tme",
                    "06-lr_cal",
                    "01-signatures",
                    "02-tme",
                    "03-lr_cal",
                    "deconvo_merged",
                )
            )
            or _path_has_part_tokens(relpath, "tme", "profile")
        ):
            return True
    if stage_ids.intersection({"salmon_quant", "salmon_merge", "prepare_salmon", "tpm_matrix"}):
        if any(token in lower for token in ("02-salmon", "03-tpm", "prepare_salmon", "merge_salmon")) or _path_has_part_token(relpath, "salmon"):
            return True
    if stage_ids.intersection({"star_quant", "star_merge", "count2tpm", "tpm_matrix"}):
        if any(token in lower for token in ("02-star", "03-tpm", "count2tpm", "merge_star_count", "readspergene")) or _path_has_part_token(relpath, "star"):
            return True
    if "fastq_qc" in stage_ids and any(token in lower for token in ("01-qc", "03-fastp", "fastp")):
        return True
    return False


def _preview_excerpt(preview: str, *, max_lines: int = 3, max_chars: int = 240) -> str:
    excerpt = " | ".join(_first_nonempty_lines(preview, limit=max_lines))
    if len(excerpt) > max_chars:
        return excerpt[: max_chars - 3] + "..."
    return excerpt


def _investigation_target_path_matches(relpath: str, target: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    lower = relpath.lower()
    basename = Path(relpath).name.lower()
    if target.get("id") == "raw_data":
        if not _looks_like_raw_fastq_path(relpath):
            return [], []
        matched_patterns = [
            pattern
            for pattern in target.get("filename_patterns", [])
            if fnmatch.fnmatch(lower, pattern.lower()) or fnmatch.fnmatch(basename, pattern.lower())
        ]
        return ["raw FASTQ-like filename"], _dedupe_terms(matched_patterns)
    matched_terms = [
        term
        for term in target.get("search_terms", [])
        if term.lower() not in _GENERIC_INVESTIGATION_TERMS and term.lower() in lower
    ]
    matched_patterns = [
        pattern
        for pattern in target.get("filename_patterns", [])
        if fnmatch.fnmatch(lower, pattern.lower()) or fnmatch.fnmatch(basename, pattern.lower())
    ]

    if not matched_terms and target.get("id") == "tpm_matrix_ready" and any(token in basename for token in ("matrix", "expr", "expression", "tpm", "fpkm", "rpkm", "abundance")):
        matched_terms.append("matrix-like filename")
    if not matched_terms and target.get("id") == "feature_scoring" and any(token in basename for token in ("signature", "score", "ssgsea", "hallmark", "reactome", "kegg")):
        matched_terms.append("score-like filename")

    return _dedupe_terms(matched_terms), _dedupe_terms(matched_patterns)


def _investigation_content_target_ids(relpath: str, preview: str) -> List[str]:
    target_ids: List[str] = []
    lower = relpath.lower()
    excluded_tpm_keywords = ("signature", "ssgsea", "score", "cibersort", "estimate", "mcpcounter", "quantiseq", "epic", "ips", "bayesprism", "ligand", "receptor", "hla", "trust4", "tcr", "bcr", "cluster", "nmf")
    tpm_like_path_keywords = ("tpm", "fpkm", "rpkm", "expression", "expr", "abundance")

    if _looks_like_fastp_report(preview):
        target_ids.append("quality_control")
    if _looks_like_salmon_quant(preview):
        target_ids.append("salmon_quant_merge")
    if _looks_like_star_reads_per_gene(preview) or (
        _looks_like_count_matrix(preview) and not any(keyword in lower for keyword in tpm_like_path_keywords)
    ):
        target_ids.append("star_quant_merge")
    if _looks_like_tpm_like_expression_matrix(preview) and not any(keyword in lower for keyword in excluded_tpm_keywords):
        target_ids.append("tpm_matrix_ready")
    if _looks_like_signature_scoring(preview):
        target_ids.append("feature_scoring")
    if any(
        checker(preview)
        for checker in (
            _looks_like_cibersort_output,
            _looks_like_ips_output,
            _looks_like_estimate_output,
            _looks_like_mcpcounter_output,
            _looks_like_quantiseq_output,
            _looks_like_epic_output,
            _looks_like_bayesprism_output,
        )
    ):
        target_ids.append("immune_deconvolution")
    if _looks_like_ligand_receptor_output(preview):
        target_ids.append("ligand_receptor_analysis")
    if _looks_like_clustering_output(preview):
        target_ids.append("clustering")
    if _looks_like_hla_merged(preview) or _looks_like_spechla_result(preview):
        target_ids.append("hla_typing_summary")
    if _looks_like_trust4_output(preview):
        target_ids.append("tcr_bcr_summary")

    return _dedupe_terms(target_ids, limit=20)


def _candidate_external_hint_ids(relpath: str, preview: str) -> List[str]:
    preview_texts = {relpath: preview} if preview else {}
    matched_hint_ids: List[str] = []
    for hint in _EXTERNAL_ANALYSIS_HINTS:
        if _collect_stage_matches([relpath], hint.get("patterns", [])):
            matched_hint_ids.append(hint["id"])
            continue
        if preview_texts and _external_hint_content_matches(preview_texts, hint["id"]):
            matched_hint_ids.append(hint["id"])
    return _dedupe_terms(matched_hint_ids)


def _build_existing_result_finding(
    relpath: str,
    preview: str,
    target: Dict[str, Any],
    matched_terms: List[str],
    matched_patterns: List[str],
    matched_by_content: bool,
) -> Dict[str, Any]:
    spec_index = {spec["id"]: spec for spec in _FUNCTION_DETECTION_SPECS}
    matched_function_ids: List[str] = []
    confirmed_related_function_ids: List[str] = []

    preview_texts = {relpath: preview} if preview else {}
    for function_id in target.get("related_function_ids", []):
        spec = spec_index.get(function_id)
        if not spec:
            continue
        name_matches = _collect_function_name_matches([relpath], spec, ignore_context=True)
        content_matches = _collect_function_content_matches(preview_texts, spec, ignore_context=True) if preview_texts else []
        if name_matches or content_matches:
            matched_function_ids.append(function_id)
        if content_matches and _candidate_has_iobrpy_context(relpath, spec):
            confirmed_related_function_ids.append(function_id)

    external_hint_ids = _candidate_external_hint_ids(relpath, preview)
    external_hint_details = [
        detail
        for detail in (
            _external_hint_profile(relpath, preview, hint_id)
            for hint_id in external_hint_ids
        )
        if detail
    ]
    if target.get("bucket_hint") == "external_tool_results" or external_hint_ids:
        bucket_id = "external_tool_results"
        rationale = "Matched external-tool-like path or content patterns outside a native IOBRpy context."
    elif confirmed_related_function_ids:
        bucket_id = "iobrpy_confirmed_results"
        rationale = f"Matched native IOBRpy-like content/context for: {', '.join(confirmed_related_function_ids)}."
    else:
        bucket_id = "agent_inferred_existing_results"
        rationale = "Matched targeted search terms, filename patterns, or content signatures, but native IOBRpy confirmation is still incomplete."

    return {
        "target_id": target["id"],
        "target_label": target.get("label"),
        "target_label_zh": target.get("label_zh", target.get("label")),
        "bucket_id": bucket_id,
        "relpath": relpath,
        "matched_terms": matched_terms,
        "matched_patterns": matched_patterns,
        "matched_by_content": matched_by_content,
        "matched_function_ids": _dedupe_terms(matched_function_ids),
        "confirmed_related_function_ids": _dedupe_terms(confirmed_related_function_ids),
        "external_hint_ids": external_hint_ids,
        "external_hint_details": external_hint_details[:2],
        "preview_excerpt": _preview_excerpt(preview) if preview else "",
        "rationale": rationale,
    }


def _target_summary_rows(
    targets: List[Dict[str, Any]],
    findings: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for target in targets:
        target_findings = [
            item
            for item in findings
            if item["target_id"] == target["id"]
            and _filter_item_specific_evidence(target["id"], [str(item.get("relpath", ""))])
        ]
        bucket_counts: Dict[str, int] = {}
        external_families: List[str] = []
        external_result_kind_labels: List[str] = []
        external_likely_tool_families: List[str] = []
        external_signal_labels: List[str] = []
        external_evidence_profiles: List[Dict[str, Any]] = []
        for item in target_findings:
            bucket_counts[item["bucket_id"]] = bucket_counts.get(item["bucket_id"], 0) + 1
            for detail in item.get("external_hint_details", []):
                if not isinstance(detail, dict):
                    continue
                if detail.get("analysis_family_label"):
                    external_families.append(str(detail.get("analysis_family_label")))
                if detail.get("result_kind_label"):
                    external_result_kind_labels.append(str(detail.get("result_kind_label")))
                external_likely_tool_families.extend(detail.get("likely_tool_families", []))
                external_signal_labels.extend(detail.get("matched_signal_labels", []))
                external_evidence_profiles.append(
                    {
                        "path": item.get("relpath"),
                        "analysis_family_label": detail.get("analysis_family_label"),
                        "confidence": detail.get("confidence"),
                        "result_kind_label": detail.get("result_kind_label"),
                        "preview_excerpt": item.get("preview_excerpt", ""),
                    }
                )
        rows.append(
            {
                "id": target["id"],
                "label": target.get("label"),
                "label_zh": target.get("label_zh", target.get("label")),
                "finding_count": len(target_findings),
                "bucket_counts": bucket_counts,
                "sample_paths": [item["relpath"] for item in target_findings[:3]],
                "external_result_families": _dedupe_terms(external_families, limit=4),
                "external_result_kind_labels": _dedupe_terms(external_result_kind_labels, limit=6),
                "external_likely_tool_families": _dedupe_terms(external_likely_tool_families, limit=6),
                "external_signal_labels": _dedupe_terms(external_signal_labels, limit=8),
                "external_evidence_profiles": external_evidence_profiles[:3],
            }
        )
    return rows


def investigate_existing_results(
    path: Path,
    payload: Dict[str, Any],
    *,
    auto_scan_retry_performed: bool = False,
    auto_triggered: bool = False,
    initial_limits: Optional[Dict[str, int]] = None,
    focus: Optional[Iterable[str]] = None,
    max_findings_per_target: int = 8,
    max_preview_files: int = 2500,
    max_files_scanned: int = 50000,
) -> Dict[str, Any]:
    fallback = payload.get("agent_fallback", {})
    targets = list(fallback.get("investigation_targets", []))
    result_buckets = list(fallback.get("result_buckets", _AGENT_FALLBACK_RESULT_BUCKETS))
    bucket_ids = [item["id"] for item in result_buckets]
    bucketed_findings: Dict[str, List[Dict[str, Any]]] = {bucket_id: [] for bucket_id in bucket_ids}
    if not payload.get("success"):
        return {
            "enabled": True,
            "recommended": False,
            "auto_scan_retry_performed": auto_scan_retry_performed,
            "result_buckets": result_buckets,
            "findings": [],
            "bucketed_findings": bucketed_findings,
            "target_summaries": [],
            "search_limits": {},
        }

    if auto_triggered:
        max_findings_per_target = min(max_findings_per_target, 5)
        max_preview_files = min(max_preview_files, 600 if (payload.get("scan", {}).get("focus_roots") or focus) else 900)
        max_files_scanned = min(max_files_scanned, 20000 if (payload.get("scan", {}).get("focus_roots") or focus) else 30000)

    search_max_depth = max(int(payload.get("scan", {}).get("max_depth", 0)), 16)
    findings: List[Dict[str, Any]] = []
    target_counts = {target["id"]: 0 for target in targets}
    finding_keys: Set[Tuple[str, str]] = set()
    files_scanned = 0
    preview_files = 0
    truncated = False

    root = path.resolve()
    focus_roots = list(payload.get("scan", {}).get("focus_roots") or list(focus or []))
    queue_entries: List[Tuple[Path, Path, int]] = []
    for focus_root in focus_roots:
        rel_dir = Path(str(focus_root))
        absolute_dir = root / rel_dir
        if absolute_dir.exists() and absolute_dir.is_dir():
            queue_entries.append((absolute_dir, rel_dir, len(rel_dir.parts)))
    queue = deque(queue_entries or [(root, Path(), 0)])
    while queue:
        current_path, rel_dir, depth = queue.popleft()
        try:
            with os.scandir(current_path) as iterator:
                entries = sorted(list(iterator), key=lambda item: item.name)
        except OSError:
            continue

        for entry in entries:
            relpath = Path(entry.name) if not rel_dir.parts else rel_dir / entry.name
            normalized_relpath = _normalize_relpath(relpath)
            if entry.is_dir(follow_symlinks=False):
                if depth < search_max_depth:
                    queue.append((current_path / entry.name, relpath, depth + 1))
                continue
            if _helper_or_nonresult_file(normalized_relpath):
                continue

            files_scanned += 1
            if files_scanned > max_files_scanned:
                truncated = True
                break

            path_matches: Dict[str, Tuple[List[str], List[str]]] = {}
            for target in targets:
                if target_counts[target["id"]] >= max_findings_per_target:
                    continue
                matched_terms, matched_patterns = _investigation_target_path_matches(normalized_relpath, target)
                if matched_terms or matched_patterns:
                    path_matches[target["id"]] = (matched_terms, matched_patterns)

            preview = ""
            content_target_ids: List[str] = []
            if _is_text_like_result(normalized_relpath) and (path_matches or preview_files < max_preview_files):
                preview = _read_text_preview(root / normalized_relpath, max_chars=64000)
                if preview:
                    preview_files += 1
                    content_target_ids = [
                        target_id
                        for target_id in _investigation_content_target_ids(normalized_relpath, preview)
                        if target_id in target_counts and target_counts[target_id] < max_findings_per_target
                    ]

            for target_id in _dedupe_terms(list(path_matches.keys()) + content_target_ids, limit=100):
                if target_counts.get(target_id, 0) >= max_findings_per_target:
                    continue
                if not _filter_item_specific_evidence(target_id, [normalized_relpath]):
                    continue
                key = (target_id, normalized_relpath)
                if key in finding_keys:
                    continue
                target = next((item for item in targets if item["id"] == target_id), None)
                if not target:
                    continue
                matched_terms, matched_patterns = path_matches.get(target_id, ([], []))
                finding = _build_existing_result_finding(
                    normalized_relpath,
                    preview,
                    target,
                    matched_terms,
                    matched_patterns,
                    target_id in content_target_ids,
                )
                bucketed_findings[finding["bucket_id"]].append(finding)
                findings.append(finding)
                target_counts[target_id] += 1
                finding_keys.add(key)

        if truncated:
            break

    target_summaries = _target_summary_rows(targets, findings)
    summary = {
        "total_findings": len(findings),
        "files_scanned": files_scanned,
        "preview_files": preview_files,
        "truncated": truncated,
        "bucket_counts": {bucket_id: len(bucketed_findings[bucket_id]) for bucket_id in bucket_ids},
    }
    return {
        "enabled": True,
        "recommended": bool(fallback.get("recommended")),
        "auto_scan_retry_performed": auto_scan_retry_performed,
        "initial_limits": initial_limits or {},
        "effective_scan_limits": {
            "max_depth": payload.get("scan", {}).get("max_depth"),
            "max_entries": payload.get("scan", {}).get("max_entries"),
        },
        "search_limits": {
            "search_max_depth": search_max_depth,
            "max_findings_per_target": max_findings_per_target,
            "max_preview_files": max_preview_files,
            "max_files_scanned": max_files_scanned,
            "focus_roots": focus_roots,
        },
        "result_buckets": result_buckets,
        "findings": findings,
        "bucketed_findings": bucketed_findings,
        "target_summaries": target_summaries,
        "summary": summary,
    }

def _ordered_labels(stage_ids: Iterable[str], stage_index: Dict[str, Dict[str, Any]], language: str = "en") -> List[str]:
    return [
        _localized_copy_value(_stage_copy(stage_id, stage_index), "label", language)
        for stage_id in _ordered_stage_ids(stage_ids)
        if stage_id in stage_index
    ]


def _aggregate_external_hint_profiles(
    hint: Dict[str, Any],
    matches: List[str],
    preview_texts: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    preview_texts = preview_texts or {}
    profiles: List[Dict[str, Any]] = []
    for relpath in matches[:8]:
        preview = preview_texts.get(relpath, "")
        profile = _external_hint_profile(
            relpath,
            preview,
            str(hint.get("id")),
            respect_native_context=True,
        )
        if profile:
            profiles.append(profile)

    inference_basis: List[str] = []
    matched_signal_ids: List[str] = []
    matched_signal_labels: List[str] = []
    likely_tool_families: List[str] = []
    result_kind_ids: List[str] = []
    result_kind_labels: List[str] = []
    confidence = "low"
    best_profile: Optional[Dict[str, Any]] = None

    for profile in profiles:
        inference_basis.extend(profile.get("inference_basis", []))
        matched_signal_ids.extend(profile.get("matched_signal_ids", []))
        matched_signal_labels.extend(profile.get("matched_signal_labels", []))
        likely_tool_families.extend(profile.get("likely_tool_families", []))
        if profile.get("result_kind_id"):
            result_kind_ids.append(str(profile.get("result_kind_id")))
        if profile.get("result_kind_label"):
            result_kind_labels.append(str(profile.get("result_kind_label")))
        profile_confidence = str(profile.get("confidence", "low"))
        if _EXTERNAL_HINT_CONFIDENCE_RANK.get(profile_confidence, 0) > _EXTERNAL_HINT_CONFIDENCE_RANK.get(confidence, 0):
            confidence = profile_confidence
            best_profile = profile
        elif best_profile is None:
            best_profile = profile

    if not best_profile and matches:
        generic_profile = _external_hint_profile(
            matches[0],
            "",
            str(hint.get("id")),
            respect_native_context=False,
        )
        if generic_profile:
            best_profile = generic_profile
            confidence = str(generic_profile.get("confidence", confidence))
            inference_basis.extend(generic_profile.get("inference_basis", []))
            likely_tool_families.extend(generic_profile.get("likely_tool_families", []))
            if generic_profile.get("result_kind_id"):
                result_kind_ids.append(str(generic_profile.get("result_kind_id")))
            if generic_profile.get("result_kind_label"):
                result_kind_labels.append(str(generic_profile.get("result_kind_label")))

    evidence_profiles = [
        {
            "path": item.get("path"),
            "confidence": item.get("confidence"),
            "result_kind_id": item.get("result_kind_id"),
            "result_kind_label": item.get("result_kind_label"),
            "matched_signal_labels": _dedupe_terms(item.get("matched_signal_labels", []), limit=4),
            "likely_tool_families": _dedupe_terms(item.get("likely_tool_families", []), limit=3),
            "preview_excerpt": item.get("preview_excerpt", ""),
        }
        for item in profiles[:4]
    ]

    return {
        "analysis_family_id": (best_profile or {}).get("analysis_family_id"),
        "analysis_family_label": (best_profile or {}).get("analysis_family_label"),
        "confidence": confidence,
        "inference_basis": _dedupe_terms(inference_basis, limit=8),
        "matched_signal_ids": _dedupe_terms(matched_signal_ids, limit=10),
        "matched_signal_labels": _dedupe_terms(matched_signal_labels, limit=10),
        "result_kind_ids": _dedupe_terms(result_kind_ids, limit=6),
        "result_kind_labels": _dedupe_terms(result_kind_labels, limit=6),
        "likely_tool_families": _dedupe_terms(likely_tool_families, limit=8),
        "family_summary": (
            (best_profile or {}).get("family_summary")
            or hint.get("description")
            or ""
        ),
        "evidence_profiles": evidence_profiles,
        "native_equivalence": "external_family_only",
    }


def _detect_external_analysis_hints(
    entries: Iterable[str],
    detected_stage_ids: Set[str],
    preview_texts: Optional[Dict[str, str]] = None,
    entry_matcher: Optional[Any] = None,
) -> List[Dict[str, Any]]:
    hints: List[Dict[str, Any]] = []
    skip_hla = _has_native_hla_result_evidence(entries, preview_texts)
    skip_tcr = "trust4" in detected_stage_ids

    for hint in _EXTERNAL_ANALYSIS_HINTS:
        if hint["id"] == "external_hla" and skip_hla:
            continue
        if hint["id"] == "external_tcr_bcr" and skip_tcr:
            continue
        matches = _collect_stage_matches(entries, hint["patterns"], matcher=entry_matcher)
        if hint["id"] == "external_tcr_bcr":
            matches = [
                match
                for match in matches
                if Path(match).name.lower().rstrip("/\\") not in {"07-tcrbcr", "tcrbcr", "trust4"}
            ]
        if preview_texts:
            matches = _dedupe_matches(matches + _external_hint_content_matches(preview_texts, hint["id"]))
        matches = _filter_external_hint_matches(matches, str(hint["id"]), preview_texts)
        if matches:
            aggregated = _aggregate_external_hint_profiles(hint, matches, preview_texts)
            hints.append(
                {
                    "id": hint["id"],
                    "label": hint["label"],
                    "description": hint["description"],
                    "evidence": matches,
                    "assume_iobrpy": False,
                    **aggregated,
                }
            )
    return hints

def analyze_pipeline_path(
    path: Path,
    *,
    max_depth: int = 8,
    max_entries: int = 8000,
    quick: bool = False,
    focus: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    path = path.expanduser()
    if quick:
        if max_depth == 8:
            max_depth = max(max_depth, _QUICK_MIN_MAX_DEPTH)
        if max_entries == 8000:
            max_entries = max(max_entries, _QUICK_MIN_MAX_ENTRIES)
    if not path.exists():
        return {
            "success": False,
            "path": str(path),
            "error": f"Path does not exist: {path}",
        }

    rules = load_pipeline_map_rules()
    stage_index = {stage["id"]: stage for stage in rules.get("stages", [])}
    priority_patterns: List[str] = []
    for stage in rules.get("stages", []):
        priority_patterns.extend(stage.get("detect_any", []))
    for hint in _EXTERNAL_ANALYSIS_HINTS:
        priority_patterns.extend(hint.get("patterns", []))
    for spec in _FUNCTION_DETECTION_SPECS:
        priority_patterns.extend(spec.get("patterns", []))
    scan = _scan_path_entries(
        path,
        max_depth=max_depth,
        max_entries=max_entries,
        priority_patterns=priority_patterns,
        quick=quick,
        focus=focus,
    )
    entries = scan["entries"]
    entry_matcher = _build_entry_matcher(entries)

    strong_detected: Set[str] = set()
    weak_detected: Set[str] = set()
    strong_evidence: Dict[str, List[str]] = {}
    weak_evidence: Dict[str, List[str]] = {}
    for stage in rules.get("stages", []):
        stage_id = stage["id"]
        stage_matches = _collect_stage_matches(entries, stage.get("detect_any", []), matcher=entry_matcher)
        weak_matches = _collect_stage_matches(entries, stage.get("detect_weak", []), matcher=entry_matcher)
        if stage_id == "raw_fastq":
            stage_matches = [match for match in stage_matches if _looks_like_raw_fastq_path(match)]
            weak_matches = [match for match in weak_matches if _looks_like_raw_fastq_path(match)]
        elif stage_id == "salmon_quant":
            stage_matches = [match for match in stage_matches if Path(match).name.lower() == "quant.sf"]
            weak_matches = [match for match in weak_matches if Path(match).name.lower() == "quant.sf"]
        elif stage_id == "tpm_matrix":
            stage_matches = [match for match in stage_matches if _looks_like_confident_tpm_path_match(match)]
        if stage_matches:
            strong_detected.add(stage_id)
            strong_evidence[stage_id] = stage_matches
        if weak_matches:
            weak_detected.add(stage_id)
            weak_evidence[stage_id] = weak_matches

    content_preview_limit = _content_signature_preview_limit(quick=quick, focus=focus)
    content_signatures = _detect_content_signatures(
        path.resolve(),
        entries,
        priority_patterns=priority_patterns if quick else None,
        max_preview_files=content_preview_limit,
        entry_matcher=entry_matcher,
    )
    for stage_id, matches in content_signatures.get("stage_matches", {}).items():
        strong_detected.add(stage_id)
        existing = strong_evidence.setdefault(stage_id, [])
        if stage_id in _STRICT_IOBRPY_STAGE_IDS:
            existing.clear()
        for match in matches:
            if match not in existing:
                existing.append(match)

    function_detections = content_signatures.get("function_detections", [])
    for stage_id in list(_STRICT_IOBRPY_STAGE_IDS):
        if stage_id in strong_detected and stage_id not in content_signatures.get("stage_matches", {}):
            strong_detected.discard(stage_id)
            strong_evidence.pop(stage_id, None)

    supported_weak = _resolve_supported_weak_stages(strong_detected, weak_detected, stage_index)
    completed = _closure_with_implies(strong_detected | supported_weak, stage_index)
    ordered_completed = _ordered_stage_ids(completed)
    current_stage = _find_current_stage(completed, stage_index)
    current_label = _localized_copy_value(_stage_copy(current_stage, stage_index), "label") if current_stage else "Unclassified"
    current_label_zh = _localized_copy_value(_stage_copy(current_stage, stage_index), "label", "zh") if current_stage else "未分类"
    next_stage_ids = _infer_next_stage_ids(completed)
    roadmap = _roadmap_links(rules)
    stage_rows = _stage_completion_rows(rules, stage_index, strong_detected, supported_weak, completed, strong_evidence, weak_evidence)

    if any(stage in completed for stage in ("prepare_salmon", "count2tpm")) and "tpm_matrix" not in completed:
        completed.add("tpm_matrix")
        ordered_completed = _ordered_stage_ids(completed)
        stage_rows = _stage_completion_rows(rules, stage_index, strong_detected, supported_weak, completed, strong_evidence, weak_evidence)
        current_stage = _find_current_stage(completed, stage_index)
        current_label = _localized_copy_value(_stage_copy(current_stage, stage_index), "label") if current_stage else current_label
        current_label_zh = _localized_copy_value(_stage_copy(current_stage, stage_index), "label", "zh") if current_stage else current_label_zh
        next_stage_ids = _infer_next_stage_ids(completed)

    deconvolution_methods = _deconvolution_methods_from_function_detections(function_detections)
    missing_downstream_analysis_suggestions = _missing_downstream_analysis_suggestions(completed, deconvolution_methods)
    suggested_commands = _extend_suggested_commands(
        _infer_suggested_commands(completed, next_stage_ids),
        missing_downstream_analysis_suggestions,
    )
    suggested_command_details = _suggested_command_details(suggested_commands)
    current_confidence = _stage_confidence(current_stage, strong_detected, supported_weak, completed) if current_stage else "none"
    supplementary_stages = _supplementary_stage_ids(completed, stage_index, current_stage)
    choices = _RESUME_CHOICES if _infer_resume_recommendation(completed) else ["start_pipeline", "inspect_commands"]
    scenario_id = _classify_scenario(completed, strong_detected | supported_weak, current_stage)
    roadmap_progress = _roadmap_progress_text(
        completed,
        current_stage,
        next_stage_ids,
        ordered_stage_ids=_ordered_stage_ids,
        stage_to_node=_STAGE_TO_NODE,
    )
    detected_stage_ids = strong_detected | supported_weak
    preview_texts = content_signatures.get("preview_texts", {})
    external_analysis_hints = _detect_external_analysis_hints(entries, detected_stage_ids, preview_texts, entry_matcher=entry_matcher)
    scan_limit_messages = _scan_limit_messages(path, scan)
    checklist_display_hints = _scan_checklist_display_hints(entries)
    scenario = _scenario_payload(
        scenario_id,
        completed_labels=_ordered_labels(completed, stage_index),
        next_labels=_ordered_labels(next_stage_ids, stage_index),
        choices=choices,
        current_stage=current_stage,
        suggested_commands=suggested_commands,
        scenario_labels=_SCENARIO_LABELS,
        choice_labels=_CHOICE_LABELS,
        next_commands=_NEXT_COMMANDS,
    )
    scenario = _enrich_scenario_with_command_details(scenario)
    recommended_action = _recommended_action(
        scenario,
        scenario.get("choice_details", []),
        choice_labels=_CHOICE_LABELS,
    )
    content_verified_stage_ids = _ordered_stage_ids(content_signatures.get("content_verified_stage_ids", []))
    detected_functions = _dedupe_matches(content_signatures.get("detected_function_ids", []), limit=len(_FUNCTION_DETECTION_SPECS))
    reusable_result_functions = _dedupe_matches(content_signatures.get("reusable_function_ids", []), limit=len(_FUNCTION_DETECTION_SPECS))
    likely_iobrpy_functions = _dedupe_matches(content_signatures.get("likely_iobrpy_function_ids", []), limit=len(_FUNCTION_DETECTION_SPECS))
    content_verified_functions = _dedupe_matches(content_signatures.get("content_verified_function_ids", []), limit=len(_FUNCTION_DETECTION_SPECS))
    workflow_checklist = _workflow_checklist(
        completed,
        strong_evidence,
        external_analysis_hints,
        deconvolution_methods,
        stage_index,
    )
    agent_fallback = _agent_fallback_payload(
        workflow_checklist,
        function_detections,
        external_analysis_hints,
        scan_limit_messages,
        likely_iobrpy_functions,
        reusable_result_functions,
    )
    roadmap_payload = {
        **roadmap,
        "current_node": _STAGE_TO_NODE.get(current_stage, "Unclassified"),
        "next_nodes": _next_nodes(next_stage_ids, _STAGE_TO_NODE),
        "completed_nodes_text": roadmap_progress["completed_text"],
        "next_nodes_text": roadmap_progress["next_text"],
        "position_summary": roadmap_progress["sentence"],
    }

    payload = {
        "success": True,
        "path": str(path),
        "path_kind": "file" if path.is_file() else "directory",
        "scan": {
            "entry_count": scan["entry_count"],
            "truncated": scan["truncated"],
            "max_depth": scan["max_depth"],
            "max_entries": scan["max_entries"],
            "entry_limit_enabled": scan.get("entry_limit_enabled", True),
            "entry_limit_hit": scan.get("entry_limit_hit", False),
            "backlog_limit_hit": scan.get("backlog_limit_hit", False),
            "unbounded_entry_scan": scan.get("unbounded_entry_scan", False),
            "truncated_reason_ids": scan.get("truncated_reason_ids", []),
            "depth_limited_dir_count": scan["depth_limited_dir_count"],
            "depth_limited_dirs": scan["depth_limited_dirs"],
            "strategy": scan["strategy"],
            "quick_skipped_generic_dir_count": scan["quick_skipped_generic_dir_count"],
            "sampled_bulk_stage_dir_count": scan["sampled_bulk_stage_dir_count"],
            "sampled_bulk_child_dir_count": scan["sampled_bulk_child_dir_count"],
            "sampled_bulk_file_count": scan["sampled_bulk_file_count"],
            "full_skipped_generic_dir_count": scan["full_skipped_generic_dir_count"],
            "focus_mode": scan["focus_mode"],
            "focus_roots": scan["focus_roots"],
            "explicit_focus_roots": scan["explicit_focus_roots"],
            "auto_focus_roots": scan["auto_focus_roots"],
            "promoted_focus_root_count": scan["promoted_focus_root_count"],
            "generic_deep_scan_limit": scan["generic_deep_scan_limit"],
        },
        **scan_limit_messages,
        "strong_detected_stages": _ordered_stage_ids(strong_detected),
        "supported_weak_stages": _ordered_stage_ids(supported_weak),
        "weak_hint_stages": _ordered_stage_ids(weak_detected - supported_weak - strong_detected),
        "detected_stages": _ordered_stage_ids(detected_stage_ids),
        "completed_stages": ordered_completed,
        "completed_stage_labels": _ordered_labels(ordered_completed, stage_index),
        "completed_stage_labels_zh": _ordered_labels(ordered_completed, stage_index, language="zh"),
        "completed_stage_summaries": _stage_summary_list(ordered_completed, stage_index, completed=True),
        "current_stage": current_stage,
        "current_label": current_label,
        "current_label_zh": current_label_zh,
        "current_stage_confidence": current_confidence,
        "supplementary_stages": supplementary_stages,
        "supplementary_stage_labels": _ordered_labels(supplementary_stages, stage_index),
        "supplementary_stage_labels_zh": _ordered_labels(supplementary_stages, stage_index, language="zh"),
        "supplementary_stage_summaries": _stage_summary_list(supplementary_stages, stage_index, completed=True),
        "next_stages": next_stage_ids,
        "next_stage_labels": _ordered_labels(next_stage_ids, stage_index),
        "next_stage_labels_zh": _ordered_labels(next_stage_ids, stage_index, language="zh"),
        "next_step_summaries": _stage_summary_list(next_stage_ids, stage_index, completed=False),
        "resume_candidates": next_stage_ids,
        "rerun_candidates": [current_stage, "runall"] if current_stage else ["runall"],
        "resume_recommended": _infer_resume_recommendation(completed),
        "rerun_recommended": not _infer_resume_recommendation(completed),
        "decision_prompt": _infer_decision_prompt(completed, next_stage_ids, current_label),
        "choices": choices,
        "suggested_commands": suggested_commands,
        "suggested_command_details": suggested_command_details,
        "missing_downstream_analysis_suggestions": missing_downstream_analysis_suggestions,
        "scenario": scenario,
        "recommended_action": recommended_action,
        "external_analysis_hints": external_analysis_hints,
        "deconvolution_methods": deconvolution_methods,
        "content_preview_file_count": content_signatures.get("preview_file_count", 0),
        "content_preview_file_limit": content_signatures.get("preview_file_limit"),
        "content_verified_stages": content_verified_stage_ids,
        "function_detections": function_detections,
        "detected_functions": detected_functions,
        "reusable_result_functions": reusable_result_functions,
        "likely_iobrpy_functions": likely_iobrpy_functions,
        "content_verified_functions": content_verified_functions,
        "agent_fallback": agent_fallback,
        "checklist_display_hints": checklist_display_hints,
        "workflow_checklist": workflow_checklist,
        "evidence": {
            stage_id: {
                "strong": strong_evidence.get(stage_id, []),
                "weak": weak_evidence.get(stage_id, []),
            }
            for stage_id in _ordered_stage_ids(strong_detected | weak_detected)
        },
        "roadmap": roadmap_payload,
        "stage_completion": stage_rows,
        "terminal_stage_map": _terminal_stage_map(workflow_checklist),
        "terminal_stage_map_zh": _terminal_stage_map(workflow_checklist, language="zh"),
    }
    return _enrich_workflow_checklist_payload(
        payload,
        strict_iobrpy_investigation_item_ids=_STRICT_IOBRPY_INVESTIGATION_ITEM_IDS,
    )


def _investigation_note_text(
    bucket_counts: Dict[str, int],
    sample_paths: List[str],
) -> Tuple[str, str]:
    sample_text = ", ".join(sample_paths[:3]) if sample_paths else "candidate files under the scanned path"
    sample_text_zh = "、".join(sample_paths[:3]) if sample_paths else "扫描路径下的候选文件"
    if bucket_counts.get("iobrpy_confirmed_results"):
        return (
            f"CLI-native existing-result investigation confirmed additional result files such as: {sample_text}.",
            f"CLI 原生已有结果调查进一步确认了额外结果文件，例如：{sample_text_zh}。",
        )
    if bucket_counts.get("external_tool_results"):
        return (
            f"CLI-native existing-result investigation found external-tool outputs such as: {sample_text}.",
            f"CLI 原生已有结果调查发现了外部工具输出，例如：{sample_text_zh}。",
        )
    return (
        f"CLI-native existing-result investigation found likely existing result files such as: {sample_text}.",
        f"CLI 原生已有结果调查发现了疑似已有结果文件，例如：{sample_text_zh}。",
    )


def _apply_existing_result_investigation(
    payload: Dict[str, Any],
    investigation: Dict[str, Any],
) -> Dict[str, Any]:
    target_summaries = {
        item["id"]: item
        for item in investigation.get("target_summaries", [])
        if item.get("finding_count")
    }
    if not target_summaries:
        return payload

    updated_checklist: List[Dict[str, Any]] = []
    applied = False
    for item in payload.get("workflow_checklist", []):
        summary = target_summaries.get(item.get("id"))
        if not summary:
            updated_checklist.append(item)
            continue

        bucket_counts = summary.get("bucket_counts", {})
        sample_paths = summary.get("sample_paths", [])
        updated_item = dict(item)
        updated_item["investigation_applied"] = True
        updated_item["investigation_finding_count"] = summary.get("finding_count", 0)
        updated_item["investigation_bucket_counts"] = bucket_counts
        updated_item["investigation_sample_paths"] = sample_paths
        strict_iobrpy_only = updated_item.get("id") in _STRICT_IOBRPY_INVESTIGATION_ITEM_IDS

        allows_external_status = updated_item.get("id") in _EXTERNAL_STATUS_CHECKLIST_ITEM_IDS
        if bucket_counts.get("external_tool_results") and not strict_iobrpy_only:
            updated_item["checked"] = True
            updated_item["status"] = "external" if allows_external_status else "completed"
        elif bucket_counts.get("iobrpy_confirmed_results") and updated_item.get("status") == "pending":
            updated_item["checked"] = True
            updated_item["status"] = "completed"
        elif (
            not strict_iobrpy_only
            and bucket_counts.get("agent_inferred_existing_results")
            and updated_item.get("status") == "pending"
        ):
            updated_item["checked"] = True
            updated_item["status"] = "partial"

        updated_checklist.append(updated_item)
        applied = True

    if not applied:
        return payload

    payload = dict(payload)
    payload["workflow_checklist"] = updated_checklist
    payload["terminal_stage_map"] = _terminal_stage_map(updated_checklist)
    payload["terminal_stage_map_zh"] = _terminal_stage_map(updated_checklist, language="zh")
    payload["existing_result_investigation_applied"] = True
    return payload


def map_pipeline_path(
    path: Path,
    *,
    max_depth: int = 8,
    max_entries: int = 8000,
    quick: bool = False,
    focus: Optional[Iterable[str]] = None,
    investigate_existing: bool = False,
    auto_scan_retry: bool = True,
    auto_investigate_existing: bool = True,
) -> Dict[str, Any]:
    requested_focus = list(focus or [])
    requested_limits = {
        "max_depth": max_depth,
        "max_entries": max_entries,
        "quick": quick,
        "focus": requested_focus,
    }
    planning_payload: Optional[Dict[str, Any]] = None
    planning_focus_roots: List[str] = []
    deep_scan_plan_performed = False
    if not quick and not requested_focus and path.exists() and path.is_dir():
        planning_adaptive_limits = _suggest_initial_scan_limits(
            path,
            max_depth=max_depth,
            max_entries=max_entries,
            quick=True,
            focus=None,
        )
        planning_effective_limits = planning_adaptive_limits["effective"]
        planning_payload = analyze_pipeline_path(
            path,
            max_depth=int(planning_effective_limits["max_depth"]),
            max_entries=int(planning_effective_limits["max_entries"]),
            quick=True,
            focus=None,
        )
        planning_focus_roots = _recommended_deep_focus_roots(planning_payload)
        if planning_focus_roots:
            focus = planning_focus_roots
        deep_scan_plan_performed = True
    use_unbounded_entries_for_focused_deep_scan = (
        not quick
        and bool(focus)
        and int(requested_limits.get("max_entries", 0)) == _DEFAULT_MAP_MAX_ENTRIES
    )
    adaptive_limits = _suggest_initial_scan_limits(
        path,
        max_depth=max_depth,
        max_entries=max_entries,
        quick=quick,
        focus=focus,
    )
    effective_limits = adaptive_limits["effective"]
    max_depth = int(effective_limits["max_depth"])
    max_entries = int(effective_limits["max_entries"])
    if use_unbounded_entries_for_focused_deep_scan:
        max_entries = 0
    initial_limits = {
        "max_depth": max_depth,
        "max_entries": max_entries,
        "quick": quick,
        "focus": list(focus or []),
    }
    auto_scan_retry_performed = False
    auto_scan_retry_reason = ""
    auto_investigate_reason_ids: List[str] = []
    payload = analyze_pipeline_path(path, max_depth=max_depth, max_entries=max_entries, quick=quick, focus=focus)

    if auto_scan_retry and payload.get("success") and payload.get("scan_retry_recommended"):
        retry_limits = payload.get("scan_retry_limits", {})
        retry_max_depth = int(retry_limits.get("max_depth", max_depth))
        retry_max_entries = int(retry_limits.get("max_entries", max_entries))
        if retry_max_depth != max_depth or retry_max_entries != max_entries:
            auto_scan_retry_reason = payload.get("scan_warning", "")
            payload = analyze_pipeline_path(path, max_depth=retry_max_depth, max_entries=retry_max_entries, quick=False, focus=focus)
            auto_scan_retry_performed = True

    if payload.get("success"):
        if auto_investigate_existing and payload.get("agent_fallback", {}).get("recommended"):
            auto_investigate_reason_ids = _auto_investigation_reason_ids(payload)
        should_investigate = investigate_existing or bool(auto_investigate_reason_ids)
        payload["requested_scan_limits"] = requested_limits
        payload["initial_scan_limits"] = initial_limits
        payload["effective_initial_scan_limits"] = initial_limits
        payload["adaptive_scan_limits"] = adaptive_limits
        payload["adaptive_scan_limits_applied"] = bool(adaptive_limits.get("applied"))
        payload["adaptive_scan_limit_reasons"] = list(adaptive_limits.get("reason_ids", []))
        payload["auto_scan_retry_performed"] = auto_scan_retry_performed
        payload["auto_investigate_existing_performed"] = should_investigate and not investigate_existing
        payload["auto_investigate_existing_reason_ids"] = auto_investigate_reason_ids
        payload["automatic_scan_behaviors"] = {
            "auto_scan_retry_enabled": auto_scan_retry,
            "auto_investigate_existing_enabled": auto_investigate_existing,
            "quick_requested": quick,
            "adaptive_scan_limits_enabled": True,
            "focused_deep_scan_uses_unbounded_entries": use_unbounded_entries_for_focused_deep_scan,
        }
        payload["deep_scan_plan_performed"] = deep_scan_plan_performed
        payload["deep_scan_plan_focus_roots"] = planning_focus_roots
        if planning_payload is not None:
            payload["deep_scan_plan"] = {
                "strategy": "quick_plan_then_focus",
                "planning_scan_limits": planning_payload.get("initial_scan_limits", planning_payload.get("scan", {})),
                "planning_detected_focus_roots": planning_focus_roots,
            }
        if auto_scan_retry_reason:
            payload["auto_scan_retry_reason"] = auto_scan_retry_reason
        if auto_scan_retry_performed and payload.get("scan_retry_recommended"):
            payload["scan_retry_recommended"] = False
            payload["manual_scan_retry_not_recommended"] = True
            payload["retry_exhausted"] = True
            payload["suppressed_scan_retry_command"] = payload.get("scan_retry_command")
            payload["suppressed_scan_retry_limits"] = payload.get("scan_retry_limits", {})
            payload["scan_retry_command"] = None
            payload["scan_retry_limits"] = {}
        if should_investigate:
            payload["existing_result_investigation"] = investigate_existing_results(
                path,
                payload,
                auto_scan_retry_performed=auto_scan_retry_performed,
                auto_triggered=bool(auto_investigate_reason_ids) and not investigate_existing,
                initial_limits=initial_limits,
                focus=focus,
            )
            payload = _enrich_function_detections_with_investigation(payload, payload["existing_result_investigation"])
            payload = _apply_existing_result_investigation(payload, payload["existing_result_investigation"])
            payload = _enrich_workflow_checklist_payload(
                payload,
                strict_iobrpy_investigation_item_ids=_STRICT_IOBRPY_INVESTIGATION_ITEM_IDS,
            )
        payload["agent_rendering_hints"] = _agent_rendering_hints_payload(payload)

    return payload

