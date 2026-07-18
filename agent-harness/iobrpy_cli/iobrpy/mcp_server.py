"""
MCP server for native IOBRpy commands.

This module exposes the native ``iobrpy`` subcommands as MCP tools so coding
agents can call the CLI directly instead of re-discovering workflows from
source code. The command list is derived from the real argparse parser in
``iobrpy.main`` and intentionally excludes the non-result ``ai`` orchestrator.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import platform
import shlex
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass, replace
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from iobrpy_cli.iobrpy import __version__
from iobrpy_cli.iobrpy.command_profiles import COMMAND_PROFILES
from iobrpy_cli.iobrpy.pipeline_map import map_pipeline_path


SUPPORTED_PROTOCOLS = ["2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"]
SERVER_INFO = {
    "name": "iobrpy-cli-mcp",
    "title": "iobrpy-cli MCP",
    "version": __version__,
    "description": "Expose native iobrpy commands as MCP tools.",
}
CAPABILITIES = {"tools": {"listChanged": False}}
EXCLUDED_COMMANDS = {"ai"}
FORWARDED_PARSER_RELATIVE_PATHS = {
    "runall": Path("workflow") / "runall.py",
    "tme_profile": Path("workflow") / "tme_profile.py",
    "trust4": Path("workflow") / "trust4.py",
    "extract_hla_read": Path("SpecHLA") / "extract_hla_read.py",
    "hla_typing": Path("workflow") / "hla_typing.py",
}
COMMAND_SEMANTIC_NOTES: Dict[str, Dict[str, List[str]]] = {
    "runall": {
        "default_behaviors": [
            "If --threads is omitted, runall uses 8 as the unified concurrency.",
            "If --batch_size is omitted, runall uses 1 for fastq_qc and quantification batching.",
            "If --signature is omitted, runall forwards signature scoring with --signature all.",
            "If --method is omitted for signature scoring, runall uses --method integration.",
            "If --mini_gene_count is omitted, runall uses --mini_gene_count 2.",
            "If --adjust_eset is omitted, runall enables --adjust_eset.",
            "If --remove_version is omitted, salmon mode still adds it to prepare_salmon and star mode still adds it to count2tpm.",
            "If quanTIseq-specific flags are omitted, runall enables --arrays, --tumor, and --scale_mrna.",
            "If EPIC reference is omitted, runall uses --reference TRef.",
            "If LR_cal flags are omitted, runall uses --data_type tpm, --id_type symbol, --cancer_type pancan, and --verbose.",
            "TRUST4 is invoked from the runall output FASTQ/QC tree and defaults its input via --fqdir unless overridden.",
        ],
        "notes": [
            "runall accepts more than the top-level main.py help shows because main.py forwards unknown args into workflow/runall.py.",
            "Flags such as --index, --suffix1, --signature, --perm, --platform, --reference, --cancer_type, and TRUST4 flags are routed to sub-steps.",
            "For salmon mode, --index routes to batch_salmon. For star mode, --index routes to batch_star_count.",
        ],
    },
    "tme_profile": {
        "default_behaviors": [
            "If --signature is omitted, tme_profile forwards signature scoring with --signature all.",
            "If --method is omitted for signature scoring, tme_profile uses --method integration.",
            "If --mini_gene_count is omitted, tme_profile uses --mini_gene_count 2.",
            "If --adjust_eset is omitted, tme_profile enables --adjust_eset.",
            "If ESTIMATE platform is omitted, tme_profile uses --platform affymetrix.",
            "If MCPcounter features are omitted, tme_profile uses --features HUGO_symbols.",
            "If quanTIseq-specific flags are omitted, tme_profile enables --arrays, --tumor, and --scale_mrna.",
            "If EPIC reference is omitted, tme_profile uses --reference TRef.",
            "If LR_cal flags are omitted, tme_profile uses --data_type tpm, --id_type symbol, --cancer_type pancan, and --verbose.",
        ],
        "notes": [
            "tme_profile accepts routed long flags beyond the top-level main.py help because main.py forwards unknown args into workflow/tme_profile.py.",
            "Flags are auto-routed to calculate_sig_score, cibersort, estimate, mcpcounter, quantiseq, epic, and LR_cal.",
        ],
    },
}
SECOND_STAGE_REASON_SUMMARIES: Dict[str, str] = {
    "low_current_stage_confidence": (
        "The deterministic scan could not place the directory on one workflow stage with high confidence."
    ),
    "scan_warning_present": (
        "The scan reported limits or warnings, so absence in the checklist should be treated cautiously."
    ),
    "agent_fallback_recommended": (
        "The CLI itself recommends fallback investigation because the directory looks mixed, nonstandard, or partially external."
    ),
    "external_analysis_hints_present": (
        "The directory contains external-analysis hints that may matter for interpretation but do not prove native iobrpy execution."
    ),
    "likely_iobrpy_results_detected": (
        "Some outputs look like iobrpy-style results, but they were not confirmed strongly enough to treat as ground truth."
    ),
    "reusable_existing_results_detected": (
        "Some existing results look reusable, so the agent should reason about resume versus rerun instead of assuming a fresh start."
    ),
    "existing_result_investigation_available": (
        "The CLI-native existing-result investigation found targeted candidate files that deserve human-style interpretation."
    ),
    "task_mentions_specialized_or_nondefault_functions": (
        "The task text mentions specialized or non-default functions, so the agent should widen its candidate-function judgment."
    ),
}
EXTERNAL_HINT_COMMAND_MAP: Dict[str, str] = {
    "external_hla": "hla_typing",
    "external_tcr_bcr": "trust4",
}
CLUSTER_DECONVOLUTION_INPUT_PREFERENCES: List[Dict[str, str]] = [
    {
        "id": "cibersort",
        "label": "CIBERSORT result matrix",
        "why": (
            "The clustering tutorials use CIBERSORT outputs as the default input, so this is the safest first-choice matrix when it is available."
        ),
    },
    {
        "id": "quantiseq",
        "label": "quanTIseq result matrix",
        "why": (
            "quanTIseq outputs can also support clustering, but they are usually a fallback behind CIBERSORT in the tutorial-style workflow."
        ),
    },
    {
        "id": "epic",
        "label": "EPIC result matrix",
        "why": (
            "EPIC outputs can support TME clustering when you want EPIC-based cell-fraction features instead of CIBERSORT."
        ),
    },
    {
        "id": "mcpcounter",
        "label": "MCPcounter result matrix",
        "why": (
            "MCPcounter can still be clustered, but it should not outrank a detected CIBERSORT matrix when the user asked for the default tutorial-style clustering route."
        ),
    },
    {
        "id": "estimate",
        "label": "ESTIMATE score matrix",
        "why": (
            "ESTIMATE scores can be used as a smaller feature matrix, though they are a narrower fallback than CIBERSORT-style deconvolution outputs."
        ),
    },
    {
        "id": "ips",
        "label": "IPS score matrix",
        "why": (
            "IPS outputs can support clustering as a limited score matrix, but they are usually less informative than a dedicated deconvolution matrix."
        ),
    },
]
RESOURCE_SENSITIVE_PARAMETER_NAMES = {
    "threads",
    "num_threads",
    "parallel_size",
    "num_processes",
    "batch_size",
    "t",
    "j",
}
RESOURCE_PROBE_COMMANDS = [
    "nproc",
    "free -h",
    "uptime",
    "ps -eo pcpu,pmem,comm --sort=-pcpu | head",
]
RESOURCE_SENSITIVE_DESCRIPTION = (
    "Resource-sensitive concurrency parameter: if the user did not specify it, inspect CPU cores, "
    "current load, and available memory on the execution host before choosing a value. Do not blindly "
    "use the native default or hard-code 8."
)
TASK_COMMAND_SIGNAL_SPECS: Dict[str, Dict[str, Any]] = {
    "runall": {
        "keywords": [
            "runall",
            "full pipeline",
            "full workflow",
            "fastq to tme",
            "from fastq",
            "whole pipeline",
            "完整流程",
            "全流程",
            "从fastq开始",
            "从原始reads开始",
        ],
        "reason": (
            "The task text suggests an end-to-end FASTQ-to-TME workflow rather than an isolated step, so runall should stay in the candidate-function set."
        ),
        "question_id": "confirm_runall_scope",
        "question": (
            "Do you want the full FASTQ-to-TME pipeline, or are you only trying to run one upstream or downstream step?"
        ),
        "question_why": (
            "This separates a true runall request from narrower QC, quantification, or downstream-matrix requests."
        ),
    },
    "fastq_qc": {
        "keywords": [
            "fastq_qc",
            "fastp",
            "quality control",
            "qc",
            "trim",
            "trimming",
            "read cleaning",
            "质控",
            "清洗",
            "去接头",
        ],
        "reason": (
            "The task text suggests FASTQ cleaning or QC rather than a complete quantification or downstream analysis workflow."
        ),
        "question_id": "confirm_fastq_qc_scope",
        "question": (
            "Do you want only FASTQ QC and cleaning, or should the analysis continue into quantification and downstream TME steps?"
        ),
        "question_why": (
            "This distinguishes a standalone QC request from runall or quantification-stage execution."
        ),
    },
    "batch_salmon": {
        "keywords": [
            "batch_salmon",
            "salmon quant",
            "salmon quantify",
            "salmon expression",
            "salmon定量",
            "salmon表达定量",
        ],
        "reason": (
            "The task text suggests the Salmon quantification stage specifically, so batch_salmon should stay visible instead of being collapsed into runall."
        ),
        "question_id": "confirm_batch_salmon_scope",
        "question": (
            "Do you want only the Salmon quantification step from FASTQ, or should the workflow continue into merge, TPM preparation, and downstream analysis?"
        ),
        "question_why": (
            "This distinguishes an isolated Salmon stage request from a full end-to-end pipeline."
        ),
    },
    "merge_salmon": {
        "keywords": [
            "merge_salmon",
            "merge salmon",
            "merge quant",
            "quant.sf",
            "合并salmon",
            "合并quant",
        ],
        "reason": (
            "The task text suggests merging existing Salmon quantification outputs rather than rerunning upstream FASTQ processing."
        ),
        "question_id": "confirm_merge_salmon_inputs",
        "question": (
            "Do you already have per-sample Salmon quantification directories and only need the merged TPM/NumReads tables?"
        ),
        "question_why": (
            "This confirms a merge-only task instead of a fresh Salmon quantification run."
        ),
    },
    "prepare_salmon": {
        "keywords": [
            "prepare_salmon",
            "salmon to tpm",
            "prepare salmon matrix",
            "salmon matrix",
            "salmon转tpm",
            "salmon结果转tpm",
        ],
        "reason": (
            "The task text suggests converting merged Salmon outputs into the cleaned downstream TPM matrix."
        ),
        "question_id": "confirm_prepare_salmon_scope",
        "question": (
            "Do you already have the merged Salmon table and need the cleaned TPM-style matrix for downstream analysis?"
        ),
        "question_why": (
            "This separates prepare_salmon from raw FASTQ quantification or already-complete downstream TPM workflows."
        ),
    },
    "batch_star_count": {
        "keywords": [
            "batch_star_count",
            "star count",
            "star quant",
            "readspergene",
            "genecounts",
            "star定量",
            "star计数",
        ],
        "reason": (
            "The task text suggests the STAR counting branch specifically, so batch_star_count should remain visible as the dedicated quantification step."
        ),
        "question_id": "confirm_batch_star_scope",
        "question": (
            "Do you want only the STAR alignment and counting step from FASTQ, or should the workflow continue into merged counts, TPM conversion, and downstream analysis?"
        ),
        "question_why": (
            "This distinguishes an isolated STAR stage request from a full pipeline run."
        ),
    },
    "merge_star_count": {
        "keywords": [
            "merge_star_count",
            "merge star count",
            "merge counts",
            "readspergene",
            "合并star计数",
        ],
        "reason": (
            "The task text suggests combining existing STAR GeneCounts outputs into one cohort-level matrix."
        ),
        "question_id": "confirm_merge_star_scope",
        "question": (
            "Do you already have per-sample STAR GeneCounts outputs and only need the merged cohort count matrix?"
        ),
        "question_why": (
            "This confirms a merge-only count-matrix task instead of rerunning STAR quantification."
        ),
    },
    "count2tpm": {
        "keywords": [
            "count2tpm",
            "counts to tpm",
            "count matrix to tpm",
            "计数转tpm",
            "counts转tpm",
        ],
        "reason": (
            "The task text suggests a STAR/count-matrix-to-TPM conversion step rather than a Salmon-based TPM workflow."
        ),
        "question_id": "confirm_count2tpm_scope",
        "question": (
            "Do you already have the merged count matrix and only need TPM conversion before downstream analysis?"
        ),
        "question_why": (
            "This distinguishes count2tpm from prepare_salmon and from broader FASTQ workflows."
        ),
    },
    "anno_eset": {
        "keywords": [
            "anno_eset",
            "annotation",
            "annotate",
            "probe",
            "probe id",
            "symbol mapping",
            "duplicate",
            "collapse duplicates",
            "annotated matrix",
            "gene symbol",
            "annot",
            "probeset",
            "probe set",
            "注释",
            "探针",
            "探针注释",
            "重复合并",
        ],
        "reason": (
            "The task text suggests identifier annotation or duplicate collapsing, which makes anno_eset a plausible matrix-preparation step."
        ),
        "question_id": "confirm_annotation_need",
        "question": (
            "Does the matrix use probe or non-symbol identifiers that need annotation or duplicate collapsing before downstream analysis?"
        ),
        "question_why": (
            "This distinguishes ordinary TPM-style downstream analysis from a preprocessing step such as anno_eset."
        ),
    },
    "mouse2human_eset": {
        "keywords": [
            "mouse2human",
            "mouse2human_eset",
            "mouse to human",
            "mouse-human",
            "mouse matrix",
            "mouse symbols",
            "mgi",
            "murine",
        ],
        "reason": (
            "The task text suggests cross-species preprocessing, which makes mouse2human_eset relevant before downstream human-focused methods."
        ),
        "question_id": "confirm_mouse_to_human_conversion",
        "question": (
            "Is the current matrix still in mouse gene symbols, and do you want a human-symbol matrix for downstream IOBRpy analysis?"
        ),
        "question_why": (
            "This confirms whether the next step is species conversion rather than direct TME analysis."
        ),
    },
    "log2_eset": {
        "keywords": [
            "log2_eset",
            "log2",
            "log transform",
            "log-transform",
            "log transformed",
            "log transformation",
            "log normalize",
            "log-normalize",
            "log转换",
            "对数转换",
            "对数矩阵",
        ],
        "reason": (
            "The task text explicitly mentions log transformation, which makes log2_eset a plausible preprocessing step."
        ),
        "question_id": "confirm_log2_transform",
        "question": (
            "Do you specifically need a log2(x+1)-style transformed matrix, or are you trying to continue with downstream analysis directly from the current matrix?"
        ),
        "question_why": (
            "This separates a matrix-transformation request from a standard downstream analysis request."
        ),
    },
    "calculate_sig_score": {
        "keywords": [
            "signature",
            "signature score",
            "pathway score",
            "ssgsea",
            "pathway",
            "gene set",
            "score matrix",
            "特征评分",
            "通路评分",
            "基因集评分",
            "signature评分",
        ],
        "reason": (
            "The task text emphasizes signatures or pathways, so calculate_sig_score may be relevant even if it is not the current narrow next step."
        ),
        "question_id": "confirm_signature_scope",
        "question": (
            "Do you want pathway or signature scores only, or do you also want the broader downstream TME workflow?"
        ),
        "question_why": (
            "This decides whether to stay with calculate_sig_score or move to a broader bundle such as tme_profile."
        ),
    },
    "tme_cluster": {
        "keywords": [
            "cluster",
            "clustering",
            "subtype",
            "subtyping",
            "cluster samples",
            "sample clusters",
            "聚类",
            "亚型",
            "分型",
        ],
        "reason": (
            "The task text suggests downstream subtype or clustering analysis, so tme_cluster is relevant."
        ),
        "question_id": "confirm_clustering_goal",
        "question": (
            "Are you trying to assign samples into TME-related clusters from existing TME or signature features?"
        ),
        "question_why": (
            "This confirms a downstream clustering goal instead of a request for TPM generation or deconvolution."
        ),
    },
    "nmf": {
        "keywords": [
            "nmf",
            "matrix factorization",
            "non-negative matrix factorization",
            "latent program",
            "decomposition",
            "分解",
            "潜在程序",
        ],
        "reason": (
            "The task text suggests factorization-style clustering or latent-program discovery, so nmf is relevant."
        ),
        "question_id": "confirm_nmf_goal",
        "question": (
            "Do you want an NMF-style decomposition and cluster selection workflow rather than the direct tme_cluster route?"
        ),
        "question_why": (
            "This separates NMF-specific downstream analysis from simpler clustering."
        ),
    },
    "bayesprism": {
        "keywords": [
            "bayesprism",
            "single-cell reference",
            "single cell reference",
            "scrna",
            "single-cell",
            "single cell",
            "单细胞参考",
            "单细胞ref",
            "贝叶斯反卷积",
        ],
        "reason": (
            "The task text suggests a BayesPrism-style Bayesian deconvolution branch, which can use the bundled reference by default or a custom single-cell reference when the user wants to override it."
        ),
        "question_id": "confirm_bayesprism_reference",
        "question": (
            "Do you want BayesPrism specifically, and should it use the bundled reference or a custom single-cell reference?"
        ),
        "question_why": (
            "BayesPrism is an optional branch outside the default bundle, and the agent should distinguish bundled-reference use from a custom-reference override."
        ),
    },
    "hla_typing": {
        "keywords": [
            "hla",
            "hla typing",
            "typing",
            "spechla",
            "hla分型",
            "分型",
        ],
        "reason": (
            "The task text suggests HLA analysis, so hla_typing should stay in the candidate-function set when relevant."
        ),
        "question_id": "confirm_hla_input_type",
        "question": (
            "Are you trying to run HLA typing from BAM or CRAM inputs, or are you only interpreting existing HLA results?"
        ),
        "question_why": (
            "This distinguishes a native hla_typing run from reuse of external HLA outputs."
        ),
    },
    "trust4": {
        "keywords": [
            "trust4",
            "tcr",
            "bcr",
            "vdj",
            "repertoire",
            "clonotype",
            "airr",
            "组库",
            "克隆型",
            "免疫组库",
        ],
        "reason": (
            "The task text suggests immune-repertoire analysis, so trust4 should stay in the candidate-function set when relevant."
        ),
        "question_id": "confirm_trust4_input_type",
        "question": (
            "Are you trying to run repertoire reconstruction from FASTQ or BAM inputs, or are you only interpreting existing TCR/BCR outputs?"
        ),
        "question_why": (
            "This distinguishes a native trust4 run from reuse of external repertoire outputs."
        ),
    },
    "tme_profile": {
        "keywords": [
            "tme_profile",
            "tme profile",
            "tpm to tme",
            "downstream tme",
            "microenvironment profiling",
            "下游tme",
            "微环境分析",
            "从tpm开始",
        ],
        "reason": (
            "The task text suggests a TPM-to-TME downstream bundle rather than raw FASTQ quantification or a single downstream method."
        ),
        "question_id": "confirm_tme_profile_scope",
        "question": (
            "Do you want the broader TPM-to-TME bundle from an existing TPM matrix, or only one downstream method such as CIBERSORT, signature scoring, or LR_cal?"
        ),
        "question_why": (
            "This distinguishes tme_profile from single-method downstream commands and from runall."
        ),
    },
    "cibersort": {
        "keywords": [
            "cibersort",
            "lm22",
            "cell fraction",
            "fraction deconvolution",
            "cibersort反卷积",
        ],
        "reason": (
            "The task text suggests CIBERSORT specifically, so the agent should keep the single-method command visible instead of jumping straight to the broader tme_profile bundle."
        ),
        "question_id": "confirm_cibersort_only",
        "question": (
            "Do you want only CIBERSORT output from the current TPM-like matrix, or would you rather run the broader tme_profile bundle?"
        ),
        "question_why": (
            "This separates a method-specific CIBERSORT request from a multi-method downstream TME analysis."
        ),
    },
    "IPS": {
        "keywords": [
            "ips",
            "immunophenoscore",
            "immunopheno",
            "免疫评分",
            "免疫表型评分",
        ],
        "reason": (
            "The task text suggests Immunophenoscore-specific scoring, so the agent should keep the standalone IPS command available."
        ),
        "question_id": "confirm_ips_only",
        "question": (
            "Do you want only IPS scoring from the current matrix, or should the analysis expand to the broader tme_profile bundle?"
        ),
        "question_why": (
            "This distinguishes a standalone IPS request from the default multi-method TME workflow."
        ),
    },
    "estimate": {
        "keywords": [
            "estimate",
            "stromal score",
            "immune score",
            "tumor purity",
            "estimate评分",
        ],
        "reason": (
            "The task text suggests ESTIMATE-style stromal or immune scoring, so the single-method estimate command should remain available."
        ),
        "question_id": "confirm_estimate_only",
        "question": (
            "Do you want only ESTIMATE scores, or do you want the broader tme_profile bundle from the current TPM-like matrix?"
        ),
        "question_why": (
            "This separates a standalone ESTIMATE request from a multi-method TME workflow."
        ),
    },
    "mcpcounter": {
        "keywords": [
            "mcpcounter",
            "mcp counter",
            "abundance score",
            "丰度评分",
        ],
        "reason": (
            "The task text suggests MCPcounter specifically, so the agent should keep the standalone method visible and confirm identifier handling."
        ),
        "question_id": "confirm_mcpcounter_scope",
        "question": (
            "Do you want only MCPcounter scores from the current matrix, and do the matrix rows use the identifier namespace that MCPcounter expects?"
        ),
        "question_why": (
            "MCPcounter is method-specific and its identifier namespace matters, so the agent should not silently substitute another downstream route."
        ),
    },
    "quantiseq": {
        "keywords": [
            "quantiseq",
            "quan tiseq",
            "quanTIseq",
            "quantiseq反卷积",
        ],
        "reason": (
            "The task text suggests quanTIseq specifically, so the agent should keep the standalone command visible instead of assuming the whole bundle is required."
        ),
        "question_id": "confirm_quantiseq_only",
        "question": (
            "Do you want only quanTIseq output from the current TPM-like matrix, or do you want the broader tme_profile bundle?"
        ),
        "question_why": (
            "This distinguishes a single-method quanTIseq request from the default downstream bundle."
        ),
    },
    "epic": {
        "keywords": [
            "epic",
            "epic deconvolution",
            "epic反卷积",
        ],
        "reason": (
            "The task text suggests EPIC specifically, so the standalone EPIC command should remain visible."
        ),
        "question_id": "confirm_epic_only",
        "question": (
            "Do you want only EPIC output from the current TPM-like matrix, or should the analysis expand to the broader tme_profile bundle?"
        ),
        "question_why": (
            "This separates a method-specific EPIC request from a multi-method downstream TME run."
        ),
    },
    "LR_cal": {
        "keywords": [
            "lr_cal",
            "ligand receptor",
            "ligand-receptor",
            "cell communication",
            "配体受体",
            "细胞通讯",
        ],
        "reason": (
            "The task text suggests ligand-receptor analysis specifically, so LR_cal should remain visible as the direct command."
        ),
        "question_id": "confirm_lr_scope",
        "question": (
            "Do you want only ligand-receptor scoring from the current TPM-like matrix, or are you looking for the broader downstream TME bundle?"
        ),
        "question_why": (
            "This distinguishes a direct LR_cal request from a broader downstream workflow."
        ),
    },
    "mouse2human_eset": {
        "keywords": [
            "mouse2human",
            "mouse2human_eset",
            "mouse to human",
            "mouse-human",
            "mouse matrix",
            "mouse symbols",
            "mgi",
            "murine",
            "小鼠转人",
            "鼠转人",
        ],
        "reason": (
            "The task text suggests cross-species preprocessing, which makes mouse2human_eset relevant before downstream human-focused methods."
        ),
        "question_id": "confirm_mouse_to_human_conversion",
        "question": (
            "Is the current matrix still in mouse gene symbols, and do you want a human-symbol matrix for downstream IOBRpy analysis?"
        ),
        "question_why": (
            "This confirms whether the next step is species conversion rather than direct TME analysis."
        ),
    },
    "extract_hla_read": {
        "keywords": [
            "extract_hla_read",
            "extract hla",
            "hla read extraction",
            "提取hla",
            "提取hla reads",
        ],
        "reason": (
            "The task text suggests the BAM-to-HLA-read extraction step specifically, so extract_hla_read should remain visible."
        ),
        "question_id": "confirm_extract_hla_read_scope",
        "question": (
            "Are you starting from sorted and indexed BAM/CRAM files and only need the HLA-read extraction step before SpecHLA?"
        ),
        "question_why": (
            "This distinguishes extract_hla_read from the batch hla_typing wrapper and from direct single-sample SpecHLA runs."
        ),
    },
    "spechla": {
        "keywords": [
            "spechla",
            "spec hla",
            "single-sample hla",
            "单样本hla",
        ],
        "reason": (
            "The task text suggests a direct per-sample SpecHLA run rather than the batch hla_typing wrapper."
        ),
        "question_id": "confirm_spechla_scope",
        "question": (
            "Are you trying to run SpecHLA for one sample from paired FASTQ inputs, rather than the batch hla_typing workflow over a BAM directory?"
        ),
        "question_why": (
            "This distinguishes single-sample SpecHLA execution from batch HLA typing and BAM-side extraction."
        ),
    },
}

FULL_PIPELINE_COMMANDS = {"runall"}
RAW_FASTQ_STAGE_COMMANDS = {"fastq_qc", "batch_salmon", "batch_star_count"}
SALMON_STAGE_COMMANDS = {"merge_salmon", "prepare_salmon"}
STAR_STAGE_COMMANDS = {"merge_star_count", "count2tpm"}
MATRIX_PREPARATION_COMMANDS = {"anno_eset", "mouse2human_eset", "log2_eset"}
DOWNSTREAM_BUNDLE_COMMANDS = {"tme_profile"}
DOWNSTREAM_SINGLE_METHOD_COMMANDS = {
    "calculate_sig_score",
    "cibersort",
    "IPS",
    "estimate",
    "mcpcounter",
    "quantiseq",
    "epic",
    "LR_cal",
    "bayesprism",
}
CLUSTERING_COMMANDS = {"tme_cluster", "nmf"}
HLA_COMMANDS = {"extract_hla_read", "spechla", "hla_typing"}
REPERTOIRE_COMMANDS = {"trust4"}
DECONVOLUTION_METHOD_COMMANDS = {"cibersort", "IPS", "estimate", "mcpcounter", "quantiseq", "epic", "bayesprism"}
DOWNSTREAM_MATRIX_COMMANDS = DOWNSTREAM_BUNDLE_COMMANDS | DOWNSTREAM_SINGLE_METHOD_COMMANDS
DECONVOLUTION_METHOD_EVIDENCE_IDS = {
    "cibersort": "cibersort",
    "IPS": "ips",
    "estimate": "estimate",
    "mcpcounter": "mcpcounter",
    "quantiseq": "quantiseq",
    "epic": "epic",
    "bayesprism": "bayesprism",
}
SALMON_STAGE_DIR_NAMES = ("11-salmon", "02-salmon")
STAR_STAGE_DIR_NAMES = ("04-star", "02-star")
PIPELINE_ROOT_MARKERS = (
    "11-salmon",
    "04-star",
    "03-tpm",
    "26-tme-profile-iobrpy",
    "22-tme",
    "13-hla",
    "19-mixcr",
    "02-fastq",
    "03-fastp",
    "02-salmon",
    "02-star",
    "01-qc",
    "07-TCRBCR",
)


@dataclass(frozen=True)
class ArgumentSpec:
    """Metadata for one exposed tool argument."""

    property_name: str
    option_strings: List[str]
    required: bool
    help_text: str
    schema: Dict[str, Any]
    is_flag: bool
    is_boolean_value: bool
    is_array: bool
    is_positional: bool
    origin: str = "native"
    help_group: str = "declared"

    @property
    def preferred_option(self) -> Optional[str]:
        for option in self.option_strings:
            if option.startswith("--"):
                return option
        return self.option_strings[0] if self.option_strings else None


_PARSER_CACHE: Optional[argparse.ArgumentParser] = None
_NATIVE_TOOL_CACHE: Optional[List[Dict[str, Any]]] = None
_HELPER_TOOL_CACHE: Optional[List[Dict[str, Any]]] = None
_TOOL_CACHE: Optional[List[Dict[str, Any]]] = None
_ARGUMENT_CACHE: Optional[Dict[str, List[ArgumentSpec]]] = None
_COMMAND_META_CACHE: Optional[Dict[str, Dict[str, Any]]] = None
_SCRIPT_PARSER_CACHE: Dict[str, argparse.ArgumentParser] = {}
_RULES_CACHE: Optional[Dict[str, Dict[str, Any]]] = None

_HELPER_TOOL_NAMES = {
    "list_native_commands",
    "map_path",
    "recommend_workflow",
    "doctor_environment",
    "analyze_path_task",
}
_INTERNAL_NATIVE_ARGUMENT_PREFIX = "_"


def _command_profile(command_name: str) -> Dict[str, Any]:
    return COMMAND_PROFILES.get(
        command_name,
        {
            "summary": "",
            "detailed_description": "",
            "input_expectations": [],
            "outputs": [],
            "when_to_use": [],
            "use_another_command_when": [],
            "important_optional_parameters": [],
        },
    )


def _pick_rule_text(value: Any) -> str:
    if isinstance(value, dict):
        for key in ("en", "zh"):
            text = value.get(key)
            if isinstance(text, str) and text.strip():
                return text.strip()
        return ""
    if isinstance(value, str):
        return value.strip()
    return ""


def _find_required_params_path() -> Path:
    package_candidate = _find_package_root() / "RAG_MCP" / "iobrpy_required_params.json"
    if package_candidate.exists():
        return package_candidate

    repo_candidate = Path(__file__).resolve().parents[3] / "src" / "iobrpy" / "RAG_MCP" / "iobrpy_required_params.json"
    if repo_candidate.exists():
        return repo_candidate

    raise RuntimeError("Could not locate RAG_MCP/iobrpy_required_params.json.")


def _load_required_param_rules() -> Dict[str, Dict[str, Any]]:
    global _RULES_CACHE
    if _RULES_CACHE is not None:
        return _RULES_CACHE

    raw = json.loads(_find_required_params_path().read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise RuntimeError("iobrpy_required_params.json must contain an object at the top level.")

    rules: Dict[str, Dict[str, Any]] = {}
    for command_name, rule in raw.items():
        if command_name in EXCLUDED_COMMANDS or not isinstance(rule, dict):
            continue

        required = [item for item in rule.get("required", []) if isinstance(item, str)]
        optional = [item for item in rule.get("optional", []) if isinstance(item, str)]
        confirm = [item for item in rule.get("confirm", []) if isinstance(item, str)]
        required_one_of = [
            [item for item in group if isinstance(item, str)]
            for group in rule.get("required_one_of", [])
            if isinstance(group, list)
        ]
        choices = {
            key: [item for item in value if isinstance(item, str)]
            for key, value in (rule.get("choices", {}) or {}).items()
            if isinstance(key, str) and isinstance(value, list)
        }
        param_hints = {
            key: value
            for key, value in (rule.get("param_hints", {}) or {}).items()
            if isinstance(key, str)
        }
        cli_flags = {
            key: value
            for key, value in (rule.get("cli_flags", {}) or {}).items()
            if isinstance(key, str) and isinstance(value, str) and value.strip()
        }
        param_types = {
            key: value
            for key, value in (rule.get("param_types", {}) or {}).items()
            if isinstance(key, str) and isinstance(value, str) and value.strip()
        }
        optional_defaults = {
            key: value
            for key, value in (rule.get("optional_defaults", {}) or {}).items()
            if isinstance(key, str)
        }
        bool_mode = {
            key: value
            for key, value in (rule.get("bool_mode", {}) or {}).items()
            if isinstance(key, str) and isinstance(value, str) and value.strip()
        }
        param_aliases = {
            key: [item for item in value if isinstance(item, str)]
            for key, value in (rule.get("param_aliases", {}) or {}).items()
            if isinstance(key, str) and isinstance(value, list)
        }

        parameter_order: List[str] = []
        for item in required + optional:
            if item not in parameter_order:
                parameter_order.append(item)

        rules[command_name] = {
            "required": required,
            "optional": optional,
            "confirm": confirm,
            "required_one_of": required_one_of,
            "choices": choices,
            "param_hints": param_hints,
            "function_summary": rule.get("function_summary"),
            "notes": rule.get("notes"),
            "optional_defaults": optional_defaults,
            "cli_flags": cli_flags,
            "param_types": param_types,
            "bool_mode": bool_mode,
            "param_aliases": param_aliases,
            "parameter_order": parameter_order,
        }

    _RULES_CACHE = rules
    return rules


def _default_cli_flag_for_key(key: str) -> str:
    return f"--{key}"


def _is_resource_sensitive_parameter(parameter_name: str) -> bool:
    normalized = parameter_name.strip().lstrip("-").replace("-", "_").lower()
    return normalized in RESOURCE_SENSITIVE_PARAMETER_NAMES


def _append_resource_sensitive_description(description: str) -> str:
    text = description.strip()
    if not text:
        return RESOURCE_SENSITIVE_DESCRIPTION
    if RESOURCE_SENSITIVE_DESCRIPTION in text:
        return text
    return f"{text} {RESOURCE_SENSITIVE_DESCRIPTION}"


def _decorate_resource_sensitive_schema(param_name: str, schema: Dict[str, Any]) -> None:
    if not _is_resource_sensitive_parameter(param_name):
        return
    schema["x-iobrpy-resource-sensitive"] = True
    schema["x-iobrpy-resource-probe-recommended"] = True
    schema["x-iobrpy-resource-probe-commands"] = list(RESOURCE_PROBE_COMMANDS)


def _decorate_resource_sensitive_parameter_card(card: Dict[str, Any]) -> Dict[str, Any]:
    parameter = str(card.get("parameter") or "").strip()
    if not _is_resource_sensitive_parameter(parameter):
        return card
    updated = dict(card)
    ask_when = str(updated.get("ask_when") or "").strip()
    resource_ask_when = (
        "Before execution, inspect CPU cores, current CPU load, and available memory on the same host "
        "that will run the command; then ask a concise natural-language question if the user did not "
        "already specify this value."
    )
    updated["decision_type"] = "compute_resource"
    updated["resource_sensitive"] = True
    updated["resource_probe_recommended"] = True
    updated["resource_probe_commands"] = list(RESOURCE_PROBE_COMMANDS)
    updated["resource_probe_scope"] = "execution_host"
    updated["ask_when"] = f"{ask_when} {resource_ask_when}".strip()
    updated["natural_language_confirmation"] = (
        "First inspect CPU and memory usage on the execution host, then recommend a thread or parallelism value; "
        "do not blindly use the default or hard-code 8."
    )
    updated["natural_language_confirmation_zh"] = (
        "先看运行机器的 CPU/内存占用，再建议线程或并行数；不要直接套默认值或硬编码 8。"
    )
    return updated


def _resolve_rule_cli_flag(command_name: str, param_name: str, rule: Dict[str, Any]) -> str:
    cli_flags = rule.get("cli_flags", {}) if isinstance(rule.get("cli_flags", {}), dict) else {}
    flag = cli_flags.get(param_name)
    if isinstance(flag, str) and flag.strip():
        return flag.strip()
    return _default_cli_flag_for_key(param_name)


def _schema_for_param(rule: Dict[str, Any], param_name: str) -> Dict[str, Any]:
    param_types = rule.get("param_types", {}) if isinstance(rule.get("param_types", {}), dict) else {}
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}
    optional_defaults = rule.get("optional_defaults", {}) if isinstance(rule.get("optional_defaults", {}), dict) else {}

    param_type = str(param_types.get(param_name, "string")).strip().lower()
    if param_type == "list":
        schema: Dict[str, Any] = {"type": "array", "items": {"type": "string"}}
    elif param_type == "int":
        schema = {"type": "integer"}
    elif param_type == "float":
        schema = {"type": "number"}
    elif param_type == "bool":
        schema = {"type": "boolean"}
    else:
        schema = {"type": "string"}

    if param_name in choices and isinstance(choices[param_name], list):
        if schema.get("type") == "array":
            schema["items"]["enum"] = list(choices[param_name])
        elif schema.get("type") == "string":
            schema["enum"] = list(choices[param_name])

    if param_name in optional_defaults and optional_defaults[param_name] is not None:
        schema["default"] = optional_defaults[param_name]

    _decorate_resource_sensitive_schema(param_name, schema)
    return schema


def _build_rule_argument_specs(command_name: str, rule: Dict[str, Any]) -> List[ArgumentSpec]:
    specs: List[ArgumentSpec] = []
    required = set(rule.get("required", []))
    param_hints = rule.get("param_hints", {}) if isinstance(rule.get("param_hints", {}), dict) else {}
    param_types = rule.get("param_types", {}) if isinstance(rule.get("param_types", {}), dict) else {}
    bool_mode = rule.get("bool_mode", {}) if isinstance(rule.get("bool_mode", {}), dict) else {}

    for param_name in rule.get("parameter_order", []):
        param_type = str(param_types.get(param_name, "string")).strip().lower()
        mode = str(bool_mode.get(param_name, "")).strip().lower()
        option = _resolve_rule_cli_flag(command_name, param_name, rule)
        schema = _schema_for_param(rule, param_name)
        is_array = schema.get("type") == "array"
        is_bool = param_type == "bool"
        is_boolean_value = is_bool and mode == "value"
        is_flag = is_bool and not is_boolean_value

        help_text = _pick_rule_text(param_hints.get(param_name))
        if _is_resource_sensitive_parameter(param_name):
            help_text = _append_resource_sensitive_description(help_text)

        specs.append(
            ArgumentSpec(
                property_name=param_name,
                option_strings=[option] if option else [],
                required=param_name in required,
                help_text=help_text,
                schema=schema,
                is_flag=is_flag,
                is_boolean_value=is_boolean_value,
                is_array=is_array,
                is_positional=False,
                origin="rules",
                help_group="declared",
            )
        )

    return specs


def _build_confirm_parameters(rule: Dict[str, Any], specs: List[ArgumentSpec]) -> List[Dict[str, Any]]:
    spec_lookup = {spec.property_name: spec for spec in specs}
    confirm_parameters: List[Dict[str, Any]] = []

    for param_name in rule.get("confirm", []):
        spec = spec_lookup.get(param_name)
        if spec is None:
            continue
        confirm_parameters.append(
            _decorate_resource_sensitive_parameter_card(
                {
                    "parameter": param_name,
                    "option_strings": list(spec.option_strings),
                    "display_name": ", ".join(spec.option_strings) if spec.option_strings else param_name,
                    "preferred_option": spec.preferred_option,
                    "required": spec.required,
                    "description": spec.help_text,
                }
            )
        )

    return confirm_parameters


def _build_rule_notes(rule: Dict[str, Any]) -> List[str]:
    notes_text = _pick_rule_text(rule.get("notes"))
    return [notes_text] if notes_text else []


def _schema_choices(spec: ArgumentSpec) -> List[str]:
    if "enum" in spec.schema and isinstance(spec.schema["enum"], list):
        return [str(item) for item in spec.schema["enum"]]
    items = spec.schema.get("items", {})
    if isinstance(items, dict) and isinstance(items.get("enum"), list):
        return [str(item) for item in items["enum"]]
    return []


def _send(obj: Dict[str, Any]) -> None:
    sys.stdout.write(json.dumps(obj, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def _jsonrpc_error(_id: Any, code: int, message: str, data: Optional[Dict[str, Any]] = None) -> None:
    err: Dict[str, Any] = {"code": code, "message": message}
    if data is not None:
        err["data"] = data
    _send({"jsonrpc": "2.0", "id": _id, "error": err})


def _find_native_main_path() -> Path:
    package_spec = importlib.util.find_spec("iobrpy")
    if package_spec and package_spec.submodule_search_locations:
        for location in package_spec.submodule_search_locations:
            candidate = Path(location) / "main.py"
            if candidate.exists():
                return candidate

    module_spec = importlib.util.find_spec("iobrpy.main")
    if module_spec and module_spec.origin:
        candidate = Path(module_spec.origin)
        if candidate.exists():
            return candidate

    raise RuntimeError("Could not locate iobrpy/main.py for MCP tool generation.")


def _find_package_root() -> Path:
    return _find_native_main_path().parent


def _find_command_source_path(command_name: str) -> Optional[Path]:
    relative = FORWARDED_PARSER_RELATIVE_PATHS.get(command_name)
    if relative is None:
        return None
    candidate = _find_package_root() / relative
    return candidate if candidate.exists() else None


def _ast_value(node: ast.AST) -> Any:
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.List):
        return [_ast_value(item) for item in node.elts]
    if isinstance(node, ast.Tuple):
        return tuple(_ast_value(item) for item in node.elts)
    if isinstance(node, ast.Set):
        return {_ast_value(item) for item in node.elts}
    if isinstance(node, ast.Dict):
        return {_ast_value(key): _ast_value(value) for key, value in zip(node.keys, node.values)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_ast_value(node.operand)
    if isinstance(node, ast.Name):
        builtin_values = {
            "True": True,
            "False": False,
            "None": None,
            "int": int,
            "float": float,
            "str": str,
            "bool": bool,
        }
        if node.id in builtin_values:
            return builtin_values[node.id]
    if isinstance(node, ast.Attribute):
        return ast.unparse(node)
    if isinstance(node, ast.Lambda):
        return bool
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "set":
        if not node.args:
            return set()
        if len(node.args) == 1:
            arg_value = _ast_value(node.args[0])
            return set(arg_value)
    raise ValueError(f"Unsupported AST node in parser extraction: {ast.dump(node)}")


def _build_parser_from_source(main_path: Path) -> argparse.ArgumentParser:
    source = main_path.read_text(encoding="utf-8")
    module = ast.parse(source, filename=str(main_path))
    main_function = next(
        (node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == "main"),
        None,
    )
    if main_function is None:
        raise RuntimeError(f"Could not find main() in {main_path}")

    parser = argparse.ArgumentParser(prog="iobrpy")
    subparsers = parser.add_subparsers(dest="command", required=True)
    parser_vars: Dict[str, argparse.ArgumentParser] = {}

    for statement in main_function.body:
        if isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call):
            call = statement.value
            target_names = [target.id for target in statement.targets if isinstance(target, ast.Name)]
        elif isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Call):
            call = statement.value
            target_names = []
        else:
            continue

        if not isinstance(call.func, ast.Attribute) or not isinstance(call.func.value, ast.Name):
            continue

        if call.func.value.id == "subparsers" and call.func.attr == "add_parser":
            if not target_names:
                continue
            command_name = _ast_value(call.args[0])
            parser_kwargs: Dict[str, Any] = {}
            for keyword in call.keywords:
                if keyword.arg in {"help", "description"}:
                    parser_kwargs[keyword.arg] = _ast_value(keyword.value)
            parser_vars[target_names[0]] = subparsers.add_parser(command_name, **parser_kwargs)
            continue

        owner = parser_vars.get(call.func.value.id)
        if owner is None or call.func.attr != "add_argument":
            continue

        arg_strings = [_ast_value(arg) for arg in call.args]
        arg_kwargs: Dict[str, Any] = {}
        for keyword in call.keywords:
            if keyword.arg is None:
                continue
            arg_kwargs[keyword.arg] = _ast_value(keyword.value)
        owner.add_argument(*arg_strings, **arg_kwargs)

    return parser


def _build_standalone_parser_from_source(source_path: Path) -> argparse.ArgumentParser:
    source = source_path.read_text(encoding="utf-8")
    module = ast.parse(source, filename=str(source_path))
    target_function = next(
        (
            node
            for node in module.body
            if isinstance(node, ast.FunctionDef) and node.name in {"build_arg_parser", "main"}
        ),
        None,
    )
    if target_function is None:
        raise RuntimeError(f"Could not find a parser-building function in {source_path}")

    parser_vars: Dict[str, argparse.ArgumentParser] = {}
    parser_order: List[str] = []

    for statement in target_function.body:
        if not isinstance(statement, ast.Assign) or not isinstance(statement.value, ast.Call):
            continue

        call = statement.value
        if not isinstance(call.func, ast.Attribute):
            continue
        if not isinstance(call.func.value, ast.Name):
            continue
        if call.func.value.id != "argparse" or call.func.attr != "ArgumentParser":
            continue

        parser_kwargs: Dict[str, Any] = {}
        for keyword in call.keywords:
            if keyword.arg in {"prog", "description"}:
                parser_kwargs[keyword.arg] = _ast_value(keyword.value)
        parser = argparse.ArgumentParser(**parser_kwargs)

        for target in statement.targets:
            if isinstance(target, ast.Name):
                parser_vars[target.id] = parser
                parser_order.append(target.id)

    if not parser_order:
        raise RuntimeError(f"Could not construct parser from {source_path}")

    for statement in target_function.body:
        call: Optional[ast.Call] = None
        if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Call):
            call = statement.value
        elif isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Call):
            call = statement.value
        if call is None:
            continue

        if not isinstance(call.func, ast.Attribute) or not isinstance(call.func.value, ast.Name):
            continue

        owner = parser_vars.get(call.func.value.id)
        if owner is None or call.func.attr != "add_argument":
            continue

        arg_strings = [_ast_value(arg) for arg in call.args]
        arg_kwargs: Dict[str, Any] = {}
        for keyword in call.keywords:
            if keyword.arg is None:
                continue
            arg_kwargs[keyword.arg] = _ast_value(keyword.value)
        owner.add_argument(*arg_strings, **arg_kwargs)

    return parser_vars[parser_order[0]]


def _capture_native_parser() -> argparse.ArgumentParser:
    global _PARSER_CACHE
    if _PARSER_CACHE is not None:
        return _PARSER_CACHE

    parser = _build_parser_from_source(_find_native_main_path())
    _PARSER_CACHE = parser
    return parser


def _capture_script_parser(command_name: str) -> Optional[argparse.ArgumentParser]:
    if command_name in _SCRIPT_PARSER_CACHE:
        return _SCRIPT_PARSER_CACHE[command_name]

    source_path = _find_command_source_path(command_name)
    if source_path is None:
        return None

    parser = _build_standalone_parser_from_source(source_path)
    _SCRIPT_PARSER_CACHE[command_name] = parser
    return parser


def _extract_top_level_value_from_source(source_path: Path, name: str) -> Any:
    source = source_path.read_text(encoding="utf-8")
    module = ast.parse(source, filename=str(source_path))
    for statement in module.body:
        value: Optional[ast.AST] = None
        if isinstance(statement, ast.Assign):
            targets = statement.targets
            if any(isinstance(target, ast.Name) and target.id == name for target in targets):
                value = statement.value
        elif isinstance(statement, ast.AnnAssign):
            if isinstance(statement.target, ast.Name) and statement.target.id == name:
                value = statement.value
        if value is not None:
            return _ast_value(value)
    return None


def _get_subparsers(parser: argparse.ArgumentParser) -> Dict[str, argparse.ArgumentParser]:
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    raise RuntimeError("Native iobrpy parser does not expose subcommands.")


def _schema_type_for_action(action: argparse.Action) -> str:
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "boolean"
    if getattr(action, "type", None) is int:
        return "integer"
    if getattr(action, "type", None) is float:
        return "number"
    if isinstance(getattr(action, "default", None), bool):
        return "boolean"
    return "string"


def _default_is_meaningful(action: argparse.Action) -> bool:
    return (
        action.default is not argparse.SUPPRESS
        and action.default is not None
        and not isinstance(action, argparse._HelpAction)
    )


def _action_property_name(action: argparse.Action) -> str:
    long_opts = [opt for opt in action.option_strings if opt.startswith("--")]
    if long_opts:
        return long_opts[0].lstrip("-")
    if action.option_strings:
        return action.dest
    return action.dest


def _action_to_argument_spec(action: argparse.Action) -> ArgumentSpec:
    is_flag = isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction))
    is_positional = not action.option_strings
    property_name = _action_property_name(action)
    value_type = _schema_type_for_action(action)
    is_array = action.nargs in ("+", "*") or (isinstance(action.nargs, int) and action.nargs > 1)

    schema: Dict[str, Any]
    if is_array:
        schema = {
            "type": "array",
            "items": {"type": value_type},
        }
    else:
        schema = {"type": value_type}

    if getattr(action, "choices", None):
        schema["enum"] = list(action.choices)

    if is_flag:
        schema["default"] = False
    elif _default_is_meaningful(action):
        schema["default"] = action.default

    help_text = action.help or ""
    if action.option_strings:
        help_text = f"Flags: {', '.join(action.option_strings)}. {help_text}".strip()

    return ArgumentSpec(
        property_name=property_name,
        option_strings=list(action.option_strings),
        required=bool(action.required),
        help_text=help_text,
        schema=schema,
        is_flag=is_flag,
        is_boolean_value=(value_type == "boolean" and not is_flag),
        is_array=is_array,
        is_positional=is_positional,
    )


def _specs_from_parser(parser: argparse.ArgumentParser, *, origin: str, help_group: str) -> List[ArgumentSpec]:
    specs: List[ArgumentSpec] = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        specs.append(replace(_action_to_argument_spec(action), origin=origin, help_group=help_group))
    return specs


def _dedupe_specs(specs: List[ArgumentSpec]) -> List[ArgumentSpec]:
    deduped: List[ArgumentSpec] = []
    seen: set[str] = set()
    for spec in specs:
        if spec.property_name in seen:
            continue
        deduped.append(spec)
        seen.add(spec.property_name)
    return deduped


def _build_routed_specs(
    command_name: str,
    subparser_lookup: Dict[str, argparse.ArgumentParser],
    command_specs: Dict[str, List[ArgumentSpec]],
) -> List[ArgumentSpec]:
    source_path = _find_command_source_path(command_name)
    if source_path is None:
        return []

    flag_buckets = _extract_top_level_value_from_source(source_path, "FLAG_BUCKETS") or {}
    if not isinstance(flag_buckets, dict):
        return []

    routed_specs: List[ArgumentSpec] = []
    forwarded_parsers = {
        name: _capture_script_parser(name)
        for name in FORWARDED_PARSER_RELATIVE_PATHS
        if _capture_script_parser(name) is not None
    }

    for target_command, flags in flag_buckets.items():
        if not isinstance(flags, (set, list, tuple)):
            continue

        target_specs = command_specs.get(target_command)
        if target_specs is None:
            target_parser = forwarded_parsers.get(target_command) or subparser_lookup.get(target_command)
            if target_parser is None:
                continue
            target_specs = _specs_from_parser(
                target_parser,
                origin=f"routed:{target_command}",
                help_group="routed",
            )
            command_specs[target_command] = target_specs

        for flag_name in flags:
            if not isinstance(flag_name, str):
                continue
            match = next(
                (
                    spec
                    for spec in target_specs
                    if spec.property_name == flag_name or spec.property_name.lower() == flag_name.lower()
                ),
                None,
            )
            if match is None:
                continue
            routed_specs.append(
                replace(
                    match,
                    required=False,
                    origin=f"routed:{target_command}",
                    help_group="routed",
                    help_text=f"Routed to `{target_command}`. {match.help_text}".strip(),
                )
            )

    return routed_specs


def _normalize_important_optional_parameters(
    command_name: str,
    merged_specs: List[ArgumentSpec],
    profile: Dict[str, Any],
) -> List[Dict[str, Any]]:
    spec_lookup = {spec.property_name: spec for spec in merged_specs}
    normalized: List[Dict[str, Any]] = []

    for item in profile.get("important_optional_parameters", []):
        if not isinstance(item, dict):
            continue

        parameter = item.get("parameter")
        if not isinstance(parameter, str):
            continue

        spec = spec_lookup.get(parameter)
        if spec is None:
            continue
        if spec.required:
            continue

        default_behavior = str(item.get("native_default_behavior", "")).strip()
        if not default_behavior:
            if spec.is_flag and spec.schema.get("default") is False:
                default_behavior = "Native default is false/off."
            elif "default" in spec.schema:
                default_behavior = f"Native default: {spec.schema['default']}"

        normalized.append(
            _decorate_resource_sensitive_parameter_card(
                {
                    "parameter": parameter,
                    "option_strings": list(spec.option_strings),
                    "display_name": ", ".join(spec.option_strings) if spec.option_strings else parameter,
                    "preferred_option": spec.preferred_option,
                    "origin": spec.origin,
                    "help_group": spec.help_group,
                    "why_it_matters": str(item.get("why_it_matters", "")).strip(),
                    "native_default_behavior": default_behavior,
                    "ask_when": str(item.get("ask_when", "")).strip(),
                    "help_text": spec.help_text,
                    "description": spec.help_text,
                }
            )
        )

    return normalized


def _merge_parameter_cards_by_parameter(*groups: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    merged: List[Dict[str, Any]] = []
    seen: set[str] = set()
    for group in groups:
        for item in group:
            if not isinstance(item, dict):
                continue
            parameter = str(item.get("parameter") or "").strip()
            if not parameter or parameter in seen:
                continue
            merged.append(item)
            seen.add(parameter)
    return merged


def _build_command_metadata(parser: argparse.ArgumentParser) -> Dict[str, Dict[str, Any]]:
    global _ARGUMENT_CACHE, _COMMAND_META_CACHE
    if _ARGUMENT_CACHE is not None and _COMMAND_META_CACHE is not None:
        return _COMMAND_META_CACHE

    metadata: Dict[str, Dict[str, Any]] = {}
    rules = _load_required_param_rules()
    required_params_path = _find_required_params_path()

    for command_name, rule in rules.items():
        profile = _command_profile(command_name)
        declared_specs = _build_rule_argument_specs(command_name, rule)
        summary_text = _pick_rule_text(rule.get("function_summary")) or profile["summary"] or command_name
        source_paths = [str(required_params_path)]
        command_source_path = _find_command_source_path(command_name)
        if command_source_path is not None:
            source_paths.append(str(command_source_path))
        else:
            source_paths.append(str(_find_native_main_path()))

        confirm_parameters = _build_confirm_parameters(rule, declared_specs)
        important_optional_parameters = _merge_parameter_cards_by_parameter(
            confirm_parameters,
            _normalize_important_optional_parameters(command_name, declared_specs, profile),
        )

        metadata[command_name] = {
            "description": summary_text,
            "argument_specs": declared_specs,
            "declared_specs": declared_specs,
            "routed_specs": [],
            "source_paths": source_paths,
            "summary": summary_text,
            "detailed_description": profile.get("detailed_description", ""),
            "input_expectations": profile["input_expectations"],
            "outputs": profile["outputs"],
            "when_to_use": profile["when_to_use"],
            "use_another_command_when": profile["use_another_command_when"],
            "confirm_parameters": confirm_parameters,
            "important_optional_parameters": important_optional_parameters,
            "project_context_notes": list(profile.get("project_context_notes", [])),
            "default_behaviors": [],
            "notes": _build_rule_notes(rule),
            "required_parameter_names": list(rule.get("required", [])),
            "optional_parameter_names": list(rule.get("optional", [])),
            "confirm_parameter_names": list(rule.get("confirm", [])),
            "required_one_of": list(rule.get("required_one_of", [])),
            "cli_flags": dict(rule.get("cli_flags", {})),
            "param_hints": dict(rule.get("param_hints", {})),
            "parameter_source": str(required_params_path),
        }

    _ARGUMENT_CACHE = {name: data["argument_specs"] for name, data in metadata.items()}
    _COMMAND_META_CACHE = metadata
    return metadata


def _build_command_argument_specs(parser: argparse.ArgumentParser) -> Dict[str, List[ArgumentSpec]]:
    return {name: data["argument_specs"] for name, data in _build_command_metadata(parser).items()}


def build_tools_from_parser(parser: argparse.ArgumentParser) -> List[Dict[str, Any]]:
    command_metadata = _build_command_metadata(parser)
    tools: List[Dict[str, Any]] = []

    for command_name, meta in command_metadata.items():
        properties: Dict[str, Any] = {}
        required: List[str] = list(meta.get("required_parameter_names", []))

        for spec in meta["argument_specs"]:
            properties[spec.property_name] = dict(spec.schema)
            if spec.help_text:
                properties[spec.property_name]["description"] = spec.help_text
            if spec.option_strings:
                properties[spec.property_name]["x-cli-flags"] = list(spec.option_strings)

        input_schema: Dict[str, Any] = {
            "type": "object",
            "properties": properties,
            "required": required,
            "additionalProperties": False,
        }
        if meta.get("required_one_of"):
            input_schema["anyOf"] = [{"required": group} for group in meta["required_one_of"] if group]

        tools.append(
            {
                "name": command_name,
                "title": f"iobrpy {command_name}",
                "description": " ".join(
                    part
                    for part in [
                        f"Native iobrpy command `{command_name}`.",
                        str(meta["summary"]).strip(),
                        str(meta["detailed_description"]).strip(),
                    ]
                    if part
                ),
                "inputSchema": {
                    **input_schema,
                },
                "annotations": {
                    "sourcePaths": meta["source_paths"],
                    "parameterSource": meta.get("parameter_source"),
                    "summary": meta["summary"],
                    "detailedDescription": meta["detailed_description"],
                    "inputExpectations": meta["input_expectations"],
                    "outputs": meta["outputs"],
                    "whenToUse": meta["when_to_use"],
                    "useAnotherCommandWhen": meta["use_another_command_when"],
                    "confirmParameters": meta["confirm_parameters"],
                    "importantOptionalParameters": meta["important_optional_parameters"],
                    "resourceSensitiveParameters": [
                        item for item in meta["important_optional_parameters"]
                        if isinstance(item, dict) and item.get("resource_sensitive")
                    ],
                    "defaultBehaviors": meta["default_behaviors"],
                    "notes": meta["notes"],
                    "requiredParameters": meta.get("required_parameter_names", []),
                    "optionalParameters": meta.get("optional_parameter_names", []),
                    "requiredOneOf": meta.get("required_one_of", []),
                    "routedArguments": [],
                },
            }
        )

    return tools


def _build_helper_tools() -> List[Dict[str, Any]]:
    return [
        {
            "name": "list_native_commands",
            "title": "List native iobrpy commands",
            "description": (
                "Return the native iobrpy command catalog with semantic summaries, inputs, outputs, "
                "defaults, and confirm-parameter hints."
            ),
            "inputSchema": {
                "type": "object",
                "properties": {
                    "details": {
                        "type": "boolean",
                        "default": False,
                        "description": "When true, include full native MCP tool schemas in the payload.",
                    },
                },
                "additionalProperties": False,
            },
            "annotations": {
                "toolType": "agent_helper",
                "returnsStructuredPayload": True,
                "recommendedUse": "Use before execution when the agent needs the native command surface.",
            },
        },
        {
            "name": "map_path",
            "title": "Map an analysis path to the IOBRpy roadmap",
            "description": (
                "Inspect a file or directory, detect the current IOBRpy workflow stage, and return the "
                "structured pipeline-map payload used for continue-versus-rerun reasoning. The default "
                "mode is the deeper scan path, which is more reliable for structurally complex directories."
            ),
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "File or directory to place on the IOBRpy roadmap.",
                    },
                    "max_depth": {
                        "type": "integer",
                        "default": 8,
                        "description": "Maximum directory depth to inspect.",
                    },
                    "max_entries": {
                        "type": "integer",
                        "default": 8000,
                        "description": "Maximum file or directory names to scan. Use 0 to disable the entry cap for a focused deep scan.",
                    },
                    "quick": {
                        "type": "boolean",
                        "default": False,
                        "description": (
                            "Deprecated for agent helper calls. Older clients may still pass this, "
                            "but map_path ignores it and uses the full/deep scan path by default."
                        ),
                    },
                    "focus": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "Deep-scan one or more likely result subdirectories while only shallowly "
                            "censusing other large generic branches. This is especially useful for "
                            "structurally complex directories when you want deeper inspection without "
                            "treating the whole tree as a flat full scan."
                        ),
                    },
                    "investigate_existing": {
                        "type": "boolean",
                        "default": False,
                        "description": "Run CLI-native fallback investigation for nonstandard existing-result layouts.",
                    },
                    "response_detail": {
                        "type": "string",
                        "enum": ["summary", "structured"],
                        "default": "summary",
                        "description": (
                            "Return summary for normal agent display. Use structured only when you need "
                            "the full structured checklist and internal detection fields."
                        ),
                    },
                    "language": {
                        "type": "string",
                        "enum": ["en", "zh", "auto"],
                        "default": "en",
                        "description": (
                            "User-facing language for compact summary rendering. Use en for English requests "
                            "and zh for Chinese requests; auto falls back to English unless CJK task context is available."
                        ),
                    },
                },
                "required": ["path"],
                "additionalProperties": False,
            },
            "annotations": {
                "toolType": "agent_helper",
                "returnsStructuredPayload": True,
                "recommendedUse": "Use first for path-driven tasks before choosing a native workflow.",
            },
        },
        {
            "name": "recommend_workflow",
            "title": "Recommend the next native workflow",
            "description": (
                "Classify an input path and task, then recommend the most appropriate native iobrpy "
                "workflow along with confirm-parameter hints."
            ),
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Input path to inspect.",
                    },
                    "task": {
                        "type": "string",
                        "default": "",
                        "description": "Natural-language task description.",
                    },
                },
                "anyOf": [{"required": ["path"]}, {"required": ["task"]}],
                "additionalProperties": False,
            },
            "annotations": {
                "toolType": "agent_helper",
                "returnsStructuredPayload": True,
                "recommendedUse": "Use after map_path when the agent needs a workflow recommendation.",
            },
        },
        {
            "name": "doctor_environment",
            "title": "Inspect the environment and preferred entrypoints",
            "description": (
                "Return environment diagnostics, installed entrypoints, optional external-tool version info, "
                "and the preferred workflow recommendation for the current path or task."
            ),
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Optional input path to inspect.",
                    },
                    "task": {
                        "type": "string",
                        "default": "",
                        "description": "Optional natural-language task description.",
                    },
                    "external": {
                        "type": "boolean",
                        "default": False,
                        "description": "When true, inspect external tools such as salmon, star, and trust4.",
                    },
                },
                "additionalProperties": False,
            },
            "annotations": {
                "toolType": "agent_helper",
                "returnsStructuredPayload": True,
                "recommendedUse": "Use only when environment diagnostics or entrypoint discovery are actually needed.",
            },
        },
        {
            "name": "analyze_path_task",
            "title": "Analyze a path and task for agent planning",
            "description": (
                "Combine path mapping and workflow recommendation into one structured payload so an LLM can "
                "reason over current stage, next choices, confirm-parameters, missing inputs, richer "
                "decision paths, contextual command hints, and a second-stage ambiguity-resolution surface "
                "for directory recognition. Use this after map_path when richer planning is still needed, "
                "not as the first scan on a large directory and not just to render a checklist. By default "
                "the returned payload is a compact summary so agents should not read client tool-result files "
                "or run JSON parser shell commands."
            ),
            "inputSchema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "File or directory to analyze.",
                    },
                    "task": {
                        "type": "string",
                        "default": "",
                        "description": "Natural-language task description to pair with the detected path state.",
                    },
                    "quick": {
                        "type": "boolean",
                        "default": False,
                        "description": (
                            "Use quick scan mode for the embedded path map only when the directory layout "
                            "is simple and a fast first pass is more important than completeness. Quick "
                            "scan can miss or misjudge information in complex or mixed-result directories."
                        ),
                    },
                    "max_depth": {
                        "type": "integer",
                        "default": 8,
                        "description": "Maximum directory depth to inspect during the embedded path map.",
                    },
                    "max_entries": {
                        "type": "integer",
                        "default": 8000,
                        "description": "Maximum file or directory names to scan during the embedded path map. Use 0 to disable the entry cap for a focused deep scan.",
                    },
                    "focus": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "Optional focus roots for the embedded path map. Prefer these when a directory "
                            "is structurally complex and you want deeper inspection of likely result roots."
                        ),
                    },
                    "investigate_existing": {
                        "type": "boolean",
                        "default": False,
                        "description": "Run existing-result investigation during the embedded path map.",
                    },
                    "response_detail": {
                        "type": "string",
                        "enum": ["summary", "structured"],
                        "default": "summary",
                        "description": (
                            "Return summary for normal agent display. Use structured only when richer "
                            "planning fields are genuinely needed."
                        ),
                    },
                    "language": {
                        "type": "string",
                        "enum": ["en", "zh", "auto"],
                        "default": "auto",
                        "description": (
                            "User-facing language for compact summary rendering. Use en for English requests, "
                            "zh for Chinese requests, or auto to infer from the task text."
                        ),
                    },
                },
                "required": ["path"],
                "additionalProperties": False,
            },
            "annotations": {
                "toolType": "agent_helper",
                "returnsStructuredPayload": True,
                "recommendedUse": "Use only after map_path when the scan result still needs richer planning; do not call it just to render workflow_checklist.",
            },
        },
    ]


def get_native_tools() -> List[Dict[str, Any]]:
    global _NATIVE_TOOL_CACHE
    if _NATIVE_TOOL_CACHE is None:
        _NATIVE_TOOL_CACHE = build_tools_from_parser(_capture_native_parser())
    return _NATIVE_TOOL_CACHE


def get_helper_tools() -> List[Dict[str, Any]]:
    global _HELPER_TOOL_CACHE
    if _HELPER_TOOL_CACHE is None:
        _HELPER_TOOL_CACHE = _build_helper_tools()
    return _HELPER_TOOL_CACHE


def get_tools() -> List[Dict[str, Any]]:
    global _TOOL_CACHE
    if _TOOL_CACHE is None:
        _TOOL_CACHE = [*get_helper_tools(), *get_native_tools()]
    return _TOOL_CACHE


def get_command_metadata(command_name: str) -> Dict[str, Any]:
    metadata = _build_command_metadata(_capture_native_parser())
    if command_name not in metadata:
        raise KeyError(command_name)
    return metadata[command_name]


def render_command_help(command_name: str) -> str:
    meta = get_command_metadata(command_name)
    declared_specs: List[ArgumentSpec] = meta["declared_specs"]
    required_specs = [spec for spec in declared_specs if spec.required]
    optional_specs = [spec for spec in declared_specs if not spec.required]

    usage_parts = [f"iobrpy-cli {command_name}"]
    for spec in declared_specs:
        label = spec.preferred_option or spec.property_name
        if spec.is_positional:
            label = spec.property_name.upper()
        if spec.is_flag:
            token = label
        elif spec.is_array:
            token = f"{label} ..."
        else:
            token = f"{label} VALUE"
        usage_parts.append(token if spec.required else f"[{token}]")
    usage = " ".join(usage_parts)

    lines = [f"Usage: {usage}", ""]
    lines.append(str(meta["summary"]).strip())
    if meta["detailed_description"]:
        lines.append("")
        lines.append(str(meta["detailed_description"]).strip())
    lines.append("")

    if meta["input_expectations"]:
        lines.append("Input Expectations:")
        for item in meta["input_expectations"]:
            lines.append(f"  - {item}")
        lines.append("")

    if meta["outputs"]:
        lines.append("Produces:")
        for item in meta["outputs"]:
            lines.append(f"  - {item}")
        lines.append("")

    if meta["when_to_use"]:
        lines.append("Use This Command When:")
        for item in meta["when_to_use"]:
            lines.append(f"  - {item}")
        lines.append("")

    if meta["use_another_command_when"]:
        lines.append("Prefer Another Command When:")
        for item in meta["use_another_command_when"]:
            lines.append(f"  - {item}")
        lines.append("")

    if meta.get("required_one_of"):
        lines.append("One Of These Parameter Groups Is Required:")
        for group in meta["required_one_of"]:
            rendered = ", ".join(group)
            lines.append(f"  - {rendered}")
        lines.append("")

    if meta["confirm_parameters"]:
        lines.append("Parameters To Confirm:")
        for item in meta["confirm_parameters"]:
            lines.append(f"  {item['display_name']}")
            if item["description"]:
                lines.append(f"    {item['description']}")
        lines.append("")

    if required_specs:
        lines.append("Required Parameters:")
        for spec in required_specs:
            flag_display = ", ".join(spec.option_strings) if spec.option_strings else spec.property_name
            lines.append(f"  {flag_display}")
            if spec.help_text:
                lines.append(f"    {spec.help_text}")
            choices = _schema_choices(spec)
            if choices:
                lines.append(f"    Choices: {', '.join(choices)}")
            if "default" in spec.schema:
                lines.append(f"    Default: {spec.schema['default']}")
        lines.append("")

    if optional_specs:
        lines.append("Optional Parameters:")
        for spec in optional_specs:
            flag_display = ", ".join(spec.option_strings) if spec.option_strings else spec.property_name
            lines.append(f"  {flag_display}")
            if spec.help_text:
                lines.append(f"    {spec.help_text}")
            choices = _schema_choices(spec)
            if choices:
                lines.append(f"    Choices: {', '.join(choices)}")
            if "default" in spec.schema:
                lines.append(f"    Default: {spec.schema['default']}")
        lines.append("")

    if meta["notes"]:
        lines.append("Notes:")
        for item in meta["notes"]:
            lines.append(f"  - {item}")
        lines.append("")

    if meta["source_paths"]:
        lines.append("Semantic Sources:")
        for path in meta["source_paths"]:
            lines.append(f"  - {path}")

    return "\n".join(lines).rstrip() + "\n"


def _serialize_value(value: Any, spec: ArgumentSpec) -> str:
    if spec.is_boolean_value:
        return "true" if bool(value) else "false"
    return str(value)


def validate_native_cli_argv(command_name: str, argv: Sequence[str]) -> None:
    meta = get_command_metadata(command_name)
    option_lookup: Dict[str, ArgumentSpec] = {}
    for spec in meta["argument_specs"]:
        for option in spec.option_strings:
            option_lookup[option] = spec

    provided: set[str] = set()
    index = 0
    args = list(argv)
    while index < len(args):
        token = args[index]
        if token in {"-h", "--help"}:
            return
        if not token.startswith("-"):
            raise ValueError(f"Unexpected positional argument for {command_name}: {token}")

        spec = option_lookup.get(token)
        if spec is None:
            allowed = ", ".join(sorted(option_lookup))
            raise ValueError(f"Unsupported option for {command_name}: {token}. Allowed options: {allowed}")

        provided.add(spec.property_name)
        if spec.is_flag:
            index += 1
            continue

        if spec.is_array:
            values: List[str] = []
            index += 1
            while index < len(args) and not args[index].startswith("-"):
                values.append(args[index])
                index += 1
            if not values:
                raise ValueError(f"Option {token} requires at least one value.")
            choices = _schema_choices(spec)
            if choices:
                invalid = [value for value in values if value not in choices]
                if invalid:
                    raise ValueError(
                        f"Unsupported value(s) for {token}: {', '.join(invalid)}. Choices: {', '.join(choices)}"
                    )
            continue

        index += 1
        if index >= len(args):
            raise ValueError(f"Option {token} requires a value.")
        value = args[index]
        choices = _schema_choices(spec)
        if choices and value not in choices:
            raise ValueError(f"Unsupported value for {token}: {value}. Choices: {', '.join(choices)}")
        index += 1

    missing_required = [
        name for name in meta.get("required_parameter_names", [])
        if name not in provided
    ]
    if missing_required:
        raise ValueError(f"Missing required parameters for {command_name}: {', '.join(missing_required)}")

    required_one_of = meta.get("required_one_of", [])
    if required_one_of:
        group_satisfied = any(all(name in provided for name in group) for group in required_one_of if group)
        if not group_satisfied:
            rendered = " | ".join(",".join(group) for group in required_one_of if group)
            raise ValueError(f"{command_name} requires one of these parameter groups: {rendered}")


def _build_cli_argv(command_name: str, arguments: Dict[str, Any]) -> List[str]:
    command_specs = _build_command_argument_specs(_capture_native_parser())
    if command_name not in command_specs:
        raise ValueError(f"Unknown or unsupported iobrpy command: {command_name}")

    allowed_keys = {spec.property_name for spec in command_specs[command_name]}
    unknown_keys = sorted(set(arguments) - allowed_keys)
    if unknown_keys:
        raise ValueError(f"Unsupported arguments for {command_name}: {', '.join(unknown_keys)}")

    meta = get_command_metadata(command_name)
    missing_required = [
        name for name in meta.get("required_parameter_names", [])
        if arguments.get(name) in (None, "", [])
    ]
    if missing_required:
        raise ValueError(f"Missing required arguments for {command_name}: {', '.join(missing_required)}")

    required_one_of = meta.get("required_one_of", [])
    if required_one_of:
        group_satisfied = any(
            all(arguments.get(name) not in (None, "", []) for name in group)
            for group in required_one_of
            if group
        )
        if not group_satisfied:
            rendered = " | ".join(",".join(group) for group in required_one_of if group)
            raise ValueError(f"{command_name} requires one of these parameter groups: {rendered}")

    argv: List[str] = [command_name]
    for spec in command_specs[command_name]:
        if spec.property_name not in arguments:
            continue

        value = arguments[spec.property_name]
        if value is None:
            continue

        if spec.is_flag:
            if value:
                if spec.preferred_option is None:
                    raise ValueError(f"Flag argument {spec.property_name} has no CLI option string.")
                argv.append(spec.preferred_option)
            continue

        if spec.is_positional:
            if spec.is_array:
                if not isinstance(value, list):
                    raise ValueError(f"Argument {spec.property_name} must be a list.")
                argv.extend(_serialize_value(item, spec) for item in value)
            else:
                argv.append(_serialize_value(value, spec))
            continue

        if spec.preferred_option is None:
            raise ValueError(f"Argument {spec.property_name} has no CLI option string.")
        argv.append(spec.preferred_option)
        if spec.is_array:
            if not isinstance(value, list):
                raise ValueError(f"Argument {spec.property_name} must be a list.")
            argv.extend(_serialize_value(item, spec) for item in value)
        else:
            argv.append(_serialize_value(value, spec))
    return argv


def _normalize_suffixes(path: Path) -> List[str]:
    return [suffix.lower() for suffix in path.suffixes]


def _is_fastq_file(path: Path) -> bool:
    suffixes = _normalize_suffixes(path)
    if suffixes[-2:] in ([".fastq", ".gz"], [".fq", ".gz"]):
        return True
    if suffixes[-1:] in ([".fastq"], [".fq"]):
        return True
    return False


def _is_bam_like(path: Path) -> bool:
    return path.suffix.lower() in {".bam", ".cram"}


def _looks_like_matrix_file(path: Path) -> bool:
    return path.suffix.lower() in {".csv", ".tsv", ".txt", ".gz"}


def _summarize_fastq_dir(path: Path) -> Dict[str, Any]:
    fastq_files = sorted(p for p in path.iterdir() if p.is_file() and _is_fastq_file(p))
    r1_files = [
        p for p in fastq_files
        if any(tag in p.name for tag in ("_R1", "_1.", "_1_", ".R1.", ".1.fastq", ".1.fq"))
    ]
    r2_files = [
        p for p in fastq_files
        if any(tag in p.name for tag in ("_R2", "_2.", "_2_", ".R2.", ".2.fastq", ".2.fq"))
    ]
    return {
        "kind": "fastq_directory",
        "fastq_count": len(fastq_files),
        "r1_count": len(r1_files),
        "r2_count": len(r2_files),
        "paired": bool(r1_files and len(r1_files) == len(r2_files)),
        "examples": [p.name for p in fastq_files[:5]],
    }


def _summarize_bam_dir(path: Path) -> Dict[str, Any]:
    bam_files = sorted(p for p in path.iterdir() if p.is_file() and _is_bam_like(p))
    return {
        "kind": "bam_directory",
        "bam_count": len(bam_files),
        "examples": [p.name for p in bam_files[:5]],
    }


def _summarize_matrix_file(path: Path) -> Dict[str, Any]:
    summary: Dict[str, Any] = {
        "kind": "expression_matrix_file",
        "path": str(path),
        "size_bytes": path.stat().st_size,
    }
    try:
        import pandas as pd

        read_kwargs: Dict[str, Any] = {"nrows": 5}
        lower_name = path.name.lower()
        if lower_name.endswith(".csv"):
            read_kwargs["sep"] = ","
        elif lower_name.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")):
            read_kwargs["sep"] = "\t"
        else:
            read_kwargs["sep"] = None
            read_kwargs["engine"] = "python"
        preview = pd.read_csv(path, **read_kwargs)
        summary["columns"] = list(preview.columns[:8])
        summary["preview_rows"] = int(preview.shape[0])
    except Exception as exc:
        summary["preview_error"] = str(exc)
    return summary


def _detect_input_summary(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {"kind": "missing", "path": str(path)}

    if path.is_dir():
        fastq_summary = _summarize_fastq_dir(path)
        if fastq_summary["fastq_count"] > 0:
            return fastq_summary

        bam_summary = _summarize_bam_dir(path)
        if bam_summary["bam_count"] > 0:
            return bam_summary

        return {
            "kind": "directory",
            "entry_count": len(list(path.iterdir())),
            "path": str(path),
        }

    if path.is_file() and _looks_like_matrix_file(path):
        return _summarize_matrix_file(path)

    return {
        "kind": "file",
        "path": str(path),
        "size_bytes": path.stat().st_size if path.exists() else None,
    }


def _safe_command_metadata_lookup(command_name: str) -> Dict[str, Any]:
    try:
        return get_command_metadata(command_name)
    except Exception:
        return {"important_optional_parameters": [], "confirm_parameters": [], "required_parameter_names": []}


def _native_command_summary_payload() -> List[Dict[str, Any]]:
    summaries: List[Dict[str, Any]] = []
    for tool in get_native_tools():
        meta = get_command_metadata(tool["name"])
        schema = tool.get("inputSchema", {})
        properties = schema.get("properties", {})
        required = list(schema.get("required", []))
        optional = [name for name in properties.keys() if name not in required]
        summaries.append(
            {
                "name": tool["name"],
                "title": tool.get("title"),
                "description": tool.get("description"),
                "summary": meta.get("summary", ""),
                "detailed_description": meta.get("detailed_description", ""),
                "input_expectations": meta.get("input_expectations", []),
                "outputs": meta.get("outputs", []),
                "when_to_use": meta.get("when_to_use", []),
                "use_another_command_when": meta.get("use_another_command_when", []),
                "confirm_parameters": meta.get("confirm_parameters", []),
                "important_optional_parameters": meta.get("important_optional_parameters", []),
                "required_args": required,
                "optional_args": optional,
                "required_one_of": meta.get("required_one_of", []),
                "routed_optional_args": [],
                "default_behaviors": meta.get("default_behaviors", []),
                "notes": meta.get("notes", []),
                "supports_extra_args": False,
            }
        )
    return summaries


def _recommend_for_summary(summary: Dict[str, Any], task: str) -> Dict[str, Any]:
    task_text = (task or "").lower()
    kind = summary.get("kind")

    if kind == "fastq_directory":
        meta = _safe_command_metadata_lookup("runall")
        known = {"fastq"}
        missing = [
            name for name in meta.get("required_parameter_names", [])
            if name not in known
        ]
        if "star" in task_text:
            mode = "star"
        elif "salmon" in task_text:
            mode = "salmon"
        else:
            mode = None
        if mode and "mode" in missing:
            missing.remove("mode")
        return {
            "recommended_command": "runall",
            "why": (
                "Detected a FASTQ directory. According to RAG_MCP/iobrpy_required_params.json, "
                "`runall` is the workflow entrypoint for this kind of input."
            ),
            "suggested_cli": "iobrpy runall --fastq <path_fastq_dir> --outdir <path_output_dir> --mode <salmon|star> --index <path_salmon_index|path_star_index> --threads <threads> --batch_size <batch_size> --project <project_name>",
            "suggested_harness_cli": "iobrpy-cli runall --fastq <path_fastq_dir> --outdir <path_output_dir> --mode <salmon|star> --index <path_salmon_index|path_star_index> --threads <threads> --batch_size <batch_size> --project <project_name>",
            "assumptions": [f"Detected {summary.get('fastq_count', 0)} FASTQ files"],
            "missing_parameters": missing,
            "confirm_parameters": meta.get("confirm_parameters", []),
            "important_optional_parameters": meta.get("important_optional_parameters", []),
            "confirmation_guidance": (
                "Confirm the parameters listed in recommendation.confirm_parameters. The command surface "
                "is intentionally constrained to RAG_MCP/iobrpy_required_params.json."
            ),
            "notes": [
                "The FASTQ path itself is already known from the detected input.",
            ],
        }

    if kind == "expression_matrix_file":
        meta = _safe_command_metadata_lookup("tme_profile")
        return {
            "recommended_command": "tme_profile",
            "why": (
                "Detected an expression matrix file. According to RAG_MCP/iobrpy_required_params.json, "
                "`tme_profile` is the matching workflow for an existing expression matrix."
            ),
            "suggested_cli": "iobrpy tme_profile --input <path_tpm_matrix> --output <path_tme_profile_outdir> --threads <threads>",
            "suggested_harness_cli": "iobrpy-cli tme_profile --input <path_tpm_matrix> --output <path_tme_profile_outdir> --threads <threads>",
            "assumptions": ["The matrix is genes x samples and TPM-like."],
            "missing_parameters": [
                name for name in meta.get("required_parameter_names", [])
                if name not in {"input"}
            ],
            "confirm_parameters": meta.get("confirm_parameters", []),
            "important_optional_parameters": meta.get("important_optional_parameters", []),
            "confirmation_guidance": (
                "Use the JSON-defined parameter set only. For tme_profile, the matrix path is already known from the detected input."
            ),
            "notes": [
                "If the matrix is counts instead of TPM-like values, use count2tpm first.",
            ],
        }

    if kind == "bam_directory":
        meta = _safe_command_metadata_lookup("hla_typing")
        return {
            "recommended_command": "hla_typing",
            "why": (
                "Detected BAM/CRAM files. According to RAG_MCP/iobrpy_required_params.json, "
                "`hla_typing` is the batch BAM-directory workflow for this input."
            ),
            "suggested_cli": "iobrpy hla_typing -b <path_bam_dir> -r <hg19|hg38> -o <path_hla_outdir> -j <threads>",
            "suggested_harness_cli": "iobrpy-cli hla_typing -b <path_bam_dir> -r <hg19|hg38> -o <path_hla_outdir> -j <threads>",
            "assumptions": [f"Detected {summary.get('bam_count', 0)} BAM/CRAM files"],
            "missing_parameters": [
                name for name in meta.get("required_parameter_names", [])
                if name not in {"bam_dir"}
            ],
            "confirm_parameters": meta.get("confirm_parameters", []),
            "important_optional_parameters": meta.get("important_optional_parameters", []),
            "confirmation_guidance": (
                "Use only the parameters declared in RAG_MCP/iobrpy_required_params.json for hla_typing."
            ),
            "notes": [
                "If the user wants repertoire analysis instead, trust4 may be the right command.",
            ],
        }

    return {
        "recommended_command": None,
        "why": "The input could not be classified confidently.",
        "suggested_cli": "iobrpy-cli commands --json",
        "suggested_harness_cli": "iobrpy-cli commands --json",
        "assumptions": [],
        "missing_parameters": [],
        "confirm_parameters": [],
        "important_optional_parameters": [],
        "confirmation_guidance": (
            "Classify the input more precisely before deciding which parameters need confirmation."
        ),
        "notes": [
            "Inspect the path with doctor_environment, map_path, or list_native_commands.",
        ],
    }


def _console_script_names() -> List[str]:
    names = sorted({
        entry.name
        for entry in importlib_metadata.entry_points(group="console_scripts")
        if entry.name in {"iobrpy", "iobrpy-cli", "iobrpy-cli-mcp"}
    })
    return list(names)


def _command_on_path(name: str) -> Optional[str]:
    return shutil.which(name)


def _distribution_version(name: str) -> Optional[str]:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return None


def _workflow_version_info(external: bool) -> Dict[str, Any]:
    if external:
        from iobrpy_cli.iobrpy.core import QuantificationWorkflow

        result = QuantificationWorkflow().get_version_info(check_external=True)
        return {
            "iobrpy_version": result.iobrpy_version,
            "iobrpy_installed": result.iobrpy_installed,
            "external_tools": result.external_tools or {},
            "python_version": result.python_version,
            "message": result.message,
        }

    try:
        return {
            "iobrpy_version": importlib_metadata.version("iobrpy"),
            "iobrpy_installed": True,
            "python_version": platform.python_version(),
        }
    except importlib_metadata.PackageNotFoundError:
        return {
            "iobrpy_version": None,
            "iobrpy_installed": False,
            "python_version": platform.python_version(),
        }


def _helper_next_steps() -> List[str]:
    return [
        "For a user-provided dataset path, start with map_path using its default full/deep scan behavior rather than asking the user to choose a scan mode.",
        "Do not ask quick-versus-full scan questions; map_path is optimized for complex directories by default.",
        "Do not use analyze_path_task just to render a checklist from a successful map_path response.",
        "Use analyze_path_task only after map_path when you still need richer planning beyond the scan payload.",
        "Never read client tool-result files or run parser shell commands just to extract fields from map_path or analyze_path_task JSON.",
        "Use recommend_workflow after mapping only when you still need a narrower command recommendation than analyze_path_task or map_path already provided.",
        "When a native command has resource-sensitive parameters such as threads, parallel_size, num_processes, or batch_size, inspect CPU cores, current load, and free memory on the execution host before choosing them.",
        "Before recommending, writing, or executing a native command invocation, validate the command and flags against RAG_MCP/iobrpy_required_params.json using list_native_commands or the native command help. Do not invent flags from general workflow concepts. In particular, tme_profile has no --mode; its JSON-defined parameters are --input, --output, and optional --threads.",
        "Run a native command tool only after the required parameters and user decisions are clear.",
    ]


def _decorate_helper_payload(tool_name: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(payload)
    out["tool"] = tool_name
    out["tool_type"] = "agent_helper"
    out.setdefault("entrypoint", "iobrpy-cli-mcp")
    out.setdefault("next_steps", _helper_next_steps())
    return out


def _limited_sequence(values: Any, *, limit: int) -> List[Any]:
    if not isinstance(values, list):
        return []
    return list(values[:limit])


def _trim_text(value: Any, *, limit: int = 280) -> Any:
    if not isinstance(value, str):
        return value
    text = value.strip()
    if len(text) <= limit:
        return text
    return f"{text[: max(limit - 3, 0)].rstrip()}..."


def _contains_cjk_text(value: Any) -> bool:
    return any("\u4e00" <= ch <= "\u9fff" for ch in str(value or ""))


def _response_language(value: Any, *, context: Any = "") -> str:
    raw = str(value or "").strip().lower().replace("_", "-")
    if raw in {"zh", "zh-cn", "zh-hans", "chinese", "cn"}:
        return "zh"
    if raw in {"en", "en-us", "english"}:
        return "en"
    if raw == "auto" and _contains_cjk_text(context):
        return "zh"
    return "en"


def _compact_parameter_cards_for_transport(items: Any, *, limit: int = 3) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if isinstance(item, dict):
            card = {
                "parameter": item.get("parameter"),
                "display_name": item.get("display_name"),
                "description": _trim_text(item.get("description"), limit=160),
                "decision_type": item.get("decision_type"),
                "resource_sensitive": item.get("resource_sensitive"),
                "resource_probe_recommended": item.get("resource_probe_recommended"),
                "resource_probe_scope": item.get("resource_probe_scope"),
                "resource_probe_commands": _limited_sequence(item.get("resource_probe_commands"), limit=4),
                "natural_language_confirmation": _trim_text(
                    item.get("natural_language_confirmation"),
                    limit=120,
                ),
                "natural_language_confirmation_zh": _trim_text(
                    item.get("natural_language_confirmation_zh"),
                    limit=120,
                ),
            }
            compacted.append({key: value for key, value in card.items() if value not in (None, [], {})})
        elif isinstance(item, str):
            compacted.append({"parameter": item, "display_name": item})
    return compacted


def _compact_numbered_options_for_transport(items: Any, *, limit: int = 4) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "recommended": item.get("recommended"),
            }
        )
    return compacted


def _compact_concise_confirmation_prompts_for_transport(items: Any, *, limit: int = 2) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "style": item.get("style"),
                "brief_prompt": _trim_text(item.get("brief_prompt"), limit=140),
                "brief_prompt_zh": _trim_text(item.get("brief_prompt_zh"), limit=80),
                "numbered_options": _compact_numbered_options_for_transport(item.get("numbered_options"), limit=4),
                "recommended_default": _trim_text(item.get("recommended_default"), limit=120),
                "avoid_raw_flag_names": item.get("avoid_raw_flag_names"),
                "resource_probe_recommended": item.get("resource_probe_recommended"),
                "resource_probe_scope": item.get("resource_probe_scope"),
                "resource_probe_commands": _limited_sequence(item.get("resource_probe_commands"), limit=4),
            }
        )
    return compacted


def _compact_existing_input_candidates_for_transport(items: Any, *, limit: int = 4) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "source_id": item.get("source_id"),
                "source_label": item.get("source_label"),
                "label": item.get("label"),
                "why_it_matches": _trim_text(item.get("why_it_matches"), limit=180),
                "paths": _limited_sequence(item.get("paths"), limit=2),
                "preferred": item.get("preferred"),
                "needs_user_choice": item.get("needs_user_choice"),
                "selection_guidance": _trim_text(item.get("selection_guidance"), limit=260),
                "features_indexing_note": _trim_text(item.get("features_indexing_note"), limit=220),
                "suggested_id_column": item.get("suggested_id_column"),
            }
        )
    return compacted


def _compact_recovery_paths_for_transport(items: Any, *, limit: int = 1) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        steps: List[Dict[str, Any]] = []
        for step in item.get("steps", [])[:4]:
            if not isinstance(step, dict):
                continue
            steps.append(
                {
                    "command": step.get("command"),
                    "why": _trim_text(step.get("why"), limit=150),
                    "suggested_arguments": step.get("suggested_arguments"),
                    "expected_output": step.get("expected_output"),
                }
            )
        compacted.append(
            {
                "id": item.get("id"),
                "title": item.get("title"),
                "preferred": item.get("preferred"),
                "why_this_path_exists": _trim_text(item.get("why_this_path_exists"), limit=180),
                "detected_upstream_paths": _limited_sequence(item.get("detected_upstream_paths"), limit=2),
                "standardized_tpm_output": item.get("standardized_tpm_output"),
                "steps": steps,
            }
        )
    return compacted


def _compact_pre_execution_questions_for_transport(items: Any, *, limit: int = 2) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "question": _trim_text(item.get("question"), limit=160),
                "why_ask_now": _trim_text(item.get("why_ask_now"), limit=150),
                "recommended_default": _trim_text(item.get("recommended_default"), limit=120),
            }
        )
    return compacted


def _compact_command_card_for_transport(card: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "command": card.get("command"),
        "category": card.get("category"),
        "priority": card.get("priority"),
        "summary": _trim_text(card.get("summary"), limit=160),
        "why_recommended": _trim_text(card.get("why_recommended"), limit=180),
        "required_inputs": _compact_parameter_cards_for_transport(card.get("required_inputs"), limit=3),
        "input_expectations": [_trim_text(item, limit=140) for item in _limited_sequence(card.get("input_expectations"), limit=2)],
        "expected_outputs": [_trim_text(item, limit=140) for item in _limited_sequence(card.get("expected_outputs"), limit=2)],
        "confirm_parameters": _compact_parameter_cards_for_transport(card.get("confirm_parameters"), limit=2),
        "important_optional_parameters": _compact_parameter_cards_for_transport(card.get("important_optional_parameters"), limit=2),
        "project_context_notes": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("project_context_notes"), limit=2)],
        "existing_input_candidates": _compact_existing_input_candidates_for_transport(card.get("existing_input_candidates"), limit=4),
        "upstream_recovery_paths": _compact_recovery_paths_for_transport(card.get("upstream_recovery_paths"), limit=1),
        "detected_input_story": _trim_text(card.get("detected_input_story"), limit=180),
        "agent_reasoning_rules": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("agent_reasoning_rules"), limit=4)],
        "execution_guardrails": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("execution_guardrails"), limit=4)],
        "pre_execution_questions": _compact_pre_execution_questions_for_transport(card.get("pre_execution_questions"), limit=2),
        "concise_confirmation_prompts": _compact_concise_confirmation_prompts_for_transport(
            card.get("concise_confirmation_prompts"),
            limit=2,
        ),
        "mention_only_if_needed": card.get("mention_only_if_needed"),
        "relevance_trigger": card.get("relevance_trigger"),
    }


def _compact_command_guidance_bucket_for_transport(items: Any, *, limit: int) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(_compact_command_card_for_transport(item))
    return compacted


def _compact_reason_cards_for_transport(items: Any, *, limit: int = 2) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "detail": _trim_text(item.get("detail"), limit=150),
                "effect_on_recommendation": _trim_text(item.get("effect_on_recommendation"), limit=120),
            }
        )
    return compacted


def _compact_next_step_card_for_transport(card: Dict[str, Any]) -> Dict[str, Any]:
    base = {
        "id": card.get("id"),
        "kind": card.get("kind"),
        "priority": card.get("priority"),
        "posture": card.get("posture"),
        "title": card.get("title"),
        "why_this_is_recommended": _trim_text(card.get("why_this_is_recommended"), limit=180),
        "ranking_reasons": _compact_reason_cards_for_transport(card.get("ranking_reasons"), limit=2),
        "required_inputs": _compact_parameter_cards_for_transport(card.get("required_inputs"), limit=3),
        "input_expectations": [_trim_text(item, limit=140) for item in _limited_sequence(card.get("input_expectations"), limit=2)],
        "expected_result": _trim_text(card.get("expected_result"), limit=160),
        "confirm_parameters": _compact_parameter_cards_for_transport(card.get("confirm_parameters"), limit=2),
        "important_optional_parameters": _compact_parameter_cards_for_transport(card.get("important_optional_parameters"), limit=2),
        "existing_input_candidates": _compact_existing_input_candidates_for_transport(card.get("existing_input_candidates"), limit=4),
        "tradeoffs": [_trim_text(item, limit=140) for item in _limited_sequence(card.get("tradeoffs"), limit=2)],
        "when_not_to_choose": _trim_text(card.get("when_not_to_choose"), limit=150),
    }
    if card.get("kind") == "command":
        base.update(
            {
                "command": card.get("command"),
                "upstream_recovery_paths": _compact_recovery_paths_for_transport(card.get("upstream_recovery_paths"), limit=1),
                "detected_input_story": _trim_text(card.get("detected_input_story"), limit=180),
                "agent_reasoning_rules": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("agent_reasoning_rules"), limit=4)],
                "execution_guardrails": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("execution_guardrails"), limit=4)],
                "pre_execution_questions": _compact_pre_execution_questions_for_transport(card.get("pre_execution_questions"), limit=2),
                "concise_confirmation_prompts": _compact_concise_confirmation_prompts_for_transport(
                    card.get("concise_confirmation_prompts"),
                    limit=2,
                ),
                "project_context_notes": [_trim_text(item, limit=150) for item in _limited_sequence(card.get("project_context_notes"), limit=2)],
                "category": card.get("category"),
            }
        )
    return {key: value for key, value in base.items() if value not in (None, [], {})}


def _compact_next_step_cards_for_transport(items: Any, *, limit: int = 6) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(_compact_next_step_card_for_transport(item))
    return compacted


def _compact_decision_paths_for_transport(items: Any, *, limit: int = 4) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "why_this_path_exists": _trim_text(item.get("why_this_path_exists"), limit=160),
                "recommended_for_current_state": item.get("recommended_for_current_state"),
                "commands": _compact_command_guidance_bucket_for_transport(item.get("commands"), limit=2),
            }
        )
    return compacted


def _compact_directory_recognition_payload_for_transport(payload: Dict[str, Any]) -> Dict[str, Any]:
    layer1 = payload.get("layer_1_cli_scan", {})
    layer2 = payload.get("layer_2_llm_reasoning_surface", {})
    compacted = {
        "layer_1_cli_scan": {
            "owner": layer1.get("owner"),
            "current_stage": layer1.get("current_stage"),
            "current_label": layer1.get("current_label"),
            "current_stage_confidence": layer1.get("current_stage_confidence"),
            "scan_warning_present": layer1.get("scan_warning_present"),
            "agent_fallback_recommended": layer1.get("agent_fallback_recommended"),
            "existing_result_investigation_available": layer1.get("existing_result_investigation_available"),
        },
        "layer_2_llm_reasoning_surface": {
            "enabled": layer2.get("enabled"),
            "should_expand_reasoning": layer2.get("should_expand_reasoning"),
            "reason_details": _limited_sequence(layer2.get("reason_details"), limit=8),
            "evidence_samples": _limited_sequence(layer2.get("evidence_samples"), limit=4),
            "candidate_function_hypotheses": _compact_command_guidance_bucket_for_transport(
                layer2.get("candidate_function_hypotheses"),
                limit=8,
            ),
            "external_result_interpretation_hints": _limited_sequence(
                layer2.get("external_result_interpretation_hints"),
                limit=4,
            ),
            "suggested_user_questions": _limited_sequence(layer2.get("suggested_user_questions"), limit=8),
            "response_guidance": layer2.get("response_guidance"),
        },
    }
    return compacted


def _compact_current_state_assessment_for_transport(assessment: Dict[str, Any]) -> Dict[str, Any]:
    compacted = dict(assessment)
    compacted["completed_labels"] = _limited_sequence(compacted.get("completed_labels"), limit=5)
    compacted["partial_labels"] = _limited_sequence(compacted.get("partial_labels"), limit=4)
    compacted["pending_labels"] = _limited_sequence(compacted.get("pending_labels"), limit=4)
    compacted["external_labels"] = _limited_sequence(compacted.get("external_labels"), limit=3)
    compacted["reusable_path_examples"] = _limited_sequence(compacted.get("reusable_path_examples"), limit=3)
    compacted["profile_tags"] = _limited_sequence(compacted.get("profile_tags"), limit=6)
    compacted["default_posture_reason_ids"] = _limited_sequence(compacted.get("default_posture_reason_ids"), limit=6)
    compacted["summary_signals"] = _limited_sequence(compacted.get("summary_signals"), limit=8)
    return compacted


def _compact_next_step_reasoning_surface_for_transport(surface: Dict[str, Any]) -> Dict[str, Any]:
    compacted = dict(surface)
    compacted["ranking_criteria"] = _limited_sequence(compacted.get("ranking_criteria"), limit=6)
    compacted["option_ids_in_priority_order"] = _limited_sequence(compacted.get("option_ids_in_priority_order"), limit=8)
    return compacted


def _compact_concise_confirmation_plan_for_transport(plan: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "prefer_brief_numbered_questions": plan.get("prefer_brief_numbered_questions"),
        "ask_at_most_this_many_items_first": plan.get("ask_at_most_this_many_items_first"),
        "opening": _trim_text(plan.get("opening"), limit=120),
        "opening_zh": _trim_text(plan.get("opening_zh"), limit=80),
        "items": _limited_sequence(plan.get("items"), limit=2),
    }


def _compact_recommendation_for_transport(recommendation_payload: Dict[str, Any]) -> Dict[str, Any]:
    recommendation = recommendation_payload.get("recommendation", {})
    if not isinstance(recommendation, dict):
        recommendation = {}
    return {
        "success": recommendation_payload.get("success"),
        "task": recommendation_payload.get("task"),
        "detected_input": recommendation_payload.get("detected_input"),
        "path_state": recommendation_payload.get("path_state"),
        "recommendation": {
            "recommended_command": recommendation.get("recommended_command"),
            "why": _trim_text(recommendation.get("why"), limit=180),
            "suggested_harness_cli": recommendation.get("suggested_harness_cli"),
            "missing_parameters": _limited_sequence(recommendation.get("missing_parameters"), limit=6),
            "confirm_parameters": _compact_parameter_cards_for_transport(recommendation.get("confirm_parameters"), limit=2),
            "important_optional_parameters": _compact_parameter_cards_for_transport(
                recommendation.get("important_optional_parameters"),
                limit=2,
            ),
            "notes": [_trim_text(item, limit=140) for item in _limited_sequence(recommendation.get("notes"), limit=3)],
        },
    }


def _compact_scan_payload(scan: Any) -> Dict[str, Any]:
    if not isinstance(scan, dict):
        return {}

    return {
        "entry_count": scan.get("entry_count"),
        "truncated": scan.get("truncated"),
        "max_depth": scan.get("max_depth"),
        "max_entries": scan.get("max_entries"),
        "entry_limit_enabled": scan.get("entry_limit_enabled"),
        "entry_limit_hit": scan.get("entry_limit_hit"),
        "backlog_limit_hit": scan.get("backlog_limit_hit"),
        "unbounded_entry_scan": scan.get("unbounded_entry_scan"),
        "truncated_reason_ids": _limited_sequence(scan.get("truncated_reason_ids"), limit=8),
        "depth_limited_dir_count": scan.get("depth_limited_dir_count"),
        "depth_limited_dirs": _limited_sequence(scan.get("depth_limited_dirs"), limit=8),
        "strategy": scan.get("strategy"),
        "quick_skipped_generic_dir_count": scan.get("quick_skipped_generic_dir_count"),
        "sampled_bulk_stage_dir_count": scan.get("sampled_bulk_stage_dir_count"),
        "sampled_bulk_child_dir_count": scan.get("sampled_bulk_child_dir_count"),
        "sampled_bulk_file_count": scan.get("sampled_bulk_file_count"),
        "full_skipped_generic_dir_count": scan.get("full_skipped_generic_dir_count"),
        "focus_mode": scan.get("focus_mode"),
        "focus_roots": _limited_sequence(scan.get("focus_roots"), limit=16),
        "explicit_focus_roots": _limited_sequence(scan.get("explicit_focus_roots"), limit=16),
        "auto_focus_roots": _limited_sequence(scan.get("auto_focus_roots"), limit=16),
        "promoted_focus_root_count": scan.get("promoted_focus_root_count"),
        "generic_deep_scan_limit": scan.get("generic_deep_scan_limit"),
    }


def _compact_external_analysis_hints(hints: Any) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(hints, list):
        return compacted
    for item in hints[:8]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "description": item.get("description"),
                "description_zh": item.get("description_zh"),
                "evidence": _limited_sequence(item.get("evidence"), limit=4),
                "analysis_family_id": item.get("analysis_family_id"),
                "analysis_family_label": item.get("analysis_family_label"),
                "confidence": item.get("confidence"),
                "inference_basis": _limited_sequence(item.get("inference_basis"), limit=4),
                "matched_signal_labels": _limited_sequence(item.get("matched_signal_labels"), limit=4),
                "result_kind_labels": _limited_sequence(item.get("result_kind_labels"), limit=3),
                "likely_tool_families": _limited_sequence(item.get("likely_tool_families"), limit=4),
                "family_summary": item.get("family_summary"),
                "evidence_profiles": [
                    {
                        "path": profile.get("path"),
                        "confidence": profile.get("confidence"),
                        "result_kind_label": profile.get("result_kind_label"),
                        "matched_signal_labels": _limited_sequence(profile.get("matched_signal_labels"), limit=3),
                    }
                    for profile in _limited_sequence(item.get("evidence_profiles"), limit=3)
                    if isinstance(profile, dict)
                ],
            }
        )
    return compacted


def _compact_function_detections(detections: Any) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(detections, list):
        return compacted
    for item in detections[:24]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "status": item.get("status"),
                "status_label": item.get("status_label"),
                "status_label_zh": item.get("status_label_zh"),
                "stage_ids": _limited_sequence(item.get("stage_ids"), limit=8),
                "note": item.get("note"),
                "note_zh": item.get("note_zh"),
                "evidence": _limited_sequence(item.get("evidence"), limit=4),
            }
        )
    return compacted


def _compact_existing_result_investigation(investigation: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(investigation, dict):
        return None

    target_summaries: List[Dict[str, Any]] = []
    for item in investigation.get("target_summaries", [])[:12]:
        if not isinstance(item, dict):
            continue
        target_summaries.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "finding_count": item.get("finding_count"),
                "bucket_counts": item.get("bucket_counts", {}),
                "sample_paths": _limited_sequence(item.get("sample_paths"), limit=4),
                "external_result_families": _limited_sequence(item.get("external_result_families"), limit=4),
                "external_result_kind_labels": _limited_sequence(item.get("external_result_kind_labels"), limit=4),
                "external_likely_tool_families": _limited_sequence(item.get("external_likely_tool_families"), limit=4),
                "external_signal_labels": _limited_sequence(item.get("external_signal_labels"), limit=5),
                "external_evidence_profiles": [
                    {
                        "path": profile.get("path"),
                        "analysis_family_label": profile.get("analysis_family_label"),
                        "confidence": profile.get("confidence"),
                        "result_kind_label": profile.get("result_kind_label"),
                    }
                    for profile in _limited_sequence(item.get("external_evidence_profiles"), limit=3)
                    if isinstance(profile, dict)
                ],
            }
        )

    compacted = {"target_summaries": target_summaries}
    if "summary" in investigation:
        compacted["summary"] = investigation.get("summary")
    if "summary_zh" in investigation:
        compacted["summary_zh"] = investigation.get("summary_zh")
    return compacted


def _compact_agent_fallback(agent_fallback: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(agent_fallback, dict):
        return None

    investigation_targets: List[Dict[str, Any]] = []
    for item in agent_fallback.get("investigation_targets", [])[:12]:
        if not isinstance(item, dict):
            continue
        investigation_targets.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "why": item.get("why"),
                "why_zh": item.get("why_zh"),
                "sample_paths": _limited_sequence(item.get("sample_paths"), limit=4),
            }
        )

    return {
        "enabled": agent_fallback.get("enabled"),
        "recommended": agent_fallback.get("recommended"),
        "reason_ids": _limited_sequence(agent_fallback.get("reason_ids"), limit=8),
        "summary": agent_fallback.get("summary"),
        "summary_zh": agent_fallback.get("summary_zh"),
        "trigger_when": _limited_sequence(agent_fallback.get("trigger_when"), limit=6),
        "trigger_when_zh": _limited_sequence(agent_fallback.get("trigger_when_zh"), limit=6),
        "actions": _limited_sequence(agent_fallback.get("actions"), limit=8),
        "actions_zh": _limited_sequence(agent_fallback.get("actions_zh"), limit=8),
        "investigation_targets": investigation_targets,
    }


def _compact_workflow_checklist(checklist: Any) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(checklist, list):
        return compacted

    for item in checklist:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                "id": item.get("id"),
                "marker": item.get("marker"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "text": item.get("text"),
                "text_zh": item.get("text_zh"),
                "details": item.get("details"),
                "details_zh": item.get("details_zh"),
                "checked": item.get("checked"),
                "status": item.get("status"),
                "status_label": item.get("status_label"),
                "status_label_zh": item.get("status_label_zh"),
                "result_source": item.get("result_source"),
                "result_source_label": item.get("result_source_label"),
                "result_source_label_zh": item.get("result_source_label_zh"),
                "investigation_applied": item.get("investigation_applied"),
                "evidence_paths": _limited_sequence(item.get("evidence_paths"), limit=6),
                "evidence_display_paths": _limited_sequence(item.get("evidence_display_paths"), limit=6),
                "raw_evidence_paths": _limited_sequence(item.get("raw_evidence_paths"), limit=6),
                "detected_column_value": item.get("detected_column_value"),
                "detected_column_value_zh": item.get("detected_column_value_zh"),
                "detected_column_source": item.get("detected_column_source"),
                "missing_items": _limited_sequence(item.get("missing_items"), limit=6),
                "missing_items_zh": _limited_sequence(item.get("missing_items_zh"), limit=6),
                "missing_column_value": item.get("missing_column_value"),
                "missing_column_value_zh": item.get("missing_column_value_zh"),
                "missing_column_source": item.get("missing_column_source"),
            }
        )
    return compacted


def _summary_workflow_checklist(checklist: Any, *, language: str = "en") -> List[Dict[str, Any]]:
    language = _response_language(language)
    compacted: List[Dict[str, Any]] = []
    if not isinstance(checklist, list):
        return compacted

    for item in checklist:
        if not isinstance(item, dict):
            continue
        if language == "zh":
            stage = item.get("label_zh") or item.get("label") or item.get("id")
            status_label = item.get("status_label_zh") or item.get("status_label")
            result_source_label = item.get("result_source_label_zh") or item.get("result_source_label")
            missing_items = _limited_sequence(item.get("missing_items_zh") or item.get("missing_items"), limit=4)
        else:
            stage = item.get("label") or item.get("label_zh") or item.get("id")
            status_label = item.get("status_label") or item.get("status_label_zh")
            result_source_label = item.get("result_source_label") or item.get("result_source_label_zh")
            missing_items = _limited_sequence(item.get("missing_items") or item.get("missing_items_zh"), limit=4)
        evidence_paths = _limited_sequence(
            item.get("evidence_display_paths")
            or item.get("evidence_paths")
            or item.get("raw_evidence_paths"),
            limit=2,
        )
        compacted.append(
            {
                key: value
                for key, value in {
                    "id": item.get("id"),
                    "marker": item.get("marker"),
                    "status": item.get("status"),
                    "status_label": status_label,
                    "stage": stage,
                    "result_source": item.get("result_source"),
                    "result_source_label": result_source_label,
                    "evidence_display_paths": evidence_paths,
                    "missing_items": missing_items,
                }.items()
                if value not in (None, [], {})
            }
        )
    return compacted


def _compact_roadmap_payload(roadmap: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(roadmap, dict):
        return None
    return {
        "title": roadmap.get("title"),
        "url": roadmap.get("url"),
        "local_doc": roadmap.get("local_doc"),
        "current_node": roadmap.get("current_node"),
        "next_nodes": _limited_sequence(roadmap.get("next_nodes"), limit=8),
        "position_summary": roadmap.get("position_summary"),
        "position_summary_zh": roadmap.get("position_summary_zh"),
    }


def _compact_agent_rendering_hints(hints: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(hints, dict):
        return None
    compacted = {
        "confirm_scan_mode_before_first_map": hints.get("confirm_scan_mode_before_first_map"),
        "default_scan_mode_for_path_requests": hints.get("default_scan_mode_for_path_requests"),
        "always_use_full_scan_for_path_requests": hints.get("always_use_full_scan_for_path_requests"),
        "do_not_ask_user_to_choose_scan_mode": hints.get("do_not_ask_user_to_choose_scan_mode"),
        "quick_scan_deprecated_for_agent_use": hints.get("quick_scan_deprecated_for_agent_use"),
        "wait_for_scan_mode_answer_before_first_map": hints.get("wait_for_scan_mode_answer_before_first_map"),
        "do_not_start_quick_scan_while_awaiting_scan_mode_choice": hints.get(
            "do_not_start_quick_scan_while_awaiting_scan_mode_choice"
        ),
        "quick_planning_scan_only_after_user_selects_full_scan": hints.get(
            "quick_planning_scan_only_after_user_selects_full_scan"
        ),
        "run_quick_planning_scan_before_full_scan": hints.get("run_quick_planning_scan_before_full_scan"),
        "use_focus_roots_from_quick_scan_when_available": hints.get(
            "use_focus_roots_from_quick_scan_when_available"
        ),
        "preferred_map_command": hints.get("preferred_map_command"),
        "path_specific_recommended_deep_map_command": hints.get("path_specific_recommended_deep_map_command"),
        "path_specific_recommended_deep_map_command_uses_unbounded_entries": hints.get(
            "path_specific_recommended_deep_map_command_uses_unbounded_entries"
        ),
        "scan_retry_command": hints.get("scan_retry_command"),
        "scan_retry_limits": hints.get("scan_retry_limits"),
        "preferred_checklist_table_field": hints.get("preferred_checklist_table_field"),
        "preferred_checklist_table_field_zh": hints.get("preferred_checklist_table_field_zh"),
        "preferred_checklist_report_field": hints.get("preferred_checklist_report_field"),
        "preferred_checklist_report_field_zh": hints.get("preferred_checklist_report_field_zh"),
        "preferred_checklist_lines_field": hints.get("preferred_checklist_lines_field"),
        "preferred_checklist_lines_field_zh": hints.get("preferred_checklist_lines_field_zh"),
        "preferred_rendering_field_order": hints.get("preferred_rendering_field_order"),
        "preferred_rendering_field_order_zh": hints.get("preferred_rendering_field_order_zh"),
        "use_compact_default_user_facing_checklist_when_present": hints.get(
            "use_compact_default_user_facing_checklist_when_present"
        ),
        "reuse_prebuilt_checklist_table_verbatim": hints.get("reuse_prebuilt_checklist_table_verbatim"),
        "do_not_replace_prebuilt_table_with_custom_summary_table": hints.get(
            "do_not_replace_prebuilt_table_with_custom_summary_table"
        ),
        "do_not_rebuild_checklist_from_deconvolution_methods": hints.get(
            "do_not_rebuild_checklist_from_deconvolution_methods"
        ),
        "workflow_checklist_missing_items_are_authoritative": hints.get(
            "workflow_checklist_missing_items_are_authoritative"
        ),
        "do_not_use_default_bundle_missing_labels_for_checklist_missing_column": hints.get(
            "do_not_use_default_bundle_missing_labels_for_checklist_missing_column"
        ),
        "immune_deconvolution_missing_column_source": hints.get("immune_deconvolution_missing_column_source"),
        "deconvolution_methods_missing_labels_include_bayesprism": hints.get(
            "deconvolution_methods_missing_labels_include_bayesprism"
        ),
        "default_bundle_missing_labels_are_only_for_tme_profile_or_runall_scope": hints.get(
            "default_bundle_missing_labels_are_only_for_tme_profile_or_runall_scope"
        ),
        "bayesprism_is_deconvolution_but_standalone_optional": hints.get(
            "bayesprism_is_deconvolution_but_standalone_optional"
        ),
        "bayesprism_missing_display_en": hints.get("bayesprism_missing_display_en"),
        "bayesprism_missing_display_zh": hints.get("bayesprism_missing_display_zh"),
        "prefer_missing_downstream_analysis_suggestions": hints.get(
            "prefer_missing_downstream_analysis_suggestions"
        ),
        "use_missing_downstream_analysis_suggestions_after_checklist": hints.get(
            "use_missing_downstream_analysis_suggestions_after_checklist"
        ),
        "do_not_omit_bayesprism_from_recommendations_when_missing": hints.get(
            "do_not_omit_bayesprism_from_recommendations_when_missing"
        ),
        "use_suggested_command_details_for_parameter_prompts": hints.get(
            "use_suggested_command_details_for_parameter_prompts"
        ),
        "suggested_command_details_include_prompt_en_and_prompt_zh": hints.get(
            "suggested_command_details_include_prompt_en_and_prompt_zh"
        ),
        "missing_downstream_suggestions_include_reason_en_and_reason_zh": hints.get(
            "missing_downstream_suggestions_include_reason_en_and_reason_zh"
        ),
        "show_parameter_placeholders_in_natural_language": hints.get(
            "show_parameter_placeholders_in_natural_language"
        ),
        "parameter_hints_source": hints.get("parameter_hints_source"),
        "parameter_source_of_truth": hints.get("parameter_source_of_truth"),
        "do_not_drop_pending_rows_or_missing_columns": hints.get("do_not_drop_pending_rows_or_missing_columns"),
        "do_not_shell_parse_mcp_tool_results": hints.get("do_not_shell_parse_mcp_tool_results"),
        "do_not_read_client_tool_result_files_for_rendering": hints.get(
            "do_not_read_client_tool_result_files_for_rendering"
        ),
        "do_not_request_user_approval_for_json_extraction_commands": hints.get(
            "do_not_request_user_approval_for_json_extraction_commands"
        ),
        "no_extra_confirmation_for_reading_map_payload": hints.get("no_extra_confirmation_for_reading_map_payload"),
        "confirmation_prompts_should_describe_intent_not_raw_command": hints.get(
            "confirmation_prompts_should_describe_intent_not_raw_command"
        ),
        "approval_prompt_style": hints.get("approval_prompt_style"),
        "approval_prompt_example_en": hints.get("approval_prompt_example_en"),
        "approval_prompt_example_zh": hints.get("approval_prompt_example_zh"),
    }
    return {key: value for key, value in compacted.items() if value is not None}


def _compact_recommended_action(action: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(action, dict):
        return None
    return {
        "id": action.get("id"),
        "label": action.get("label"),
        "label_zh": action.get("label_zh"),
        "why": action.get("why"),
        "explanation": action.get("explanation"),
        "explanation_zh": action.get("explanation_zh"),
        "suggested_commands": _limited_sequence(action.get("suggested_commands"), limit=12),
        "suggested_command_details": _compact_suggested_command_details(
            action.get("suggested_command_details"),
            limit=12,
        ),
    }


def _summary_recommended_action(action: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(action, dict):
        return None
    return {
        key: value
        for key, value in {
            "id": action.get("id"),
            "label": action.get("label"),
            "label_zh": action.get("label_zh"),
            "explanation": _trim_text(action.get("explanation") or action.get("why"), limit=220),
            "explanation_zh": _trim_text(action.get("explanation_zh"), limit=220),
            "suggested_commands": _limited_sequence(action.get("suggested_commands"), limit=3),
        }.items()
        if value not in (None, [], {})
    }


def _compact_deconvolution_methods(methods: Any) -> Any:
    if not isinstance(methods, dict):
        return methods
    evidence: Dict[str, List[Any]] = {}
    for key, value in methods.get("evidence", {}).items():
        evidence[str(key)] = _limited_sequence(value, limit=3)
    compacted = dict(methods)
    if evidence:
        compacted["evidence"] = evidence
    for key in (
        "detected_ids",
        "detected_labels",
        "missing_ids",
        "missing_labels",
        "default_bundle_ids",
        "default_bundle_detected_ids",
        "default_bundle_detected_labels",
        "default_bundle_missing_ids",
        "default_bundle_missing_labels",
        "optional_ids",
        "optional_labels",
        "optional_detected_ids",
        "optional_detected_labels",
        "optional_missing_ids",
        "optional_missing_labels",
        "inferred_from_wrapper_ids",
    ):
        compacted[key] = _limited_sequence(compacted.get(key), limit=12)
    return compacted


def _compact_agent_parameter_hints(items: Any, *, limit: int = 10) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                key: value
                for key, value in {
                    "name": item.get("name"),
                    "flag": item.get("flag"),
                    "required": item.get("required"),
                    "type": item.get("type"),
                    "placeholder": item.get("placeholder"),
                    "prompt_zh": _trim_text(item.get("prompt_zh"), limit=120),
                    "prompt_en": _trim_text(item.get("prompt_en"), limit=120),
                    "choices": _limited_sequence(item.get("choices"), limit=8),
                    "default": item.get("default"),
                    "resource_sensitive": item.get("resource_sensitive"),
                    "contextual_placeholders": _limited_sequence(item.get("contextual_placeholders"), limit=4),
                }.items()
                if value not in (None, [], {})
            }
        )
    return compacted


def _compact_suggested_command_details(
    items: Any,
    *,
    limit: int = 12,
    include_command_template: bool = True,
    include_parameter_hints: bool = True,
    parameter_hint_limit: int = 10,
) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                key: value
                for key, value in {
                    "command_id": item.get("command_id"),
                    "command": item.get("command"),
                    "command_template": item.get("command_template") if include_command_template else None,
                    "summary": _trim_text(item.get("summary"), limit=180),
                    "summary_zh": _trim_text(item.get("summary_zh"), limit=180),
                    "parameter_hints": (
                        _compact_agent_parameter_hints(item.get("parameter_hints"), limit=parameter_hint_limit)
                        if include_parameter_hints
                        else None
                    ),
                    "parameter_groups": (
                        _limited_sequence(item.get("parameter_groups"), limit=4)
                        if include_parameter_hints
                        else None
                    ),
                    "parameter_source": item.get("parameter_source") if include_parameter_hints else None,
                    "agent_hint_source": item.get("agent_hint_source") if include_parameter_hints else None,
                }.items()
                if value not in (None, [], {})
            }
        )
    return compacted


def _compact_missing_downstream_suggestions(
    items: Any,
    *,
    limit: int = 12,
    include_command_details: bool = True,
) -> List[Dict[str, Any]]:
    compacted: List[Dict[str, Any]] = []
    if not isinstance(items, list):
        return compacted
    for item in items[:limit]:
        if not isinstance(item, dict):
            continue
        compacted.append(
            {
                key: value
                for key, value in {
                    "id": item.get("id"),
                    "command_id": item.get("command_id"),
                    "label": item.get("label"),
                    "label_zh": item.get("label_zh"),
                    "reason": _trim_text(item.get("reason"), limit=220),
                    "reason_zh": _trim_text(item.get("reason_zh"), limit=220),
                    "missing_analysis_ids": _limited_sequence(item.get("missing_analysis_ids"), limit=12),
                    "missing_analysis_labels": _limited_sequence(item.get("missing_analysis_labels"), limit=12),
                    "missing_analysis_labels_zh": _limited_sequence(
                        item.get("missing_analysis_labels_zh"),
                        limit=12,
                    ),
                    "suggested_commands": _limited_sequence(item.get("suggested_commands"), limit=4),
                    "suggested_command_details": (
                        _compact_suggested_command_details(
                            item.get("suggested_command_details"),
                            limit=2,
                            include_command_template=False,
                            include_parameter_hints=False,
                        )
                        if include_command_details
                        else None
                    ),
                    "notes": _limited_sequence(item.get("notes"), limit=4),
                    "notes_zh": _limited_sequence(item.get("notes_zh"), limit=4),
                }.items()
                if value not in (None, [], {})
            }
        )
    return compacted


def _compact_scenario_payload(scenario: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(scenario, dict):
        return None

    choice_details: List[Dict[str, Any]] = []
    for item in scenario.get("choice_details", [])[:6]:
        if not isinstance(item, dict):
            continue
        choice_details.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "label_zh": item.get("label_zh"),
                "when_to_choose": item.get("when_to_choose"),
                "when_to_choose_zh": item.get("when_to_choose_zh"),
                "suggested_commands": _limited_sequence(item.get("suggested_commands"), limit=12),
                "suggested_command_details": _compact_suggested_command_details(
                    item.get("suggested_command_details"),
                    limit=12,
                ),
            }
        )

    return {
        "id": scenario.get("id"),
        "summary": scenario.get("summary"),
        "summary_zh": scenario.get("summary_zh"),
        "communication_goals": _limited_sequence(scenario.get("communication_goals"), limit=8),
        "required_user_decisions": _limited_sequence(scenario.get("required_user_decisions"), limit=8),
        "recommended_choice": scenario.get("recommended_choice"),
        "choice_details": choice_details,
    }


def _summary_scenario_payload(scenario: Any) -> Optional[Dict[str, Any]]:
    if not isinstance(scenario, dict):
        return None

    choice_details: List[Dict[str, Any]] = []
    for item in scenario.get("choice_details", [])[:4]:
        if not isinstance(item, dict):
            continue
        choice_details.append(
            {
                key: value
                for key, value in {
                    "id": item.get("id"),
                    "label": item.get("label"),
                    "label_zh": item.get("label_zh"),
                    "when_to_choose": _trim_text(item.get("when_to_choose"), limit=180),
                    "when_to_choose_zh": _trim_text(item.get("when_to_choose_zh"), limit=180),
                }.items()
                if value not in (None, [], {})
            }
        )

    return {
        key: value
        for key, value in {
            "id": scenario.get("id"),
            "summary": _trim_text(scenario.get("summary"), limit=260),
            "summary_zh": _trim_text(scenario.get("summary_zh"), limit=260),
            "recommended_choice": scenario.get("recommended_choice"),
            "choice_details": choice_details,
        }.items()
        if value not in (None, [], {})
    }


def _default_user_facing_checklist(payload: Dict[str, Any], *, zh: bool) -> Tuple[Optional[str], Optional[str]]:
    text_fields = (
        ["workflow_checklist_table_zh", "workflow_checklist_report_zh"]
        if zh
        else ["workflow_checklist_table", "workflow_checklist_report"]
    )
    for field in text_fields:
        value = payload.get(field)
        if isinstance(value, str) and value.strip():
            return field, value

    lines_field = "workflow_checklist_lines_zh" if zh else "workflow_checklist_lines"
    lines = payload.get(lines_field)
    if isinstance(lines, list) and lines:
        rendered = "\n".join(str(line) for line in lines if line is not None).strip()
        if rendered:
            return lines_field, rendered
    return None, None


def _summary_map_payload_for_transport(payload: Dict[str, Any], *, language: str = "en") -> Dict[str, Any]:
    language = _response_language(language)
    zh_table = payload.get("workflow_checklist_table_zh") or payload.get("default_user_facing_checklist_zh")
    en_table = payload.get("workflow_checklist_table") or payload.get("default_user_facing_checklist")
    if language == "zh":
        default_checklist = zh_table or en_table
        default_checklist_field = "workflow_checklist_table_zh" if zh_table else "workflow_checklist_table"
        default_checklist_language = "zh" if zh_table else "en"
    else:
        default_checklist = en_table or zh_table
        default_checklist_field = "workflow_checklist_table" if en_table else "workflow_checklist_table_zh"
        default_checklist_language = "en" if en_table else "zh"
    scenario = _summary_scenario_payload(payload.get("scenario"))
    compacted: Dict[str, Any] = {
        "success": payload.get("success"),
        "path": payload.get("path"),
        "path_kind": payload.get("path_kind"),
        "current_stage": payload.get("current_stage"),
        "current_label": payload.get("current_label"),
        "current_label_zh": payload.get("current_label_zh"),
        "current_stage_confidence": payload.get("current_stage_confidence"),
        "scan": _compact_scan_payload(payload.get("scan")),
        "scan_warning": payload.get("scan_warning"),
        "scan_warning_zh": payload.get("scan_warning_zh"),
        "scan_retry_recommended": payload.get("scan_retry_recommended"),
        "scan_retry_command": payload.get("scan_retry_command"),
        "default_user_facing_checklist": default_checklist,
        "default_user_facing_checklist_field": default_checklist_field,
        "default_user_facing_checklist_language": default_checklist_language,
        "default_user_facing_checklist_format": "markdown_table",
        "available_checklist_languages": [
            lang
            for lang, table in (("en", en_table), ("zh", zh_table))
            if isinstance(table, str) and table.strip()
        ],
        "workflow_checklist_title": payload.get("workflow_checklist_title"),
        "workflow_checklist_title_zh": payload.get("workflow_checklist_title_zh"),
        "workflow_checklist_authoritative_for_missing_items": True,
        "do_not_rebuild_checklist_from_deconvolution_methods": True,
        "recommended_action": _summary_recommended_action(payload.get("recommended_action")),
        "missing_downstream_analysis_suggestions": _compact_missing_downstream_suggestions(
            payload.get("missing_downstream_analysis_suggestions"),
            limit=8,
            include_command_details=False,
        ),
        "scenario": scenario,
        "response_guidance": {
            "render_default_user_facing_checklist_directly": True,
            "do_not_read_tool_result_files_for_checklist": True,
            "do_not_shell_parse_mcp_json": True,
            "request_structured_detail_only_for_debug_or_command_planning": True,
            "keep_prebuilt_table_layout_even_for_tall_rows": True,
            "do_not_convert_prebuilt_table_to_field_list": True,
            "default_user_facing_checklist_format": "markdown_table",
            "do_not_convert_prebuilt_markdown_table_to_ascii_box_table": True,
            "preserve_absolute_evidence_paths_in_checklist": True,
            "prefer_line_breaks_between_multiple_evidence_paths": True,
            "immune_deconvolution_multiline_spacing": "double_break",
            "do_not_collapse_multiline_evidence_paths_inside_table_cells": True,
            "multi_value_cell_visible_fallback_separator": "•",
            "do_not_shorten_detected_paths_in_tables": True,
            "do_not_use_zh_checklist_fields_in_english_rendering": True,
            "immune_deconvolution_detected_column_source": "workflow_checklist[].evidence_display_paths",
            "immune_deconvolution_missing_column_source_en": "workflow_checklist[].missing_column_value",
            "immune_deconvolution_missing_column_source_zh": "workflow_checklist[].missing_column_value_zh",
            "tpm_matrix_detected_column_source": "workflow_checklist[].evidence_display_paths",
            "do_not_borrow_feature_scoring_outputs_for_tpm_matrix_row": True,
            "if_rebuilding_table_use_evidence_display_paths_not_deconvolution_method_labels": True,
            "use_detected_and_missing_column_value_fields_when_rebuilding_table": True,
            "if_html_line_breaks_are_stripped_keep_visible_cell_separator": True,
        },
        "next_steps": [
            "Render default_user_facing_checklist directly; do not read client tool-results files.",
            "Request response_detail=structured only for debug or command-planning details.",
        ],
    }
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _compact_map_payload_for_transport(
    payload: Dict[str, Any],
    *,
    detail: str = "summary",
    language: str = "en",
) -> Dict[str, Any]:
    if detail == "structured":
        return _structured_map_payload_for_transport(payload)
    return _summary_map_payload_for_transport(payload, language=language)


def _structured_map_payload_for_transport(payload: Dict[str, Any]) -> Dict[str, Any]:
    default_checklist_field, default_checklist = _default_user_facing_checklist(payload, zh=False)
    default_checklist_field_zh, default_checklist_zh = _default_user_facing_checklist(payload, zh=True)
    compacted: Dict[str, Any] = {
        "success": payload.get("success"),
        "path": payload.get("path"),
        "path_kind": payload.get("path_kind"),
        "default_user_facing_checklist": default_checklist,
        "default_user_facing_checklist_field": default_checklist_field,
        "default_user_facing_checklist_zh": default_checklist_zh,
        "default_user_facing_checklist_field_zh": default_checklist_field_zh,
        "scan": _compact_scan_payload(payload.get("scan")),
        "scan_warning": payload.get("scan_warning"),
        "scan_warning_zh": payload.get("scan_warning_zh"),
        "scan_retry_recommended": payload.get("scan_retry_recommended"),
        "scan_retry_command": payload.get("scan_retry_command"),
        "scan_retry_limits": payload.get("scan_retry_limits"),
        "auto_scan_retry_performed": payload.get("auto_scan_retry_performed"),
        "current_stage": payload.get("current_stage"),
        "current_label": payload.get("current_label"),
        "current_label_zh": payload.get("current_label_zh"),
        "current_stage_confidence": payload.get("current_stage_confidence"),
        "completed_stages": _limited_sequence(payload.get("completed_stages"), limit=24),
        "completed_stage_labels": _limited_sequence(payload.get("completed_stage_labels"), limit=24),
        "completed_stage_labels_zh": _limited_sequence(payload.get("completed_stage_labels_zh"), limit=24),
        "completed_stage_summaries": _limited_sequence(payload.get("completed_stage_summaries"), limit=16),
        "next_stages": _limited_sequence(payload.get("next_stages"), limit=12),
        "next_stage_labels": _limited_sequence(payload.get("next_stage_labels"), limit=12),
        "next_stage_labels_zh": _limited_sequence(payload.get("next_stage_labels_zh"), limit=12),
        "next_step_summaries": _limited_sequence(payload.get("next_step_summaries"), limit=12),
        "resume_recommended": payload.get("resume_recommended"),
        "rerun_recommended": payload.get("rerun_recommended"),
        "decision_prompt": payload.get("decision_prompt"),
        "recommended_action": _compact_recommended_action(payload.get("recommended_action")),
        "missing_downstream_analysis_suggestions": _compact_missing_downstream_suggestions(
            payload.get("missing_downstream_analysis_suggestions")
        ),
        "suggested_command_details": _compact_suggested_command_details(payload.get("suggested_command_details")),
        "scenario": _compact_scenario_payload(payload.get("scenario")),
        "external_analysis_hints": _compact_external_analysis_hints(payload.get("external_analysis_hints")),
        "deconvolution_methods": _compact_deconvolution_methods(payload.get("deconvolution_methods")),
        "function_detections": _compact_function_detections(payload.get("function_detections")),
        "detected_functions": _limited_sequence(payload.get("detected_functions"), limit=24),
        "reusable_result_functions": _limited_sequence(payload.get("reusable_result_functions"), limit=24),
        "likely_iobrpy_functions": _limited_sequence(payload.get("likely_iobrpy_functions"), limit=24),
        "content_verified_functions": _limited_sequence(payload.get("content_verified_functions"), limit=24),
        "agent_fallback": _compact_agent_fallback(payload.get("agent_fallback")),
        "existing_result_investigation": _compact_existing_result_investigation(
            payload.get("existing_result_investigation")
        ),
        "existing_result_investigation_applied": payload.get("existing_result_investigation_applied"),
        "checklist_display_hints": payload.get("checklist_display_hints"),
        "workflow_checklist_title": payload.get("workflow_checklist_title"),
        "workflow_checklist_title_zh": payload.get("workflow_checklist_title_zh"),
        "workflow_checklist": _compact_workflow_checklist(payload.get("workflow_checklist")),
        "workflow_checklist_table": payload.get("workflow_checklist_table"),
        "workflow_checklist_table_zh": payload.get("workflow_checklist_table_zh"),
        "roadmap": _compact_roadmap_payload(payload.get("roadmap")),
        "agent_rendering_hints": _compact_agent_rendering_hints(payload.get("agent_rendering_hints")),
    }
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _compact_map_payload_for_analyze_summary(pipeline_map: Any) -> Dict[str, Any]:
    if not isinstance(pipeline_map, dict):
        return {}
    keep_keys = [
        "success",
        "path",
        "path_kind",
        "current_stage",
        "current_label",
        "current_label_zh",
        "current_stage_confidence",
        "scan",
        "scan_warning",
        "scan_warning_zh",
        "scan_retry_recommended",
        "default_user_facing_checklist",
        "default_user_facing_checklist_field",
        "default_user_facing_checklist_language",
        "default_user_facing_checklist_zh",
        "default_user_facing_checklist_field_zh",
        "workflow_checklist_title",
        "workflow_checklist_title_zh",
        "workflow_checklist_table",
        "workflow_checklist_table_zh",
        "workflow_checklist",
        "recommended_action",
        "scenario",
        "deconvolution_methods",
        "agent_fallback",
        "agent_rendering_hints",
    ]
    compacted = {key: pipeline_map.get(key) for key in keep_keys}
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _slim_next_step_card_for_transport(card: Any) -> Dict[str, Any]:
    if not isinstance(card, dict):
        return {}
    compacted = {
        "id": card.get("id"),
        "kind": card.get("kind"),
        "priority": card.get("priority"),
        "title": card.get("title"),
        "command": card.get("command"),
        "why_this_is_recommended": _trim_text(card.get("why_this_is_recommended"), limit=160),
        "expected_result": _trim_text(card.get("expected_result"), limit=140),
        "existing_input_candidates": _compact_existing_input_candidates_for_transport(
            card.get("existing_input_candidates"),
            limit=2,
        ),
        "concise_confirmation_prompts": _compact_concise_confirmation_prompts_for_transport(
            card.get("concise_confirmation_prompts"),
            limit=1,
        ),
    }
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _compact_analyze_guidance_for_default_transport(guidance: Any, pipeline_map: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(guidance, dict):
        guidance = {}
    ranked = [
        card
        for card in (
            _slim_next_step_card_for_transport(item)
            for item in _limited_sequence(guidance.get("ranked_next_step_options"), limit=4)
        )
        if card
    ]
    concise_confirmation_plan = _compact_concise_confirmation_plan_for_transport(
        guidance.get("concise_confirmation_plan", {})
        if isinstance(guidance.get("concise_confirmation_plan"), dict)
        else {}
    )
    response_guidance = guidance.get("response_guidance", {})
    if not isinstance(response_guidance, dict):
        response_guidance = {}
    compacted = {
        "current_stage": guidance.get("current_stage") or pipeline_map.get("current_stage"),
        "current_label": guidance.get("current_label") or pipeline_map.get("current_label"),
        "current_label_zh": pipeline_map.get("current_label_zh"),
        "current_stage_confidence": guidance.get("current_stage_confidence")
        or pipeline_map.get("current_stage_confidence"),
        "recommended_action": guidance.get("recommended_action") or pipeline_map.get("recommended_action"),
        "resume_recommended": guidance.get("resume_recommended"),
        "rerun_recommended": guidance.get("rerun_recommended"),
        "recommended_command": guidance.get("recommended_command"),
        "missing_parameters": _limited_sequence(guidance.get("missing_parameters"), limit=4),
        "confirm_parameters": _compact_parameter_cards_for_transport(guidance.get("confirm_parameters"), limit=2),
        "important_optional_parameters": _compact_parameter_cards_for_transport(
            guidance.get("important_optional_parameters"),
            limit=2,
        ),
        "ranked_next_step_options": ranked,
        "concise_confirmation_plan": concise_confirmation_plan,
        "response_guidance": {
            "render_from_default_user_facing_checklist_first": True,
            "do_not_read_tool_result_files": True,
            "do_not_shell_parse_mcp_json": True,
            "do_not_call_analyze_path_task_again_for_checklist_rendering": True,
            "prefer_concise_natural_language_confirmation_questions": response_guidance.get(
                "prefer_concise_natural_language_confirmation_questions",
                True,
            ),
            "ask_only_the_highest_impact_one_or_two_confirmation_items_first": response_guidance.get(
                "ask_only_the_highest_impact_one_or_two_confirmation_items_first",
                True,
            ),
        },
    }
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _compact_analyze_path_task_payload_for_default_transport(payload: Dict[str, Any]) -> Dict[str, Any]:
    pipeline_map = payload.get("pipeline_map", {})
    compact_pipeline_map = _compact_map_payload_for_analyze_summary(pipeline_map)
    guidance = _compact_analyze_guidance_for_default_transport(payload.get("agent_guidance"), compact_pipeline_map)
    compacted = {
        "success": payload.get("success"),
        "entrypoint": payload.get("entrypoint"),
        "path": payload.get("path"),
        "task": payload.get("task"),
        "current_stage": compact_pipeline_map.get("current_stage"),
        "current_label": compact_pipeline_map.get("current_label"),
        "current_label_zh": compact_pipeline_map.get("current_label_zh"),
        "current_stage_confidence": compact_pipeline_map.get("current_stage_confidence"),
        "default_user_facing_checklist": compact_pipeline_map.get("default_user_facing_checklist"),
        "default_user_facing_checklist_field": compact_pipeline_map.get("default_user_facing_checklist_field"),
        "default_user_facing_checklist_language": compact_pipeline_map.get("default_user_facing_checklist_language"),
        "pipeline_map": compact_pipeline_map,
        "recommendation": payload.get("recommendation"),
        "agent_guidance": guidance,
        "next_steps": [
            "Render default_user_facing_checklist directly; it is already localized when a localized checklist is available.",
            "Do not run python3 -c, jq, or similar parsers just to extract fields from this payload.",
            "Do not call analyze_path_task again just to render the checklist; use map_path or this compact payload.",
        ],
        "transport_compaction": {
            **(
                payload.get("transport_compaction")
                if isinstance(payload.get("transport_compaction"), dict)
                else {}
            ),
            "applied": True,
            "default_analyze_payload": "compact_summary",
            "rich_guidance_omitted_by_default": True,
            "prevents_tool_result_file_roundtrip": True,
        },
    }
    return {key: value for key, value in compacted.items() if value not in (None, [], {})}


def _summarize_map_payload_for_recommendation(payload: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not isinstance(payload, dict):
        return None

    scenario = payload.get("scenario", {})
    scenario_id = scenario.get("id") if isinstance(scenario, dict) else None
    summary: Dict[str, Any] = {
        "current_stage": payload.get("current_stage"),
        "current_label": payload.get("current_label"),
        "current_label_zh": payload.get("current_label_zh"),
        "current_stage_confidence": payload.get("current_stage_confidence"),
        "resume_recommended": payload.get("resume_recommended"),
        "rerun_recommended": payload.get("rerun_recommended"),
        "scan_warning": payload.get("scan_warning"),
        "completed_stage_labels": _limited_sequence(payload.get("completed_stage_labels"), limit=12),
        "completed_stage_labels_zh": _limited_sequence(payload.get("completed_stage_labels_zh"), limit=12),
        "next_stage_labels": _limited_sequence(payload.get("next_stage_labels"), limit=8),
        "next_stage_labels_zh": _limited_sequence(payload.get("next_stage_labels_zh"), limit=8),
    }
    if scenario_id:
        summary["scenario"] = {
            "id": scenario_id,
            "summary": scenario.get("summary"),
            "summary_zh": scenario.get("summary_zh"),
        }
    return {key: value for key, value in summary.items() if value not in (None, [], {})}


def _call_list_native_commands_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    details = bool(arguments.get("details", False))
    tool_summaries = _native_command_summary_payload()
    payload: Dict[str, Any] = {
        "success": True,
        "entrypoint": "iobrpy-cli-mcp",
        "native_invocation_pattern": "iobrpy-cli <native_command> [native_args...]",
        "path_shorthand": "iobrpy-cli <existing_path> -> map --path <existing_path>",
        "helper_tools": [
            {"name": tool["name"], "description": tool.get("description", "")}
            for tool in get_helper_tools()
        ],
        "native_command_count": len(tool_summaries),
        "commands": tool_summaries,
        "policy": {
            "excluded_native_commands": ["ai"],
            "preferred_first_step": "map_path",
            "compatibility_namespaces": ["analyze", "quantify", "workflow", "immune", "hla_tcr"],
        },
    }
    if details:
        payload["tool_schemas"] = get_native_tools()
    return _decorate_helper_payload("list_native_commands", payload)


def _call_map_path_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    path = arguments.get("path")
    if not isinstance(path, str) or not path:
        raise ValueError("map_path requires a non-empty string path.")

    compact_for_transport = bool(arguments.get("_compact_for_transport", True))
    max_depth = int(arguments.get("max_depth", 8))
    max_entries = int(arguments.get("max_entries", 8000))
    requested_quick = bool(arguments.get("quick", False))
    quick = False
    focus = arguments.get("focus") or ()
    investigate_existing = bool(arguments.get("investigate_existing", False))
    response_detail = str(arguments.get("response_detail") or "summary").strip().lower()
    if response_detail not in {"summary", "structured"}:
        response_detail = "summary"
    language = _response_language(arguments.get("language"))

    payload = map_pipeline_path(
        Path(path),
        max_depth=max_depth,
        max_entries=max_entries,
        quick=quick,
        focus=focus,
        investigate_existing=investigate_existing,
        auto_scan_retry=True,
        auto_investigate_existing=True,
    )
    if requested_quick:
        payload["quick_request_ignored"] = True
        payload["quick_mode_deprecated"] = True
        payload["quick_mode_deprecation_note"] = (
            "map_path now uses the full/deep scan path for agent calls even when older clients pass quick=true."
        )
    payload["entrypoint"] = "iobrpy-cli-mcp"
    payload["next_steps"] = [
        "Continue from the current stage when resume_recommended is true instead of rerunning the full pipeline by default.",
        "Run recommend_workflow after map_path only when you still need a narrower command recommendation for the next biological task.",
        "Confirm environment-specific required parameters such as index and output paths with the user before execution.",
    ]
    if compact_for_transport:
        payload = _compact_map_payload_for_transport(payload, detail=response_detail, language=language)
        if requested_quick:
            payload["quick_request_ignored"] = True
            payload["quick_mode_deprecated"] = True
        payload["transport_compaction"] = {
            "applied": True,
            "response_detail": response_detail,
            "language": language,
            "removed_large_duplicate_renderings": True,
            "removed_pre_rendered_line_duplicates": True,
            "kept_structured_checklist": bool(payload.get("workflow_checklist")),
            "single_language_checklist_in_summary": response_detail == "summary",
            "prevents_tool_result_file_roundtrip": response_detail == "summary",
        }
    return _decorate_helper_payload("map_path", payload)


def _call_recommend_workflow_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    path = arguments.get("path")
    task = str(arguments.get("task") or "")
    if not path and not task:
        raise ValueError("recommend_workflow requires at least one of path or task.")
    if path is not None and not isinstance(path, str):
        raise ValueError("recommend_workflow path must be a string when provided.")

    detected_override = arguments.get("_detected_input")
    pipeline_map_override = arguments.get("_pipeline_map_payload")
    include_pipeline_map = bool(arguments.get("_include_pipeline_map", False))
    compact_for_transport = bool(arguments.get("_compact_for_transport", True))

    detected = (
        detected_override
        if isinstance(detected_override, dict)
        else (_detect_input_summary(Path(path)) if path else {"kind": "unspecified", "path": None})
    )
    pipeline_map = pipeline_map_override if isinstance(pipeline_map_override, dict) else None
    recommendation = _recommend_for_summary(detected, task)
    payload = {
        "success": True,
        "entrypoint": "iobrpy-cli-mcp",
        "task": task,
        "detected_input": detected,
        "path_state": _summarize_map_payload_for_recommendation(pipeline_map),
        "recommendation": recommendation,
        "next_steps": [
            "Use map_path first when you need to decide between continuing an existing analysis and rerunning it.",
            "Use list_native_commands to inspect the available native commands.",
            "Review recommendation.confirm_parameters before execution so high-impact settings are not guessed.",
            "Prefer the top-level native command tools rather than legacy wrapper namespaces.",
        ],
    }
    if include_pipeline_map and isinstance(pipeline_map, dict):
        payload["pipeline_map"] = (
            _compact_map_payload_for_transport(pipeline_map) if compact_for_transport else pipeline_map
        )
    return _decorate_helper_payload("recommend_workflow", payload)


def _call_doctor_environment_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    path = arguments.get("path")
    task = str(arguments.get("task") or "")
    external = bool(arguments.get("external", False))
    if path is not None and not isinstance(path, str):
        raise ValueError("doctor_environment path must be a string when provided.")

    detected_override = arguments.get("_detected_input")
    pipeline_map_override = arguments.get("_pipeline_map_payload")
    include_pipeline_map = bool(arguments.get("_include_pipeline_map", False))
    compact_for_transport = bool(arguments.get("_compact_for_transport", True))
    version_info = _workflow_version_info(external)
    detected = detected_override if isinstance(detected_override, dict) else (_detect_input_summary(Path(path)) if path else None)
    pipeline_map = pipeline_map_override if isinstance(pipeline_map_override, dict) else None
    recommendation = _recommend_for_summary(detected, task) if detected else None
    payload = {
        "success": True,
        "entrypoint": "iobrpy-cli-mcp",
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "cwd": str(Path.cwd()),
        "console_scripts": _console_script_names(),
        "commands_on_path": {
            name: _command_on_path(name)
            for name in ("iobrpy", "iobrpy-cli", "iobrpy-cli-mcp")
        },
        "installed_distributions": {
            "iobrpy": _distribution_version("iobrpy"),
            "iobrpy-cli": _distribution_version("iobrpy-cli"),
        },
        "package_conflicts": (
            [
                "A standalone iobrpy-cli distribution is also installed. If its code is stale, it can shadow the harness bundled with iobrpy."
            ]
            if _distribution_version("iobrpy-cli")
            else []
        ),
        "version_info": version_info,
        "detected_input": detected,
        "path_state": _summarize_map_payload_for_recommendation(pipeline_map),
        "recommendation": recommendation,
        "guidance": [
            "Prefer native iobrpy tools or agent helper tools before scanning source files.",
            "Do not substitute the R IOBR package when Python iobrpy is available.",
            "For path-driven requests, run map_path with its default full/deep scan behavior before analyze_path_task; do not ask a quick-versus-deep scan-mode question.",
            "When recommend_workflow returns confirm_parameters, confirm those before execution.",
        ],
    }
    if include_pipeline_map and isinstance(pipeline_map, dict):
        payload["pipeline_map"] = (
            _compact_map_payload_for_transport(pipeline_map) if compact_for_transport else pipeline_map
        )
    return _decorate_helper_payload("doctor_environment", payload)


def _extract_command_name_from_invocation(text: str) -> Optional[str]:
    if not isinstance(text, str) or not text.strip():
        return None
    for token in text.strip().split():
        if token in {"iobrpy-cli", "iobrpy"}:
            continue
    parts = text.strip().split()
    for idx, token in enumerate(parts[:-1]):
        if token in {"iobrpy-cli", "iobrpy"}:
            candidate = parts[idx + 1].strip()
            if candidate in {tool["name"] for tool in get_native_tools()}:
                return candidate
    return None


def _required_parameter_cards(meta: Dict[str, Any]) -> List[Dict[str, Any]]:
    cards: List[Dict[str, Any]] = []
    for spec in meta.get("argument_specs", []):
        if not spec.required:
            continue
        cards.append(
            {
                "parameter": spec.property_name,
                "display_name": spec.preferred_option or spec.property_name,
                "description": spec.help_text,
            }
        )

    required_one_of = meta.get("required_one_of", [])
    for group in required_one_of:
        if not group:
            continue
        cards.append(
            {
                "parameter_group": list(group),
                "display_name": " / ".join(group),
                "description": "One of these parameter groups must be provided together.",
            }
        )
    return cards


def _dedupe_text_values(values: Sequence[Any], *, limit: int = 6) -> List[str]:
    deduped: List[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if not text or text in seen:
            continue
        deduped.append(text)
        seen.add(text)
        if len(deduped) >= limit:
            break
    return deduped


def _guidance_checklist_index(map_payload: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {
        str(item.get("id")): item
        for item in map_payload.get("workflow_checklist", [])
        if isinstance(item, dict) and item.get("id")
    }


def _checklist_item_paths_for_guidance(item: Dict[str, Any], *, limit: int = 4) -> List[str]:
    for field in ("evidence_display_paths", "evidence_paths", "raw_evidence_paths"):
        value = item.get(field)
        if isinstance(value, list) and value:
            return _dedupe_text_values(value, limit=limit)
    return []


def _candidate_from_checklist_item(
    checklist_index: Dict[str, Dict[str, Any]],
    item_id: str,
    *,
    label: str,
    why_it_matches: str,
    limit: int = 4,
) -> Optional[Dict[str, Any]]:
    item = checklist_index.get(item_id, {})
    if not isinstance(item, dict):
        return None
    paths = _checklist_item_paths_for_guidance(item, limit=limit)
    if not paths:
        return None
    return {
        "source_id": item_id,
        "source_label": str(item.get("label") or item.get("label_zh") or item_id),
        "label": label,
        "why_it_matches": why_it_matches,
        "paths": paths,
    }


def _cluster_features_indexing_note(command_name: str) -> str:
    label = "tme_cluster" if command_name == "tme_cluster" else "nmf"
    return (
        f"For `{label}`, the first column is the sample ID and is excluded before counting `--features`, "
        "so a range such as `1:22` refers to the second through twenty-third actual CSV columns. The left bound is commonly `1`."
    )


def _candidate_from_paths(
    *,
    source_id: str,
    source_label: str,
    label: str,
    why_it_matches: str,
    paths: Sequence[Any],
    preferred: bool = False,
    needs_user_choice: bool = False,
    selection_guidance: str = "",
    features_indexing_note: str = "",
    suggested_id_column: str = "",
    limit: int = 4,
) -> Optional[Dict[str, Any]]:
    deduped_paths = _dedupe_text_values(paths, limit=limit)
    if not deduped_paths:
        return None
    candidate: Dict[str, Any] = {
        "source_id": source_id,
        "source_label": source_label,
        "label": label,
        "why_it_matches": why_it_matches,
        "paths": deduped_paths,
    }
    if preferred:
        candidate["preferred"] = True
    if needs_user_choice:
        candidate["needs_user_choice"] = True
    if selection_guidance:
        candidate["selection_guidance"] = selection_guidance
    if features_indexing_note:
        candidate["features_indexing_note"] = features_indexing_note
    if suggested_id_column:
        candidate["suggested_id_column"] = suggested_id_column
    return candidate


def _cluster_deconvolution_input_candidates(command_name: str, map_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    methods = map_payload.get("deconvolution_methods", {})
    if not isinstance(methods, dict):
        return []
    evidence = methods.get("evidence", {})
    if not isinstance(evidence, dict):
        return []

    candidates: List[Dict[str, Any]] = []
    for spec in CLUSTER_DECONVOLUTION_INPUT_PREFERENCES:
        paths = evidence.get(spec["id"])
        if not isinstance(paths, list) or not paths:
            continue
        candidate = _candidate_from_paths(
            source_id=f"deconvolution:{spec['id']}",
            source_label=spec["label"],
            label=f"Detected {spec['label']}",
            why_it_matches=spec["why"],
            paths=paths,
            preferred=not candidates,
            selection_guidance=(
                "This is usually easier to execute directly than a very wide signature-score matrix."
            ),
            features_indexing_note=_cluster_features_indexing_note(command_name),
            suggested_id_column="ID",
        )
        if candidate:
            candidates.append(candidate)
    return candidates


def _cluster_signature_input_candidate(command_name: str, checklist_index: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    item = checklist_index.get("feature_scoring", {})
    if not isinstance(item, dict):
        return None
    return _candidate_from_paths(
        source_id="feature_scoring",
        source_label=str(item.get("label") or item.get("label_zh") or "feature_scoring"),
        label="Detected signature-score outputs",
        why_it_matches=(
            "These signature-score outputs can be clustered, but they should usually be used only after the user confirms which signature subset is biologically meaningful."
        ),
        paths=_checklist_item_paths_for_guidance(item, limit=4),
        needs_user_choice=True,
        selection_guidance=(
            "Do not auto-run clustering on the whole signature matrix when it may contain many hundreds or thousands of columns. Ask the user which signature subset or family they want to cluster on first."
        ),
        features_indexing_note=_cluster_features_indexing_note(command_name),
    )


def _command_existing_input_candidates(command_name: str, map_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    checklist_index = _guidance_checklist_index(map_payload)
    candidates: List[Dict[str, Any]] = []

    def add(item_id: str, label: str, why_it_matches: str, **extras: Any) -> None:
        candidate = _candidate_from_checklist_item(
            checklist_index,
            item_id,
            label=label,
            why_it_matches=why_it_matches,
        )
        if candidate:
            if extras:
                candidate.update(extras)
            candidates.append(candidate)

    if command_name in FULL_PIPELINE_COMMANDS | RAW_FASTQ_STAGE_COMMANDS | REPERTOIRE_COMMANDS:
        add(
            "raw_data",
            "Detected raw sequencing inputs",
            "These detected raw sequencing files are the natural starting point for a fresh upstream or end-to-end run.",
        )

    if command_name in {"runall", "batch_salmon", "batch_star_count", "trust4"}:
        add(
            "quality_control",
            "Detected cleaned FASTQ outputs",
            "If cleaned FASTQ outputs already exist, they may be the practical sequencing handoff for the next upstream step.",
        )

    if command_name == "merge_salmon":
        add(
            "salmon_quant_merge",
            "Detected Salmon quantification outputs",
            "These detected Salmon quantification outputs are the direct merge inputs for merge_salmon.",
        )

    if command_name == "prepare_salmon":
        add(
            "salmon_quant_merge",
            "Detected merged Salmon outputs",
            "These merged Salmon outputs are the most likely direct inputs for prepare_salmon before downstream TPM use.",
        )

    if command_name in {
        "calculate_sig_score",
        "cibersort",
        "quantiseq",
        "epic",
        "mcpcounter",
        "estimate",
        "IPS",
        "LR_cal",
        "bayesprism",
        "tme_profile",
    } | MATRIX_PREPARATION_COMMANDS:
        add(
            "tpm_matrix_ready",
            "Detected TPM-like matrix candidates",
            "These are the most likely bulk-expression inputs that can feed this downstream command directly.",
        )

    if command_name in {"tme_cluster", "nmf"}:
        candidates.extend(_cluster_deconvolution_input_candidates(command_name, map_payload))
        add(
            "immune_deconvolution",
            "Detected deconvolution result directories",
            "These detected deconvolution directories are valid clustering inputs, and they should usually be preferred before a very wide signature-score matrix.",
        )
        signature_candidate = _cluster_signature_input_candidate(command_name, checklist_index)
        if signature_candidate:
            candidates.append(signature_candidate)
        add(
            "tpm_matrix_ready",
            "Detected TPM-like matrix fallback inputs",
            "If you want to rebuild or curate the feature matrix first, these TPM-like inputs are the upstream fallback starting point.",
        )

    if command_name == "merge_star_count":
        add(
            "star_quant_merge",
            "Detected STAR result directories",
            "These directories already contain STAR-side outputs that can be merged into a cohort-level count matrix.",
        )

    if command_name == "count2tpm":
        add(
            "star_quant_merge",
            "Detected merged STAR counts",
            "These STAR-side outputs are the most likely direct inputs for count2tpm before downstream TPM analysis.",
        )

    if command_name in HLA_COMMANDS:
        add(
            "raw_data",
            "Detected upstream sequencing inputs",
            "These paths show the dataset origin, but HLA commands still need the right BAM/CRAM or paired-FASTQ handoff depending on the chosen workflow.",
        )
        add(
            "hla_typing_summary",
            "Detected HLA result outputs",
            "These existing HLA outputs help decide whether to interpret current results or rerun the native HLA workflow.",
        )

    if command_name in REPERTOIRE_COMMANDS:
        add(
            "tcr_bcr_summary",
            "Detected repertoire-analysis outputs",
            "These existing repertoire outputs help decide whether to reuse current results or rerun native TRUST4.",
        )

    if command_name == "LR_cal":
        add(
            "ligand_receptor_analysis",
            "Detected ligand-receptor outputs",
            "These existing LR_cal-style outputs help decide whether to interpret current results or rerun ligand-receptor scoring.",
        )

    if command_name in {"tme_profile"} | DECONVOLUTION_METHOD_COMMANDS:
        add(
            "immune_deconvolution",
            "Detected deconvolution output directories",
            "These existing downstream TME outputs indicate that part of the bundle may already exist and should be reused or summarized before rerunning.",
        )

    return candidates[:6]


def _first_existing_path(candidates: Sequence[Path]) -> Optional[Path]:
    for candidate in candidates:
        try:
            if candidate.exists():
                return candidate
        except OSError:
            continue
    return None


def _find_stage_dir(root: Path, names: Sequence[str]) -> Optional[Path]:
    return _first_existing_path([root / name for name in names])


def _find_first_matching_file(search_roots: Sequence[Path], patterns: Sequence[str]) -> Optional[Path]:
    for search_root in search_roots:
        try:
            if not search_root.exists() or not search_root.is_dir():
                continue
        except OSError:
            continue
        for pattern in patterns:
            try:
                for match in sorted(search_root.glob(pattern)):
                    if match.is_file():
                        return match
            except OSError:
                continue
    return None


def _guess_project_name_from_path(path: Optional[Path]) -> str:
    if path is None:
        return "project"
    name = path.name
    for suffix in (
        "_merged_salmon_tpm.tsv.gz",
        "_merged_salmon_tpm.tsv",
        "_salmon_tpm.tsv.gz",
        "_salmon_tpm.tsv",
        "_tpm_gene_symbols.tsv.gz",
        "_tpm_gene_symbols.tsv",
        ".STAR.count.tsv.gz",
        ".STAR.count.tsv",
        ".csv.gz",
        ".tsv.gz",
        ".csv",
        ".tsv",
    ):
        if name.endswith(suffix):
            trimmed = name[: -len(suffix)].strip("_-.")
            if trimmed:
                return trimmed
    stem = path.stem.strip("_-.")
    return stem or path.parent.name or "project"


def _standardized_tpm_paths(project_root: Path) -> Dict[str, Path]:
    tpm_dir = project_root / "03-tpm"
    return {
        "dir": tpm_dir,
        "salmon_intermediate": tpm_dir / "prepare_salmon.csv",
        "star_intermediate": tpm_dir / "count2tpm.csv",
        "standardized_tpm": tpm_dir / "tpm_matrix.csv",
    }


def _build_target_retry_arguments(
    command_name: str,
    standardized_tpm_path: Path,
    original_arguments: Optional[Dict[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    if command_name not in DOWNSTREAM_MATRIX_COMMANDS:
        return None
    if original_arguments:
        retry_args = dict(original_arguments)
        retry_args["input"] = str(standardized_tpm_path)
        return retry_args
    return {"input": str(standardized_tpm_path)}


def _downstream_recovery_paths_from_project_root(
    command_name: str,
    project_root: Optional[Path],
    *,
    original_arguments: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    if command_name not in DOWNSTREAM_MATRIX_COMMANDS or project_root is None:
        return []
    try:
        if not project_root.exists() or not project_root.is_dir():
            return []
    except OSError:
        return []

    tpm_paths = _standardized_tpm_paths(project_root)
    salmon_dir = _find_stage_dir(project_root, SALMON_STAGE_DIR_NAMES)
    star_dir = _find_stage_dir(project_root, STAR_STAGE_DIR_NAMES)
    salmon_search_roots = [project_root]
    if salmon_dir is not None:
        salmon_search_roots.append(salmon_dir)
    star_search_roots = [project_root]
    if star_dir is not None:
        star_search_roots.append(star_dir)

    salmon_merged_matrix = _find_first_matching_file(
        salmon_search_roots,
        (
            "*_merged_salmon_tpm.tsv.gz",
            "*_merged_salmon_tpm.tsv",
            "*_salmon_tpm.tsv.gz",
            "*_salmon_tpm.tsv",
        ),
    )
    star_merged_matrix = _find_first_matching_file(
        star_search_roots,
        (
            "*.STAR.count.tsv.gz",
            "*.STAR.count.tsv",
            "*_merged_star_count.tsv.gz",
            "*_merged_star_count.tsv",
        ),
    )

    recovery_paths: List[Dict[str, Any]] = []
    target_retry_arguments = _build_target_retry_arguments(
        command_name,
        tpm_paths["standardized_tpm"],
        original_arguments=original_arguments,
    )

    if salmon_merged_matrix is not None or salmon_dir is not None:
        project_name = _guess_project_name_from_path(salmon_merged_matrix or salmon_dir)
        salmon_steps: List[Dict[str, Any]] = []
        detected_paths: List[str] = []
        if salmon_merged_matrix is not None:
            detected_paths.append(str(salmon_merged_matrix))
        elif salmon_dir is not None:
            detected_paths.append(str(salmon_dir))

        if salmon_merged_matrix is None and salmon_dir is not None:
            inferred_merged = salmon_dir / f"{project_name}_merged_salmon_tpm.tsv.gz"
            salmon_steps.append(
                {
                    "command": "merge_salmon",
                    "why": "The current downstream matrix looks unusable, so first rebuild the merged Salmon TPM table from existing Salmon quantification outputs.",
                    "suggested_arguments": {
                        "path_salmon": str(salmon_dir),
                        "project": project_name,
                    },
                    "expected_output": str(inferred_merged),
                }
            )
            salmon_merged_matrix = inferred_merged

        if salmon_merged_matrix is not None:
            salmon_steps.append(
                {
                    "command": "prepare_salmon",
                    "why": "Regenerate a clean TPM matrix from the merged Salmon table using IOBRpy's standard Salmon-to-TPM path.",
                    "suggested_arguments": {
                        "input": str(salmon_merged_matrix),
                        "output": str(tpm_paths["salmon_intermediate"]),
                        "return_feature": "symbol",
                        "remove_version": True,
                    },
                    "expected_output": str(tpm_paths["salmon_intermediate"]),
                }
            )
            salmon_steps.append(
                {
                    "command": "log2_eset",
                    "why": "Standardize the regenerated TPM with log2(x+1), because downstream IOBRpy functions should consume the log2_eset output rather than the raw regenerated TPM.",
                    "suggested_arguments": {
                        "input": str(tpm_paths["salmon_intermediate"]),
                        "output": str(tpm_paths["standardized_tpm"]),
                    },
                    "expected_output": str(tpm_paths["standardized_tpm"]),
                }
            )
            retry_step: Dict[str, Any] = {
                "command": command_name,
                "why": "Retry the downstream command on the standardized log2_eset TPM output instead of the failing current matrix.",
            }
            if target_retry_arguments:
                retry_step["suggested_arguments"] = target_retry_arguments
            salmon_steps.append(retry_step)
            recovery_paths.append(
                {
                    "id": f"{command_name}:salmon_regeneration",
                    "title": f"Rebuild standardized TPM from Salmon outputs, then retry {command_name}",
                    "preferred": True,
                    "why_this_path_exists": (
                        "Existing Salmon outputs were detected, so the safest rollback path is to regenerate TPM through prepare_salmon and then standardize it with log2_eset before retrying the downstream command."
                    ),
                    "detected_upstream_paths": detected_paths,
                    "standardized_tpm_output": str(tpm_paths["standardized_tpm"]),
                    "steps": salmon_steps,
                }
            )

    if star_merged_matrix is not None or star_dir is not None:
        project_name = _guess_project_name_from_path(star_merged_matrix or star_dir)
        star_steps: List[Dict[str, Any]] = []
        detected_paths: List[str] = []
        if star_merged_matrix is not None:
            detected_paths.append(str(star_merged_matrix))
        elif star_dir is not None:
            detected_paths.append(str(star_dir))

        if star_merged_matrix is None and star_dir is not None:
            inferred_merged = star_dir / f"{project_name}.STAR.count.tsv.gz"
            star_steps.append(
                {
                    "command": "merge_star_count",
                    "why": "The current downstream matrix looks unusable, so first merge the existing STAR outputs into one cohort count matrix.",
                    "suggested_arguments": {
                        "path": str(star_dir),
                        "project": project_name,
                    },
                    "expected_output": str(inferred_merged),
                }
            )
            star_merged_matrix = inferred_merged

        if star_merged_matrix is not None:
            star_steps.append(
                {
                    "command": "count2tpm",
                    "why": "Convert the merged STAR count matrix back into TPM before retrying the downstream command.",
                    "suggested_arguments": {
                        "input": str(star_merged_matrix),
                        "output": str(tpm_paths["star_intermediate"]),
                        "idtype": "ensembl",
                        "org": "hsa",
                        "remove_version": True,
                    },
                    "expected_output": str(tpm_paths["star_intermediate"]),
                }
            )
            star_steps.append(
                {
                    "command": "log2_eset",
                    "why": "Standardize the regenerated TPM with log2(x+1), because downstream IOBRpy functions should consume the log2_eset output rather than the raw regenerated TPM.",
                    "suggested_arguments": {
                        "input": str(tpm_paths["star_intermediate"]),
                        "output": str(tpm_paths["standardized_tpm"]),
                    },
                    "expected_output": str(tpm_paths["standardized_tpm"]),
                }
            )
            retry_step = {
                "command": command_name,
                "why": "Retry the downstream command on the standardized log2_eset TPM output instead of the failing current matrix.",
            }
            if target_retry_arguments:
                retry_step["suggested_arguments"] = target_retry_arguments
            star_steps.append(retry_step)
            recovery_paths.append(
                {
                    "id": f"{command_name}:star_regeneration",
                    "title": f"Rebuild standardized TPM from STAR outputs, then retry {command_name}",
                    "preferred": salmon_merged_matrix is None and salmon_dir is None,
                    "why_this_path_exists": (
                        "Existing STAR outputs were detected, so IOBRpy can roll back through count merging and TPM regeneration, then standardize the new TPM with log2_eset before retrying the downstream command."
                    ),
                    "detected_upstream_paths": detected_paths,
                    "standardized_tpm_output": str(tpm_paths["standardized_tpm"]),
                    "steps": star_steps,
                }
            )

    return recovery_paths


def _project_root_from_map_payload(map_payload: Dict[str, Any]) -> Optional[Path]:
    path_value = map_payload.get("path")
    if not isinstance(path_value, str) or not path_value:
        return None
    candidate = Path(path_value)
    if candidate.is_dir():
        return candidate
    return candidate.parent if candidate.parent else None


def _pathlike_argument_paths(arguments: Dict[str, Any]) -> List[Path]:
    paths: List[Path] = []
    pathish_keys = {
        "input",
        "output",
        "path",
        "path_salmon",
        "fastq",
        "fqdir",
        "bam",
        "bam_dir",
        "read1",
        "read2",
        "sc_dat",
        "cell_state_labels",
        "cell_type_labels",
    }
    for key, value in arguments.items():
        if key not in pathish_keys:
            continue
        values = value if isinstance(value, list) else [value]
        for item in values:
            if not isinstance(item, str) or not item:
                continue
            paths.append(Path(item))
    return paths


def _score_pipeline_root_candidate(candidate: Path) -> int:
    score = 0
    try:
        if not candidate.exists() or not candidate.is_dir():
            return score
    except OSError:
        return score
    for marker in PIPELINE_ROOT_MARKERS:
        try:
            if (candidate / marker).exists():
                score += 3
        except OSError:
            continue
    try:
        if any(candidate.glob("*_merged_salmon_tpm.tsv*")):
            score += 4
        if any(candidate.glob("*.STAR.count.tsv*")):
            score += 4
        if any(candidate.glob("*tpm*.tsv*")) or any(candidate.glob("*tpm*.csv*")):
            score += 2
    except OSError:
        return score
    return score


def _guess_project_root_from_arguments(arguments: Dict[str, Any]) -> Optional[Path]:
    best_candidate: Optional[Path] = None
    best_score = -1
    for argument_path in _pathlike_argument_paths(arguments):
        start = argument_path if argument_path.is_dir() else argument_path.parent
        lineage = [start, *list(start.parents)[:6]]
        for candidate in lineage:
            score = _score_pipeline_root_candidate(candidate)
            if score > best_score:
                best_candidate = candidate
                best_score = score
    if best_candidate is not None and best_score > 0:
        return best_candidate
    fallback_paths = _pathlike_argument_paths(arguments)
    if fallback_paths:
        fallback = fallback_paths[0]
        return fallback if fallback.is_dir() else fallback.parent
    return None


def _native_failure_classification(command_name: str, stdout: str, stderr: str) -> Dict[str, Any]:
    combined = f"{stdout}\n{stderr}".lower()
    if "cannot save file into a non-existent directory" in combined:
        return {
            "id": "missing_output_directory",
            "summary": "The command failed because the parent output directory does not exist yet.",
            "retryable": True,
        }
    if command_name == "bayesprism":
        if "gene names of reference and mixture do not match" in combined:
            return {
                "id": "bayesprism_gene_name_mismatch",
                "summary": "BayesPrism rejected the current matrix because its gene names do not match the reference gene space closely enough.",
                "retryable": True,
            }
        if "zero-size array to reduction operation maximum" in combined or "a total of 0 gene" in combined:
            return {
                "id": "bayesprism_empty_or_filtered_matrix",
                "summary": "BayesPrism ended up with an effectively unusable mixture matrix after preprocessing and filtering.",
                "retryable": True,
            }
    if command_name in DOWNSTREAM_MATRIX_COMMANDS and (
        "gene names" in combined
        or "index" in combined
        or "empty" in combined
        or "no identity" in combined
        or "zero-size array" in combined
    ):
        return {
            "id": "downstream_matrix_input_incompatible",
            "summary": "The downstream command appears to have failed on the current matrix representation rather than on a missing executable.",
            "retryable": True,
        }
    return {
        "id": "native_command_failed",
        "summary": "The native command failed and may need parameter review or an upstream input fix.",
        "retryable": False,
    }


def _native_failure_guidance(
    command_name: str,
    arguments: Dict[str, Any],
    stdout: str,
    stderr: str,
) -> Dict[str, Any]:
    classification = _native_failure_classification(command_name, stdout, stderr)
    project_root = _guess_project_root_from_arguments(arguments)
    recovery_paths = _downstream_recovery_paths_from_project_root(
        command_name,
        project_root,
        original_arguments=arguments,
    )
    host_instructions: List[str] = []
    if classification["id"] == "missing_output_directory":
        host_instructions.append("Create the missing parent output directory, then retry the same native command.")
    if recovery_paths:
        host_instructions.append(
            "If the current TPM-like matrix appears unusable, follow the preferred upstream recovery path instead of stopping at the first downstream failure."
        )
        host_instructions.append(
            "When regenerating TPM from Salmon or STAR outputs, run log2_eset on the regenerated TPM before retrying downstream functions."
        )
        host_instructions.append(
            "After an upstream recovery path runs, report the generated_output_paths from auto_recovery so users can see the concrete regenerated files."
        )
    if command_name == "bayesprism":
        host_instructions.append(
            "Do not jump straight to demanding a custom single-cell reference if a standardized TPM regeneration path is still available from existing Salmon or STAR outputs."
        )
    return {
        "classification": classification,
        "project_root_guess": str(project_root) if project_root is not None else None,
        "upstream_recovery_paths": recovery_paths,
        "host_agent_instructions": host_instructions,
    }


def _command_agent_reasoning_rules(
    command_name: str,
    existing_input_candidates: List[Dict[str, Any]],
) -> List[str]:
    if command_name == "runall":
        return [
            "Prefer `runall` only when the user truly wants the full FASTQ-to-TME workflow from raw sequencing data.",
            "If the scanned path already contains TPM or downstream TME outputs, summarize resume-versus-rerun before recommending a fresh full-pipeline run.",
            "Treat `runall` as a broad wrapper that subsumes QC, quantification, TPM generation, signature scoring, six default TME methods, LR_cal, and TRUST4.",
        ]
    if command_name == "fastq_qc":
        return [
            "Keep `fastq_qc` as the QC-only step when the user wants cleaning or report generation without immediate downstream quantification.",
            "Do not suggest `fastq_qc` as the main next step if the user already has a ready TPM matrix or clearly wants downstream TME analysis.",
        ]
    if command_name in {"batch_salmon", "batch_star_count"}:
        quantifier = "Salmon" if command_name == "batch_salmon" else "STAR"
        return [
            f"Keep `{command_name}` as the quantification-only {quantifier} step when the user wants {quantifier} outputs without the rest of the full pipeline.",
            "Do not silently upgrade a quantification-only request into `runall` unless the user explicitly wants the whole workflow.",
            "If merged matrices already exist, prefer the merge or TPM-conversion step instead of rerunning quantification by default.",
        ]
    if command_name in {"merge_salmon", "merge_star_count"}:
        return [
            f"Use `{command_name}` only after per-sample quantification outputs already exist; it is a merge stage, not a raw-FASTQ entrypoint.",
            "Explain that this step creates a cohort-level matrix but does not by itself finish downstream TPM profiling or TME analysis.",
        ]
    if command_name in {"prepare_salmon", "count2tpm"}:
        return [
            f"Use `{command_name}` as the matrix-conversion step after merged quantification outputs already exist, not as a raw-FASTQ or downstream-TME entrypoint.",
            "If a ready TPM-like matrix already exists in the scanned path, default to reuse or inspection before suggesting a regeneration step.",
            "Identifier handling matters here, so the agent should keep return-feature/idtype/version-stripping choices visible instead of assuming they are irrelevant.",
        ]
    if command_name == "anno_eset":
        return [
            "Treat `anno_eset` as a preprocessing repair step for probe annotation, identifier harmonization, or duplicate-collapsing, not as a default downstream TME stage.",
            "Only mention it when the matrix identifiers or annotation source genuinely need clarification.",
        ]
    if command_name == "mouse2human_eset":
        return [
            "Treat `mouse2human_eset` as an explicit cross-species preprocessing step, not as a default requirement for human matrices.",
            "Only recommend it when the matrix appears to use mouse gene symbols or the user explicitly asks for mouse-to-human conversion.",
        ]
    if command_name == "log2_eset":
        return [
            "Treat `log2_eset` as an explicit matrix-transformation request, not as a default downstream step.",
            "Only recommend it when the user asks for a log-transformed matrix or a downstream task clearly requires that representation.",
        ]
    if command_name == "calculate_sig_score":
        return [
            "Keep `calculate_sig_score` visible when the user wants signatures or pathway scores only, without the rest of the downstream TME bundle.",
            "If the user really wants several downstream methods plus LR_cal, `tme_profile` is usually the broader and more efficient route.",
        ]
    if command_name == "tme_profile":
        return [
            "Prefer `tme_profile` when the user has a TPM-like matrix and wants the standard downstream bundle, not just one method.",
            "Do not suggest `tme_profile` if the user explicitly wants only one downstream method such as CIBERSORT, EPIC, or LR_cal.",
            "Treat it as the TPM-to-TME wrapper, not as a path-parsing or FASTQ-entrypoint command.",
            "If the current TPM-like matrix looks unusable, roll back to existing Salmon or STAR outputs, regenerate TPM, and run log2_eset before retrying downstream analysis.",
        ]
    if command_name in {"cibersort", "quantiseq", "epic", "mcpcounter", "estimate", "IPS"}:
        return [
            f"Keep `{command_name}` visible when the user explicitly wants this method rather than the broader `tme_profile` bundle.",
            "Use a detected TPM-like matrix as input rather than a directory path or raw FASTQ root.",
            "If the user wants several TME methods together, `tme_profile` is usually the cleaner recommendation.",
            "If the current TPM-like matrix fails, prefer regenerating standardized TPM from Salmon or STAR outputs and then applying log2_eset before retrying this downstream method.",
        ]
    if command_name == "LR_cal":
        return [
            "Treat `LR_cal` as the ligand-receptor branch on top of a ready TPM-like matrix, not as a quantification or deconvolution command.",
            "If ligand-receptor outputs already exist, prefer interpretation or explicit rerun reasoning instead of blindly recommending another LR_cal execution.",
            "If the current TPM-like matrix fails, prefer regenerating standardized TPM from Salmon or STAR outputs and then applying log2_eset before retrying LR_cal.",
        ]
    if command_name == "bayesprism":
        return [
            "Treat `bayesprism` as the optional Bayesian deconvolution branch outside the default six-method bundle.",
            "Default to IOBRpy's bundled reference behavior unless the user explicitly wants to override it with a custom single-cell reference.",
            "Use a detected TPM-like matrix as the bulk input and explain that BayesPrism writes its own dedicated output directory.",
            "If the current TPM-like matrix fails, first try regenerating standardized TPM from detected Salmon or STAR outputs and run log2_eset before considering a custom single-cell reference override.",
        ]
    if command_name in CLUSTERING_COMMANDS:
        return [
            "Prefer detected deconvolution matrices before very wide signature-score matrices for clustering-style requests.",
            "Keep `tme_cluster` ahead of `nmf` for ordinary subtype discovery unless the user explicitly asks for factorization-style latent programs.",
        ]
    if command_name == "extract_hla_read":
        return [
            "Use `extract_hla_read` only when the user is starting from sorted and indexed BAM/CRAM files and needs the HLA-read extraction step specifically.",
            "Do not substitute it for `spechla` or `hla_typing` when the user actually wants a direct SpecHLA run or the batch wrapper.",
        ]
    if command_name == "spechla":
        return [
            "Keep `spechla` as the per-sample paired-FASTQ SpecHLA route.",
            "Prefer the batch `hla_typing` wrapper when the user has a BAM directory and wants a cohort-scale HLA workflow.",
        ]
    if command_name == "hla_typing":
        return [
            "Prefer `hla_typing` for batch HLA typing over a BAM directory when the user wants the whole extract-plus-SpecHLA workflow.",
            "Do not collapse it into `spechla` unless the user is clearly working on one sample with paired FASTQ inputs.",
        ]
    if command_name == "trust4":
        return [
            "Treat `trust4` as the native repertoire-reconstruction command for FASTQ, BAM, or directory inputs.",
            "Do not treat external MixCR-style outputs as proof that native TRUST4 has already been run.",
            "When repertoire outputs already exist, explain reuse-versus-rerun instead of recommending a blind native rerun.",
        ]
    return []


def _command_detected_input_story(command_name: str, existing_input_candidates: List[Dict[str, Any]]) -> str:
    if not existing_input_candidates:
        return ""
    if command_name == "runall":
        return (
            "Use detected raw or cleaned FASTQ roots as the real entrypoints for `runall`. If the path already contains substantial downstream outputs, explain reuse-versus-rerun before proposing a fresh full-pipeline run."
        )
    if command_name in RAW_FASTQ_STAGE_COMMANDS:
        return (
            "Use detected raw or cleaned FASTQ roots as the sequencing entrypoint for this upstream step, and keep it separate from full-pipeline reruns unless the user explicitly wants the whole workflow."
        )
    if command_name in {"merge_salmon", "prepare_salmon"}:
        return (
            "Treat the detected Salmon-side outputs below as the direct upstream handoff. These are mid-pipeline matrix-building steps, not raw-FASTQ entrypoints."
        )
    if command_name in {"merge_star_count", "count2tpm"}:
        return (
            "Treat the detected STAR-side outputs below as the direct upstream handoff. These are count merging or TPM-conversion steps rather than full quantification reruns."
        )
    if command_name in MATRIX_PREPARATION_COMMANDS:
        return (
            "The most likely inputs are ready TPM-like matrices. These commands are preprocessing repairs, so they should only be used when annotation, species conversion, or log transformation is actually needed."
        )
    if command_name == "bayesprism":
        return (
            "The most likely bulk input is the detected TPM-like matrix below; BayesPrism can use IOBRpy's bundled reference by default, and a custom single-cell reference is only needed when overriding it."
        )
    if command_name in {"tme_profile"} | DOWNSTREAM_SINGLE_METHOD_COMMANDS:
        return (
            "The most likely direct input is the detected TPM-like matrix below. Use it as the concrete bulk-expression handoff instead of referring vaguely to a directory."
        )
    if command_name in {"tme_cluster", "nmf"}:
        return (
            "Prefer a detected deconvolution matrix first, especially CIBERSORT when it is available. Use signature-score outputs only after the user confirms which subset to cluster on, because those files can be very wide. Remember that `--features` counts only feature columns and excludes the first sample-ID column."
        )
    if command_name in HLA_COMMANDS:
        return (
            "Use the detected sequencing or HLA-result paths below to explain whether the user is resuming from existing HLA outputs or starting a fresh native HLA workflow. The exact required handoff still differs across extract_hla_read, spechla, and hla_typing."
        )
    if command_name in REPERTOIRE_COMMANDS:
        return (
            "Use the detected sequencing or repertoire-result paths below to explain whether TRUST4 should rerun natively or whether existing repertoire outputs should be reused first."
        )
    return "Use the detected inputs below when explaining what this command can run on right now."


def _command_execution_guardrails(
    command_name: str,
    existing_input_candidates: List[Dict[str, Any]],
) -> List[str]:
    if command_name == "runall":
        return [
            "Use `runall` only when the user truly wants the full FASTQ-to-TME workflow from raw sequencing inputs.",
            "Do not default to `runall` if the scanned path already contains reusable TPM or downstream TME outputs; explain resume-versus-rerun first.",
            "Confirm `--mode`, `--index`, `--project`, sequencing layout (`--se` or paired-end), and FASTQ suffix handling before execution.",
            "Create or choose a clean output root instead of pointing `runall` at an existing downstream-results directory.",
        ]
    if command_name == "fastq_qc":
        return [
            "Use `fastq_qc` only when the user wants QC and cleaning specifically, rather than the full quantification pipeline.",
            "Confirm FASTQ pairing or single-end layout and the read suffix pattern before execution.",
            "Create or verify the QC output directory before writing fastp and MultiQC results.",
        ]
    if command_name in {"batch_salmon", "batch_star_count"}:
        quantifier = "Salmon" if command_name == "batch_salmon" else "STAR"
        return [
            f"Keep `{command_name}` as the quantification-only {quantifier} step instead of silently upgrading the request to `runall`.",
            "Confirm the index path, FASTQ pairing or single-end layout, and read suffix handling before execution.",
            "If merged quantification outputs already exist, prefer merge or TPM-conversion steps before rerunning quantification by default.",
        ]
    if command_name in {"merge_salmon", "merge_star_count"}:
        return [
            f"Use `{command_name}` only after per-sample quantification outputs already exist.",
            "Explain that this step produces a cohort-level matrix, not the final downstream TME bundle by itself.",
            "Create or verify the output destination before writing the merged matrix.",
        ]
    if command_name in {"prepare_salmon", "count2tpm"}:
        return [
            f"Use `{command_name}` only after the correct merged upstream matrix has been chosen; do not treat it as a raw-FASTQ command.",
            "If a ready TPM-like matrix already exists, do not regenerate it by default unless the user explicitly wants a clean rerun.",
            "Keep identifier-handling choices visible, because feature space and version stripping can change downstream compatibility.",
        ]
    if command_name == "anno_eset":
        return [
            "Use `anno_eset` only when the matrix identifiers truly need annotation or duplicate collapsing.",
            "Confirm the annotation source and duplicate-resolution method before execution.",
            "Do not substitute `anno_eset` for TPM-conversion commands such as prepare_salmon or count2tpm.",
        ]
    if command_name == "mouse2human_eset":
        return [
            "Use `mouse2human_eset` only when the matrix genuinely uses mouse gene symbols and downstream human-space analysis is intended.",
            "Confirm whether the input is a matrix or a table with a dedicated symbol column before execution.",
            "Do not suggest this step for already-human expression matrices.",
        ]
    if command_name == "log2_eset":
        return [
            "Use `log2_eset` only on a matrix that is not already log-transformed.",
            "Do not auto-apply log2 transformation unless the user explicitly wants it or a downstream task clearly requires it.",
            "Create or verify the output path before writing the transformed matrix.",
        ]
    if command_name == "calculate_sig_score":
        return [
            "Use a TPM-like matrix as input rather than a directory path.",
            "Confirm the signature family or subset and scoring method before execution when the user does not want the full built-in collection.",
            "If the user wants deconvolution or LR_cal too, prefer `tme_profile` instead of a standalone signature-only run.",
        ]
    if command_name == "tme_profile":
        return [
            "Use `tme_profile` only with a concrete TPM-like matrix file, not with a directory path or raw FASTQ root.",
            "Prefer standalone methods when the user explicitly asks for only one downstream output.",
            "Create or verify the output directory before execution, because the bundled workflow writes several downstream subdirectories.",
        ]
    if command_name == "cibersort":
        return [
            "Use a TPM-like matrix as input, not a raw FASTQ directory.",
            "Do not default to this standalone command if the user really wants the broader `tme_profile` bundle.",
            "Change `--QN`, absolute scoring, or related options only when the user intentionally needs those CIBERSORT-specific modes.",
        ]
    if command_name == "mcpcounter":
        return [
            "Use a TPM-like matrix as input and confirm the row-identifier namespace before execution.",
            "Do not substitute MCPcounter for the broader `tme_profile` bundle when the user wants multiple methods together.",
            "Treat MCPcounter outputs as abundance-style scores rather than fractions when explaining results.",
        ]
    if command_name == "estimate":
        return [
            "Use a TPM-like matrix as input and keep platform handling visible when it matters.",
            "Do not default to the standalone command if the user wants the full downstream TME bundle.",
            "Explain that ESTIMATE produces stromal/immune-style scores rather than cell-fraction outputs.",
        ]
    if command_name == "quantiseq":
        return [
            "Use a TPM-like matrix as input and do not silently toggle quanTIseq preprocessing flags without a reason.",
            "Prefer `tme_profile` when the user wants several downstream methods together.",
            "Treat quanTIseq as a deconvolution branch, not as a clustering or signature-scoring command.",
        ]
    if command_name == "epic":
        return [
            "Use a TPM-like matrix as input and keep EPIC reference choice visible when it matters.",
            "Prefer `tme_profile` when the user wants the broader downstream bundle instead of EPIC alone.",
            "Explain that EPIC output depends on the selected reference mode rather than treating all EPIC runs as identical.",
        ]
    if command_name == "IPS":
        return [
            "Use a TPM-like matrix as input and explain that IPS is score-based rather than a fraction-style deconvolution matrix.",
            "Prefer `tme_profile` if the user wants multiple downstream methods, not only IPS.",
        ]
    if command_name == "LR_cal":
        return [
            "Use a TPM-like matrix as input, not a deconvolution matrix or clustering table.",
            "If LR outputs already exist, prefer interpretation or explicit rerun reasoning before launching another LR_cal run.",
            "Keep cancer-type context visible when it materially changes ligand-receptor scoring interpretation.",
        ]
    if command_name == "bayesprism":
        return [
            "Use a TPM-like matrix as the bulk input rather than a directory path.",
            "Default to IOBRpy's bundled reference unless the user explicitly wants to override it with a custom single-cell reference.",
            "Create or verify the dedicated output directory before execution.",
        ]
    if command_name == "tme_cluster":
        return [
            "Prefer `tme_cluster` over `nmf` for ordinary TME subtype or clustering requests unless the user explicitly asks for NMF or latent-factor decomposition.",
            "When a CIBERSORT result matrix is detected, treat it as the default clustering input before MCPcounter or a very wide signature-score matrix.",
            "For `--features`, count columns after excluding the sample-ID column. Tutorial-style ranges often start from `1`, so `1:22` means the second through twenty-third actual CSV columns.",
            "If the only obvious input is a signature-score matrix, ask the user which signature subset to cluster on before execution instead of auto-using all columns.",
            "Create or verify the parent directory of `--output` before execution, because the workflow does not create missing output parents automatically.",
        ]
    if command_name == "nmf":
        return [
            "Keep `nmf` as the factorization-oriented alternative behind `tme_cluster` unless the user explicitly wants NMF or latent-program discovery.",
            "Tutorial-style NMF examples also use CIBERSORT outputs first, with `--features` counted after excluding the sample-ID column.",
            "If the input is a signature-score matrix, ask the user which subset to factorize before execution instead of auto-using an extremely wide matrix.",
            "Create or verify the output directory path before execution.",
        ]
    if command_name == "extract_hla_read":
        return [
            "Use `extract_hla_read` only when the input is a sorted and indexed BAM/CRAM file.",
            "Confirm the reference setting before execution, because extraction behavior depends on the chosen genome context.",
            "Do not substitute this extraction step for a direct SpecHLA run or the batch hla_typing wrapper.",
        ]
    if command_name == "spechla":
        return [
            "Use `spechla` only for a single sample with paired FASTQ inputs.",
            "Confirm whether the run is exon/RNA-style or WGS-style before execution.",
            "Prefer `hla_typing` when the user has a BAM directory and wants the whole batch workflow.",
        ]
    if command_name == "hla_typing":
        return [
            "Use `hla_typing` as the batch BAM-directory wrapper rather than a single-sample FASTQ command.",
            "Do not collapse it into `spechla` or `extract_hla_read` unless the user explicitly wants to control those lower-level steps.",
            "Create or verify the output directory before execution, because the workflow writes extraction outputs, SpecHLA results, and merged summaries.",
        ]
    if command_name == "trust4":
        return [
            "Confirm whether TRUST4 should run from BAM, a BAM directory, paired FASTQ, single-end FASTQ, or a FASTQ directory before execution.",
            "Do not treat external repertoire outputs such as MixCR as proof that native TRUST4 has already been run.",
            "Create or verify the output directory before execution.",
        ]
    return []


def _command_pre_execution_questions(
    command_name: str,
    existing_input_candidates: List[Dict[str, Any]],
) -> List[Dict[str, str]]:
    if command_name == "runall":
        return [
            {
                "id": "confirm_runall_rerun_scope",
                "question": (
                    "Should this be a true fresh full-pipeline rerun from raw FASTQ, or would you rather continue from the existing intermediate results that were already detected?"
                ),
                "why_ask_now": (
                    "runall is expensive and broad, so the agent should not launch it when resume-or-reuse is the real goal."
                ),
                "recommended_default": "Reuse existing intermediate results before any full rerun unless the user explicitly wants a clean restart.",
            },
            {
                "id": "confirm_runall_mode_and_layout",
                "question": (
                    "Which quantification branch should runall use (Salmon or STAR), and are the sequencing files paired-end or single-end with any non-default suffix pattern?"
                ),
                "why_ask_now": (
                    "These choices change the quantification branch and FASTQ matching behavior across the whole pipeline."
                ),
                "recommended_default": "Use the branch and FASTQ layout that match the actual sequencing files instead of guessing.",
            },
        ]
    if command_name == "fastq_qc":
        return [
            {
                "id": "confirm_fastq_layout_for_qc",
                "question": (
                    "Are these FASTQ files paired-end or single-end, and do they use the default read suffix pattern?"
                ),
                "why_ask_now": (
                    "fastq_qc needs the correct sequencing layout and suffix interpretation before cleaning starts."
                ),
                "recommended_default": "Assume paired-end only when the file naming clearly supports it.",
            }
        ]
    if command_name in {"batch_salmon", "batch_star_count"}:
        quantifier = "Salmon" if command_name == "batch_salmon" else "STAR"
        return [
            {
                "id": f"confirm_{command_name}_scope",
                "question": (
                    f"Do you want only the {quantifier} quantification stage, or should I plan for the later merge, TPM, and downstream analysis steps too?"
                ),
                "why_ask_now": (
                    "This distinguishes a quantification-only request from a broader pipeline request."
                ),
                "recommended_default": f"Stay with `{command_name}` only when the user explicitly wants the isolated {quantifier} stage.",
            },
            {
                "id": f"confirm_{command_name}_layout",
                "question": (
                    "Are the sequencing files paired-end or single-end, and do they use a non-default read suffix pattern?"
                ),
                "why_ask_now": (
                    "Quantification commands depend on correct FASTQ pairing and naming."
                ),
                "recommended_default": "Assume the default suffix only when the filenames clearly match it.",
            },
        ]
    if command_name in {"prepare_salmon", "count2tpm"}:
        upstream_label = "merged Salmon matrix" if command_name == "prepare_salmon" else "merged count matrix"
        return [
            {
                "id": f"confirm_{command_name}_input_matrix",
                "question": (
                    f"Which detected {upstream_label} should be used as the direct input for `{command_name}`?"
                ),
                "why_ask_now": (
                    "Matrix-conversion steps should operate on the correct merged upstream table rather than a guessed file."
                ),
                "recommended_default": f"Use the most obvious detected {upstream_label} unless the user points to a different matrix.",
            }
        ]
    if command_name == "anno_eset":
        return [
            {
                "id": "confirm_annotation_source",
                "question": (
                    "Which annotation source or platform should anno_eset use, and how should duplicate probes be collapsed?"
                ),
                "why_ask_now": (
                    "Annotation source and duplicate handling materially change the resulting matrix."
                ),
                "recommended_default": "Do not guess the annotation source if the platform is not obvious from the matrix.",
            }
        ]
    if command_name == "mouse2human_eset":
        return [
            {
                "id": "confirm_mouse2human_mode",
                "question": (
                    "Is the input a matrix or a general table, and which column contains the mouse symbols if it is not a matrix?"
                ),
                "why_ask_now": (
                    "mouse2human_eset changes behavior between matrix mode and table mode."
                ),
                "recommended_default": "Use matrix mode only when the file is clearly a genes-by-samples matrix.",
            }
        ]
    if command_name == "log2_eset":
        return [
            {
                "id": "confirm_log2_needed_now",
                "question": (
                    "Do you explicitly want a log2(x+1) matrix now, or are you trying to continue directly with downstream analysis from the current TPM-like matrix?"
                ),
                "why_ask_now": (
                    "The agent should not apply log transformation by default if the current matrix is already suitable for the next intended step."
                ),
                "recommended_default": "Skip log2_eset unless the user explicitly needs the transformed matrix.",
            }
        ]
    if command_name == "calculate_sig_score":
        return [
            {
                "id": "confirm_signature_family",
                "question": (
                    "Which signature family or subset should be scored, and do you want only signature output or the broader downstream TME workflow?"
                ),
                "why_ask_now": (
                    "Signature scope changes the biological interpretation and whether calculate_sig_score should stay standalone."
                ),
                "recommended_default": "Stay with a narrower signature-only run only when the user explicitly wants scores without the rest of tme_profile.",
            }
        ]
    if command_name == "tme_profile":
        return [
            {
                "id": "confirm_tme_profile_matrix",
                "question": (
                    "Which detected TPM-like matrix should feed tme_profile, and do you want the full downstream bundle rather than only one method?"
                ),
                "why_ask_now": (
                    "tme_profile needs a concrete matrix file and should not be used when the user really wants only one downstream method."
                ),
                "recommended_default": "Use the clearest detected TPM-like matrix when the user wants the standard downstream bundle.",
            }
        ]
    if command_name == "bayesprism":
        return [
            {
                "id": "confirm_bayesprism_reference_mode",
                "question": (
                    "Should BayesPrism use IOBRpy's bundled reference behavior, or do you want to provide a custom single-cell reference?"
                ),
                "why_ask_now": (
                    "Reference mode changes whether BayesPrism stays on the simplest built-in path or becomes a custom-reference run."
                ),
                "recommended_default": "Use the bundled reference path unless the user explicitly wants a custom single-cell reference.",
            }
        ]
    if command_name == "extract_hla_read":
        return [
            {
                "id": "confirm_extract_hla_read_input",
                "question": (
                    "Are you starting from sorted and indexed BAM/CRAM inputs for HLA-read extraction, or are you actually trying to run a later HLA step?"
                ),
                "why_ask_now": (
                    "extract_hla_read is only the BAM-side preparation step and should not be confused with SpecHLA or batch hla_typing."
                ),
                "recommended_default": "Use extract_hla_read only when the user explicitly wants the BAM-to-HLA-read handoff.",
            }
        ]
    if command_name == "spechla":
        return [
            {
                "id": "confirm_spechla_single_sample_inputs",
                "question": (
                    "Are you running SpecHLA for one sample from paired FASTQ inputs, and should it use exon/RNA mode or WGS mode?"
                ),
                "why_ask_now": (
                    "SpecHLA execution changes with sample scope and exon-versus-WGS interpretation."
                ),
                "recommended_default": "Prefer spechla only for single-sample paired FASTQ workflows.",
            }
        ]
    if command_name == "hla_typing":
        return [
            {
                "id": "confirm_hla_typing_batch_scope",
                "question": (
                    "Do you want the batch hla_typing workflow over a BAM directory, or are you only trying to interpret existing HLA outputs?"
                ),
                "why_ask_now": (
                    "The batch wrapper is expensive and should not be launched if the user only wants to inspect existing HLA results."
                ),
                "recommended_default": "Prefer hla_typing only when the user truly wants a native batch rerun from BAM inputs.",
            }
        ]
    if command_name == "trust4":
        return [
            {
                "id": "confirm_trust4_mode",
                "question": (
                    "Should TRUST4 run from BAM, a BAM directory, paired FASTQ, single-end FASTQ, or a FASTQ directory?"
                ),
                "why_ask_now": (
                    "TRUST4 supports multiple input modes, so the agent should not guess the wrong entrypoint."
                ),
                "recommended_default": "Use the input mode that matches the detected or user-specified sequencing assets.",
            }
        ]
    if command_name not in {"tme_cluster", "nmf"}:
        return []

    questions: List[Dict[str, str]] = []
    deconvolution_candidates = [
        item
        for item in existing_input_candidates
        if isinstance(item, dict) and str(item.get("source_id", "")).startswith("deconvolution:")
    ]
    preferred_candidate = next(
        (item for item in deconvolution_candidates if bool(item.get("preferred"))),
        deconvolution_candidates[0] if deconvolution_candidates else None,
    )
    signature_candidate = next(
        (
            item
            for item in existing_input_candidates
            if isinstance(item, dict) and str(item.get("source_id")) == "feature_scoring"
        ),
        None,
    )

    if preferred_candidate:
        preferred_label = str(
            preferred_candidate.get("source_label")
            or preferred_candidate.get("label")
            or "the preferred detected deconvolution matrix"
        )
        other_labels = [
            str(item.get("source_label") or item.get("label") or item.get("source_id") or "another detected matrix")
            for item in deconvolution_candidates
            if item is not preferred_candidate
        ]
        matrix_phrase = preferred_label
        if other_labels:
            matrix_phrase = f"{preferred_label} instead of {', '.join(other_labels[:2])}"
        questions.append(
            {
                "id": f"confirm_{command_name}_input_matrix",
                "question": (
                    f"I detected {preferred_label} as the best current clustering matrix. "
                    f"Should I use {matrix_phrase}, or do you want a different deconvolution or curated signature matrix?"
                ),
                "why_ask_now": (
                    "Matrix choice changes the clustering feature space. The tutorials prefer CIBERSORT first when it is available, and the agent should not silently switch to another matrix without asking."
                ),
                "recommended_default": (
                    f"Use {preferred_label} as the tutorial-style default when the user has no special preference."
                ),
            }
        )

    if signature_candidate:
        questions.append(
            {
                "id": f"confirm_{command_name}_signature_subset",
                "question": (
                    "If you want to cluster on signature-score outputs instead, which signature subset or column family should be used rather than the whole matrix?"
                ),
                "why_ask_now": (
                    "Signature-score files can contain hundreds or thousands of columns, so the agent should not auto-run clustering on the whole matrix."
                ),
                "recommended_default": (
                    "Stay with the preferred detected deconvolution matrix unless the user explicitly wants a curated signature subset."
                ),
            }
        )

    return questions


def _candidate_choice_label(candidate: Dict[str, Any]) -> str:
    label = str(
        candidate.get("source_label")
        or candidate.get("label")
        or candidate.get("source_id")
        or "detected input"
    )
    paths = candidate.get("paths", [])
    if isinstance(paths, list) and paths:
        first_path = str(paths[0])
        basename = Path(first_path.rstrip("/\\")).name
        if basename and basename not in label:
            return f"{label} ({basename})"
    return label


def _numbered_confirmation_option(
    option_id: str,
    *,
    label: str,
    label_zh: str,
    recommended: bool = False,
) -> Dict[str, Any]:
    return {
        "id": option_id,
        "label": label,
        "label_zh": label_zh,
        "recommended": recommended,
    }


def _matrix_choice_options(existing_input_candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    options: List[Dict[str, Any]] = []
    for candidate in existing_input_candidates[:3]:
        if not isinstance(candidate, dict):
            continue
        label = _candidate_choice_label(candidate)
        source_id = str(candidate.get("source_id") or label)
        options.append(
            _numbered_confirmation_option(
                f"use:{source_id}",
                label=f"Use {label}",
                label_zh=f"用 {label}",
                recommended=bool(candidate.get("preferred")),
            )
        )
    options.append(
        _numbered_confirmation_option(
            "specify_other_input",
            label="Use another file I specify",
            label_zh="我来指定其他文件",
        )
    )
    return options[:4]


def _brief_confirmation_prompt(
    command_name: str,
    question: Dict[str, str],
    existing_input_candidates: List[Dict[str, Any]],
) -> Dict[str, Any]:
    question_id = str(question.get("id") or "")
    prompt = {
        "id": question_id,
        "style": "short_free_text",
        "brief_prompt": str(question.get("question") or ""),
        "brief_prompt_zh": str(question.get("question") or ""),
        "numbered_options": [],
        "recommended_default": str(question.get("recommended_default") or ""),
        "avoid_raw_flag_names": True,
    }

    if question_id == "confirm_runall_rerun_scope":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "Before we run anything broad, which path do you want?",
                "brief_prompt_zh": "继续前先确认 1 件事：你想走哪条路？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "reuse_existing_results",
                        label="Reuse existing results first",
                        label_zh="先复用现有结果（推荐）",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "rerun_current_stage",
                        label="Rerun only the current stage",
                        label_zh="只重跑当前阶段",
                    ),
                    _numbered_confirmation_option(
                        "rerun_full_pipeline",
                        label="Rerun the whole pipeline",
                        label_zh="从头重跑全流程",
                    ),
                ],
            }
        )
    elif question_id in {
        "confirm_fastq_layout_for_qc",
        "confirm_batch_salmon_layout",
        "confirm_batch_star_count_layout",
    }:
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "How should I treat the sequencing layout?",
                "brief_prompt_zh": "继续前先确认测序布局：",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "paired_end",
                        label="Paired-end",
                        label_zh="双端测序（推荐）",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "single_end",
                        label="Single-end",
                        label_zh="单端测序",
                    ),
                    _numbered_confirmation_option(
                        "inspect_suffix_first",
                        label="Check filenames first",
                        label_zh="先检查文件名再决定",
                    ),
                ],
            }
        )
    elif question_id in {
        "confirm_prepare_salmon_input_matrix",
        "confirm_count2tpm_input_matrix",
        "confirm_tme_profile_matrix",
        "confirm_tme_cluster_input_matrix",
        "confirm_nmf_input_matrix",
    }:
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "Which detected matrix should I use?",
                "brief_prompt_zh": "继续前只确认 1 件事：用哪张矩阵？",
                "numbered_options": _matrix_choice_options(existing_input_candidates),
            }
        )
    elif question_id in {"confirm_tme_cluster_signature_subset", "confirm_nmf_signature_subset"}:
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "If you want the signature-score file instead, which route should I take?",
                "brief_prompt_zh": "如果你想用 signature 文件聚类，请选：",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "stay_with_detected_deconvolution_matrix",
                        label="Stay with the detected deconvolution matrix",
                        label_zh="继续用已检测到的反卷积矩阵（推荐）",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "specify_signature_subset",
                        label="I will specify the signature subset",
                        label_zh="我来指定 signature 子集",
                    ),
                ],
            }
        )
    elif question_id == "confirm_bayesprism_reference_mode":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "How should BayesPrism choose its reference?",
                "brief_prompt_zh": "BayesPrism 这一步只确认 1 件事：参考数据怎么选？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "use_bundled_reference",
                        label="Use IOBRpy bundled reference",
                        label_zh="用 IOBRpy 内置 reference（推荐）",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "use_custom_single_cell_reference",
                        label="Use my custom single-cell reference",
                        label_zh="我提供自定义单细胞 reference",
                    ),
                    _numbered_confirmation_option(
                        "skip_bayesprism_for_now",
                        label="Skip BayesPrism for now",
                        label_zh="先跳过 BayesPrism",
                    ),
                ],
            }
        )
    elif question_id == "confirm_log2_needed_now":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "Do you want a log2 matrix now?",
                "brief_prompt_zh": "现在要不要先做 log2_eset？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "run_log2_now",
                        label="Yes, generate the log2 matrix now",
                        label_zh="要，现在就生成 log2 矩阵",
                    ),
                    _numbered_confirmation_option(
                        "continue_downstream_without_log2_now",
                        label="No, continue downstream first",
                        label_zh="先不做，继续后续分析（推荐）",
                        recommended=True,
                    ),
                ],
            }
        )
    elif question_id == "confirm_hla_typing_batch_scope":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "What do you want from HLA typing?",
                "brief_prompt_zh": "HLA 这一步你更想做哪种？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "run_batch_hla_typing",
                        label="Run the batch HLA workflow",
                        label_zh="跑批量 HLA 分型流程",
                    ),
                    _numbered_confirmation_option(
                        "inspect_existing_hla_results",
                        label="Only inspect existing HLA results",
                        label_zh="只看现有 HLA 结果（推荐）",
                        recommended=True,
                    ),
                ],
            }
        )
    elif question_id == "confirm_extract_hla_read_input":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "Which HLA route do you mean?",
                "brief_prompt_zh": "你现在想走哪条 HLA 路线？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "extract_from_bam",
                        label="Extract HLA reads from BAM/CRAM",
                        label_zh="先从 BAM/CRAM 提取 HLA reads",
                    ),
                    _numbered_confirmation_option(
                        "run_later_hla_step",
                        label="I actually want the later HLA step",
                        label_zh="我其实是想跑后面的 HLA 步骤",
                    ),
                ],
            }
        )
    elif question_id == "confirm_spechla_single_sample_inputs":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "How should I run SpecHLA?",
                "brief_prompt_zh": "SpecHLA 这一步怎么跑？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "paired_fastq_exon_rna",
                        label="Single sample paired FASTQ in exon/RNA mode",
                        label_zh="单样本双端 FASTQ，exon/RNA 模式",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "paired_fastq_wgs",
                        label="Single sample paired FASTQ in WGS mode",
                        label_zh="单样本双端 FASTQ，WGS 模式",
                    ),
                ],
            }
        )
    elif question_id == "confirm_trust4_mode":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "Which TRUST4 input mode should I use?",
                "brief_prompt_zh": "TRUST4 用哪种输入模式？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "bam_directory",
                        label="BAM directory",
                        label_zh="BAM 目录",
                    ),
                    _numbered_confirmation_option(
                        "paired_fastq",
                        label="Paired FASTQ",
                        label_zh="双端 FASTQ",
                        recommended=True,
                    ),
                    _numbered_confirmation_option(
                        "single_fastq",
                        label="Single-end FASTQ",
                        label_zh="单端 FASTQ",
                    ),
                    _numbered_confirmation_option(
                        "fastq_directory",
                        label="FASTQ directory",
                        label_zh="FASTQ 目录",
                    ),
                ],
            }
        )
    elif question_id == "confirm_signature_family":
        prompt.update(
            {
                "style": "numbered_options",
                "brief_prompt": "What do you want from signature scoring?",
                "brief_prompt_zh": "特征评分这一步你想怎么做？",
                "numbered_options": [
                    _numbered_confirmation_option(
                        "signature_only",
                        label="Only generate signature scores",
                        label_zh="只生成特征评分",
                    ),
                    _numbered_confirmation_option(
                        "full_tme_bundle",
                        label="Run the full downstream TME bundle instead",
                        label_zh="改为跑完整 TME 下游流程",
                        recommended=True,
                    ),
                ],
            }
        )

    return prompt


def _parameter_confirmation_prompt(parameter: Dict[str, Any]) -> Dict[str, Any]:
    parameter_name = str(parameter.get("parameter") or parameter.get("display_name") or "parameter")
    display_name = str(parameter.get("display_name") or parameter_name)
    description = str(parameter.get("description") or "")
    return {
        "id": f"confirm_parameter:{parameter_name}",
        "style": "short_free_text",
        "brief_prompt": f"Before execution, please confirm {display_name}.",
        "brief_prompt_zh": f"执行前再确认一下：{display_name}。",
        "numbered_options": [],
        "recommended_default": description,
        "avoid_raw_flag_names": False,
    }


def _resource_confirmation_prompt(parameter: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    parameter_name = str(parameter.get("parameter") or "").strip()
    if not _is_resource_sensitive_parameter(parameter_name):
        return None
    display_name = str(parameter.get("display_name") or parameter_name)
    return {
        "id": f"confirm_resource:{parameter_name}",
        "style": "numbered_options",
        "brief_prompt": (
            f"{display_name} controls runtime and resource use. First inspect CPU cores, load, and free "
            "memory on the execution host, then choose a sensible value instead of blindly using the default."
        ),
        "brief_prompt_zh": "执行前先看运行机器的 CPU/内存占用，再确认线程或并行数；不要直接套默认值。",
        "numbered_options": [
            _numbered_confirmation_option(
                "probe_then_recommend",
                label="Check resources, then recommend a value",
                label_zh="先查资源，再建议数值（推荐）",
                recommended=True,
            ),
            _numbered_confirmation_option(
                "use_conservative_default",
                label="Use a conservative default",
                label_zh="用保守默认值",
            ),
            _numbered_confirmation_option(
                "user_will_specify",
                label="I will specify the value",
                label_zh="我来指定数值",
            ),
        ],
        "recommended_default": "Probe CPU cores, current load, and free memory on the same host before execution.",
        "avoid_raw_flag_names": True,
        "resource_probe_recommended": True,
        "resource_probe_scope": "execution_host",
        "resource_probe_commands": list(RESOURCE_PROBE_COMMANDS),
    }


def _command_concise_confirmation_prompts(
    command_name: str,
    existing_input_candidates: List[Dict[str, Any]],
    confirm_parameters: List[Dict[str, Any]],
    important_optional_parameters: List[Dict[str, Any]],
    pre_execution_questions: List[Dict[str, str]],
) -> List[Dict[str, Any]]:
    prompts: List[Dict[str, Any]] = []
    seen_ids: set[str] = set()

    def add_prompt(prompt: Dict[str, Any]) -> None:
        prompt_id = str(prompt.get("id") or "")
        if prompt_id and prompt_id not in seen_ids:
            prompts.append(prompt)
            seen_ids.add(prompt_id)

    pre_prompts = [
        _brief_confirmation_prompt(command_name, question, existing_input_candidates)
        for question in pre_execution_questions
    ]
    resource_prompts: List[Dict[str, Any]] = []
    resource_seen: set[str] = set()
    for parameter in [*confirm_parameters, *important_optional_parameters]:
        if not isinstance(parameter, dict):
            continue
        prompt = _resource_confirmation_prompt(parameter)
        if not prompt:
            continue
        prompt_id = str(prompt.get("id") or "")
        if prompt_id and prompt_id not in resource_seen:
            resource_prompts.append(prompt)
            resource_seen.add(prompt_id)

    if pre_prompts:
        add_prompt(pre_prompts[0])
        if resource_prompts:
            add_prompt(resource_prompts[0])
        for prompt in pre_prompts[1:]:
            add_prompt(prompt)
    else:
        for prompt in resource_prompts:
            add_prompt(prompt)

    if not prompts:
        for parameter in confirm_parameters[:2]:
            if not isinstance(parameter, dict):
                continue
            prompt = _parameter_confirmation_prompt(parameter)
            add_prompt(prompt)

    return prompts


def _command_guidance_card(
    command_name: str,
    *,
    map_payload: Optional[Dict[str, Any]] = None,
    reason: str,
    category: str,
    priority: int,
    mention_only_if_needed: bool = False,
    relevance_trigger: Optional[str] = None,
) -> Dict[str, Any]:
    meta = get_command_metadata(command_name)
    profile = _command_profile(command_name)
    effective_map_payload = map_payload or {}
    existing_input_candidates = _command_existing_input_candidates(command_name, effective_map_payload)
    confirm_parameters = list(meta.get("confirm_parameters", []))
    important_optional_parameters = list(meta.get("important_optional_parameters", []))
    pre_execution_questions = _command_pre_execution_questions(command_name, existing_input_candidates)
    upstream_recovery_paths = _downstream_recovery_paths_from_project_root(
        command_name,
        _project_root_from_map_payload(effective_map_payload),
    )
    return {
        "command": command_name,
        "category": category,
        "priority": priority,
        "summary": meta.get("summary", ""),
        "why_recommended": reason,
        "required_inputs": _required_parameter_cards(meta),
        "input_expectations": list(meta.get("input_expectations", [])),
        "expected_outputs": list(meta.get("outputs", [])),
        "confirm_parameters": confirm_parameters,
        "important_optional_parameters": important_optional_parameters,
        "project_context_notes": list(profile.get("project_context_notes", [])),
        "existing_input_candidates": existing_input_candidates,
        "upstream_recovery_paths": upstream_recovery_paths,
        "detected_input_story": _command_detected_input_story(command_name, existing_input_candidates),
        "agent_reasoning_rules": _command_agent_reasoning_rules(command_name, existing_input_candidates),
        "execution_guardrails": _command_execution_guardrails(command_name, existing_input_candidates),
        "pre_execution_questions": pre_execution_questions,
        "concise_confirmation_prompts": _command_concise_confirmation_prompts(
            command_name,
            existing_input_candidates,
            confirm_parameters,
            important_optional_parameters,
            pre_execution_questions,
        ),
        "mention_only_if_needed": mention_only_if_needed,
        "relevance_trigger": relevance_trigger,
    }


def _append_command_guidance(
    bucket: List[Dict[str, Any]],
    seen: set[str],
    command_name: Optional[str],
    *,
    map_payload: Optional[Dict[str, Any]] = None,
    reason: str,
    category: str,
    priority: int,
    mention_only_if_needed: bool = False,
    relevance_trigger: Optional[str] = None,
) -> None:
    if not command_name or command_name in seen:
        return
    try:
        card = _command_guidance_card(
            command_name,
            map_payload=map_payload,
            reason=reason,
            category=category,
            priority=priority,
            mention_only_if_needed=mention_only_if_needed,
            relevance_trigger=relevance_trigger,
        )
    except KeyError:
        return
    bucket.append(card)
    seen.add(command_name)


def _append_enriched_command_guidance(
    bucket: List[Dict[str, Any]],
    seen: set[str],
    command_name: Optional[str],
    *,
    map_payload: Optional[Dict[str, Any]] = None,
    reason: str,
    category: str,
    priority: int,
    extras: Optional[Dict[str, Any]] = None,
    mention_only_if_needed: bool = False,
    relevance_trigger: Optional[str] = None,
) -> None:
    if not command_name or command_name in seen:
        return
    try:
        card = _command_guidance_card(
            command_name,
            map_payload=map_payload,
            reason=reason,
            category=category,
            priority=priority,
            mention_only_if_needed=mention_only_if_needed,
            relevance_trigger=relevance_trigger,
        )
    except KeyError:
        return
    if extras:
        card.update(extras)
    bucket.append(card)
    seen.add(command_name)


def _sorted_guidance_cards(cards: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return sorted(cards, key=lambda item: (int(item.get("priority", 999)), str(item.get("command", ""))))


def _task_text_tokens(text: str) -> str:
    return (text or "").strip().lower()


def _task_mentions_any(text: str, keywords: Sequence[str]) -> bool:
    haystack = _task_text_tokens(text)
    return any(keyword.lower() in haystack for keyword in keywords)


def _has_matrix_context(map_payload: Dict[str, Any], recommendation: Dict[str, Any]) -> bool:
    completed = set(map_payload.get("completed_stages", []))
    current_stage = map_payload.get("current_stage")
    recommended_command = recommendation.get("recommended_command")
    return (
        "tpm_matrix" in completed
        or current_stage in {"tpm_matrix", "signature_scoring", "deconvolution", "ligand_receptor", "clustering"}
        or recommended_command in {
            "tme_profile",
            "calculate_sig_score",
            "cibersort",
            "quantiseq",
            "epic",
            "mcpcounter",
            "estimate",
            "IPS",
            "LR_cal",
            "tme_cluster",
            "nmf",
            "bayesprism",
        }
    )


def _task_signal_matches(command_name: str, task_text: str) -> bool:
    spec = TASK_COMMAND_SIGNAL_SPECS.get(command_name, {})
    return bool(spec) and _task_mentions_any(task_text, spec.get("keywords", []))


def _second_stage_reason_details(map_payload: Dict[str, Any], task: str) -> List[Dict[str, str]]:
    details: List[Dict[str, str]] = []
    seen: set[str] = set()

    def add(reason_id: str) -> None:
        if reason_id in seen:
            return
        summary = SECOND_STAGE_REASON_SUMMARIES.get(reason_id)
        if not summary:
            return
        details.append({"id": reason_id, "summary": summary})
        seen.add(reason_id)

    current_stage_confidence = str(map_payload.get("current_stage_confidence") or "none")
    if current_stage_confidence in {"low", "none"}:
        add("low_current_stage_confidence")
    if map_payload.get("scan_warning"):
        add("scan_warning_present")
    if bool(map_payload.get("agent_fallback", {}).get("recommended")):
        add("agent_fallback_recommended")
    if map_payload.get("external_analysis_hints"):
        add("external_analysis_hints_present")
    if map_payload.get("likely_iobrpy_functions"):
        add("likely_iobrpy_results_detected")
    if map_payload.get("reusable_result_functions"):
        add("reusable_existing_results_detected")
    investigation = map_payload.get("existing_result_investigation", {})
    if isinstance(investigation, dict) and investigation.get("target_summaries"):
        add("existing_result_investigation_available")

    task_text = _task_text_tokens(task)
    if any(_task_signal_matches(command_name, task_text) for command_name in TASK_COMMAND_SIGNAL_SPECS):
        add("task_mentions_specialized_or_nondefault_functions")

    return details


def _directory_evidence_samples(map_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    samples: List[Dict[str, Any]] = []
    detection_priority = {
        "confirmed_by_content": 0,
        "reusable_result": 1,
        "likely_iobrpy_result": 2,
    }
    detections = [
        item
        for item in map_payload.get("function_detections", [])
        if isinstance(item, dict)
        and item.get("status") in detection_priority
        and item.get("evidence")
    ]
    detections.sort(key=lambda item: (detection_priority.get(str(item.get("status")), 99), str(item.get("id", ""))))
    for item in detections[:4]:
        samples.append(
            {
                "kind": "function_detection",
                "id": item.get("id"),
                "label": item.get("label", item.get("id")),
                "recognition_status": item.get("status"),
                "status_label": item.get("status_label", item.get("status")),
                "paths": list(item.get("evidence", []))[:3],
                "note": item.get("note", ""),
            }
        )

    for hint in map_payload.get("external_analysis_hints", [])[:2]:
        if not isinstance(hint, dict):
            continue
        samples.append(
            {
                "kind": "external_analysis_hint",
                "id": hint.get("id"),
                "label": hint.get("label", hint.get("id")),
                "description": hint.get("description", ""),
                "paths": list(hint.get("evidence", []))[:3],
                "analysis_family_label": hint.get("analysis_family_label"),
                "confidence": hint.get("confidence"),
                "matched_signal_labels": list(hint.get("matched_signal_labels", []))[:4],
                "result_kind_labels": list(hint.get("result_kind_labels", []))[:3],
                "likely_tool_families": list(hint.get("likely_tool_families", []))[:3],
            }
        )

    investigation = map_payload.get("existing_result_investigation", {})
    if isinstance(investigation, dict):
        for summary in investigation.get("target_summaries", [])[:2]:
            if not isinstance(summary, dict) or not summary.get("finding_count"):
                continue
            samples.append(
                {
                    "kind": "existing_result_investigation",
                    "id": summary.get("id"),
                    "label": summary.get("label", summary.get("id")),
                    "finding_count": summary.get("finding_count", 0),
                    "bucket_counts": dict(summary.get("bucket_counts", {})),
                    "paths": list(summary.get("sample_paths", []))[:3],
                    "external_result_families": list(summary.get("external_result_families", []))[:3],
                    "external_result_kind_labels": list(summary.get("external_result_kind_labels", []))[:3],
                    "external_likely_tool_families": list(summary.get("external_likely_tool_families", []))[:3],
                }
            )

    return samples[:8]


def _external_result_interpretation_hints(map_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    hints: List[Dict[str, Any]] = []
    for item in map_payload.get("external_analysis_hints", [])[:4]:
        if not isinstance(item, dict):
            continue
        hints.append(
            {
                "source": "external_analysis_hint",
                "id": item.get("id"),
                "analysis_family_label": item.get("analysis_family_label", item.get("label", item.get("id"))),
                "confidence": item.get("confidence"),
                "family_summary": item.get("family_summary", item.get("description", "")),
                "matched_signal_labels": list(item.get("matched_signal_labels", []))[:5],
                "result_kind_labels": list(item.get("result_kind_labels", []))[:4],
                "likely_tool_families": list(item.get("likely_tool_families", []))[:4],
                "evidence_paths": list(item.get("evidence", []))[:3],
                "native_equivalence": item.get("native_equivalence", "external_family_only"),
                "reuse_guidance": (
                    "Treat this as external analysis-family evidence that may be reusable for interpretation, but do not claim the equivalent native iobrpy immune command already ran."
                ),
            }
        )

    investigation = map_payload.get("existing_result_investigation", {})
    if isinstance(investigation, dict):
        for summary in investigation.get("target_summaries", [])[:3]:
            if not isinstance(summary, dict):
                continue
            families = list(summary.get("external_result_families", []))
            if not families:
                continue
            hints.append(
                {
                    "source": "existing_result_investigation",
                    "id": summary.get("id"),
                    "analysis_family_label": ", ".join(families[:2]),
                    "confidence": "medium",
                    "family_summary": (
                        "Targeted existing-result investigation found external immune-side result families that deserve interpretation before deciding whether to reuse or rerun."
                    ),
                    "matched_signal_labels": list(summary.get("external_signal_labels", []))[:5],
                    "result_kind_labels": list(summary.get("external_result_kind_labels", []))[:4],
                    "likely_tool_families": list(summary.get("external_likely_tool_families", []))[:4],
                    "evidence_paths": list(summary.get("sample_paths", []))[:3],
                    "native_equivalence": "external_family_only",
                    "reuse_guidance": (
                        "Use these findings to explain what family of HLA/TCR result is already present, but keep native trust4/hla_typing provenance separate."
                    ),
                }
            )

    deduped: List[Dict[str, Any]] = []
    seen_keys: set[tuple[str, str]] = set()
    for item in hints:
        key = (str(item.get("source")), str(item.get("analysis_family_label")))
        if key in seen_keys:
            continue
        deduped.append(item)
        seen_keys.add(key)
    return deduped[:6]


def _candidate_function_hypotheses(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    task: str,
) -> List[Dict[str, Any]]:
    hypotheses: List[Dict[str, Any]] = []
    seen: set[str] = set()
    task_text = _task_text_tokens(task)
    detection_by_id = {
        str(item.get("id")): item
        for item in map_payload.get("function_detections", [])
        if isinstance(item, dict) and item.get("id")
    }
    status_priority = {
        "confirmed_by_content": 10,
        "reusable_result": 20,
        "likely_iobrpy_result": 30,
    }
    status_reason = {
        "confirmed_by_content": (
            "The scanned directory contains direct filename or content evidence for outputs associated with this command."
        ),
        "reusable_result": (
            "The directory appears to contain reusable results associated with this command, even if they were not promoted to fully confirmed native outputs."
        ),
        "likely_iobrpy_result": (
            "The directory looks like it may contain iobrpy-style results for this command, but the evidence is weaker and should be explained as a hypothesis."
        ),
    }

    detections = [
        item
        for item in detection_by_id.values()
        if item.get("status") in status_priority
    ]
    detections.sort(key=lambda item: (status_priority.get(str(item.get("status")), 99), str(item.get("id", ""))))
    for item in detections:
        _append_enriched_command_guidance(
            hypotheses,
            seen,
            str(item.get("id")),
            map_payload=map_payload,
            reason=status_reason[str(item.get("status"))],
            category="directory_recognition_hypothesis",
            priority=status_priority[str(item.get("status"))],
            extras={
                "recognition_source": "function_detection",
                "recognition_status": item.get("status"),
                "detection_note": item.get("note", ""),
                "evidence_paths": list(item.get("evidence", []))[:3],
            },
        )

    for idx, hint in enumerate(map_payload.get("external_analysis_hints", [])):
        if not isinstance(hint, dict):
            continue
        command_name = EXTERNAL_HINT_COMMAND_MAP.get(str(hint.get("id")))
        if not command_name:
            continue
        _append_enriched_command_guidance(
            hypotheses,
            seen,
            command_name,
            map_payload=map_payload,
            reason=(
                "The directory naming or content suggests this analysis family may matter here, but it should be framed as an external or family-level hint rather than proof of native iobrpy execution."
            ),
            category="external_hint_hypothesis",
            priority=40 + idx,
            extras={
                "recognition_source": "external_analysis_hint",
                "hint_id": hint.get("id"),
                "hint_label": hint.get("label", hint.get("id")),
                "evidence_paths": list(hint.get("evidence", []))[:3],
                "analysis_family_label": hint.get("analysis_family_label"),
                "external_hint_confidence": hint.get("confidence"),
                "external_result_kind_labels": list(hint.get("result_kind_labels", []))[:3],
                "external_signal_labels": list(hint.get("matched_signal_labels", []))[:4],
                "likely_tool_families": list(hint.get("likely_tool_families", []))[:4],
                "family_summary": hint.get("family_summary"),
            },
        )

    for idx, (command_name, spec) in enumerate(TASK_COMMAND_SIGNAL_SPECS.items()):
        if not _task_signal_matches(command_name, task_text):
            continue
        detection = detection_by_id.get(command_name, {})
        _append_enriched_command_guidance(
            hypotheses,
            seen,
            command_name,
            map_payload=map_payload,
            reason=str(spec.get("reason", "")),
            category="task_signal_hypothesis",
            priority=60 + idx,
            extras={
                "recognition_source": "task_signal",
                "task_signal_keywords": list(spec.get("keywords", [])),
                "evidence_paths": list(detection.get("evidence", []))[:3],
                "detection_note": detection.get("note", ""),
            },
            mention_only_if_needed=command_name in {"anno_eset", "mouse2human_eset", "log2_eset"} and not detection,
            relevance_trigger=(
                str(spec.get("question_why", ""))
                if command_name in {"anno_eset", "mouse2human_eset", "log2_eset"}
                else None
            ),
        )

    if _has_matrix_context(map_payload, recommendation):
        for offset, command_name in enumerate(["anno_eset", "mouse2human_eset", "log2_eset"]):
            if command_name in seen:
                continue
            spec = TASK_COMMAND_SIGNAL_SPECS.get(command_name, {})
            _append_enriched_command_guidance(
                hypotheses,
                seen,
                command_name,
                map_payload=map_payload,
                reason=(
                    "This specialized matrix-preparation step is not part of the narrow default next-step list, but it becomes relevant when the user needs matrix cleanup before downstream analysis."
                ),
                category="latent_matrix_preparation_hypothesis",
                priority=120 + offset,
                extras={
                    "recognition_source": "matrix_context",
                    "evidence_paths": [],
                },
                mention_only_if_needed=True,
                relevance_trigger=str(spec.get("question_why", "")) or "Mention only if the matrix truly needs this preprocessing step.",
            )

    return _sorted_guidance_cards(hypotheses)


def _suggested_user_questions(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    task: str,
) -> List[Dict[str, str]]:
    questions: List[Dict[str, str]] = []
    seen: set[str] = set()

    def add(question_id: str, question: str, why_ask_now: str) -> None:
        if question_id in seen:
            return
        questions.append(
            {
                "id": question_id,
                "question": question,
                "why_ask_now": why_ask_now,
            }
        )
        seen.add(question_id)

    current_stage_confidence = str(map_payload.get("current_stage_confidence") or "none")
    if current_stage_confidence in {"low", "none"} or map_payload.get("scan_warning"):
        add(
            "confirm_expected_existing_outputs",
            "What outputs do you expect this directory to already contain, so I can check the ambiguous parts against your expectation?",
            "The deterministic scan still carries uncertainty, so expected outputs are the fastest way to disambiguate the directory state.",
        )

    if map_payload.get("external_analysis_hints"):
        hinted_labels = ", ".join(
            str(item.get("label", item.get("id")))
            for item in map_payload.get("external_analysis_hints", [])[:2]
            if isinstance(item, dict)
        )
        add(
            "confirm_external_reuse_or_rerun",
            f"Should the detected external-style results ({hinted_labels}) be reused as existing outputs, or do you want to rerun the equivalent native iobrpy step?",
            "External hints are informative, but they should not be silently treated as native iobrpy results.",
        )

    if map_payload.get("reusable_result_functions") or map_payload.get("likely_iobrpy_functions"):
        add(
            "resume_or_rerun_existing_results",
            "Do you want to continue from the detected intermediate results, or rerun from an earlier stage for a cleaner pipeline state?",
            "The directory appears to contain reusable or likely-iobrpy results, so resume-versus-rerun changes the best next command.",
        )

    task_text = _task_text_tokens(task)
    for command_name, spec in TASK_COMMAND_SIGNAL_SPECS.items():
        if _task_signal_matches(command_name, task_text):
            add(
                str(spec.get("question_id", f"confirm_{command_name}")),
                str(spec.get("question", "")),
                str(spec.get("question_why", "")),
            )

    current_stage = str(map_payload.get("current_stage") or "")
    next_stages = set(map_payload.get("next_stages", []))
    if (
        "clustering" in next_stages
        or current_stage in {"deconvolution", "ligand_receptor", "clustering"}
        or _task_signal_matches("tme_cluster", task_text)
        or _task_signal_matches("nmf", task_text)
    ):
        clustering_candidates = _command_existing_input_candidates("tme_cluster", map_payload)
        for item in _command_pre_execution_questions("tme_cluster", clustering_candidates):
            add(
                str(item.get("id") or "confirm_clustering_matrix"),
                str(item.get("question") or ""),
                str(item.get("why_ask_now") or ""),
            )

    confirm_parameters = recommendation.get("confirm_parameters", [])
    if confirm_parameters:
        parameter_labels: List[str] = []
        for item in confirm_parameters:
            if isinstance(item, dict):
                label = (
                    item.get("display_name")
                    or item.get("parameter")
                    or item.get("name")
                    or item.get("id")
                )
            else:
                label = str(item)
            if label:
                parameter_labels.append(str(label))
        if parameter_labels:
            add(
                "confirm_recommendation_parameters",
                f"Before execution, can you confirm these high-impact parameters: {', '.join(parameter_labels[:3])}?",
                "These parameters materially change the chosen workflow or output interpretation.",
            )

    return questions


def _directory_recognition_payload(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    task: str,
) -> Dict[str, Any]:
    reason_details = _second_stage_reason_details(map_payload, task)
    evidence_samples = _directory_evidence_samples(map_payload)
    candidate_function_hypotheses = _candidate_function_hypotheses(map_payload, recommendation, task)
    suggested_user_questions = _suggested_user_questions(map_payload, recommendation, task)
    external_result_interpretation_hints = _external_result_interpretation_hints(map_payload)
    return {
        "layer_1_cli_scan": {
            "owner": "iobrpy-cli pipeline_map",
            "purpose": "Deterministically detect current stage, reusable outputs, and known workflow evidence.",
            "current_stage": map_payload.get("current_stage"),
            "current_label": map_payload.get("current_label"),
            "current_stage_confidence": map_payload.get("current_stage_confidence"),
            "scan_warning_present": bool(map_payload.get("scan_warning")),
            "agent_fallback_recommended": bool(map_payload.get("agent_fallback", {}).get("recommended")),
            "existing_result_investigation_available": bool(
                isinstance(map_payload.get("existing_result_investigation"), dict)
                and map_payload.get("existing_result_investigation", {}).get("target_summaries")
            ),
        },
        "layer_2_llm_reasoning_surface": {
            "enabled": True,
            "purpose": (
                "Use the host LLM to resolve directory ambiguity, explain uncertainty, and widen candidate-function judgment without replacing the deterministic pipeline_map scan."
            ),
            "should_expand_reasoning": bool(reason_details or candidate_function_hypotheses or suggested_user_questions),
            "reason_details": reason_details,
            "evidence_samples": evidence_samples,
            "candidate_function_hypotheses": candidate_function_hypotheses,
            "external_result_interpretation_hints": external_result_interpretation_hints,
            "suggested_user_questions": suggested_user_questions,
            "response_guidance": {
                "treat_pipeline_map_as_ground_truth_for_confirmed_results": True,
                "use_candidate_function_hypotheses_for_ambiguous_or_nondefault_capabilities": True,
                "use_external_result_interpretation_hints_for_hla_and_tcr_bcr_families": True,
                "explain_uncertainty_when_using_hypotheses": True,
                "do_not_claim_external_hints_are_native_iobrpy_results": True,
                "do_not_limit_function_judgment_to_the_current_stage_table": True,
            },
        },
    }


def _decision_paths_from_map_payload(map_payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    scenario = map_payload.get("scenario", {}) if isinstance(map_payload.get("scenario", {}), dict) else {}
    recommended_choice = scenario.get("recommended_choice")
    details = scenario.get("choice_details", []) if isinstance(scenario.get("choice_details", []), list) else []
    paths: List[Dict[str, Any]] = []
    for item in details:
        if not isinstance(item, dict):
            continue
        command_cards: List[Dict[str, Any]] = []
        local_seen: set[str] = set()
        for idx, invocation in enumerate(item.get("suggested_commands", []) or []):
            command_name = _extract_command_name_from_invocation(invocation)
            _append_command_guidance(
                command_cards,
                local_seen,
                command_name,
                map_payload=map_payload,
                reason=f"This command is attached to the `{item.get('label', item.get('id', 'decision'))}` decision path.",
                category="decision_path",
                priority=idx,
            )
        paths.append(
            {
                "id": item.get("id"),
                "label": item.get("label"),
                "why_this_path_exists": item.get("when_to_choose", ""),
                "recommended_for_current_state": item.get("id") == recommended_choice,
                "commands": _sorted_guidance_cards(command_cards),
            }
        )
    return paths


def _recommended_and_alternative_commands(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    task: str,
) -> Dict[str, List[Dict[str, Any]]]:
    recommended: List[Dict[str, Any]] = []
    alternatives: List[Dict[str, Any]] = []
    latent: List[Dict[str, Any]] = []
    seen: set[str] = set()

    checklist_index = _workflow_checklist_index(map_payload)
    next_stages = set(map_payload.get("next_stages", []))
    current_stage = str(map_payload.get("current_stage") or "")
    completed_stages = set(map_payload.get("completed_stages", []))
    deconv_payload = map_payload.get("deconvolution_methods", {})
    if isinstance(deconv_payload, dict):
        deconv_methods = {
            str(item)
            for item in (
                list(deconv_payload.get("detected_ids", []) or [])
                + list(deconv_payload.get("optional_detected_ids", []) or [])
            )
            if item
        }
        optional_missing_deconv_methods = {
            str(item) for item in (deconv_payload.get("optional_missing_ids", []) or []) if item
        }
    else:
        deconv_methods = {
            str(item.get("id"))
            for item in deconv_payload
            if isinstance(item, dict) and item.get("id")
        }
        optional_missing_deconv_methods = set()
    deconvolution_item = checklist_index.get("immune_deconvolution", {})
    star_merge_item = checklist_index.get("star_quant_merge", {})
    task_text = _task_text_tokens(task)
    recommended_command = recommendation.get("recommended_command")
    demote_broad_bundle_recommendation = (
        isinstance(recommended_command, str)
        and recommended_command in {"tme_profile", "runall"}
        and (
            current_stage in {"deconvolution", "ligand_receptor", "clustering"}
            or {"deconvolution", "ligand_receptor"}.issubset(completed_stages)
        )
    )

    if isinstance(recommended_command, str) and recommended_command:
        if demote_broad_bundle_recommendation:
            _append_command_guidance(
                alternatives,
                seen,
                recommended_command,
                map_payload=map_payload,
                reason=(
                    "This broad bundled workflow still matches the input type, but the directory already appears to be far downstream, so it is better treated as an optional rerun path rather than the default next step."
                ),
                category="broad_rerun_bundle",
                priority=45,
            )
        else:
            _append_command_guidance(
                recommended,
                seen,
                recommended_command,
                map_payload=map_payload,
                reason=recommendation.get("why", "This is the most likely workflow entrypoint for the detected input."),
                category="primary_recommendation",
                priority=10,
            )

    suggested_invocations = []
    for key in ("suggested_commands",):
        suggested_invocations.extend(map_payload.get(key, []) or [])
    recommended_action = map_payload.get("recommended_action", {})
    if isinstance(recommended_action, dict):
        suggested_invocations.extend(recommended_action.get("suggested_commands", []) or [])
    scenario = map_payload.get("scenario", {})
    if isinstance(scenario, dict):
        for choice in scenario.get("choice_details", []) or []:
            if isinstance(choice, dict) and choice.get("id") == scenario.get("recommended_choice"):
                suggested_invocations.extend(choice.get("suggested_commands", []) or [])
                break

    priority_seed = 20
    for invocation in suggested_invocations:
        command_name = _extract_command_name_from_invocation(invocation)
        _append_command_guidance(
            recommended,
            seen,
            command_name,
            map_payload=map_payload,
            reason="This command is directly tied to the current roadmap position and next detected step.",
            category="next_step",
            priority=priority_seed,
        )
        priority_seed += 1

    if "clustering" in next_stages or current_stage == "ligand_receptor":
        _append_command_guidance(
            recommended,
            seen,
            "tme_cluster",
            map_payload=map_payload,
            reason="Clustering is the most direct missing downstream step after the currently detected TME and ligand-receptor results.",
            category="next_step",
            priority=30,
        )
        _append_command_guidance(
            alternatives,
            seen,
            "nmf",
            map_payload=map_payload,
            reason="NMF is a good alternative when you want latent programs or decomposition-based subtypes instead of only a direct clustering run.",
            category="alternative_analysis",
            priority=40,
        )

    if "ligand_receptor" in next_stages or current_stage == "deconvolution":
        _append_command_guidance(
            recommended,
            seen,
            "LR_cal",
            map_payload=map_payload,
            reason="Ligand-receptor analysis is the next standard downstream step after deconvolution or TPM-stage analysis.",
            category="next_step",
            priority=31,
        )

    if "signature_scoring" in next_stages:
        _append_command_guidance(
            recommended,
            seen,
            "calculate_sig_score",
            map_payload=map_payload,
            reason="Signature scoring is one of the next missing downstream steps on the detected roadmap.",
            category="next_step",
            priority=32,
        )

    if "deconvolution" in next_stages or current_stage in {"tpm_matrix", "signature_scoring"}:
        for idx, command_name in enumerate(["cibersort", "quantiseq", "epic", "mcpcounter", "estimate", "IPS"]):
            _append_command_guidance(
                alternatives,
                seen,
                command_name,
                map_payload=map_payload,
                reason="This is a valid standalone downstream TME method if the user wants a narrower analysis instead of the full bundled workflow.",
                category="standalone_downstream_method",
                priority=60 + idx,
            )
        if "bayesprism" not in deconv_methods:
            _append_command_guidance(
                alternatives,
                seen,
                "bayesprism",
                map_payload=map_payload,
                reason="BayesPrism is not part of the default bundled TME methods, but it is relevant if the user wants an additional Bayesian deconvolution branch and it can use IOBRpy's bundled reference by default.",
                category="optional_downstream_method",
                priority=70,
            )

    if (
        isinstance(deconvolution_item, dict)
        and str(deconvolution_item.get("status")) == "partial"
        and "bayesprism" in optional_missing_deconv_methods
    ):
        _append_command_guidance(
            alternatives,
            seen,
            "bayesprism",
            map_payload=map_payload,
            reason="The default TME bundle is partially complete. BayesPrism is not part of that bundle, but it is still available as an optional standalone Bayesian branch that can usually reuse the existing TPM-like matrix and the bundled BayesPrism reference.",
            category="missing_optional_downstream_method",
            priority=42,
        )

    if isinstance(star_merge_item, dict) and str(star_merge_item.get("status")) == "partial":
        _append_command_guidance(
            alternatives,
            seen,
            "merge_star_count",
            map_payload=map_payload,
            reason="STAR-side outputs are already present but the merged count matrix still looks incomplete, so merge_star_count is a concrete backfill option rather than a vague upstream rerun.",
            category="backfill_upstream_gap",
            priority=43,
        )

    if _task_mentions_any(task_text, ["tcr", "bcr", "vdj", "repertoire", "clonotype"]) or any(
        hint.get("id") == "external_tcr_bcr"
        for hint in map_payload.get("external_analysis_hints", [])
        if isinstance(hint, dict)
    ):
        _append_command_guidance(
            alternatives,
            seen,
            "trust4",
            map_payload=map_payload,
            reason="The task or detected directory hints suggest TCR/BCR repertoire analysis may be relevant.",
            category="immune_branch",
            priority=80,
        )

    if _task_mentions_any(task_text, ["hla", "typing", "spechla"]) or any(
        hint.get("id") == "external_hla"
        for hint in map_payload.get("external_analysis_hints", [])
        if isinstance(hint, dict)
    ):
        _append_command_guidance(
            alternatives,
            seen,
            "hla_typing",
            map_payload=map_payload,
            reason="The task or detected directory hints suggest HLA typing may be relevant.",
            category="immune_branch",
            priority=81,
        )

    if _has_matrix_context(map_payload, recommendation):
        _append_command_guidance(
            latent,
            seen,
            "anno_eset",
            map_payload=map_payload,
            reason="This becomes relevant when the matrix needs identifier annotation or duplicate-probe collapsing before downstream scoring or deconvolution.",
            category="specialized_matrix_preparation",
            priority=100,
            mention_only_if_needed=True,
            relevance_trigger="Mention when the user needs annotation, probe-to-symbol mapping, or duplicate collapsing.",
        )
        _append_command_guidance(
            latent,
            seen,
            "mouse2human_eset",
            map_payload=map_payload,
            reason="This becomes relevant when the matrix uses mouse symbols but downstream workflows need human symbols.",
            category="specialized_matrix_preparation",
            priority=101,
            mention_only_if_needed=True,
            relevance_trigger="Mention when the user asks for mouse-to-human symbol conversion or cross-species preprocessing.",
        )
        _append_command_guidance(
            latent,
            seen,
            "log2_eset",
            map_payload=map_payload,
            reason="This becomes relevant when the user explicitly needs log2(x+1) transformed values before a downstream step or visualization.",
            category="specialized_matrix_preparation",
            priority=102,
            mention_only_if_needed=True,
            relevance_trigger="Mention when the user explicitly asks for log transformation or a downstream step requires it.",
        )

    return {
        "recommended": _sorted_guidance_cards(recommended),
        "alternatives": _sorted_guidance_cards(alternatives),
        "latent": _sorted_guidance_cards(latent),
    }


def _workflow_checklist_index(map_payload: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {
        str(item.get("id")): item
        for item in map_payload.get("workflow_checklist", [])
        if isinstance(item, dict) and item.get("id")
    }


def _user_visible_checklist_labels(map_payload: Dict[str, Any], statuses: Sequence[str], *, limit: int = 4) -> List[str]:
    labels: List[str] = []
    for item in map_payload.get("workflow_checklist", []):
        if not isinstance(item, dict):
            continue
        if str(item.get("status")) not in statuses:
            continue
        label = str(item.get("label") or item.get("label_zh") or item.get("id") or "").strip()
        if label and label not in labels:
            labels.append(label)
        if len(labels) >= limit:
            break
    return labels


def _checklist_path_examples(
    map_payload: Dict[str, Any],
    statuses: Sequence[str],
    *,
    limit_items: int = 3,
    limit_paths_per_item: int = 2,
) -> List[Dict[str, Any]]:
    examples: List[Dict[str, Any]] = []
    for item in map_payload.get("workflow_checklist", []):
        if not isinstance(item, dict):
            continue
        if str(item.get("status")) not in statuses:
            continue
        paths = _dedupe_text_values(
            item.get("evidence_display_paths", [])
            or item.get("evidence_paths", [])
            or item.get("raw_evidence_paths", []),
            limit=limit_paths_per_item,
        )
        if not paths:
            continue
        examples.append(
            {
                "id": str(item.get("id") or ""),
                "label": str(item.get("label") or item.get("label_zh") or item.get("id") or ""),
                "paths": paths,
            }
        )
        if len(examples) >= limit_items:
            break
    return examples


_STAGE_LABEL_OVERRIDES: Dict[str, str] = {
    "raw_fastq": "Raw FASTQ data",
    "fastq_qc": "Read QC and trimming",
    "salmon_quant": "Salmon quantification",
    "salmon_merge": "Salmon merge",
    "prepare_salmon": "TPM matrix preparation",
    "star_quant": "STAR quantification",
    "star_merge": "STAR count merge",
    "count2tpm": "TPM matrix generation",
    "tpm_matrix": "TPM matrix",
    "signature_scoring": "Signature scoring",
    "deconvolution": "Immune deconvolution",
    "ligand_receptor": "Ligand-receptor analysis",
    "clustering": "TME clustering or subtyping",
    "hla_typing": "HLA typing",
    "trust4": "TCR/BCR repertoire analysis",
}


def _user_visible_value_labels(
    values: Sequence[Any],
    *,
    limit: int = 4,
    treat_as_stage_ids: bool = False,
) -> List[str]:
    labels: List[str] = []
    for value in values:
        if isinstance(value, dict):
            text = (
                value.get("label")
                or value.get("label_zh")
                or value.get("title")
                or value.get("id")
                or value.get("command")
                or ""
            )
        else:
            text = value
        label = str(text or "").strip()
        if not label:
            continue
        if treat_as_stage_ids:
            label = _STAGE_LABEL_OVERRIDES.get(label, label.replace("_", " ").title())
        if label not in labels:
            labels.append(label)
        if len(labels) >= limit:
            break
    return labels


def _reason_card(
    *,
    reason_id: str,
    priority: int,
    category: str,
    detail: str,
    evidence_labels: Sequence[str] = (),
    effect_on_recommendation: str = "",
) -> Dict[str, Any]:
    return {
        "id": reason_id,
        "priority": priority,
        "category": category,
        "detail": detail,
        "evidence_labels": list(evidence_labels),
        "effect_on_recommendation": effect_on_recommendation,
    }


def _current_state_assessment(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    command_guidance: Dict[str, List[Dict[str, Any]]],
) -> Dict[str, Any]:
    scenario = map_payload.get("scenario", {}) if isinstance(map_payload.get("scenario", {}), dict) else {}
    scenario_id = str(scenario.get("id") or "")
    current_label = str(map_payload.get("current_label") or map_payload.get("current_stage") or "Unclassified")
    current_stage = str(map_payload.get("current_stage") or "")
    current_stage_confidence = str(map_payload.get("current_stage_confidence") or "none")
    completed_stages = set(map_payload.get("completed_stages", []))
    next_stages = set(map_payload.get("next_stages", []))
    if not scenario_id:
        if current_stage in {"deconvolution", "ligand_receptor", "clustering"} or {"deconvolution", "ligand_receptor"}.issubset(completed_stages):
            scenario_id = "downstream_partial" if next_stages else "downstream_complete"
        elif current_stage in {"tpm_matrix", "signature_scoring"} or "tpm_matrix" in completed_stages:
            scenario_id = "tpm_ready"
        elif current_stage in {"fastq_qc", "salmon_quant", "salmon_merge", "prepare_salmon", "star_quant", "star_merge", "count2tpm"}:
            scenario_id = "upstream_partial"
        elif "raw_fastq" in completed_stages:
            scenario_id = "raw_fastq_only"
    completed_labels = _user_visible_checklist_labels(map_payload, ["completed"], limit=5)
    partial_labels = _user_visible_checklist_labels(map_payload, ["partial"], limit=4)
    pending_labels = _user_visible_checklist_labels(map_payload, ["pending"], limit=4)
    external_labels = _user_visible_checklist_labels(map_payload, ["external"], limit=3)
    if not completed_labels:
        completed_labels = _user_visible_value_labels(map_payload.get("completed_stage_labels", []), limit=5)
    if not completed_labels and completed_stages:
        completed_labels = _user_visible_value_labels(
            list(completed_stages),
            limit=5,
            treat_as_stage_ids=True,
        )
    if not completed_labels:
        completed_labels = _user_visible_value_labels(
            map_payload.get("reusable_result_functions", []),
            limit=5,
            treat_as_stage_ids=True,
        )
    if not partial_labels:
        partial_labels = _user_visible_value_labels(map_payload.get("supplementary_stage_labels", []), limit=4)
    if not pending_labels:
        pending_labels = _user_visible_value_labels(map_payload.get("next_stage_labels", []), limit=4)
    if not pending_labels and next_stages:
        pending_labels = _user_visible_value_labels(
            list(next_stages),
            limit=4,
            treat_as_stage_ids=True,
        )
    if not external_labels:
        external_labels = _user_visible_value_labels(
            map_payload.get("external_analysis_hints", []),
            limit=3,
        )
    if not external_labels:
        investigation = map_payload.get("existing_result_investigation", {})
        if isinstance(investigation, dict):
            external_labels = _user_visible_value_labels(
                investigation.get("target_summaries", []),
                limit=3,
            )
    reusable_path_examples = _checklist_path_examples(map_payload, ["completed", "partial"], limit_items=4)
    summary_signals: List[Dict[str, Any]] = []
    profile_tags: List[str] = []

    if scenario_id in {"downstream_complete", "downstream_partial"}:
        overall_assessment = (
            "This path looks much more like an already-processed result directory than a fresh pipeline input. "
            "The main question is what to add or interpret next, not which broad workflow to start from scratch."
        )
        default_posture = "reuse_existing_results_before_any_rerun"
        default_posture_reason = (
            "Most of the expensive upstream and downstream work appears to be reusable, so broad reruns should be treated as optional rather than default."
        )
        default_posture_reason_ids = ["processed_results_dominate", "reuse_cheaper_than_restart"]
        profile_tags.extend(["processed_results", "downstream_context"])
    elif scenario_id == "tpm_ready":
        overall_assessment = (
            "This path looks ready for downstream matrix-based analysis, so the next step should usually be a focused downstream method rather than upstream quantification."
        )
        default_posture = "continue_from_matrix_stage"
        default_posture_reason = (
            "A TPM-like matrix already appears to exist, so the user can choose targeted downstream analyses directly."
        )
        default_posture_reason_ids = ["matrix_ready", "downstream_methods_preferred"]
        profile_tags.extend(["matrix_ready", "downstream_context"])
    elif scenario_id == "upstream_partial":
        overall_assessment = (
            "This path looks resumable from an intermediate upstream stage rather than empty, so continuing from the current stage is usually cheaper than restarting."
        )
        default_posture = "resume_current_pipeline_position"
        default_posture_reason = (
            "Reusable intermediate outputs appear to exist, so a full rerun would likely repeat work unnecessarily."
        )
        default_posture_reason_ids = ["intermediate_results_present", "resume_cheaper_than_restart"]
        profile_tags.extend(["resumable", "upstream_context"])
    elif scenario_id == "raw_fastq_only":
        overall_assessment = (
            "This path mainly looks like raw input data, so the default next step is to choose the right workflow entrypoint rather than interpret existing outputs."
        )
        default_posture = "choose_pipeline_entrypoint"
        default_posture_reason = (
            "There is not enough reusable processed output yet to justify a downstream-only plan."
        )
        default_posture_reason_ids = ["raw_input_only", "no_reusable_processed_outputs"]
        profile_tags.extend(["raw_input", "pipeline_entry_needed"])
    else:
        overall_assessment = (
            "The directory state is not cleanly classifiable, so the safest default is to explain uncertainty and narrow the task before executing anything broad."
        )
        default_posture = "resolve_ambiguity_first"
        default_posture_reason = (
            "The current stage or result provenance is still ambiguous enough that a one-shot workflow choice would be premature."
        )
        default_posture_reason_ids = ["directory_state_ambiguous"]
        profile_tags.append("ambiguous")

    if current_label:
        overall_assessment = f"{overall_assessment} The best deterministic roadmap position right now is `{current_label}`."

    reuse_highlights = completed_labels or _user_visible_value_labels(
        map_payload.get("completed_stage_labels", []),
        limit=5,
    )
    missing_or_optional_followups = pending_labels + partial_labels
    caution_highlights: List[str] = []
    if external_labels:
        caution_highlights.append(
            f"External-only result hints are present for: {', '.join(external_labels[:3])}."
        )
    if map_payload.get("scan_warning"):
        caution_highlights.append(str(map_payload.get("scan_warning")))
    elif map_payload.get("scan_note"):
        caution_highlights.append(str(map_payload.get("scan_note")))
    if current_stage_confidence in {"low", "none"}:
        caution_highlights.append(
            "Current-stage confidence is not strong, so suggestions should be framed with uncertainty."
        )

    next_step_bias = [
        item.get("command")
        for item in command_guidance.get("recommended", [])
        if isinstance(item, dict) and item.get("command")
    ][:4]

    if reuse_highlights:
        summary_signals.append(
            _reason_card(
                reason_id="reusable_outputs_detected",
                priority=10,
                category="positive",
                detail="Several stages already look reusable, so the answer should start from what exists rather than from a fresh workflow template.",
                evidence_labels=reuse_highlights[:4],
                effect_on_recommendation="prefer_reuse_or_fill_missing_branches",
            )
        )
    if partial_labels or pending_labels:
        summary_signals.append(
            _reason_card(
                reason_id="remaining_followup_gaps",
                priority=20,
                category="gap",
                detail="There are still unresolved or optional follow-up branches, so the next-step explanation should distinguish missing essentials from optional add-ons.",
                evidence_labels=(partial_labels + pending_labels)[:4],
                effect_on_recommendation="rank_gap_filling_above_broad_reruns",
            )
        )
    if external_labels:
        profile_tags.append("external_hints")
        summary_signals.append(
            _reason_card(
                reason_id="external_result_hints_present",
                priority=30,
                category="caution",
                detail="Some detected outputs look external rather than natively generated by iobrpy, so reuse and provenance should be discussed explicitly.",
                evidence_labels=external_labels[:3],
                effect_on_recommendation="ask_reuse_vs_rerun_for_external_results",
            )
        )
    if map_payload.get("scan_warning"):
        profile_tags.append("scan_uncertainty")
        summary_signals.append(
            _reason_card(
                reason_id="scan_warning_present",
                priority=40,
                category="caution",
                detail=str(map_payload.get("scan_warning")),
                effect_on_recommendation="keep_uncertainty_visible_in_the_summary",
            )
        )
    elif map_payload.get("scan_note"):
        summary_signals.append(
            _reason_card(
                reason_id="scan_strategy_note",
                priority=45,
                category="scan_note",
                detail=str(map_payload.get("scan_note")),
                effect_on_recommendation="explain_scan_scope_when_relevant",
            )
        )
    if current_stage_confidence in {"low", "none"}:
        profile_tags.append("low_confidence")
        summary_signals.append(
            _reason_card(
                reason_id="current_stage_confidence_low",
                priority=50,
                category="uncertainty",
                detail="Current-stage confidence is not strong, so the final summary should sound appropriately tentative rather than overly certain.",
                effect_on_recommendation="prefer_disambiguation_before_committing_to_one_workflow",
            )
        )
    if any(
        isinstance(item, dict) and item.get("category") == "broad_rerun_bundle"
        for item in command_guidance.get("alternatives", [])
    ):
        summary_signals.append(
            _reason_card(
                reason_id="broad_bundle_demoted",
                priority=60,
                category="ranking",
                detail="A broad bundled workflow still matches the input type, but it should be explained as an optional rerun path rather than the default next step.",
                effect_on_recommendation="keep_broad_reruns_below_targeted_followups",
            )
        )

    return {
        "directory_profile": scenario_id or "unknown",
        "profile_tags": profile_tags,
        "overall_assessment": overall_assessment,
        "fallback_one_line_assessment": overall_assessment,
        "default_posture": default_posture,
        "default_posture_reason": default_posture_reason,
        "default_posture_reason_ids": default_posture_reason_ids,
        "summary_signals": sorted(summary_signals, key=lambda item: (int(item.get("priority", 999)), str(item.get("id", "")))),
        "reuse_highlights": reuse_highlights,
        "reusable_path_examples": reusable_path_examples,
        "missing_or_optional_followups": missing_or_optional_followups,
        "caution_highlights": caution_highlights,
        "next_step_bias": next_step_bias,
        "recommended_conversation_order": [
            "show_workflow_checklist",
            "summarize_current_state_in_plain_language",
            "explain_why_reuse_or_resume_is_or_is_not_default",
            "present_ranked_next_step_cards",
        ],
        "do_not_default_to": [
            "broad bundled reruns when most downstream outputs already exist",
            "assuming external hints prove native iobrpy execution",
        ],
        "recommended_command_from_input_only": recommendation.get("recommended_command"),
        "llm_summary_contract": {
            "compose_summary_in_your_own_words": True,
            "prefer_summary_signals_over_fallback_one_line_assessment": True,
            "must_cover": [
                "what_state_the_directory_most_likely_represents",
                "what_can_be_reused",
                "what_is_still_missing_or_optional",
            ],
            "name_reusable_paths_when_helpful": True,
            "mention_when_present": [
                "external_result_hints",
                "scan_uncertainty",
                "low_stage_confidence",
            ],
            "avoid": [
                "reading the fallback assessment aloud verbatim",
                "jumping straight to a command before framing reuse_vs_rerun",
            ],
        },
    }


def _non_command_next_step_card(
    *,
    card_id: str,
    title: str,
    why: str,
    required_inputs: Sequence[str],
    expected_result: str,
    priority: int,
    posture: str,
    ranking_reasons: Optional[Sequence[Dict[str, Any]]] = None,
    tradeoffs: Optional[Sequence[str]] = None,
    when_not_to_choose: str = "",
    existing_input_candidates: Optional[Sequence[Dict[str, Any]]] = None,
    project_context_notes: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    return {
        "id": card_id,
        "kind": "non_command",
        "priority": priority,
        "posture": posture,
        "title": title,
        "why_this_is_recommended": why,
        "ranking_reasons": list(ranking_reasons or []),
        "required_inputs": list(required_inputs),
        "expected_result": expected_result,
        "tradeoffs": list(tradeoffs or []),
        "when_not_to_choose": when_not_to_choose,
        "existing_input_candidates": list(existing_input_candidates or []),
        "project_context_notes": list(project_context_notes or []),
    }


def _command_next_step_card(
    card: Dict[str, Any],
    *,
    priority: int,
    posture: str,
    ranking_reasons: Optional[Sequence[Dict[str, Any]]] = None,
    tradeoffs: Optional[Sequence[str]] = None,
    when_not_to_choose: str = "",
) -> Dict[str, Any]:
    return {
        "id": f"command:{card.get('command')}",
        "kind": "command",
        "priority": priority,
        "posture": posture,
        "command": card.get("command"),
        "title": str(card.get("command") or ""),
        "why_this_is_recommended": card.get("why_recommended", ""),
        "ranking_reasons": list(ranking_reasons or []),
        "required_inputs": card.get("required_inputs", []),
        "input_expectations": card.get("input_expectations", []),
        "expected_result": ", ".join(card.get("expected_outputs", [])) if card.get("expected_outputs") else "",
        "confirm_parameters": card.get("confirm_parameters", []),
        "important_optional_parameters": card.get("important_optional_parameters", []),
        "existing_input_candidates": card.get("existing_input_candidates", []),
        "upstream_recovery_paths": card.get("upstream_recovery_paths", []),
        "detected_input_story": card.get("detected_input_story", ""),
        "agent_reasoning_rules": card.get("agent_reasoning_rules", []),
        "execution_guardrails": card.get("execution_guardrails", []),
        "pre_execution_questions": card.get("pre_execution_questions", []),
        "concise_confirmation_prompts": card.get("concise_confirmation_prompts", []),
        "project_context_notes": card.get("project_context_notes", []),
        "tradeoffs": list(tradeoffs or []),
        "when_not_to_choose": when_not_to_choose,
        "category": card.get("category"),
    }


def _concise_confirmation_plan(next_step_cards: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    items: List[Dict[str, Any]] = []
    seen_prompt_ids: set[str] = set()
    for card in next_step_cards:
        if not isinstance(card, dict):
            continue
        for prompt in card.get("concise_confirmation_prompts", []):
            if not isinstance(prompt, dict):
                continue
            prompt_id = str(prompt.get("id") or "")
            if not prompt_id or prompt_id in seen_prompt_ids:
                continue
            items.append(
                {
                    **prompt,
                    "for_card_id": card.get("id"),
                    "for_command": card.get("command"),
                    "priority": card.get("priority"),
                }
            )
            seen_prompt_ids.add(prompt_id)
            if len(items) >= 4:
                break
        if len(items) >= 4:
            break

    return {
        "prefer_brief_numbered_questions": True,
        "ask_at_most_this_many_items_first": 2,
        "opening": "Before execution, ask only the one or two highest-impact confirmation questions first.",
        "opening_zh": "继续执行前，只先确认最关键的 1 到 2 件事。",
        "items": items,
    }


def _user_facing_next_step_cards(
    map_payload: Dict[str, Any],
    recommendation: Dict[str, Any],
    command_guidance: Dict[str, List[Dict[str, Any]]],
) -> List[Dict[str, Any]]:
    scenario = map_payload.get("scenario", {}) if isinstance(map_payload.get("scenario", {}), dict) else {}
    scenario_id = str(scenario.get("id") or "")
    current_stage = str(map_payload.get("current_stage") or "")
    completed_stages = set(map_payload.get("completed_stages", []))
    next_stages = set(map_payload.get("next_stages", []))
    if not scenario_id:
        if current_stage in {"deconvolution", "ligand_receptor", "clustering"} or {"deconvolution", "ligand_receptor"}.issubset(completed_stages):
            scenario_id = "downstream_partial" if next_stages else "downstream_complete"
        elif current_stage in {"tpm_matrix", "signature_scoring"} or "tpm_matrix" in completed_stages:
            scenario_id = "tpm_ready"
    checklist_index = _workflow_checklist_index(map_payload)
    tme_profile_input_candidates = _command_existing_input_candidates("tme_profile", map_payload)
    bayesprism_input_candidates = _command_existing_input_candidates("bayesprism", map_payload)
    bayesprism_context_notes = list(_command_profile("bayesprism").get("project_context_notes", []))
    cards: List[Dict[str, Any]] = []
    seen_ids: set[str] = set()

    def add_card(card: Dict[str, Any]) -> None:
        card_id = str(card.get("id") or "")
        if not card_id or card_id in seen_ids:
            return
        cards.append(card)
        seen_ids.add(card_id)

    if scenario_id in {"downstream_complete", "downstream_partial"}:
        add_card(
            _non_command_next_step_card(
                card_id="inspect_or_summarize_existing_results",
                title="Inspect or summarize the existing results before rerunning anything broad",
                why=(
                    "This directory already looks heavily processed, so the first high-value step is to decide whether you want interpretation, one missing downstream branch, or a true rerun."
                ),
                required_inputs=[
                    "The existing result directories already detected under the path",
                ],
                expected_result=(
                    "A clear decision about whether to add one focused downstream analysis, reuse external results, or rerun a broader workflow."
                ),
                priority=5,
                posture="recommended",
                ranking_reasons=[
                    _reason_card(
                        reason_id="processed_directory_first_decision",
                        priority=5,
                        category="ranking",
                        detail="The path already looks heavily processed, so the first recommendation should frame the decision before it jumps into commands.",
                        effect_on_recommendation="explain_decision_frame_before_listing_actions",
                    )
                ],
                tradeoffs=[
                    "This does not run a new analysis by itself; it improves the choice of the next real action.",
                ],
                existing_input_candidates=bayesprism_input_candidates[:2],
            )
        )

    priority_seed = 10
    for card in command_guidance.get("recommended", []):
        if not isinstance(card, dict):
            continue
        add_card(
            _command_next_step_card(
                card,
                priority=priority_seed,
                posture="recommended",
                ranking_reasons=[
                    _reason_card(
                        reason_id=f"rank_{card.get('command')}",
                        priority=priority_seed,
                        category="ranking",
                        detail=str(card.get("why_recommended", "")),
                        effect_on_recommendation="keep_this_above_optional_or_broader_paths",
                    )
                ],
            )
        )
        priority_seed += 5

    for card in command_guidance.get("alternatives", []):
        if not isinstance(card, dict):
            continue
        if str(card.get("category")) == "broad_rerun_bundle":
            add_card(
                _command_next_step_card(
                    card,
                    priority=priority_seed,
                    posture="conditional",
                    ranking_reasons=[
                        _reason_card(
                            reason_id=f"rank_{card.get('command')}_conditional",
                            priority=priority_seed,
                            category="ranking",
                            detail="This is still relevant, but it should be explained as a larger rerun path rather than the default next move.",
                            effect_on_recommendation="keep_below_targeted_gap_filling_steps",
                        )
                    ],
                    tradeoffs=[
                        "This usually repeats more work and compute than a focused downstream follow-up.",
                    ],
                    when_not_to_choose="Skip this when you only need to fill one missing downstream branch instead of regenerating a broad bundled workflow.",
                )
            )
        else:
            add_card(
                _command_next_step_card(
                    card,
                    priority=priority_seed,
                    posture="optional",
                    ranking_reasons=[
                        _reason_card(
                            reason_id=f"rank_{card.get('command')}_optional",
                            priority=priority_seed,
                            category="ranking",
                            detail=str(card.get("why_recommended", "")),
                            effect_on_recommendation="present_as_optional_breadth_or_specialization",
                        )
                    ],
                )
            )
        priority_seed += 5

    immune_item = checklist_index.get("immune_deconvolution", {})
    if isinstance(immune_item, dict) and str(immune_item.get("status")) == "partial":
        add_card(
            _non_command_next_step_card(
                card_id="consider_missing_deconvolution_methods",
                title="Decide whether to fill the missing default TME methods",
                why=(
                    "The directory already contains some deconvolution outputs, so the next decision is whether to keep the current usable result set or add the remaining default tme_profile/runall methods."
                ),
                required_inputs=[
                    "The detected TPM-like matrix that can feed the missing default TME methods",
                ],
                expected_result=(
                    "A deliberate choice between keeping the current deconvolution set as-is or adding concrete default methods such as EPIC, quanTIseq, MCPcounter, ESTIMATE, or IPS. BayesPrism should be presented separately as an optional standalone branch, not as part of tme_profile."
                ),
                priority=35,
                posture="optional",
                ranking_reasons=[
                    _reason_card(
                        reason_id="missing_default_deconvolution_methods",
                        priority=35,
                        category="gap",
                        detail="The current directory is usable already, so default-method backfilling should be described as a targeted completeness decision rather than a mandatory rerun.",
                        effect_on_recommendation="treat_as_targeted_default_method_backfill",
                    )
                ],
                existing_input_candidates=tme_profile_input_candidates or bayesprism_input_candidates,
                project_context_notes=[
                    "The default downstream bundle has six methods: CIBERSORT, IPS, ESTIMATE, MCPcounter, quanTIseq, and EPIC.",
                    "BayesPrism is standalone and should not be counted as missing from tme_profile or runall.",
                    *bayesprism_context_notes,
                ],
            )
        )

    if map_payload.get("external_analysis_hints"):
        add_card(
            _non_command_next_step_card(
                card_id="decide_external_reuse_vs_native_rerun",
                title="Decide whether external immune-side results should be reused or rerun natively",
                why=(
                    "The path contains external-analysis hints, so the next decision is whether those outputs are good enough to reuse or whether you want native iobrpy provenance."
                ),
                required_inputs=[
                    "The existing external result files and the user goal for reproducibility",
                ],
                expected_result=(
                    "A clear choice between reusing external outputs and rerunning the native counterpart such as trust4 or hla_typing."
                ),
                priority=40,
                posture="conditional",
                ranking_reasons=[
                    _reason_card(
                        reason_id="external_results_need_provenance_decision",
                        priority=40,
                        category="provenance",
                        detail="External results change how strongly the agent can claim provenance and reproducibility.",
                        effect_on_recommendation="make_reuse_vs_native_rerun_explicit",
                    )
                ],
            )
        )

    rerun_choices = {
        str(item.get("id")): item
        for item in scenario.get("choice_details", [])
        if isinstance(item, dict)
    }
    if "rerun_current_stage" in rerun_choices:
        add_card(
            _non_command_next_step_card(
                card_id="rerun_current_stage_only_if_needed",
                title="Rerun only the current stage if the existing outputs are wrong or incomplete",
                why=str(rerun_choices["rerun_current_stage"].get("when_to_choose", "")),
                required_inputs=[
                    "The same inputs required by the current stage",
                ],
                expected_result="A regenerated version of the current stage without repeating the entire workflow.",
                priority=80,
                posture="conditional",
                ranking_reasons=[
                    _reason_card(
                        reason_id="rerun_current_stage_is_fallback",
                        priority=80,
                        category="rerun",
                        detail="This belongs below reuse-oriented options unless the current outputs are actually wrong or incomplete.",
                        effect_on_recommendation="keep_below_reuse_and_gap_filling_options",
                    )
                ],
                tradeoffs=[
                    "This repeats work for the current stage even if most downstream outputs are otherwise reusable.",
                ],
                when_not_to_choose="Skip this when the existing outputs look reusable and the real need is only one downstream branch or interpretation.",
            )
        )
    if "rerun_full_pipeline" in rerun_choices:
        add_card(
            _non_command_next_step_card(
                card_id="rerun_full_pipeline_only_for_clean_rebuild",
                title="Use a full rerun only when you want a clean rebuild in a new output directory",
                why=str(rerun_choices["rerun_full_pipeline"].get("when_to_choose", "")),
                required_inputs=[
                    "Raw FASTQ inputs, a fresh output directory, and the correct reference index",
                ],
                expected_result="A clean end-to-end rebuilt analysis tree.",
                priority=90,
                posture="conditional",
                ranking_reasons=[
                    _reason_card(
                        reason_id="full_rerun_is_last_resort",
                        priority=90,
                        category="rerun",
                        detail="A full rerun should stay near the bottom unless the user explicitly wants a clean rebuild or the existing outputs are unusable.",
                        effect_on_recommendation="keep_as_last_resort",
                    )
                ],
                tradeoffs=[
                    "This is the most expensive option in time and compute and usually needs a fresh output directory.",
                ],
                when_not_to_choose="Skip this when most of the expensive intermediate or downstream outputs are already acceptable.",
            )
        )

    return sorted(cards, key=lambda item: (int(item.get("priority", 999)), str(item.get("id", ""))))


def _next_step_reasoning_surface(
    map_payload: Dict[str, Any],
    current_state_assessment: Dict[str, Any],
    next_step_cards: List[Dict[str, Any]],
) -> Dict[str, Any]:
    scenario = map_payload.get("scenario", {}) if isinstance(map_payload.get("scenario", {}), dict) else {}
    scenario_id = str(scenario.get("id") or current_state_assessment.get("directory_profile") or "unknown")
    if scenario_id in {"downstream_complete", "downstream_partial"}:
        primary_decision_frame = (
            "Explain the reuse-versus-add-one-more-step decision before presenting any rerun path."
        )
    elif scenario_id == "tpm_ready":
        primary_decision_frame = (
            "Explain that the matrix stage is already usable, then rank focused downstream analyses above upstream work."
        )
    elif scenario_id == "upstream_partial":
        primary_decision_frame = (
            "Explain that the directory is resumable, then rank resume-oriented steps above restart-oriented ones."
        )
    elif scenario_id == "raw_fastq_only":
        primary_decision_frame = (
            "Explain that the main decision is which workflow entrypoint to launch from raw input."
        )
    else:
        primary_decision_frame = (
            "Explain uncertainty first, then present the safest narrowed options before any broad workflow."
        )

    ranking_criteria = [
        {
            "id": "reuse_before_rerun_when_possible",
            "detail": "Prefer options that reuse good existing outputs before options that repeat broad computation.",
        },
        {
            "id": "fill_specific_gaps_before_broad_bundles",
            "detail": "When the directory is already far downstream, rank missing focused branches above bundled reruns.",
        },
        {
            "id": "make_external_provenance_explicit",
            "detail": "If external results are present, explain the provenance choice instead of silently treating them as native outputs.",
        },
        {
            "id": "name_detected_inputs_when_available",
            "detail": "When a card already includes detected input candidates, name those concrete files or directories instead of describing the input generically.",
        },
    ]
    if any(str(item.get("command") or "") in {"tme_cluster", "nmf"} for item in next_step_cards if isinstance(item, dict)):
        ranking_criteria.append(
            {
                "id": "confirm_clustering_matrix_before_execution",
                "detail": "For clustering, confirm which detected matrix to use before execution; prefer tutorial-style CIBERSORT inputs when available, and do not auto-use a very wide signature matrix.",
            }
        )

    return {
        "primary_decision_frame": primary_decision_frame,
        "ranking_criteria": ranking_criteria,
        "option_ids_in_priority_order": [str(item.get("id")) for item in next_step_cards if item.get("id")],
        "llm_recommendation_contract": {
            "compose_recommendations_in_your_own_words": True,
            "explain_ranking_before_listing_commands": True,
            "quote_card_titles_only_when_useful": False,
            "connect_each_option_to_inputs_and_outputs": True,
            "name_detected_files_or_directories_when_present": True,
            "use_project_context_notes_to_explain_builtin_defaults": True,
            "use_ranking_reasons_and_tradeoffs_instead_of_parroting_one_card_field": True,
        },
    }


def _call_analyze_path_task_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    path = arguments.get("path")
    if not isinstance(path, str) or not path:
        raise ValueError("analyze_path_task requires a non-empty string path.")
    task_text = str(arguments.get("task") or "")
    compact_for_transport = bool(arguments.get("_compact_for_transport", True))
    include_rich_guidance = bool(arguments.get("_include_rich_guidance", False))
    language = _response_language(arguments.get("language", "auto"), context=task_text)

    map_payload = _call_map_path_tool(
        {
            **arguments,
            "_compact_for_transport": False,
        }
    )
    detected_input = _detect_input_summary(Path(path))
    recommendation_payload = _call_recommend_workflow_tool(
        {
            "path": path,
            "task": task_text,
            "_detected_input": detected_input,
            "_pipeline_map_payload": map_payload,
            "_include_pipeline_map": False,
            "_compact_for_transport": False,
        }
    )

    scenario = map_payload.get("scenario", {}) if isinstance(map_payload.get("scenario", {}), dict) else {}
    recommendation = (
        recommendation_payload.get("recommendation", {})
        if isinstance(recommendation_payload.get("recommendation", {}), dict)
        else {}
    )
    decision_paths = _decision_paths_from_map_payload(map_payload)
    command_guidance = _recommended_and_alternative_commands(
        map_payload,
        recommendation,
        task_text,
    )
    directory_recognition = _directory_recognition_payload(
        map_payload,
        recommendation,
        task_text,
    )
    current_state_assessment = _current_state_assessment(
        map_payload,
        recommendation,
        command_guidance,
    )
    next_step_cards = _user_facing_next_step_cards(
        map_payload,
        recommendation,
        command_guidance,
    )
    concise_confirmation_plan = _concise_confirmation_plan(next_step_cards)
    next_step_reasoning_surface = _next_step_reasoning_surface(
        map_payload,
        current_state_assessment,
        next_step_cards,
    )
    compact_map_detail = "structured" if include_rich_guidance else "summary"
    compact_map_payload = _compact_map_payload_for_transport(
        map_payload,
        detail=compact_map_detail,
        language=language,
    )
    compact_recommendation = _compact_recommendation_for_transport(recommendation_payload)
    compact_decision_paths = _compact_decision_paths_for_transport(decision_paths)
    compact_command_guidance = {
        "recommended": _compact_command_guidance_bucket_for_transport(command_guidance["recommended"], limit=4),
        "alternatives": _compact_command_guidance_bucket_for_transport(command_guidance["alternatives"], limit=6),
        "latent": _compact_command_guidance_bucket_for_transport(command_guidance["latent"], limit=4),
    }
    compact_directory_recognition = _compact_directory_recognition_payload_for_transport(directory_recognition)
    compact_current_state_assessment = _compact_current_state_assessment_for_transport(current_state_assessment)
    compact_next_step_cards = _compact_next_step_cards_for_transport(next_step_cards, limit=6)
    compact_concise_confirmation_plan = _compact_concise_confirmation_plan_for_transport(concise_confirmation_plan)
    compact_next_step_reasoning_surface = _compact_next_step_reasoning_surface_for_transport(next_step_reasoning_surface)
    payload = {
        "success": bool(map_payload.get("success")) and bool(recommendation_payload.get("success")),
        "entrypoint": "iobrpy-cli-mcp",
        "path": path,
        "task": task_text,
        "pipeline_map": compact_map_payload,
        "recommendation": compact_recommendation if compact_for_transport else {
            "success": recommendation_payload.get("success"),
            "task": recommendation_payload.get("task"),
            "detected_input": recommendation_payload.get("detected_input"),
            "path_state": recommendation_payload.get("path_state"),
            "recommendation": recommendation_payload.get("recommendation"),
        },
        "agent_guidance": {
            "current_stage": map_payload.get("current_stage"),
            "current_label": map_payload.get("current_label"),
            "current_stage_confidence": map_payload.get("current_stage_confidence"),
            "recommended_action": map_payload.get("recommended_action"),
            "resume_recommended": map_payload.get("resume_recommended"),
            "rerun_recommended": map_payload.get("rerun_recommended"),
            "scenario_summary": scenario.get("summary", ""),
            "communication_goals": scenario.get("communication_goals", []),
            "required_user_decisions": scenario.get("required_user_decisions", []),
            "recommended_command": recommendation.get("recommended_command"),
            "missing_parameters": recommendation.get("missing_parameters", []),
            "confirm_parameters": recommendation.get("confirm_parameters", []),
            "important_optional_parameters": recommendation.get("important_optional_parameters", []),
            "decision_paths": compact_decision_paths if compact_for_transport else decision_paths,
            "recommended_commands": compact_command_guidance["recommended"] if compact_for_transport else command_guidance["recommended"],
            "alternative_commands": compact_command_guidance["alternatives"] if compact_for_transport else command_guidance["alternatives"],
            "latent_specialized_commands": compact_command_guidance["latent"] if compact_for_transport else command_guidance["latent"],
            "current_state_reasoning_surface": compact_current_state_assessment if compact_for_transport else current_state_assessment,
            "current_state_assessment": compact_current_state_assessment if compact_for_transport else current_state_assessment,
            "ranked_next_step_options": compact_next_step_cards if compact_for_transport else next_step_cards,
            "next_step_reasoning_surface": compact_next_step_reasoning_surface if compact_for_transport else next_step_reasoning_surface,
            "user_facing_next_step_cards": compact_next_step_cards if compact_for_transport else next_step_cards,
            "concise_confirmation_plan": compact_concise_confirmation_plan if compact_for_transport else concise_confirmation_plan,
            "directory_recognition": compact_directory_recognition if compact_for_transport else directory_recognition,
            "response_guidance": {
                "assume_beginner_friendly_explanation": True,
                "do_not_limit_recommendations_to_three": True,
                "for_each_presented_recommendation_explain": [
                    "why_it_matches_the_current_state",
                    "what_input_data_is_needed",
                    "what_result_it_will_produce",
                    "which_detected_files_or_directories_are_the_best_current_inputs",
                ],
                "mention_latent_specialized_commands_only_when_relevant": True,
                "prefer_current_state_assessment_over_bare_stage_label": True,
                "prefer_user_facing_next_step_cards_over_a_single_generic_next_step": True,
                "prefer_current_state_reasoning_surface_over_prewritten_assessment": True,
                "compose_current_state_summary_in_your_own_words": True,
                "prefer_ranked_next_step_options_and_next_step_reasoning_surface_over_flat_commands": True,
                "compose_next_step_explanations_in_your_own_words": True,
                "prefer_detected_input_candidates_over_generic_input_phrasing": True,
                "use_agent_reasoning_rules_to_avoid_command_misfires": True,
                "use_upstream_recovery_paths_when_a_downstream_matrix_input_fails": True,
                "use_project_context_notes_to_explain_builtin_defaults_and_iobrpy_specific_behavior": True,
                "use_execution_guardrails_before_constructing_clustering_commands": True,
                "ask_pre_execution_questions_when_a_card_requires_user_choice": True,
                "prefer_concise_confirmation_plan_over_raw_confirm_parameters": True,
                "prefer_concise_natural_language_confirmation_questions": True,
                "render_confirmation_choices_as_numbered_list_when_possible": True,
                "ask_only_the_highest_impact_one_or_two_confirmation_items_first": True,
                "avoid_dumping_raw_flag_names_before_execution_when_plain_language_will_do": True,
                "resource_sensitive_parameters_require_execution_host_probe": True,
                "do_not_blindly_use_thread_defaults_or_hardcode_8": True,
            },
        },
        "next_steps": [
            "Use pipeline_map.workflow_checklist as the source of truth for what is already complete.",
            "Use agent_guidance.current_state_reasoning_surface to compose the current-state summary in your own words instead of parroting a canned assessment.",
            "Use agent_guidance.ranked_next_step_options plus next_step_reasoning_surface to explain recommendation order in your own words instead of collapsing everything into one generic next step.",
            "Use agent_guidance.decision_paths plus recommended_commands and alternative_commands when explaining next options to the user.",
            "Explain why each presented recommendation matches the current state, what input data it needs, which detected files or directories are the best current inputs, and what result it will produce.",
            "When a next-step card includes existing_input_candidates or detected_input_story, use those concrete paths and explanations instead of vague phrases such as 'existing TME feature files'.",
            "When a card includes agent_reasoning_rules, use them to explain why that command fits or does not fit the current request before execution.",
            "When a downstream matrix command card includes upstream_recovery_paths, use them if the current TPM-like matrix is unusable or the command fails. Regenerated TPM should be standardized with log2_eset before retrying downstream functions.",
            "When a clustering card includes execution_guardrails, follow them before choosing the input matrix or constructing `--features` ranges.",
            "When a card includes pre_execution_questions, ask or answer those matrix-choice questions before executing the native clustering command.",
            "When concise_confirmation_plan or concise_confirmation_prompts are present, prefer them over raw confirm_parameters. Ask only the one or two highest-impact questions first, and phrase them as short natural-language numbered options when possible.",
            "When a command has resource-sensitive parameters such as threads, num_threads, parallel_size, num_processes, or batch_size, probe CPU cores, load, and free memory on the same execution host before choosing a value; do not blindly keep the native default or hard-code 8.",
            "Before recommending, writing, or executing a native command invocation, validate the command and flags against RAG_MCP/iobrpy_required_params.json using list_native_commands or the native command help. Do not invent flags from general workflow concepts. In particular, tme_profile has no --mode; its JSON-defined parameters are --input, --output, and optional --threads.",
            "When a card includes project_context_notes, use them to explain iobrpy-specific defaults such as BayesPrism's bundled reference behavior.",
            "When directory_recognition.layer_2_llm_reasoning_surface.should_expand_reasoning is true, use its reason_details, evidence_samples, candidate_function_hypotheses, external_result_interpretation_hints, and suggested_user_questions to explain ambiguity or specialized-function possibilities.",
            "Only mention latent_specialized_commands such as anno_eset, mouse2human_eset, or log2_eset when the user task actually needs those matrix-preparation steps.",
            "Once the path state and parameters are clear, call the chosen native command tool.",
        ],
        "transport_compaction": {
            "applied": compact_for_transport,
            "pipeline_map_compacted": True,
            "recommendation_compacted": compact_for_transport,
            "agent_guidance_compacted": compact_for_transport,
            "duplicate_map_scan_avoided": True,
            "language": language,
        },
    }
    if compact_for_transport and not include_rich_guidance:
        payload = _compact_analyze_path_task_payload_for_default_transport(payload)
    return _decorate_helper_payload("analyze_path_task", payload)


def _call_helper_tool(command_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    helper_handlers = {
        "list_native_commands": _call_list_native_commands_tool,
        "map_path": _call_map_path_tool,
        "recommend_workflow": _call_recommend_workflow_tool,
        "doctor_environment": _call_doctor_environment_tool,
        "analyze_path_task": _call_analyze_path_task_tool,
    }
    handler = helper_handlers.get(command_name)
    if handler is None:
        raise ValueError(f"Unknown helper tool: {command_name}")
    return handler(arguments)


def _split_internal_native_arguments(arguments: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    public_arguments: Dict[str, Any] = {}
    internal_arguments: Dict[str, Any] = {}
    for key, value in arguments.items():
        if isinstance(key, str) and key.startswith(_INTERNAL_NATIVE_ARGUMENT_PREFIX):
            internal_arguments[key] = value
        else:
            public_arguments[key] = value
    return public_arguments, internal_arguments


def _coerce_internal_bool(value: Any, default: bool = True) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"false", "0", "no", "off"}:
            return False
        if normalized in {"true", "1", "yes", "on"}:
            return True
    return bool(value)


def _run_native_process(argv: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "iobrpy.main", *argv],
        capture_output=True,
        text=True,
        check=False,
    )


def _native_process_result(
    command_name: str,
    argv: Sequence[str],
    process: subprocess.CompletedProcess[str],
) -> Dict[str, Any]:
    invocation_parts = [sys.executable, "-m", "iobrpy.main", *argv]
    return {
        "success": process.returncode == 0,
        "command": command_name,
        "argv": list(argv),
        "invocation": " ".join(shlex.quote(part) for part in invocation_parts),
        "returncode": process.returncode,
        "stdout": process.stdout,
        "stderr": process.stderr,
    }


def _output_path_looks_like_file(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith(
        (
            ".csv",
            ".csv.gz",
            ".tsv",
            ".tsv.gz",
            ".txt",
            ".txt.gz",
            ".json",
            ".json.gz",
            ".rds",
            ".pkl",
            ".pickle",
            ".xlsx",
            ".xls",
            ".gz",
        )
    )


def _directory_target_for_output_path(path: Path) -> Path:
    return path.parent if _output_path_looks_like_file(path) else path


def _materialize_directory_targets(paths: Sequence[Path]) -> List[str]:
    created: List[str] = []
    seen: set[str] = set()
    for path in paths:
        rendered = str(path)
        if not rendered or rendered in seen:
            continue
        path.mkdir(parents=True, exist_ok=True)
        seen.add(rendered)
        created.append(rendered)
    return created


def _materialize_output_targets_from_arguments(arguments: Dict[str, Any]) -> List[str]:
    output_dirs: List[Path] = []
    for key in ("output", "outdir"):
        value = arguments.get(key)
        if not isinstance(value, str) or not value:
            continue
        output_dirs.append(_directory_target_for_output_path(Path(value)))
    return _materialize_directory_targets(output_dirs)


def _materialize_output_targets_for_step(step: Dict[str, Any]) -> List[str]:
    paths: List[Path] = []
    suggested_arguments = step.get("suggested_arguments", {})
    if isinstance(suggested_arguments, dict):
        for key in ("output", "outdir"):
            value = suggested_arguments.get(key)
            if isinstance(value, str) and value:
                paths.append(_directory_target_for_output_path(Path(value)))
    expected_output = step.get("expected_output")
    if isinstance(expected_output, str) and expected_output:
        paths.append(_directory_target_for_output_path(Path(expected_output)))
    return _materialize_directory_targets(paths)


def _trim_output_excerpt(text: str, limit: int = 1200) -> str:
    if len(text) <= limit:
        return text
    return f"...\n{text[-limit:]}"


def _summarize_failure_for_transport(result: Dict[str, Any]) -> Dict[str, Any]:
    summary = {
        "command": result.get("command"),
        "argv": result.get("argv"),
        "invocation": result.get("invocation"),
        "returncode": result.get("returncode"),
    }
    stdout = result.get("stdout")
    stderr = result.get("stderr")
    if isinstance(stdout, str) and stdout:
        summary["stdout_excerpt"] = _trim_output_excerpt(stdout)
    if isinstance(stderr, str) and stderr:
        summary["stderr_excerpt"] = _trim_output_excerpt(stderr)
    if "failure_guidance" in result:
        summary["failure_guidance"] = result["failure_guidance"]
    if result.get("recovered_from_failure"):
        summary["recovered_from_failure"] = True
    return summary


def _auto_recovery_step_output_paths(step: Dict[str, Any]) -> List[str]:
    expected_output = step.get("expected_output")
    if isinstance(expected_output, str):
        return _dedupe_text_values([expected_output], limit=12)
    if isinstance(expected_output, list):
        return _dedupe_text_values(expected_output, limit=12)
    return []


def _generated_output_paths_from_step_records(step_records: Sequence[Dict[str, Any]]) -> List[str]:
    generated_paths: List[str] = []
    for record in step_records:
        if not record.get("success"):
            continue
        generated_paths.extend(record.get("generated_output_paths", []))
    return _dedupe_text_values(generated_paths, limit=24)


def _auto_recovery_step_record(
    index: int,
    step: Dict[str, Any],
    result: Dict[str, Any],
    created_output_targets: Sequence[str],
) -> Dict[str, Any]:
    record = {
        "index": index,
        "command": step.get("command"),
        "why": step.get("why"),
        "suggested_arguments": step.get("suggested_arguments", {}),
        "created_output_targets": list(created_output_targets),
        "success": result.get("success"),
        "returncode": result.get("returncode"),
        "invocation": result.get("invocation"),
    }
    expected_outputs = _auto_recovery_step_output_paths(step)
    if expected_outputs:
        record["expected_output"] = expected_outputs[0] if len(expected_outputs) == 1 else expected_outputs
        if result.get("success"):
            record["generated_output_paths"] = expected_outputs
    if result.get("recovered_from_failure"):
        record["recovered_from_failure"] = True
    if not result.get("success"):
        record["failure"] = _summarize_failure_for_transport(result)
    return record


def _attempt_retry_after_preparing_output_targets(
    command_name: str,
    arguments: Dict[str, Any],
) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    created_output_targets = _materialize_output_targets_from_arguments(arguments)
    if not created_output_targets:
        return None, None

    retry_result = _execute_native_command(
        command_name,
        arguments,
        allow_recovery_path=False,
        allow_output_dir_retry=False,
    )
    retry_step = _auto_recovery_step_record(
        1,
        {
            "command": command_name,
            "why": "Retry the original command after creating the missing output target directories.",
            "suggested_arguments": dict(arguments),
        },
        retry_result,
        created_output_targets,
    )
    payload = {
        "attempted": True,
        "mode": "retry_same_command_after_preparing_output_targets",
        "used_path_id": "retry_same_command_after_preparing_output_targets",
        "success": bool(retry_result.get("success")),
        "steps": [retry_step],
        "created_output_targets": created_output_targets,
    }
    if not retry_result.get("success"):
        payload["failed_step"] = retry_step
    return payload, retry_result if retry_result.get("success") else None


def _attempt_preferred_upstream_recovery_path(
    command_name: str,
    original_arguments: Dict[str, Any],
    recovery_paths: Sequence[Dict[str, Any]],
) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    if not recovery_paths:
        return None, None

    chosen_path = next((path for path in recovery_paths if path.get("preferred")), recovery_paths[0])
    step_records: List[Dict[str, Any]] = []
    final_result: Optional[Dict[str, Any]] = None

    for index, step in enumerate(chosen_path.get("steps", []), start=1):
        step_command = step.get("command")
        if not isinstance(step_command, str) or not step_command:
            continue
        step_arguments = step.get("suggested_arguments", {})
        if not isinstance(step_arguments, dict):
            step_arguments = {}
        if not step_arguments and step_command == command_name:
            step_arguments = dict(original_arguments)
        created_output_targets = _materialize_output_targets_for_step(step)
        step_result = _execute_native_command(
            step_command,
            step_arguments,
            allow_recovery_path=False,
            allow_output_dir_retry=True,
        )
        step_record = _auto_recovery_step_record(index, step, step_result, created_output_targets)
        step_records.append(step_record)
        if not step_result.get("success"):
            payload = {
                "attempted": True,
                "mode": "follow_preferred_upstream_recovery_path",
                "used_path_id": chosen_path.get("id"),
                "path_title": chosen_path.get("title"),
                "path_reason": chosen_path.get("why_this_path_exists"),
                "success": False,
                "steps": step_records,
                "failed_step": step_record,
                "generated_output_paths": _generated_output_paths_from_step_records(step_records),
            }
            return payload, None
        final_result = step_result

    generated_output_paths = _generated_output_paths_from_step_records(step_records)
    payload = {
        "attempted": True,
        "mode": "follow_preferred_upstream_recovery_path",
        "used_path_id": chosen_path.get("id"),
        "path_title": chosen_path.get("title"),
        "path_reason": chosen_path.get("why_this_path_exists"),
        "success": True,
        "steps": step_records,
    }
    if generated_output_paths:
        payload["generated_output_paths"] = generated_output_paths
    return payload, final_result


def _maybe_auto_recover_native_failure(
    command_name: str,
    arguments: Dict[str, Any],
    failure_guidance: Dict[str, Any],
    *,
    allow_recovery_path: bool,
    allow_output_dir_retry: bool,
) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    classification = failure_guidance.get("classification", {})
    classification_id = classification.get("id")

    if allow_output_dir_retry and classification_id == "missing_output_directory":
        auto_recovery, recovered_result = _attempt_retry_after_preparing_output_targets(command_name, arguments)
        if auto_recovery is not None:
            return auto_recovery, recovered_result

    recovery_paths = failure_guidance.get("upstream_recovery_paths", [])
    if (
        allow_recovery_path
        and classification.get("retryable")
        and isinstance(recovery_paths, list)
        and recovery_paths
    ):
        return _attempt_preferred_upstream_recovery_path(command_name, arguments, recovery_paths)

    return None, None


def _execute_native_command(
    command_name: str,
    arguments: Dict[str, Any],
    *,
    allow_recovery_path: bool,
    allow_output_dir_retry: bool,
) -> Dict[str, Any]:
    argv = _build_cli_argv(command_name, arguments)
    process = _run_native_process(argv)
    result = _native_process_result(command_name, argv, process)

    if process.returncode == 0:
        return result

    failure_guidance = _native_failure_guidance(
        command_name,
        arguments,
        process.stdout,
        process.stderr,
    )
    result["failure_guidance"] = failure_guidance

    auto_recovery, recovered_result = _maybe_auto_recover_native_failure(
        command_name,
        arguments,
        failure_guidance,
        allow_recovery_path=allow_recovery_path,
        allow_output_dir_retry=allow_output_dir_retry,
    )
    if auto_recovery is not None:
        result["auto_recovery"] = auto_recovery
        if recovered_result is not None and recovered_result.get("success"):
            recovered = dict(recovered_result)
            nested_auto_recovery = recovered.get("auto_recovery")
            nested_initial_failure = recovered.get("initial_failure")
            if nested_auto_recovery is not None:
                auto_recovery["final_step_auto_recovery"] = nested_auto_recovery
            if nested_initial_failure is not None:
                auto_recovery["final_step_initial_failure"] = nested_initial_failure
            recovered["recovered_from_failure"] = True
            recovered["initial_failure"] = _summarize_failure_for_transport(result)
            recovered["failure_guidance"] = failure_guidance
            recovered["auto_recovery"] = auto_recovery
            if auto_recovery.get("generated_output_paths"):
                recovered["recovery_generated_output_paths"] = auto_recovery["generated_output_paths"]
            return recovered

    return result


def call_tool(command_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    if command_name in _HELPER_TOOL_NAMES:
        return _call_helper_tool(command_name, arguments)

    public_arguments, internal_arguments = _split_internal_native_arguments(arguments)
    allow_auto_recovery = _coerce_internal_bool(
        internal_arguments.get("_allow_auto_recovery"),
        default=True,
    )
    return _execute_native_command(
        command_name,
        public_arguments,
        allow_recovery_path=allow_auto_recovery,
        allow_output_dir_retry=allow_auto_recovery,
    )


def _handle_initialize(_id: Any, params: Optional[Dict[str, Any]]) -> None:
    requested = (params or {}).get("protocolVersion")
    chosen = requested if requested in SUPPORTED_PROTOCOLS else SUPPORTED_PROTOCOLS[-1]
    _send(
        {
            "jsonrpc": "2.0",
            "id": _id,
            "result": {
                "protocolVersion": chosen,
                "capabilities": CAPABILITIES,
                "serverInfo": SERVER_INFO,
            },
        }
    )


def _handle_tools_list(_id: Any) -> None:
    _send({"jsonrpc": "2.0", "id": _id, "result": {"tools": get_tools()}})


def _handle_tools_call(_id: Any, params: Optional[Dict[str, Any]]) -> None:
    params = params or {}
    name = params.get("name")
    arguments = params.get("arguments") or {}

    if not isinstance(name, str) or not name:
        raise ValueError("Tool name is required.")
    if not isinstance(arguments, dict):
        raise ValueError("Tool arguments must be an object.")

    result = call_tool(name, arguments)
    _send(
        {
            "jsonrpc": "2.0",
            "id": _id,
            "result": {
                "content": [{"type": "text", "text": json.dumps(result, ensure_ascii=False, separators=(",", ":"))}],
                "isError": not result["success"],
            },
        }
    )


def main() -> None:
    for raw_line in sys.stdin:
        line = raw_line.strip()
        if not line:
            continue

        try:
            message = json.loads(line)
            method = message.get("method")
            _id = message.get("id")
            params = message.get("params")

            if method == "initialize":
                _handle_initialize(_id, params)
            elif method == "notifications/initialized":
                continue
            elif method == "tools/list":
                _handle_tools_list(_id)
            elif method == "tools/call":
                _handle_tools_call(_id, params)
            else:
                _jsonrpc_error(_id, -32601, f"Unknown method: {method}")
        except Exception as exc:
            _jsonrpc_error(
                message.get("id") if "message" in locals() else None,
                -32603,
                "Internal error",
                data={"exception": repr(exc), "traceback": traceback.format_exc()},
            )


if __name__ == "__main__":
    main()
    def _schema_choices(spec: ArgumentSpec) -> List[str]:
        if "enum" in spec.schema and isinstance(spec.schema["enum"], list):
            return [str(item) for item in spec.schema["enum"]]
        items = spec.schema.get("items", {})
        if isinstance(items, dict) and isinstance(items.get("enum"), list):
            return [str(item) for item in items["enum"]]
        return []

