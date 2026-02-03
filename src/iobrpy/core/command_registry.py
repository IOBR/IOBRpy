from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


@dataclass(frozen=True)
class CommandSpec:
    name: str
    required: List[str]
    optional: List[str] = field(default_factory=list)
    output_params: List[str] = field(default_factory=list)
    aliases: List[str] = field(default_factory=list)
    param_flags: Dict[str, str] = field(default_factory=dict)


COMMAND_SPECS: Dict[str, CommandSpec] = {
    "runall": CommandSpec(
        name="runall",
        required=["mode", "outdir", "fastq"],
        optional=["resume", "dry_run"],
        output_params=["outdir"],
        aliases=["run all", "pipeline", "end-to-end"],
    ),
    "tme_profile": CommandSpec(
        name="tme_profile",
        required=["input", "output"],
        optional=["threads"],
        output_params=["output"],
        aliases=["tme", "tme profile"],
    ),
    "fastq_qc": CommandSpec(
        name="fastq_qc",
        required=["path1_fastq", "path2_fastp"],
        optional=["num_threads", "suffix1", "batch_size", "se", "length_required"],
        output_params=["path2_fastp"],
        aliases=["qc", "fastp"],
    ),
    "trust4": CommandSpec(
        name="trust4",
        required=["fqdir", "o"],
        optional=["t", "f", "ref"],
        output_params=["o"],
        aliases=["trust4", "tcr", "bcr"],
        param_flags={"fqdir": "--fqdir", "o": "-o", "t": "-t"},
    ),
    "spechla": CommandSpec(
        name="spechla",
        required=["name", "read1", "read2", "outdir"],
        optional=["threads"],
        output_params=["outdir"],
        aliases=["spec hla", "hla typing"],
        param_flags={"name": "--name", "read1": "--read1", "read2": "--read2", "outdir": "--outdir", "threads": "--threads"},
    ),
    "prepare_salmon": CommandSpec(
        name="prepare_salmon",
        required=["input", "output"],
        optional=["return_feature", "remove_version"],
        output_params=["output"],
    ),
    "count2tpm": CommandSpec(
        name="count2tpm",
        required=["input", "output"],
        optional=[
            "effLength_csv",
            "idtype",
            "org",
            "source",
            "id",
            "length",
            "gene_symbol",
            "check_data",
            "remove_version",
        ],
        output_params=["output"],
    ),
    "anno_eset": CommandSpec(
        name="anno_eset",
        required=["input", "output", "annotation"],
        optional=["annotation_file", "annotation_key", "symbol", "probe", "method", "remove_version"],
        output_params=["output"],
    ),
    "calculate_sig_score": CommandSpec(
        name="calculate_sig_score",
        required=["input", "output", "signature"],
        optional=["method", "mini_gene_count", "adjust_eset", "parallel_size"],
        output_params=["output"],
    ),
    "cibersort": CommandSpec(
        name="cibersort",
        required=["input", "output"],
        optional=["perm", "QN", "absolute", "abs_method", "threads"],
        output_params=["output"],
    ),
    "IPS": CommandSpec(
        name="IPS",
        required=["input", "output"],
        output_params=["output"],
    ),
    "estimate": CommandSpec(
        name="estimate",
        required=["input", "output"],
        optional=["platform"],
        output_params=["output"],
    ),
    "mcpcounter": CommandSpec(
        name="mcpcounter",
        required=["input", "features", "output"],
        output_params=["output"],
    ),
    "quantiseq": CommandSpec(
        name="quantiseq",
        required=["input", "output"],
        optional=["arrays", "signame", "tumor", "mRNAscale", "method", "rmgenes"],
        output_params=["output"],
    ),
    "epic": CommandSpec(
        name="epic",
        required=["input", "output"],
        optional=["reference", "header", "filter"],
        output_params=["output"],
    ),
    "deside": CommandSpec(
        name="deside",
        required=["input", "output"],
        optional=["sctype", "sc_dat", "weight", "trans", "cancer_type", "python"],
        output_params=["output"],
    ),
    "tme_cluster": CommandSpec(
        name="tme_cluster",
        required=["input", "output"],
        optional=["k", "pca", "k_method", "seed", "features", "scale"],
        output_params=["output"],
    ),
    "LR_cal": CommandSpec(
        name="LR_cal",
        required=["input", "output"],
        optional=["cancer_type"],
        output_params=["output"],
    ),
    "nmf": CommandSpec(
        name="nmf",
        required=["input", "output"],
        optional=["k", "seed", "use", "method", "iter", "exclude", "feature"],
        output_params=["output"],
    ),
    "mouse2human_eset": CommandSpec(
        name="mouse2human_eset",
        required=["input", "output"],
        optional=["mode", "symbol"],
        output_params=["output"],
    ),
    "batch_salmon": CommandSpec(
        name="batch_salmon",
        required=["index", "path_fq", "path_out"],
        optional=["suffix1", "batch_size", "num_threads", "libtype", "gtf"],
        output_params=["path_out"],
    ),
    "merge_salmon": CommandSpec(
        name="merge_salmon",
        required=["path_salmon", "project"],
        optional=["num_processes"],
        output_params=["path_salmon"],
    ),
    "merge_star_count": CommandSpec(
        name="merge_star_count",
        required=["path", "project"],
        output_params=["path"],
    ),
    "batch_star_count": CommandSpec(
        name="batch_star_count",
        required=["index", "path_fq", "path_out"],
        optional=["suffix1", "batch_size", "num_threads"],
        output_params=["path_out"],
    ),
    "log2_eset": CommandSpec(
        name="log2_eset",
        required=["input", "output"],
        output_params=["output"],
    ),
    "extract_hla_read": CommandSpec(
        name="extract_hla_read",
        required=[],
    ),
    "hla_typing": CommandSpec(
        name="hla_typing",
        required=[],
    ),
    "bayesprism": CommandSpec(
        name="bayesprism",
        required=["input", "output"],
        optional=["threads", "sc_dat", "cell_state_labels", "cell_type_labels", "key"],
        output_params=["output"],
    ),
}


def all_command_names() -> Set[str]:
    return set(COMMAND_SPECS.keys())


def find_command(text: str) -> Optional[str]:
    lowered = text.lower()
    for name, spec in COMMAND_SPECS.items():
        if name.lower() in lowered:
            return name
        for alias in spec.aliases:
            if alias.lower() in lowered:
                return name
    return None
