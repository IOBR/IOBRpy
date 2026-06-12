#!/usr/bin/env python3
"""
iobrpy_cli.py - Main CLI for IOBRpy CLI harness.

This module provides a Click-based command-line interface with REPL support
and JSON output mode for agent consumption.
"""

import json
import io
import os
import platform
import shutil
import subprocess
import sys
import textwrap
from contextlib import redirect_stdout
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Optional, List, Dict, Any, Sequence, Tuple

import click

from iobrpy_cli.iobrpy.core import (
    Project,
    ProjectManager,
    SessionManager,
    Exporter,
    TMEAnalyzer,
    DeconvolutionMethod,
    ScoringMethod,
    QuantificationWorkflow,
    QuantificationMode,
    WorkflowExecutor,
    WorkflowResult,
)

from iobrpy_cli.iobrpy.utils import (
    to_json,
    success_response,
    error_response,
    format_table,
    format_file_info,
    format_directory_info,
)
from iobrpy_cli.iobrpy.mcp_server import (
    get_native_tools,
    get_command_metadata,
    render_command_help,
    validate_native_cli_argv,
)
from iobrpy_cli.iobrpy.agent_install import (
    SUPPORTED_AGENT_CLIENTS,
    agent_status_bundle,
    install_agent_bundle,
)
from iobrpy_cli.iobrpy import __version__
from iobrpy_cli.iobrpy.pipeline_map import (
    analyze_pipeline_path,
    format_pipeline_map_report,
    map_pipeline_path,
)


# Global JSON output mode flag
_json_output = False
_ORIGINAL_CLICK_ECHO = click.echo


def _routed_click_echo(message=None, file=None, nl=True, err=False, color=None):
    """Send non-JSON progress output to stderr when JSON mode is enabled."""
    if _json_output and file is None and not err:
        err = True
    return _ORIGINAL_CLICK_ECHO(message, file=file, nl=nl, err=err, color=color)


click.echo = _routed_click_echo


def _run_with_stdout_routed(func, *args, **kwargs):
    """Route plain stdout from native helpers to stderr when JSON mode is enabled."""
    if not _json_output:
        return func(*args, **kwargs)

    stdout_buffer = io.StringIO()
    with redirect_stdout(stdout_buffer):
        result = func(*args, **kwargs)

    buffered_output = stdout_buffer.getvalue()
    if buffered_output:
        print(buffered_output, file=sys.stderr, end="")

    return result


def print_json(data: dict) -> None:
    """Print data as JSON."""
    print(to_json(data, pretty=_json_output))


def print_result(data: dict, message: str = "") -> None:
    """Print result in appropriate format."""
    success = data.get('success') if isinstance(data, dict) else None

    if _json_output:
        if isinstance(data, dict) and 'success' in data:
            output = dict(data)
            if message and not (success is False and message == "Failed" and output.get('message')):
                output['message'] = message
        else:
            output = success_response(data, message=message) if message else {'success': True, 'data': data}
        print_json(output)
    else:
        if message:
            click.echo(click.style(message, fg='green'))
        if isinstance(data, dict) and 'data' in data:
            data = data['data']
        if data:
            print(json.dumps(data, indent=2, default=str))

    if success is False:
        raise click.exceptions.Exit(1)


def _print_structured(data: Dict[str, Any]) -> None:
    """Always print machine-readable JSON for agent-facing helper commands."""
    print_json(data)


def _format_agent_record_line(record: Dict[str, Any]) -> str:
    client = record.get("client", "unknown")
    component = record.get("component", "component")
    status = str(record.get("status", "unknown")).replace("_", " ")
    message = record.get("message", "")
    line = f"- {client} {component}: {status}"
    if message:
        line += f" - {message}"
    return line


def _print_agent_install_result(payload: Dict[str, Any]) -> None:
    if _json_output:
        _print_structured(payload)
        return

    results = payload.get("results", [])
    if results and all(record.get("status") == "dry_run" for record in results):
        click.echo("Agent setup preview")
    else:
        click.echo("Agent setup: success" if payload.get("success") else "Agent setup: failed")

    for record in results:
        click.echo(_format_agent_record_line(record))

    next_steps = payload.get("next_steps", [])
    if next_steps:
        click.echo("Next steps:")
        for step in next_steps:
            click.echo(f"- {step}")


def _print_agent_status_result(payload: Dict[str, Any]) -> None:
    if _json_output:
        _print_structured(payload)
        return

    click.echo("Agent status: healthy" if payload.get("healthy") else "Agent status: needs attention")
    for record in payload.get("results", []):
        click.echo(_format_agent_record_line(record))

    next_steps = payload.get("next_steps", [])
    if next_steps:
        click.echo("Suggested fixes:")
        for step in next_steps:
            click.echo(f"- {step}")


def _print_pipeline_map_result(payload: Dict[str, Any]) -> None:
    if _json_output:
        _print_structured(payload)
        return

    click.echo(format_pipeline_map_report(payload))


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


def _native_command_summary() -> List[Dict[str, Any]]:
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


def _safe_command_metadata(command_name: str) -> Dict[str, Any]:
    try:
        return get_command_metadata(command_name)
    except Exception:
        return {"important_optional_parameters": [], "confirm_parameters": [], "required_parameter_names": []}


def _recommend_for_summary(summary: Dict[str, Any], task: str) -> Dict[str, Any]:
    task_text = (task or "").lower()
    kind = summary.get("kind")

    if kind == "fastq_directory":
        meta = _safe_command_metadata("runall")
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
            "suggested_cli": "iobrpy runall --fastq <path> --outdir <outdir> --mode <salmon|star> --index <index> --threads <int> --batch_size <int> --project <name>",
            "suggested_harness_cli": "iobrpy-cli runall --fastq <path> --outdir <outdir> --mode <salmon|star> --index <index> --threads <int> --batch_size <int> --project <name>",
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
        meta = _safe_command_metadata("tme_profile")
        return {
            "recommended_command": "tme_profile",
            "why": (
                "Detected an expression matrix file. According to RAG_MCP/iobrpy_required_params.json, "
                "`tme_profile` is the matching workflow for an existing expression matrix."
            ),
            "suggested_cli": "iobrpy tme_profile --input <matrix> --output <outdir> --threads <int>",
            "suggested_harness_cli": "iobrpy-cli tme_profile --input <matrix> --output <outdir> --threads <int>",
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
        meta = _safe_command_metadata("hla_typing")
        return {
            "recommended_command": "hla_typing",
            "why": (
                "Detected BAM/CRAM files. According to RAG_MCP/iobrpy_required_params.json, "
                "`hla_typing` is the batch BAM-directory workflow for this input."
            ),
            "suggested_cli": "iobrpy hla_typing --bam-dir <path> --ref <reference> --outdir <outdir>",
            "suggested_harness_cli": "iobrpy-cli hla_typing --bam-dir <path> --ref <reference> --outdir <outdir>",
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
            "Inspect the path with iobrpy-cli info/ls or provide a task description.",
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


_LEGACY_NAMESPACE_ALIASES: Dict[str, Dict[str, str]] = {
    "analyze": {
        "tme-profile": "tme_profile",
        "signature-score": "calculate_sig_score",
        "estimate": "estimate",
        "ips": "IPS",
        "mcpcounter": "mcpcounter",
        "cibersort": "cibersort",
        "quantiseq": "quantiseq",
        "epic": "epic",
        "bayesprism": "bayesprism",
        "lr-analysis": "LR_cal",
    },
    "quantify": {
        "qc": "fastq_qc",
        "runall": "runall",
        "prepare-salmon": "prepare_salmon",
        "count2tpm": "count2tpm",
        "batch-salmon": "batch_salmon",
        "merge-salmon": "merge_salmon",
        "batch-star-count": "batch_star_count",
        "merge-star-count": "merge_star_count",
    },
    "workflow": {
        "prepare-salmon": "prepare_salmon",
        "count2tpm": "count2tpm",
        "anno-eset": "anno_eset",
        "calculate-sig-score": "calculate_sig_score",
        "batch-salmon": "batch_salmon",
        "merge-salmon": "merge_salmon",
        "batch-star-count": "batch_star_count",
        "merge-star-count": "merge_star_count",
        "fastq-qc": "fastq_qc",
        "log2-eset": "log2_eset",
        "tme-cluster": "tme_cluster",
        "nmf": "nmf",
        "mouse2human-eset": "mouse2human_eset",
        "lr-cal": "LR_cal",
    },
    "immune": {
        "trust4": "trust4",
        "spechla": "spechla",
        "extract-hla-read": "extract_hla_read",
        "hla": "hla_typing",
        "hla-typing": "hla_typing",
    },
    "hla-tcr": {
        "trust4": "trust4",
        "spechla": "spechla",
        "extract-hla-read": "extract_hla_read",
        "hla": "hla_typing",
        "hla-typing": "hla_typing",
    },
}
_LEGACY_NAMESPACE_HELPERS: Dict[str, List[str]] = {
    "analyze": ["combine", "deconvolute"],
    "quantify": ["scan-fastq", "clean", "validate"],
}
_ROOT_HELPERS = {
    "agent",
    "commands",
    "doctor",
    "map",
    "recommend",
    "repl",
    "info",
    "ls",
    "project",
    "export",
    "version",
}
_NATIVE_COMMAND_ALIAS_CACHE: Optional[Dict[str, str]] = None
_HELPER_DESCRIPTIONS: Dict[str, str] = {
    "agent": "Install persistent Codex and Claude Code iobrpy entrypoints, plugin assets, and MCP integrations.",
    "commands": "List native iobrpy commands with semantic summaries.",
    "doctor": "Show environment diagnostics and preferred entrypoints.",
    "map": "Detect the current IOBRpy pipeline stage for a path and suggest whether to continue or rerun.",
    "recommend": "Classify an input path, suggest the right workflow, and surface confirm_parameters from the required-params JSON.",
    "project": "Manage harness-local project metadata only; not a biological analysis command.",
    "info": "Show file or directory information.",
    "ls": "List directory contents in structured form.",
    "export": "Export or summarize existing result files.",
    "version": "Show installed iobrpy and external tool versions.",
    "repl": "Start the harness REPL.",
}


def _normalize_cli_token(token: str) -> str:
    return token.strip().lower().replace("_", "-")


def _native_command_aliases() -> Dict[str, str]:
    global _NATIVE_COMMAND_ALIAS_CACHE
    if _NATIVE_COMMAND_ALIAS_CACHE is not None:
        return _NATIVE_COMMAND_ALIAS_CACHE

    aliases: Dict[str, str] = {}
    for tool in get_native_tools():
        name = tool["name"]
        for alias in {
            name,
            name.lower(),
            name.replace("_", "-"),
            name.lower().replace("_", "-"),
            _normalize_cli_token(name),
        }:
            aliases[alias] = name

    aliases["ips"] = "IPS"
    aliases["lr-cal"] = "LR_cal"
    aliases["extract-hla-read"] = "extract_hla_read"
    aliases["tme-profile"] = "tme_profile"
    _NATIVE_COMMAND_ALIAS_CACHE = aliases
    return aliases


def _resolve_native_passthrough(args: Sequence[str]) -> Optional[List[str]]:
    if not args:
        return None

    native_aliases = _native_command_aliases()
    head = args[0]
    direct = native_aliases.get(head) or native_aliases.get(_normalize_cli_token(head))
    if direct:
        return [direct, *args[1:]]

    namespace = _normalize_cli_token(head)
    namespace_aliases = _LEGACY_NAMESPACE_ALIASES.get(namespace)
    if namespace_aliases and len(args) >= 2:
        subcommand = _normalize_cli_token(args[1])
        mapped = namespace_aliases.get(subcommand)
        if mapped:
            return [mapped, *args[2:]]

    return None


_LEGACY_CLICK_SUBCOMMAND_ALIASES: Dict[Tuple[str, str], str] = {
    ("workflow", "nmf"): "nmf-cluster",
    ("immune", "hla"): "hla-typing",
    ("hla-tcr", "hla-typing"): "hla",
}
_LEGACY_NATIVE_STYLE_OPTIONS: Dict[Tuple[str, str], set[str]] = {
    ("analyze", "tme-profile"): {"--input"},
    ("analyze", "signature-score"): {"--input"},
    ("analyze", "lr-analysis"): {"--input"},
    ("quantify", "runall"): {"--fastq", "--outdir", "--batch_size", "--project"},
}


def _looks_like_native_style_legacy_invocation(namespace: str, subcommand: str, args: Sequence[str]) -> bool:
    native_style_options = _LEGACY_NATIVE_STYLE_OPTIONS.get((namespace, subcommand), set())
    if any(token in native_style_options for token in args):
        return True
    return any(token.startswith("--") and "_" in token for token in args)


def _resolve_legacy_click_passthrough(args: Sequence[str]) -> Optional[List[str]]:
    """Route compatibility namespaces through the legacy Click wrappers.

    The preferred modern path is still the top-level native command
    (`iobrpy-cli runall ...`, `iobrpy-cli tme_profile ...`).  This helper only
    preserves the older namespace surface (`workflow`, `analyze`, etc.) where
    users and tests still expect hyphenated options and JSON wrapping.
    """
    if not args:
        return None

    output_args = list(args)
    json_prefix: List[str] = []
    namespace_index = 0
    if output_args[0] == "--json":
        if len(output_args) == 1:
            return None
        json_prefix = ["--json"]
        namespace_index = 1

    namespace = _normalize_cli_token(output_args[namespace_index])
    if namespace not in _LEGACY_NAMESPACE_ALIASES:
        return None

    output_args[namespace_index] = namespace
    group = cli.commands.get(namespace) if "cli" in globals() else None
    if not isinstance(group, click.Group):
        return None

    command_names = set(group.commands)
    subcommand_index = namespace_index + 1
    if len(output_args) <= subcommand_index or output_args[subcommand_index] in {"-h", "--help"}:
        return [*json_prefix, *output_args[namespace_index:]]

    subcommand = _normalize_cli_token(output_args[subcommand_index])
    subcommand = _LEGACY_CLICK_SUBCOMMAND_ALIASES.get((namespace, subcommand), subcommand)
    if subcommand not in command_names:
        return None
    if _looks_like_native_style_legacy_invocation(namespace, subcommand, output_args[subcommand_index + 1:]):
        return None

    output_args[subcommand_index] = subcommand
    return [*json_prefix, *output_args[namespace_index:]]


def _helper_command_names() -> List[str]:
    return sorted(_ROOT_HELPERS)


def _path_shorthand_click_args(args: Sequence[str]) -> Optional[List[str]]:
    non_json_args = [arg for arg in args if arg != "--json"]
    if len(non_json_args) != 1:
        return None

    candidate = Path(non_json_args[0])
    if not candidate.exists():
        return None

    click_args = ["map", "--path", non_json_args[0]]
    if "--json" in args:
        click_args = ["--json", *click_args]
    return click_args


def _native_help_lines() -> List[str]:
    lines: List[str] = []
    for summary in _native_command_summary():
        description = summary.get("summary") or summary.get("description") or ""
        description = description.splitlines()[0].strip()
        description = textwrap.shorten(description, width=78, placeholder="...")
        lines.append(f"  {summary['name']:<20} {description}")
    return lines


def _print_root_help() -> None:
    click.echo("Usage: iobrpy-cli [OPTIONS] COMMAND [ARGS]...")
    click.echo()
    click.echo(
        "Native-first IOBRpy harness. Agent-facing command schemas and help for all "
        "non-deside/ai commands are constrained by `src/iobrpy/RAG_MCP/iobrpy_required_params.json`."
    )
    click.echo()
    click.echo("Quick start:")
    click.echo("  iobrpy-cli doctor --json")
    click.echo("  iobrpy-cli map --path <path> --json              # stage map: continue vs rerun")
    click.echo("  iobrpy-cli recommend --path <path> --task '<task>' --json   # review confirm_parameters")
    click.echo("  iobrpy-cli runall --fastq <dir> --outdir <out> --mode salmon --index <idx> --threads <n> --batch_size <n> --project <name>")
    click.echo("  iobrpy-cli /path/to/data    # shorthand for map --path")
    click.echo()
    click.echo("Helper commands:")
    for name in _helper_command_names():
        click.echo(f"  {name:<10} {_HELPER_DESCRIPTIONS.get(name, '')}")
    click.echo()
    click.echo("Compatibility namespaces:")
    click.echo("  analyze, quantify, workflow, immune, hla_tcr")
    click.echo("  These namespaces forward raw args to the matching native command.")
    click.echo()
    click.echo("Native commands:")
    for line in _native_help_lines():
        click.echo(line)


def _print_namespace_help(namespace: str) -> None:
    mappings = _LEGACY_NAMESPACE_ALIASES[namespace]
    click.echo(f"Usage: iobrpy-cli {namespace} <subcommand> [ARGS]...")
    click.echo()
    click.echo(
        "Compatibility namespace. For native iobrpy commands, raw arguments are "
        "forwarded unchanged to the matching native subcommand."
    )
    click.echo()
    click.echo("Native mappings:")
    for alias, native in sorted(mappings.items()):
        click.echo(f"  {alias:<18} -> {native}")
    helper_subcommands = _LEGACY_NAMESPACE_HELPERS.get(namespace, [])
    if helper_subcommands:
        click.echo()
        click.echo("Helper-only subcommands:")
        for name in helper_subcommands:
            click.echo(f"  {name}")
    click.echo()
    sample_native = next(iter(sorted(set(mappings.values()))), None)
    if sample_native:
        click.echo("For real parameter help, run the native help directly. Example:")
        click.echo(f"  iobrpy-cli {sample_native} --help")


def _run_native_cli(args: Sequence[str]) -> int:
    if args:
        try:
            validate_native_cli_argv(str(args[0]), list(args[1:]))
        except ValueError as exc:
            click.echo(f"Error: {exc}", err=True)
            return 2
    result = subprocess.run([sys.executable, "-m", "iobrpy.main", *args])
    return int(result.returncode)


def _print_native_passthrough_help(args: Sequence[str]) -> int:
    passthrough = _resolve_native_passthrough(args)
    if passthrough is None:
        return 1
    command_name = passthrough[0]
    click.echo(render_command_help(command_name))
    return 0


def _run_click_cli(args: Sequence[str]) -> int:
    try:
        cli.main(args=list(args), prog_name="iobrpy-cli", standalone_mode=True)
        return 0
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        return 1


def _dispatch_cli(argv: Sequence[str]) -> int:
    args = list(argv)
    if not args:
        _print_root_help()
        return 0

    if args[0] in {"-h", "--help"}:
        _print_root_help()
        return 0

    path_shorthand = _path_shorthand_click_args(args)
    if path_shorthand is not None:
        return _run_click_cli(path_shorthand)

    legacy_click_args = _resolve_legacy_click_passthrough(args)
    if legacy_click_args is not None:
        return _run_click_cli(legacy_click_args)

    if any(token in {"-h", "--help"} for token in args[1:]):
        passthrough = _resolve_native_passthrough([token for token in args if token not in {"-h", "--help"}])
        if passthrough is not None:
            return _print_native_passthrough_help(args)

    if "--json" in args and args[0] != "--json":
        json_index = args.index("--json")
        click_args = args[:json_index] + args[json_index + 1:]
        passthrough = _resolve_native_passthrough(click_args)
        if passthrough is not None:
            click.echo(
                "Error: --json is only supported for helper commands such as doctor, "
                "map, recommend, and commands.",
                err=True,
            )
            return 2
        return _run_click_cli(["--json", *click_args])

    if args[0] == "--json":
        passthrough = _resolve_native_passthrough(args[1:])
        if passthrough is not None:
            click.echo(
                "Error: --json is only supported for helper commands such as doctor, "
                "map, recommend, and commands.",
                err=True,
            )
            return 2
        return _run_click_cli(args)

    namespace = _normalize_cli_token(args[0])
    if namespace in _LEGACY_NAMESPACE_ALIASES and len(args) == 1:
        _print_namespace_help(namespace)
        return 0

    if namespace in _LEGACY_NAMESPACE_ALIASES and args[1] in {"-h", "--help"}:
        _print_namespace_help(namespace)
        return 0

    passthrough = _resolve_native_passthrough(args)
    if passthrough is not None:
        return _run_native_cli(passthrough)

    return _run_click_cli(args)


# ============================================================================
# CLI Group
# ============================================================================

@click.group()
@click.version_option(version=__version__, prog_name="iobrpy-cli")
@click.option('--json', is_flag=True, help='Output in JSON format for agent consumption')
def cli(json: bool):
    """IOBRpy CLI harness - Stateful CLI for TME analysis."""
    global _json_output
    _json_output = json


# ============================================================================
# Agent Setup Commands
# ============================================================================


@cli.group()
def agent():
    """Install persistent agent integrations such as Codex skills and MCP configs."""
    pass


@agent.command("install")
@click.option(
    "--client",
    "clients",
    multiple=True,
    type=click.Choice([*SUPPORTED_AGENT_CLIENTS, "all"]),
    default=("codex",),
    show_default=True,
    help="Agent client(s) to configure.",
)
@click.option("--skill/--no-skill", default=True, help="Install packaged `/iobrpy` and `/iobrpy-result` guidance, plugins, and entrypoints when supported.")
@click.option("--mcp/--no-mcp", default=True, help="Configure MCP integration for each selected client.")
@click.option("--force", is_flag=True, help="Overwrite the packaged Codex skill when it already exists.")
@click.option("--dry-run", is_flag=True, help="Show what would change without writing files or invoking client CLIs.")
@click.option("--server-name", default="iobrpy", show_default=True, help="MCP server name to register.")
@click.option("--codex-home", type=click.Path(path_type=Path, file_okay=False), help="Override the Codex home directory.")
@click.option("--claude-home", type=click.Path(path_type=Path, file_okay=False), help="Override the Claude Code home directory.")
@click.option("--claude-command", default="claude", show_default=True, help="Claude Code CLI command name.")
def agent_install(
    clients: Sequence[str],
    skill: bool,
    mcp: bool,
    force: bool,
    dry_run: bool,
    server_name: str,
    codex_home: Optional[Path],
    claude_home: Optional[Path],
    claude_command: str,
):
    """Install Codex and Claude Code iobrpy entrypoints, result assets, and MCP registrations."""
    payload = install_agent_bundle(
        list(clients),
        include_skill=skill,
        include_mcp=mcp,
        force=force,
        dry_run=dry_run,
        server_name=server_name,
        codex_home=codex_home,
        claude_home=claude_home,
        claude_command=claude_command,
    )
    _print_agent_install_result(payload)
    if not payload["success"]:
        raise click.exceptions.Exit(1)


@agent.command("status")
@click.option(
    "--client",
    "clients",
    multiple=True,
    type=click.Choice([*SUPPORTED_AGENT_CLIENTS, "all"]),
    default=("all",),
    show_default=True,
    help="Agent client(s) to inspect.",
)
@click.option("--skill/--no-skill", default=True, help="Inspect packaged `/iobrpy` and `/iobrpy-result` guidance, plugins, and entrypoint status.")
@click.option("--mcp/--no-mcp", default=True, help="Inspect MCP integration status for each selected client.")
@click.option("--server-name", default="iobrpy", show_default=True, help="MCP server name to inspect.")
@click.option("--codex-home", type=click.Path(path_type=Path, file_okay=False), help="Override the Codex home directory.")
@click.option("--claude-home", type=click.Path(path_type=Path, file_okay=False), help="Override the Claude Code home directory.")
@click.option("--claude-command", default="claude", show_default=True, help="Claude Code CLI command name.")
def agent_status(
    clients: Sequence[str],
    skill: bool,
    mcp: bool,
    server_name: str,
    codex_home: Optional[Path],
    claude_home: Optional[Path],
    claude_command: str,
):
    """Show IOBRpy workflow, result, plugin, and MCP status without modifying user configuration."""
    payload = agent_status_bundle(
        list(clients),
        include_skill=skill,
        include_mcp=mcp,
        server_name=server_name,
        codex_home=codex_home,
        claude_home=claude_home,
        claude_command=claude_command,
    )
    _print_agent_status_result(payload)


# ============================================================================
# Project Commands
# ============================================================================

@cli.group()
def project():
    """Manage harness-local project metadata; this does not analyze a dataset path by itself."""
    pass


@project.command()
@click.argument('name')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
@click.option('--description', '-d', help='Project description')
@click.option('--input-type', type=click.Choice(['fastq', 'tpm', 'counts']), help='Project input type')
@click.option('--fastq-dir', type=click.Path(), help='FASTQ input directory')
@click.option('--tpm-matrix', type=click.Path(), help='TPM matrix file')
@click.option('--count-matrix', type=click.Path(), help='Count matrix file')
@click.option('--organism', default='hsa', help='Organism (hsa or mmus)')
@click.option('--mode', type=click.Choice(['salmon', 'star']), help='Quantification mode')
@click.option('--threads', default=8, help='Thread count')
def create(name: str, root: str, description: str, input_type: str, fastq_dir: str,
           tpm_matrix: str, count_matrix: str, organism: str, mode: str, threads: int):
    """Create a new IOBRpy project."""
    manager = ProjectManager(Path(root))

    try:
        kwargs = {
            'description': description,
            'organism': organism,
            'threads': threads,
        }
        if input_type:
            kwargs['input_type'] = input_type
        if fastq_dir:
            kwargs['fastq_dir'] = Path(fastq_dir)
            kwargs['input_type'] = 'fastq'
        if tpm_matrix:
            kwargs['tpm_matrix'] = Path(tpm_matrix)
            kwargs['input_type'] = 'tpm'
        if count_matrix:
            kwargs['count_matrix'] = Path(count_matrix)
            kwargs['input_type'] = 'counts'
        if mode:
            kwargs['mode'] = mode

        proj = manager.create_project(name, **kwargs)
        print_result(proj.get_status(), f"Project '{name}' created at {proj.project_dir}")
    except FileExistsError:
        print_json(error_response(f"Project '{name}' already exists"))
        sys.exit(1)
    except Exception as e:
        print_json(error_response(f"Failed to create project: {str(e)}"))
        sys.exit(1)


@project.command('list')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
def list_projects(root: str):
    """List all IOBRpy projects."""
    manager = ProjectManager(Path(root))
    projects = manager.list_projects()

    if _json_output:
        print_json(success_response(projects, message=f"Found {len(projects)} projects"))
    else:
        if not projects:
            click.echo("No projects found.")
            return

        click.echo(click.style(f"Found {len(projects)} project(s):", fg='blue'))
        for p in projects:
            status_icon = "READY" if p.get('input_files', 0) > 0 else "EMPTY"
            click.echo(f"  {status_icon} {p['name']:20} | {p['input_type'] or 'No input type':10} | {p['organism']:4} | {p['input_files']} files")


@project.command()
@click.argument('name')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
def status(name: str, root: str):
    """Show project status."""
    manager = ProjectManager(Path(root))

    try:
        proj = manager.get_project(name)
        status = proj.get_status()
        print_result(status, f"Project '{name}' status")
    except FileNotFoundError:
        print_json(error_response(f"Project '{name}' not found"))
        sys.exit(1)


@project.command()
@click.argument('name')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
def show(name: str, root: str):
    """Show project details."""
    manager = ProjectManager(Path(root))

    try:
        proj = manager.get_project(name)
        details = {
            'name': proj.name,
            'path': str(proj.project_dir),
            'config': proj.config.to_dict(),
            'input_dir': str(proj.get_input_dir()),
            'output_dir': str(proj.get_output_dir()),
            'logs_dir': str(proj.get_logs_dir()),
        }
        if _json_output:
            print_json(success_response(details))
        else:
            click.echo(click.style(f"Project: {proj.name}", fg='blue', bold=True))
            click.echo(f"Path: {proj.project_dir}")
            click.echo(f"Input: {proj.get_input_dir()}")
            click.echo(f"Output: {proj.get_output_dir()}")
            click.echo(f"Logs: {proj.get_logs_dir()}")
            click.echo(f"Input type: {proj.config.input_type or 'Not set'}")
            click.echo(f"Organism: {proj.config.organism}")
            click.echo(f"Mode: {proj.config.mode or 'Not set'}")
    except FileNotFoundError:
        print_json(error_response(f"Project '{name}' not found"))
        sys.exit(1)


@project.command()
@click.argument('name')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
@click.option('--force', '-f', is_flag=True, help='Force deletion without confirmation')
def delete(name: str, root: str, force: bool):
    """Delete a project."""
    manager = ProjectManager(Path(root))

    try:
        if not force:
            if not click.confirm(f"Are you sure you want to delete project '{name}'?"):
                click.echo("Aborted.")
                return

        if manager.delete_project(name):
            click.echo(f"Project '{name}' deleted.")
        else:
            click.echo(f"Failed to delete project '{name}'.")
    except FileNotFoundError:
        print_json(error_response(f"Project '{name}' not found"))
        sys.exit(1)


# ============================================================================
# Analysis Commands
# ============================================================================

@cli.group()
def analyze():
    """Run TME analysis workflows."""
    pass


@analyze.command('tme-profile')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--signatures', default='all', help='Signature groups')
@click.option('--methods', help='Deconvolution methods (comma-separated)')
@click.option('--no-lr', is_flag=True, help='Skip ligand-receptor analysis')
@click.option('--project', '-p', help='Project name to save results')
@click.option('--root', '-r', default='.', type=click.Path(), help='Project root directory')
def tme_profile(input_matrix: str, output: str, threads: int,
                signatures: str, methods: str, no_lr: bool,
                project: str, root: str):
    """Run complete TME profile analysis."""
    analyzer = TMEAnalyzer()

    # Handle output directory
    if project:
        manager = ProjectManager(Path(root))
        try:
            proj = manager.get_project(project)
            output = str(proj.get_output_dir() / 'tme_profile')
        except FileNotFoundError:
            print_json(error_response(f"Project '{project}' not found"))
            sys.exit(1)
    elif not output:
        output = Path(input_matrix).parent / 'tme_profile_output'

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    # Parse methods
    deconv_methods = None
    if methods:
        try:
            deconv_methods = [
                TMEAnalyzer.parse_deconvolution_method(method_name)
                for method_name in methods.split(',')
                if method_name.strip()
            ]
        except ValueError as exc:
            print_json(error_response(str(exc)))
            sys.exit(1)

    # Run analysis
    click.echo(f"Running TME profile analysis...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")
    click.echo(f"  Signatures: {signatures}")
    click.echo()

    results = analyzer.run_tme_profile(
        input_matrix=Path(input_matrix),
        output_dir=output,
        threads=threads,
        signature_groups=signatures,
        deconvolution_methods=deconv_methods,
        include_lr=not no_lr,
    )

    # Format results
    result_data = {
        'output_dir': str(output),
        'signatures': {
            'success': results['signatures'].success,
            'output': str(results['signatures'].output_file) if results['signatures'].output_file else None,
        },
        'deconvolution': {
            method: {
                'success': r.success,
                'output': str(r.output_file) if r.output_file else None,
                'message': r.message,
            }
            for method, r in results['deconvolution'].items()
        },
        'lr': {
            'success': results['lr'].success if results['lr'] else False,
            'output': str(results['lr'].output_file) if results['lr'] and results['lr'].output_file else None,
        } if results['lr'] else None,
    }

    print_result(result_data, "TME profile analysis completed")


@analyze.command('signature-score')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--signatures', default='all', help='Signature groups')
@click.option('--method', type=click.Choice(['pca', 'zscore', 'ssgsea', 'integration']),
              default='integration', help='Scoring method')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--mini-gene-count', default=3, type=int, help='Minimum genes per signature')
@click.option('--adjust-eset', is_flag=True, help='Apply additional filtering before scoring')
def signature_score(input_matrix: str, output: str, signatures: str,
                   method: str, threads: int, mini_gene_count: int, adjust_eset: bool):
    """Calculate signature scores."""
    analyzer = TMEAnalyzer()

    if not output:
        output = Path(input_matrix).parent / 'signature_scores.csv'

    click.echo(f"Calculating signature scores...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Method: {method}")
    click.echo(f"  Signatures: {signatures}")

    # Prepare method-specific parameters
    method_params = {
        'input_matrix': Path(input_matrix),
        'output_file': Path(output),
        'signature_groups': signatures,
        'method': ScoringMethod(method),
        'threads': threads,
        'mini_gene_count': mini_gene_count,
        'adjust_eset': adjust_eset,
    }

    result = analyzer.calculate_signature_scores(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'signature_count': result.signature_count,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('estimate')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file path')
@click.option('--platform', type=click.Choice(['affymetrix', 'agilent', 'illumina']),
              default='affymetrix', help='Platform type for input data')
def estimate_command(input_matrix: str, output: str, platform: str):
    """Calculate ESTIMATE score."""
    from iobrpy.workflow.estimate import main as estimate_main

    # Handle default output
    if not output:
        output = Path(input_matrix).parent / 'estimate_scores.csv'
    else:
        output = Path(output)

    click.echo(f"Calculating ESTIMATE scores...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Platform: {platform}")
    click.echo(f"  Output: {output}")

    # Prepare sys.argv to simulate command line arguments
    original_argv = sys.argv.copy()
    sys.argv = ['iobrpy', '--input', input_matrix, '--output', str(output), '--platform', platform]

    try:
        _run_with_stdout_routed(estimate_main)
        result_data = {
            'success': True,
            'input': input_matrix,
            'output': str(output),
            'platform': platform,
            'method': 'ESTIMATE'
        }
        print_result(result_data, "ESTIMATE calculation completed")
    except Exception as e:
        result_data = {
            'success': False,
            'input': input_matrix,
            'output': str(output),
            'platform': platform,
            'method': 'ESTIMATE',
            'error': str(e)
        }
        print_result(result_data, f"ESTIMATE calculation failed: {str(e)}")
    finally:
        sys.argv = original_argv


@analyze.command('ips')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file path')
def ips_command(input_matrix: str, output: str):
    """Calculate Immunophenoscore (IPS)."""
    import os
    import sys

    from iobrpy.workflow.IPS import main as ips_main

    # Handle default output
    if not output:
        output = Path(input_matrix).parent / 'IPS_results.txt'
    else:
        output = Path(output)

    click.echo(f"Calculating Immunophenoscore...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")

    # Prepare sys.argv to simulate command line arguments
    original_argv = sys.argv.copy()
    sys.argv = ['iobrpy', '--input', input_matrix, '--output', str(output)]

    try:
        _run_with_stdout_routed(ips_main)
        result_data = {
            'success': True,
            'input': input_matrix,
            'output': str(output),
            'method': 'IPS'
        }
        print_result(result_data, "IPS calculation completed")
    except Exception as e:
        result_data = {
            'success': False,
            'input': input_matrix,
            'output': str(output),
            'method': 'IPS',
            'error': str(e)
        }
        print_result(result_data, f"IPS calculation failed: {str(e)}")
    finally:
        sys.argv = original_argv


@analyze.command('mcpcounter')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file path')
@click.option('--features', required=True,
              type=click.Choice(['affy133P2_probesets', 'HUGO_symbols', 'ENTREZ_ID', 'ENSEMBL_ID']),
              help='Type of gene identifiers')
def mcpcounter_command(input_matrix: str, output: str, features: str):
    """Run MCPcounter estimation."""
    from iobrpy.workflow.mcpcounter import MCPcounter_estimate, preprocess_input
    import pandas as pd

    # Handle default output
    if not output:
        output = Path(input_matrix).parent / 'mcpcounter_results.txt'
    else:
        output = Path(output)

    click.echo(f"Running MCPcounter estimation...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Features: {features}")
    click.echo(f"  Output: {output}")

    # Load and preprocess input
    ext = Path(input_matrix).suffix.lower()
    if ext == '.csv':
        sep = ','
    elif ext in ('.tsv', '.txt'):
        sep = '\t'
    else:
        sep = None

    df = pd.read_csv(input_matrix, sep=sep, index_col=0,
                     engine='python' if sep is None else None)

    # Clean gene names
    df.index = df.index.str.replace(r'\..*', '', regex=True)
    df.index = df.index.str.upper().str.strip()

    if df.index.duplicated().any():
        dup_count = df.index.duplicated(keep=False).sum()
        click.echo(f"Warning: {dup_count} duplicate genes found; keeping max expression per gene.")
        df = df.groupby(df.index).max()

    # Run MCPcounter
    try:
        results = _run_with_stdout_routed(MCPcounter_estimate, df, features)
        results.T.to_csv(output, sep='\t', float_format='%.6f')

        result_data = {
            'success': True,
            'input': input_matrix,
            'output': str(output),
            'features': features,
            'method': 'MCPcounter'
        }
        print_result(result_data, "MCPcounter estimation completed")
    except Exception as e:
        result_data = {
            'success': False,
            'input': input_matrix,
            'output': str(output),
            'features': features,
            'method': 'MCPcounter',
            'error': str(e)
        }
        print_result(result_data, f"MCPcounter estimation failed: {str(e)}")


@analyze.command('deconvolute')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--method', type=click.Choice(['cibersort', 'ips', 'estimate', 'mcpcounter', 'quantiseq', 'epic', 'bayesprism'], case_sensitive=False),
              required=True, help='Deconvolution method')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--cibersort-perm', type=int, help='Number of permutations for CIBERSORT')
@click.option('--cibersort-qn/--no-cibersort-qn', default=True, help='Quantile normalization for CIBERSORT')
@click.option('--cibersort-absolute/--no-cibersort-absolute', default=False, help='Run CIBERSORT in absolute mode')
@click.option('--cibersort-abs-method', type=click.Choice(['sig.score', 'no.sumto1']), help='Absolute scoring method for CIBERSORT')
@click.option('--estimate-platform', type=click.Choice(['affymetrix', 'agilent', 'illumina']), help='Platform for ESTIMATE')
@click.option('--epic-reference', type=click.Choice(['TRef', 'BRef', 'both']), default='TRef', help='Reference for EPIC')
@click.option('--quantiseq-arrays', is_flag=True, help='Apply quantile normalization for array data')
@click.option('--quantiseq-signame', help='Signature name for quanTIseq')
@click.option('--quantiseq-tumor', is_flag=True, help='Calculate tumor score for quanTIseq')
@click.option('--quantiseq-scale-mrna', is_flag=True, help='Scale mRNA for quanTIseq')
@click.option('--quantiseq-method', type=click.Choice(['lsei', 'hampel', 'huber', 'bisquare']), help='Method for quanTIseq')
@click.option('--quantiseq-rmgenes', help='Genes to remove for quanTIseq')
@click.option('--mcpcounter-features', help='Feature columns for MCPcounter')
@click.option('--bayesprism-sc-dat', type=click.Path(exists=True), help='Single-cell data for BayesPrism')
@click.option('--bayesprism-cell-state-labels', type=click.Path(exists=True), help='Cell state labels for BayesPrism')
@click.option('--bayesprism-cell-type-labels', type=click.Path(exists=True), help='Cell type labels for BayesPrism')
@click.option('--bayesprism-key', help='Key for BayesPrism')
def deconvolute(input_matrix: str, output: str, method: str, threads: int,
                cibersort_perm: int, cibersort_qn: bool, cibersort_absolute: bool,
                cibersort_abs_method: str, estimate_platform: str, epic_reference: str,
                quantiseq_arrays: bool, quantiseq_signame: str, quantiseq_tumor: bool,
                quantiseq_scale_mrna: bool, quantiseq_method: str, quantiseq_rmgenes: str,
                mcpcounter_features: str, bayesprism_sc_dat: str, bayesprism_cell_state_labels: str,
                bayesprism_cell_type_labels: str, bayesprism_key: str):
    """Run a deconvolution method with advanced options."""
    analyzer = TMEAnalyzer()

    if not output:
        output = Path(input_matrix).parent / 'deconvolution_results'
    else:
        output = Path(output)

    output.mkdir(parents=True, exist_ok=True)

    click.echo(f"Running {method} deconvolution...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")

    # Prepare method-specific parameters
    method_params = {
        'method': TMEAnalyzer.parse_deconvolution_method(method),
        'input_matrix': Path(input_matrix),
        'output_file': output / f"{method}_results.csv",
        'threads': threads,
    }

    # Add method-specific parameters
    if method == 'cibersort':
        if cibersort_perm:
            method_params['perm'] = cibersort_perm
        method_params['qn'] = cibersort_qn
        method_params['absolute'] = cibersort_absolute
        if cibersort_abs_method:
            method_params['abs_method'] = cibersort_abs_method
    elif method == 'estimate':
        if estimate_platform:
            method_params['platform'] = estimate_platform
    elif method == 'epic':
        method_params['reference'] = epic_reference
    elif method == 'quantiseq':
        if quantiseq_arrays:
            method_params['arrays'] = True
        if quantiseq_signame:
            method_params['signame'] = quantiseq_signame
        if quantiseq_tumor:
            method_params['tumor'] = True
        if quantiseq_scale_mrna:
            method_params['scale_mrna'] = True
        if quantiseq_method:
            method_params['method'] = quantiseq_method
        if quantiseq_rmgenes:
            method_params['rmgenes'] = quantiseq_rmgenes
    elif method == 'mcpcounter':
        if mcpcounter_features:
            method_params['features'] = mcpcounter_features
    elif method == 'bayesprism':
        if bayesprism_sc_dat:
            method_params['sc_dat'] = Path(bayesprism_sc_dat)
        if bayesprism_cell_state_labels:
            method_params['cell_state_labels'] = Path(bayesprism_cell_state_labels)
        if bayesprism_cell_type_labels:
            method_params['cell_type_labels'] = Path(bayesprism_cell_type_labels)
        if bayesprism_key:
            method_params['key'] = bayesprism_key

    result = analyzer.run_deconvolution(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('lr-analysis')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--cancer-type', default='pancan', help='Cancer type')
@click.option('--id-type', default='ensembl', help='Gene ID type')
def lr_analysis(input_matrix: str, output: str, cancer_type: str, id_type: str):
    """Run ligand-receptor analysis."""
    analyzer = TMEAnalyzer()

    if not output:
        output = Path(input_matrix).parent / 'lr_scores.csv'

    click.echo(f"Running ligand-receptor analysis...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Cancer type: {cancer_type}")

    result = analyzer.calculate_lr_scores(
        input_matrix=Path(input_matrix),
        output_file=Path(output),
        cancer_type=cancer_type,
        id_type=id_type,
    )

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'pair_count': result.pair_count,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('combine')
@click.argument('result-dir', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output combined file')
@click.option('--methods', help='Comma-separated deconvolution methods to combine (default: all)')
@click.option('--prefix', default='combined', help='Prefix for output file')
def combine_results(result_dir: str, output: str, methods: str, prefix: str):
    """Combine multiple deconvolution results into a summary."""
    import glob

    click.echo(f"Combining deconvolution results from: {result_dir}")
    click.echo(f"  Methods: {methods if methods else 'all'}")
    click.echo(f"  Prefix: {prefix}")
    click.echo()

    result_path = Path(result_dir)
    if not output:
        output = result_path / f'{prefix}_deconvolution.csv'
    else:
        output = Path(output)

    # Find all deconvolution result files
    method_patterns = [
        '*cibersort_results*.csv',
        '*IPS_results*.csv',
        '*estimate_results*.csv',
        '*mcpcounter_results*.csv',
        '*quantiseq_results*.csv',
        '*epic_results*.csv',
        '*bayesprism_results*.csv',
    ]

    if methods:
        # Filter by specified methods
        method_names = [m.strip().lower() for m in methods.split(',')]
        method_patterns = [p for p in method_patterns
                       if any(m in p.lower() for m in method_names)]

    combined_data = None
    sample_names = None

    # Load all deconvolution results
    for pattern in method_patterns:
        for f in result_path.glob(pattern):
            try:
                import pandas as pd
                df = pd.read_csv(f, index_col=0)
                method_name = f.stem.replace('_results', '').replace('_result', '')

                # Add method prefix to columns
                df.columns = [f'{method_name}_{col}' for col in df.columns]

                if combined_data is not None:
                    # Merge with existing data
                    combined_data = pd.concat([combined_data, df], axis=1)
                else:
                    combined_data = df
                    sample_names = df.index.tolist()

                click.echo(f"  Loaded: {f.name}")
            except Exception as e:
                click.echo(f"  Warning: Failed to load {f.name}: {e}", err=True)

    if combined_data is None or combined_data.empty:
        result_data = {
            'success': False,
            'output': str(output),
            'message': 'No deconvolution results found',
        }
        print_result(result_data, result_data['message'])
        return

    # Write combined output
    try:
        combined_data.to_csv(output)
        result_data = {
            'success': True,
            'output': str(output),
            'methods': sorted({c.split('_')[0] for c in combined_data.columns if '_' in c}),
            'sample_count': len(sample_names) if sample_names else 0,
            'message': f'Combined {len(combined_data.columns)} deconvolution results',
        }
    except Exception as e:
        result_data = {
            'success': False,
            'output': str(output),
            'message': f'Failed to write combined output: {e}',
        }

    # Print summary
    if not _json_output and result_data.get('success'):
        click.echo(click.style(f"\nCombined {result_data['methods']} results", fg='green'))
        click.echo(f"  Samples: {result_data['sample_count']}")
        click.echo(f"  Output: {output}")
        if sample_names:
            click.echo(f"\nSamples: {', '.join(sample_names[:10])}")
            if len(sample_names) > 10:
                click.echo(f"  ... and {len(sample_names) - 10} more")

    print_result(result_data, result_data['message'])


# Deconvolution shortcut commands - convenient aliases for common methods

@analyze.command('cibersort')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--perm', type=int, help='Number of permutations')
@click.option('--qn/--no-qn', default=True, help='Quantile normalization')
@click.option('--absolute/--no-absolute', default=False, help='Run in absolute mode')
@click.option('--abs-method', type=click.Choice(['sig.score', 'no.sumto1']), help='Absolute scoring method')
def cibersort_cmd(input_matrix: str, output: str, threads: int, perm: int,
                  qn: bool, absolute: bool, abs_method: str):
    """Run CIBERSORT deconvolution (shortcut for 'analyze deconvolute --method cibersort')."""
    analyzer = TMEAnalyzer()

    if not output:
        output = str(Path(input_matrix).parent / 'cibersort_results.csv')

    click.echo(f"Running CIBERSORT deconvolution...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")

    method_params = {
        'method': DeconvolutionMethod.CIBERSORT,
        'input_matrix': Path(input_matrix),
        'output_file': Path(output),
        'threads': threads,
        'qn': qn,
        'absolute': absolute,
    }

    if perm:
        method_params['perm'] = perm
    if abs_method:
        method_params['abs_method'] = abs_method

    result = analyzer.run_deconvolution(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('quantiseq')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--arrays', is_flag=True, help='Apply array-mode normalization')
@click.option('--signame', help='Signature name')
@click.option('--tumor', is_flag=True, help='Calculate tumor score')
@click.option('--scale-mrna', is_flag=True, help='Scale mRNA')
@click.option('--method', type=click.Choice(['lsei', 'hampel', 'huber', 'bisquare']), help='Deconvolution method')
@click.option('--rmgenes', help='Genes to remove')
def quantiseq_cmd(input_matrix: str, output: str, threads: int,
                 arrays: bool, signame: str, tumor: bool,
                 scale_mrna: bool, method: str, rmgenes: str):
    """Run quanTIseq deconvolution (shortcut for 'analyze deconvolute --method quantiseq')."""
    analyzer = TMEAnalyzer()

    if not output:
        output = str(Path(input_matrix).parent / 'quantiseq_results.csv')

    click.echo(f"Running quanTIseq deconvolution...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")

    method_params = {
        'method': DeconvolutionMethod.QUANTISEQ,
        'input_matrix': Path(input_matrix),
        'output_file': Path(output),
        'threads': threads,
    }

    if arrays:
        method_params['arrays'] = True
    if signame:
        method_params['signame'] = signame
    if tumor:
        method_params['tumor'] = True
    if scale_mrna:
        method_params['scale_mrna'] = True
    if method:
        method_params['method'] = method
    if rmgenes:
        method_params['rmgenes'] = rmgenes

    result = analyzer.run_deconvolution(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('epic')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--threads', '-t', default=1, help='Thread count')
@click.option('--reference', type=click.Choice(['TRef', 'BRef', 'both']), default='TRef',
              help='Reference profile')
def epic_cmd(input_matrix: str, output: str, threads: int, reference: str):
    """Run EPIC deconvolution (shortcut for 'analyze deconvolute --method epic')."""
    analyzer = TMEAnalyzer()

    if not output:
        output = str(Path(input_matrix).parent / 'epic_results.csv')

    click.echo(f"Running EPIC deconvolution...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Reference: {reference}")

    method_params = {
        'method': DeconvolutionMethod.EPIC,
        'input_matrix': Path(input_matrix),
        'output_file': Path(output),
        'threads': threads,
        'reference': reference,
    }

    result = analyzer.run_deconvolution(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@analyze.command('bayesprism')
@click.argument('input-matrix', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of CPU cores')
@click.option('--sc-dat', type=click.Path(), help='Custom single-cell count matrix (genes x cells)')
@click.option('--cell-state-labels', type=click.Path(), help='Custom cell state labels file')
@click.option('--cell-type-labels', type=click.Path(), help='Custom cell type labels file')
@click.option('--key', help='Tumor key for custom reference (default: Malignant_cells)')
def bayesprism_cmd(input_matrix: str, output: str, threads: int,
                  sc_dat: str, cell_state_labels: str, cell_type_labels: str, key: str):
    """Run BayesPrism deconvolution (Bayesian with uncertainty estimates)."""
    analyzer = TMEAnalyzer()

    click.echo(f"Running BayesPrism deconvolution...")
    click.echo(f"  Input: {input_matrix}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")

    method_params = {
        'method': DeconvolutionMethod.BAYESPRISM,
        'input_matrix': Path(input_matrix),
        'output_file': Path(output),
        'threads': threads,
    }

    if sc_dat:
        method_params['sc_dat'] = Path(sc_dat)
    if cell_state_labels:
        method_params['cell_state_labels'] = Path(cell_state_labels)
    if cell_type_labels:
        method_params['cell_type_labels'] = Path(cell_type_labels)
    if key:
        method_params['key'] = key

    result = analyzer.run_deconvolution(**method_params)

    result_data = {
        'success': result.success,
        'output': str(result.output_file) if result.output_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


# ============================================================================
# Immune Analysis Commands (HLA/TCR/BCR)
# ============================================================================

@cli.group()
def immune():
    """Run immune repertoire analysis (HLA typing, TCR/BCR reconstruction)."""
    pass


@immune.command('trust4')
@click.option('--bam', '-b', type=click.Path(exists=True), help='Single BAM file')
@click.option('--bam-dir', type=click.Path(exists=True), help='Directory with BAM files')
@click.option('--fastq-dir', type=click.Path(exists=True), help='Directory with FASTQ files')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of threads')
@click.option('--batch', is_flag=True, help='Batch mode (process multiple samples)')
def trust4_cmd(bam: str, bam_dir: str, fastq_dir: str,
              output: str, threads: int, batch: bool):
    """Run TRUST4 TCR/BCR repertoire reconstruction."""
    from .core.hla_tcr import HLATCRAnalyzer, Trust4Mode

    analyzer = HLATCRAnalyzer()

    mode = Trust4Mode.BATCH if (batch or bam_dir or fastq_dir) else Trust4Mode.SINGLE

    # Validate input
    if not bam and not bam_dir and not fastq_dir:
        click.echo("Error: Must specify --bam, --bam-dir, or --fastq-dir")
        sys.exit(1)

    click.echo(f"Running TRUST4...")
    click.echo(f"  Output: {output}")
    click.echo(f"  Mode: {mode.value}")
    click.echo(f"  Threads: {threads}")

    result = analyzer.run_trust4(
        bam_file=Path(bam) if bam else None,
        bam_dir=Path(bam_dir) if bam_dir else None,
        fastq_dir=Path(fastq_dir) if fastq_dir else None,
        output_dir=Path(output),
        mode=mode,
        threads=threads,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'immune_data_file': str(result.immune_data_file) if result.immune_data_file else None,
        'indices_file': str(result.indices_file) if result.indices_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@immune.command('spechla')
@click.option('--sample-name', '-n', required=True, help='Sample name')
@click.option('--read1', '-1', type=click.Path(exists=True), required=True,
              help='Read1 FASTQ.gz')
@click.option('--read2', '-2', type=click.Path(exists=True), required=True,
              help='Read2 FASTQ.gz')
@click.option('--use-exon', '-u', type=click.Choice(['0', '1']), default='1',
              help='SpecHLA pipeline type (1=exon/RNA, 0=WGS)')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Number of threads')
def spechla_cmd(sample_name: str, read1: str, read2: str,
                use_exon: str, output: str, threads: int):
    """Run SpecHLA RNA-seq exon-level HLA typing."""
    from .core.hla_tcr import HLATCRAnalyzer

    analyzer = HLATCRAnalyzer()

    click.echo(f"Running SpecHLA HLA typing...")
    click.echo(f"  Sample: {sample_name}")
    click.echo(f"  Read1: {read1}")
    click.echo(f"  Read2: {read2}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Use exon: {use_exon}")
    click.echo(f"  Threads: {threads}")

    result = analyzer.run_spechla(
        sample_name=sample_name,
        read1=Path(read1),
        read2=Path(read2),
        output_dir=Path(output),
        threads=threads,
        use_exon=int(use_exon),
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


@immune.command('extract-hla-read')
@click.option('--sample', '-s', required=True, help='Sample name (e.g., NA12878)')
@click.option('--bam', '-b', type=click.Path(exists=True), required=True,
              help='Sorted and indexed BAM or CRAM file')
@click.option('--ref', '-r', type=click.Choice(['hg38', 'hg19']), required=True,
              help='Reference genome (hg38 or hg19)')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--no-auto-install', is_flag=True,
              help='Disable automatic installation of missing tools')
def extract_hla_read_cmd(sample: str, bam: str, ref: str,
                        output: str, no_auto_install: bool):
    """Extract HLA-related reads from BAM/CRAM files for HLA typing.

    This command extracts HLA-related reads from BAM/CRAM files and converts
    them to FASTQ format, which is a prerequisite for SpecHLA HLA typing.

    Required dependencies: samtools, bamutil, libdeflate=1.25, htslib=1.21
    These will be auto-installed via conda/mamba if not present.
    """
    from .core.hla_tcr import HLATCRAnalyzer

    analyzer = HLATCRAnalyzer()

    click.echo(f"Extracting HLA reads from BAM file...")
    click.echo(f"  Sample: {sample}")
    click.echo(f"  BAM: {bam}")
    click.echo(f"  Reference: {ref}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Auto-install: {not no_auto_install}")
    click.echo()

    result = analyzer.run_extract_hla_read(
        sample_id=sample,
        bam_path=Path(bam),
        ref=ref,
        outdir=Path(output),
        no_auto_install=no_auto_install,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_id': result.sample_id,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


@immune.command('hla-typing')
@click.option('--bam-dir', '-b', type=click.Path(exists=True), required=True,
              help='Directory containing BAM files')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--reference', '-r', type=click.Choice(['hg19', 'hg38']), default='hg38',
              help='Reference genome')
@click.option('--threads', '-t', default=8, help='Number of threads')
@click.option('--wgs', is_flag=True, help='Use WGS mode (default: RNA-seq)')
def hla_typing_cmd(bam_dir: str, output: str, reference: str,
                   threads: int, wgs: bool):
    """Run batch HLA typing workflow (extract + SpecHLA)."""
    from .core.hla_tcr import HLATCRAnalyzer

    analyzer = HLATCRAnalyzer()

    click.echo(f"Running HLA typing...")
    click.echo(f"  BAM directory: {bam_dir}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Reference: {reference}")
    click.echo(f"  Threads: {threads}")
    click.echo(f"  WGS mode: {wgs}")

    result = analyzer.run_hla_typing(
        bam_dir=Path(bam_dir),
        output_dir=Path(output),
        reference=reference,
        threads=threads,
        use_wgs=wgs,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'hla_result_file': str(result.hla_result_file) if result.hla_result_file else None,
        'message': result.message,
    }

    print_result(result_data, result.message if result.success else "Failed")


# ============================================================================
# Quantification Commands
# ============================================================================

@cli.group()
def quantify():
    """Run quantification workflows."""
    pass


@quantify.command('qc')
@click.argument('fastq-dir', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output directory')
@click.option('--threads', '-t', default=8, help='Thread count')
@click.option('--suffix', default='_R1', help='FASTQ suffix pattern')
@click.option('--single-end', is_flag=True, help='Single-end reads')
@click.option('--min-length', default=50, help='Minimum read length')
def fastq_qc(fastq_dir: str, output: str, threads: int, suffix: str,
              single_end: bool, min_length: int):
    """Run FASTQ quality control."""
    workflow = QuantificationWorkflow()

    if not output:
        output = Path(fastq_dir).parent / 'fastq_qc_output'

    click.echo(f"Running FASTQ QC...")
    click.echo(f"  Input: {fastq_dir}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")

    result = workflow.run_fastq_qc(
        fastq_dir=Path(fastq_dir),
        output_dir=Path(output),
        threads=threads,
        suffix1=suffix,
        single_end=single_end,
        length_required=min_length,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


@quantify.command('runall')
@click.argument('fastq-dir', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--mode', type=click.Choice(['salmon', 'star']), default='salmon', help='Quantification mode')
@click.option('--index', type=click.Path(), required=True, help='Path to Salmon/STAR index')
@click.option('--threads', '-t', default=8, help='Thread count')
@click.option('--batch-size', default=1, help='Batch size for parallel processing')
@click.option('--suffix', default='_R1', help='FASTQ suffix pattern')
@click.option('--min-length', default=50, help='Minimum read length')
@click.option('--return-feature', default='symbol', help='Gene feature type (ENST, ENSG, symbol)')
@click.option('--idtype', default='ensembl', help='Gene ID type for count2tpm')
@click.option('--org', default='hsa', help='Organism (hsa or mmus)')
@click.option('--deconvolution-methods', help='Comma-separated deconvolution methods')
@click.option('--no-lr', is_flag=True, help='Skip ligand-receptor analysis')
@click.option('--signatures', default='all', help='Signature groups')
@click.option('--run-trust4', is_flag=True, help='Run TRUST4 TCR/BCR analysis')
@click.option('--resume', is_flag=True, help='Resume from existing outputs')
@click.option('--dry-run', is_flag=True, help='Print commands without executing')
@click.option('--lr-cancer-type', default='pancan', help='Cancer type for LR analysis')
def runall_cmd(fastq_dir: str, output: str, mode: str, index: str,
                threads: int, batch_size: int, suffix: str,
                min_length: int, return_feature: str, idtype: str,
                org: str, deconvolution_methods: str, no_lr: bool,
                signatures: str, run_trust4: bool, resume: bool,
                dry_run: bool, lr_cancer_type: str):
    """Run complete end-to-end pipeline from FASTQ to TME results."""
    workflow = QuantificationWorkflow()

    # Parse deconvolution methods
    deconv_methods = None
    if deconvolution_methods:
        try:
            deconv_methods = [
                TMEAnalyzer.parse_deconvolution_method(method_name).value
                for method_name in deconvolution_methods.split(',')
                if method_name.strip()
            ]
        except ValueError as exc:
            print_json(error_response(str(exc)))
            sys.exit(1)

    click.echo(f"Running complete pipeline...")
    click.echo(f"  Mode: {mode}")
    click.echo(f"  FASTQ input: {fastq_dir}")
    click.echo(f"  Index: {index}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")
    click.echo(f"  Resume: {resume}")
    click.echo()

    results = workflow.runall(
        fastq_dir=Path(fastq_dir),
        output_dir=Path(output),
        mode=QuantificationMode(mode),
        index=Path(index),
        threads=threads,
        batch_size=batch_size,
        suffix1=suffix,
        length_required=min_length,
        return_feature=return_feature,
        idtype=idtype,
        org=org,
        deconvolution_methods=deconv_methods,
        include_lr=not no_lr,
        signature_groups=signatures,
        run_trust4=run_trust4,
        resume=resume,
        dry_run=dry_run,
        lr_cancer_type=lr_cancer_type,
    )

    # Format results
    result_data = {
        'output_dir': str(output),
        'steps': results['steps'],
        'final': results['final'],
    }

    # Print step-by-step summary
    if not _json_output:
        click.echo()
        click.echo(click.style("Pipeline Steps Summary:", fg='blue', bold=True))
        for step, result in results['steps'].items():
            status_icon = "OK" if result.get('success') else "FAIL"
            click.echo(f"  {status_icon} {step:15}: {result.get('message', '')}")
        click.echo()
        click.echo(click.style(results['final']['message'], fg='green' if results['final']['success'] else 'red'))

    print_result(result_data, results['final']['message'])


@quantify.command('clean')
@click.argument('output-dir', type=click.Path())
@click.option('--keep-intermediates', is_flag=True, help='Keep intermediate files, only clean failed run artifacts')
@click.option('--dry-run', is_flag=True, help='Show what would be removed without deleting')
def quantify_clean(output_dir: str, keep_intermediates: bool, dry_run: bool):
    """Clean up intermediate files from pipeline run."""
    from .core import QuantificationWorkflow

    workflow = QuantificationWorkflow()

    click.echo(f"Cleaning pipeline output: {output_dir}")
    click.echo(f"  Keep intermediates: {keep_intermediates}")
    click.echo(f"  Dry run: {dry_run}")
    click.echo()

    result = workflow.clean_pipeline(
        output_dir=Path(output_dir),
        keep_intermediates=keep_intermediates,
        dry_run=dry_run,
    )

    result_data = {
        'success': result.success,
        'files_removed': result.files_removed,
        'dirs_removed': result.dirs_removed,
        'space_freed_mb': result.space_freed_mb,
        'message': result.message,
    }

    print_result(result_data, result.message)


@quantify.command('validate')
@click.argument('output-dir', type=click.Path(exists=True))
@click.option('--mode', type=click.Choice(['salmon', 'star']), default='salmon',
              help='Quantification mode')
@click.option('--check-tpm/--no-check-tpm', default=True, help='Verify TPM matrix exists')
def quantify_validate(output_dir: str, mode: str, check_tpm: bool):
    """Validate pipeline output directory for completeness."""
    from .core import QuantificationWorkflow, QuantificationMode

    workflow = QuantificationWorkflow()

    click.echo(f"Validating pipeline output: {output_dir}")
    click.echo(f"  Mode: {mode}")
    click.echo(f"  Check TPM: {check_tpm}")
    click.echo()

    result = workflow.validate_pipeline(
        output_dir=Path(output_dir),
        mode=QuantificationMode(mode),
        check_tpm=check_tpm,
    )

    result_data = {
        'success': result.success,
        'is_valid': result.is_valid,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'completed_steps': result.completed_steps or [],
        'missing_steps': result.missing_steps or [],
        'issues': result.issues or [],
    }

    # Print detailed validation results
    if not _json_output:
        click.echo(click.style(f"{result.message}", fg='green' if result.is_valid else 'red'))
        if result.completed_steps:
            click.echo("\nCompleted steps:")
            for step in result.completed_steps:
                click.echo(f"  OK {step}")
        if result.missing_steps:
            click.echo("\nMissing steps:")
            for step in result.missing_steps:
                click.echo(f"  FAIL {step}", color='red')
        if result.issues:
            click.echo("\nIssues found:")
            for issue in result.issues:
                click.echo(f"  ! {issue}", color='yellow')

    print_result(result_data, result.message)


@quantify.command('scan-fastq')
@click.argument('fastq-dir', type=click.Path(exists=True))
@click.option('--suffix1', default='_R1', help='R1 FASTQ suffix pattern')
@click.option('--suffix2', default='_R2', help='R2 FASTQ suffix pattern')
@click.option('--single-end', is_flag=True, help='Single-end reads')
def quantify_scan_fastq(fastq_dir: str, suffix1: str, suffix2: str, single_end: bool):
    """Scan FASTQ directory and report structure."""
    from .core import QuantificationWorkflow

    workflow = QuantificationWorkflow()

    click.echo(f"Scanning FASTQ directory: {fastq_dir}")
    click.echo(f"  R1 suffix: {suffix1}")
    click.echo(f"  R2 suffix: {suffix2}")
    click.echo(f"  Single-end: {single_end}")
    click.echo()

    result = workflow.scan_fastq(
        fastq_dir=Path(fastq_dir),
        suffix1=suffix1,
        suffix2=suffix2,
        single_end=single_end,
    )

    result_data = {
        'success': result.success,
        'sample_count': result.sample_count,
        'samples': result.samples or [],
        'issues': result.issues or [],
    }

    # Print detailed scan results
    if not _json_output:
        click.echo(click.style(result.message, fg='green' if result.success else 'red'))
        if result.samples:
            click.echo(f"\n{len(result.samples)} samples found:")
            for i, sample in enumerate(result.samples[:20], 1):  # Show first 20
                click.echo(f"  {i}. {sample}")
            if len(result.samples) > 20:
                click.echo(f"  ... and {len(result.samples) - 20} more")
        if result.issues:
            click.echo("\nIssues found:")
            for issue in result.issues:
                click.echo(f"  ! {issue}", color='yellow')

    print_result(result_data, result.message)


# ============================================================================
# HLA/TCR/BCR Commands
# ============================================================================

@cli.group()
def hla_tcr():
    """Run HLA typing and TCR/BCR analysis."""
    pass


@hla_tcr.command('trust4')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Thread count')
@click.option('--mode', type=click.Choice(['single', 'batch']), default='batch', help='Execution mode')
def trust4_cmd(input_path: str, output: str, threads: int, mode: str):
    """Run TRUST4 TCR/BCR repertoire reconstruction."""
    from .hla_tcr import HLATCRAnalyzer, Trust4Mode

    analyzer = HLATCRAnalyzer()
    input_dir = Path(input_path)

    trust4_mode = Trust4Mode.SINGLE if mode == 'single' else Trust4Mode.BATCH

    click.echo(f"Running TRUST4 analysis...")
    click.echo(f"  Input: {input_path}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Mode: {mode}")
    click.echo(f"  Threads: {threads}")
    click.echo()

    result = analyzer.run_trust4(
        bam_dir=input_dir if input_dir.is_dir() else None,
        bam_file=input_dir if input_dir.is_file() else None,
        fastq_dir=input_dir,
        output_dir=Path(output),
        mode=trust4_mode,
        threads=threads,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'immune_data_file': str(result.immune_data_file) if result.immune_data_file else None,
        'indices_file': str(result.indices_file) if result.indices_file else None,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


@hla_tcr.command('spechla')
@click.option('--use-exon', '-u', type=click.Choice([0, 1]), default=1, help='SpecHLA pipeline type (0=WGS, 1=exon/RNA)')
@click.argument('bam-dir', type=click.Path(exists=True), required=True)
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--threads', '-t', default=8, help='Thread count')
def spechla_cmd(use_exon: int, bam_dir: str, output: str, threads: int):
    """Run SpecHLA HLA typing."""
    from .hla_tcr import HLATCRAnalyzer, SpecHLAResult

    analyzer = HLATCRAnalyzer()

    click.echo(f"Running SpecHLA...")
    click.echo(f"  BAM directory: {bam_dir}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Threads: {threads}")
    click.echo(f"  Mode: {'WGS' if use_exon else 'exon/RNA'}")
    click.echo()

    result = analyzer.run_spechla(
        bam_dir=Path(bam_dir),
        output_dir=Path(output),
        threads=threads,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'hla_result_file': str(result.hla_result_file) if result.hla_result_file else None,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


@hla_tcr.command('hla')
@click.option('--bam-dir', '-b', type=click.Path(exists=True), required=True, help='Directory containing BAM files')
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--reference', '-r', type=click.Choice(['hg19', 'hg38']), default='hg38', help='Reference genome')
@click.option('--threads', '-t', default=8, help='Thread count')
@click.option('--use-wgs', '-w', is_flag=True, help='Use WGS mode')
def hla_typing_cmd(bam_dir: str, output: str, reference: str, threads: int, use_wgs: bool):
    """Run complete HLA typing workflow (ExtractHLA + SpecHLA)."""
    from .hla_tcr import HLATCRAnalyzer

    analyzer = HLATCRAnalyzer()

    click.echo(f"Running HLA typing...")
    click.echo(f"  BAM directory: {bam_dir}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Reference: {reference}")
    click.echo(f"  Threads: {threads}")
    click.echo(f"  Mode: {'WGS' if use_wgs else 'exon/RNA'}")
    click.echo()

    result = analyzer.run_hla_typing(
        bam_dir=Path(bam_dir),
        output_dir=Path(output),
        reference=reference,
        threads=threads,
        use_exon=1 if use_wgs else 0,
    )

    result_data = {
        'success': result.success,
        'output_dir': str(result.output_dir) if result.output_dir else None,
        'sample_count': result.sample_count,
        'spechla_result_dir': str(result.spechla_result_dir) if result.spechla_result_dir else None,
        'message': result.message,
        'duration_ms': result.duration_ms,
    }

    print_result(result_data, result.message if result.success else "Failed")


# ============================================================================
# Workflow Commands
# ============================================================================

@cli.group()
def workflow():
    """Run individual workflow steps."""
    pass


@workflow.command('prepare-salmon')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output matrix file')
@click.option('--return-feature', default='symbol', type=click.Choice(['ENST','ENSG','symbol']),
              help='Gene feature to retain')
@click.option('--remove-version', is_flag=True, help='Remove version suffix from gene IDs')
def prepare_salmon(input_path: str, output: str, return_feature: str, remove_version: bool):
    """Prepare Salmon data matrix."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.prepare_salmon(input_path, output, return_feature, remove_version)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('count2tpm')
@click.argument('count-mat', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output TPM matrix file')
@click.option('--efflength-csv', type=click.Path(), help='Effective length CSV file')
@click.option('--idtype', default='ensembl', type=click.Choice(['ensembl','entrez','symbol','mgi']),
              help='Gene ID type')
@click.option('--org', default='hsa', type=click.Choice(['hsa','mmus']), help='Organism')
@click.option('--source', default='local', type=click.Choice(['local','biomart']),
              help='Gene length source')
@click.option('--id', 'id_col', default='id', help='ID column name in effLength CSV')
@click.option('--length', 'length_col', default='eff_length', help='Length column name')
@click.option('--gene-symbol', 'gene_symbol_col', default='symbol', help='Symbol column name')
@click.option('--check-data', is_flag=True, help='Check and remove missing values')
@click.option('--remove-version', is_flag=True, help='Remove version suffix from gene IDs')
def count2tpm(count_mat: str, output: str, efflength_csv: str, idtype: str, org: str,
              source: str, id_col: str, length_col: str, gene_symbol_col: str,
              check_data: bool, remove_version: bool):
    """Convert count matrix to TPM."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.count2tpm(
        count_mat, output, efflength_csv, idtype, org, source,
        id_col, length_col, gene_symbol_col, check_data, remove_version
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('anno-eset')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output annotated file')
@click.option('--annotation', required=True,
              type=click.Choice(['anno_hug133plus2','anno_rnaseq','anno_illumina','anno_grch38']),
              help='Annotation key')
@click.option('--annotation-file', type=click.Path(), help='External annotation file')
@click.option('--annotation-key', help='Key for external annotation file')
@click.option('--symbol', default='symbol', help='Annotation symbol column')
@click.option('--probe', default='id', help='Annotation probe column')
@click.option('--method', default='mean', type=click.Choice(['mean','sd','sum']),
              help='Duplicate handling method')
@click.option('--remove-version', is_flag=True, help='Remove version suffix from gene IDs')
def anno_eset(input_path: str, output: str, annotation: str, annotation_file: str,
               annotation_key: str, symbol: str, probe: str, method: str, remove_version: bool):
    """Annotate expression set."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.anno_eset(
        input_path, output, annotation, symbol, probe, method,
        annotation_file, annotation_key, remove_version
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('calculate-sig-score')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output scores file')
@click.option('--signature', required=True, help='Signature groups (comma-separated)')
@click.option('--method', default='pca', type=click.Choice(['pca','zscore','ssgsea','integration']),
              help='Scoring method')
@click.option('--mini-gene-count', default=3, type=int, help='Minimum genes per signature')
@click.option('--adjust-eset', is_flag=True, help='Apply additional filtering')
@click.option('--parallel-size', default=1, type=int, help='Number of threads')
def calculate_sig_score(input_path: str, output: str, signature: str, method: str,
                        mini_gene_count: int, adjust_eset: bool, parallel_size: int):
    """Calculate signature scores."""
    # Split signature by comma or space
    sig_list = [s.strip() for s in signature.replace(',', ' ').split() if s.strip()]
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.calculate_sig_score(
        input_path, output, sig_list, method, mini_gene_count,
        adjust_eset, parallel_size
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('batch-salmon')
@click.option('--index', required=True, type=click.Path(exists=True), help='Salmon index path')
@click.option('--path-fq', required=True, type=click.Path(), help='FASTQ directory')
@click.option('--path-out', required=True, type=click.Path(), help='Output directory')
@click.option('--suffix1', default='_1.fastq.gz', help='R1 suffix pattern')
@click.option('--batch-size', default=1, type=int, help='Number of concurrent samples')
@click.option('--num-threads', default=8, type=int, help='Threads per process')
@click.option('--gtf', type=click.Path(), help='Optional GTF file')
def batch_salmon(index: str, path_fq: str, path_out: str, suffix1: str,
                 batch_size: int, num_threads: int, gtf: str):
    """Batch-run Salmon quantification."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.batch_salmon(index, path_fq, path_out, suffix1, batch_size, num_threads, gtf)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('merge-salmon')
@click.option('--path-salmon', required=True, type=click.Path(exists=True),
              help='Salmon output directory')
@click.option('--project', required=True, help='Output file prefix')
@click.option('--num-processes', type=int, help='Number of processes for loading')
def merge_salmon(path_salmon: str, project: str, num_processes: int):
    """Merge Salmon quant.sf files."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.merge_salmon(path_salmon, project, num_processes)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('batch-star-count')
@click.option('--index', required=True, type=click.Path(exists=True), help='STAR index path')
@click.option('--path-fq', required=True, type=click.Path(), help='FASTQ directory')
@click.option('--path-out', required=True, type=click.Path(), help='Output directory')
@click.option('--suffix1', default='_1.fastq.gz', help='R1 suffix pattern')
@click.option('--batch-size', default=1, type=int, help='Number of samples per batch')
@click.option('--num-threads', default=8, type=int, help='Threads per process')
def batch_star_count(index: str, path_fq: str, path_out: str, suffix1: str,
                     batch_size: int, num_threads: int):
    """Batch-run STAR with GeneCounts."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.batch_star_count(index, path_fq, path_out, suffix1, batch_size, num_threads)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('merge-star-count')
@click.option('--path', required=True, type=click.Path(exists=True), help='STAR output directory')
@click.option('--project', required=True, help='Output file prefix')
def merge_star_count(path: str, project: str):
    """Merge STAR *_ReadsPerGene.out.tab files."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.merge_star_count(path, project)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('fastq-qc')
@click.option('--path1-fastq', required=True, type=click.Path(), help='FASTQ input directory')
@click.option('--path2-fastp', required=True, type=click.Path(), help='Output directory')
@click.option('--num-threads', default=8, type=int, help='Number of threads')
@click.option('--suffix1', default='_1.fastq.gz', help='R1 suffix pattern')
@click.option('--batch-size', default=1, type=int, help='Number of concurrent samples')
@click.option('--se', is_flag=True, help='Single-end mode')
@click.option('--length-required', default=50, type=int, help='Minimum read length')
def fastq_qc(path1_fastq: str, path2_fastp: str, num_threads: int, suffix1: str,
              batch_size: int, se: bool, length_required: int):
    """Run FASTQ quality control."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.fastq_qc(path1_fastq, path2_fastp, num_threads, suffix1, batch_size, se, length_required)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('log2-eset')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
def log2_eset(input_path: str, output: str):
    """Apply log2(x+1) to expression matrix."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.log2_eset(input_path, output)
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('tme-cluster')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
@click.option('--features', help='Feature columns (e.g., "1:22")')
@click.option('--pattern', help='Regex pattern for feature columns')
@click.option('--id', help='Sample ID column name')
@click.option('--scale/--no-scale', default=True, help='Z-score scaling')
@click.option('--min-nc', default=2, type=int, help='Minimum clusters')
@click.option('--max-nc', default=6, type=int, help='Maximum clusters')
@click.option('--max-iter', default=10, type=int, help='Max iterations')
@click.option('--tol', default=1e-4, type=float, help='Convergence tolerance')
@click.option('--print-intermediate', is_flag=True, help='Print intermediate results')
@click.option('--input-sep', help='Input field separator')
@click.option('--output-sep', help='Output field separator')
def tme_cluster(input_path: str, output: str, features: str, pattern: str, id: str,
                scale: bool, min_nc: int, max_nc: int, max_iter: int, tol: float,
                print_intermediate: bool, input_sep: str, output_sep: str):
    """Run TME clustering."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.tme_cluster(
        input_path, output, features, pattern, id, scale, min_nc, max_nc,
        max_iter, tol, print_intermediate, input_sep, output_sep
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('nmf-cluster')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output directory')
@click.option('--kmin', default=2, type=int, help='Minimum k value')
@click.option('--kmax', default=8, type=int, help='Maximum k value')
@click.option('--features', help='Feature columns (e.g., "2-10")')
@click.option('--log1p', is_flag=True, help='Apply log1p transform')
@click.option('--normalize', is_flag=True, help='Apply L1 normalization')
@click.option('--shift', type=float, help='Shift values to make non-negative')
@click.option('--random-state', default=42, type=int, help='Random state')
@click.option('--max-iter', default=1000, type=int, help='Max iterations')
@click.option('--skip-k-2', is_flag=True, help='Skip k=2')
def nmf_cluster(input_path: str, output: str, kmin: int, kmax: int, features: str,
                log1p: bool, normalize: bool, shift: float, random_state: int,
                max_iter: int, skip_k_2: bool):
    """Run NMF-based clustering."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.nmf_cluster(
        input_path, output, kmin, kmax, features, log1p, normalize, shift,
        random_state, max_iter, skip_k_2
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('mouse2human')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
@click.option('--is-matrix', is_flag=True, help='Input is a matrix (genes x samples)')
@click.option('--column-of-symbol', help='Column with gene symbols')
@click.option('--verbose', is_flag=True, help='Verbose output')
@click.option('--sep', help='Input separator')
@click.option('--out-sep', help='Output separator')
@click.option('--progress/--no-progress', default=True, help='Show progress bar')
def mouse2human(input_path: str, output: str, is_matrix: bool, column_of_symbol: str,
                verbose: bool, sep: str, out_sep: str, progress: bool):
    """Convert mouse gene symbols to human symbols."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.mouse2human(
        input_path, output, is_matrix, column_of_symbol, verbose, sep, out_sep, progress
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('remove-version')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
@click.option('--is-matrix', is_flag=True, help='Input is a matrix (genes x samples)')
@click.option('--column-of-symbol', help='Column with gene symbols')
@click.option('--verbose', is_flag=True, help='Verbose output')
@click.option('--sep', help='Input separator')
@click.option('--out-sep', help='Output separator')
def remove_version(input_path: str, output: str, is_matrix: bool, column_of_symbol: str,
                   verbose: bool, sep: str, out_sep: str):
    """Remove version suffixes from gene IDs."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.remove_version(
        input_path, output, is_matrix, column_of_symbol, verbose, sep, out_sep
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('mouse2human-eset')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
@click.option('--is-matrix', is_flag=True, help='Input is a matrix (genes x samples)')
@click.option('--column-of-symbol', help='Column with gene symbols')
@click.option('--verbose', is_flag=True, help='Verbose output')
@click.option('--sep', help='Input separator')
@click.option('--out-sep', help='Output separator')
@click.option('--progress/--no-progress', default=True, help='Show progress bar')
def mouse2human_eset(input_path: str, output: str, is_matrix: bool, column_of_symbol: str,
                     verbose: bool, sep: str, out_sep: str, progress: bool):
    """Convert mouse gene symbols to human symbols for expression set."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.mouse2human_eset(
        input_path, output, is_matrix, column_of_symbol, verbose, sep, out_sep, progress
    )
    print_result(result.to_dict(), result.message if result.success else "Failed")


@workflow.command('lr-cal')
@click.argument('input-path', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), required=True, help='Output file')
@click.option('--data-type', default='tpm', type=click.Choice(['count','tpm']),
              help='Input data type')
@click.option('--id-type', default='ensembl', help='Gene ID type')
@click.option('--cancer-type', default='pancan', help='Cancer type')
@click.option('--verbose', is_flag=True, help='Verbose output')
def lr_cal(input_path: str, output: str, data_type: str, id_type: str, cancer_type: str,
            verbose: bool):
    """Calculate ligand-receptor interactions."""
    executor = WorkflowExecutor(verbose=not _json_output)
    result = executor.lr_cal(input_path, output, data_type, id_type, cancer_type, verbose)
    print_result(result.to_dict(), result.message if result.success else "Failed")


# ============================================================================
# Version Command
# ============================================================================

@cli.command('version')
@click.option('--external/--no-external', default=True, help='Check external tool versions')
def version_cmd(external: bool):
    """Show version information for IOBRpy and external tools."""
    from .core import QuantificationWorkflow

    workflow = QuantificationWorkflow()

    click.echo("Checking version information...")
    click.echo(f"  External tools: {external}")
    click.echo()

    result = workflow.get_version_info(check_external=external)

    result_data = {
        'iobrpy_version': result.iobrpy_version,
        'iobrpy_installed': result.iobrpy_installed,
        'external_tools': result.external_tools or {},
        'python_version': result.python_version,
    }

    # Print detailed version results
    if not _json_output:
        click.echo(click.style(result.message, fg='blue'))
        if result.external_tools:
            click.echo("\nExternal Tools:")
            for tool, version in result.external_tools.items():
                status_color = 'green' if 'not found' not in version and 'error' not in version else 'red'
                click.echo(f"  {tool:10} {click.style(version, fg=status_color)}")

    print_result(result_data, result.message)


# ============================================================================
# Export Commands
# ============================================================================

@cli.group()
def export():
    """Export analysis results."""
    pass


@export.command('json')
@click.argument('data-file', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
@click.option('--pretty', is_flag=True, help='Pretty-print JSON')
def export_json(data_file: str, output: str, pretty: bool):
    """Export file as JSON."""
    exporter = Exporter(Path(data_file).parent)

    # Read the data file
    try:
        import pandas as pd
        df = pd.read_csv(data_file)
        data = df.to_dict(orient='records')
    except Exception:
        with open(data_file, 'r') as f:
            data = f.read()

    result = exporter.export_json(data, Path(output) if output else None, pretty=pretty)
    print_result({'output': result.output_path}, result.message if result.success else "Failed")


@export.command('summary')
@click.argument('results-dir', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), help='Output file')
def export_summary(results_dir: str, output: str):
    """Create a summary report of analysis results."""
    exporter = Exporter(Path(results_dir))

    result = exporter.create_summary_report(
        Path(results_dir),
        Path(output) if output else None,
    )

    print_result({'output': result.output_path}, result.message if result.success else "Failed")


# ============================================================================
# Utility Commands
# ============================================================================

@cli.command()
@click.argument('path', type=click.Path())
def info(path: str):
    """Show file or directory information."""
    p = Path(path)

    if p.is_file():
        result = format_file_info(p)
    elif p.is_dir():
        result = format_directory_info(p, recursive=False)
    else:
        result = error_response(f"Path not found: {path}")

    print_json(result)


@cli.command()
@click.argument('path', type=click.Path(exists=True))
def ls(path: str):
    """List directory contents."""
    p = Path(path)

    if not p.is_dir():
        print_json(error_response(f"Not a directory: {path}"))
        sys.exit(1)

    result = format_directory_info(p, recursive=False)
    print_json(result)


@cli.command()
@click.option('--details', is_flag=True, help='Include full tool schemas')
def commands(details: bool):
    """List native iobrpy commands with semantic summaries, inputs, outputs, and default behaviors."""
    tool_summaries = _native_command_summary()
    payload: Dict[str, Any] = {
        "success": True,
        "entrypoint": "iobrpy-cli",
        "native_invocation_pattern": "iobrpy-cli <native_command> [native_args...]",
        "path_shorthand": "iobrpy-cli <existing_path> -> map --path <existing_path>",
        "helper_commands": [
            {"name": name, "description": _HELPER_DESCRIPTIONS.get(name, "")}
            for name in _helper_command_names()
        ],
        "native_command_count": len(tool_summaries),
        "commands": tool_summaries,
        "policy": {
            "excluded_native_commands": ["deside", "ai"],
            "preferred_first_step": "map",
            "compatibility_namespaces": ["analyze", "quantify", "workflow", "immune", "hla_tcr"],
        },
    }
    if details:
        payload["tool_schemas"] = get_native_tools()
    _print_structured(payload)


@cli.command()
@click.option('--path', default='.', show_default=True, type=click.Path(), help='Input path to place on the IOBRpy roadmap')
@click.option('--max-depth', default=8, show_default=True, type=int, help='Maximum directory depth to inspect')
@click.option('--max-entries', default=8000, show_default=True, type=int, help='Maximum file or directory names to scan. Use 0 to disable the entry cap, which is useful for focused deep scans on complex result directories.')
@click.option('--quick', is_flag=True, hidden=True, help='Deprecated and ignored; map now always uses the default deep/full scan path.')
@click.option('--focus', multiple=True, type=click.Path(), help='Deep-scan one or more likely result subdirectories while only shallowly censusing other large generic branches')
@click.option('--investigate-existing', is_flag=True, help='Run CLI-native fallback investigation for nonstandard existing-result layouts')
def map(path: str, max_depth: int, max_entries: int, quick: bool, focus: Tuple[str, ...], investigate_existing: bool):
    """Detect the current IOBRpy pipeline stage for a path and suggest whether to continue or rerun.

    The default behavior is the deeper scan path, which is safer for structurally
    complex directories. The legacy ``--quick`` flag is accepted for backward
    compatibility but ignored. For deeper focused scans, ``--max-entries 0``
    disables the entry cap.
    """
    requested_quick = quick
    payload = map_pipeline_path(
        Path(path),
        max_depth=max_depth,
        max_entries=max_entries,
        quick=False,
        focus=focus,
        investigate_existing=investigate_existing,
        auto_scan_retry=True,
        auto_investigate_existing=True,
    )
    if requested_quick:
        payload["quick_request_ignored"] = True
        payload["quick_mode_deprecated"] = True
        payload["quick_mode_deprecation_note"] = (
            "The legacy --quick flag is ignored; map now uses the default deep/full scan path."
        )
    payload["entrypoint"] = "iobrpy-cli"
    payload["next_steps"] = [
        "Continue from the current stage when resume_recommended is true instead of rerunning the full pipeline by default.",
        "Run iobrpy-cli recommend --path <path> --task '<task>' --json after map when you want a command recommendation for the next biological task.",
        "Confirm environment-specific required parameters such as index and output paths with the user before execution.",
    ]
    _print_pipeline_map_result(payload)
    if not payload.get("success"):
        raise click.exceptions.Exit(1)


@cli.command()
@click.option('--path', type=click.Path(), help='Input path to inspect')
@click.option('--task', default='', help='Natural-language task description')
@click.option('--external/--no-external', default=False, help='Check external tools like salmon/star/trust4')
def doctor(path: str, task: str, external: bool):
    """Show environment diagnostics plus the preferred analysis entrypoints for the current task/path."""
    version_info: Dict[str, Any] = {}
    root_version = _distribution_version("iobrpy")
    standalone_harness_version = _distribution_version("iobrpy-cli")
    if external:
        workflow = QuantificationWorkflow()
        result = workflow.get_version_info(check_external=True)
        version_info = {
            "iobrpy_version": result.iobrpy_version,
            "iobrpy_installed": result.iobrpy_installed,
            "external_tools": result.external_tools or {},
            "python_version": result.python_version,
            "message": result.message,
        }
    else:
        try:
            version_info = {
                "iobrpy_version": importlib_metadata.version("iobrpy"),
                "iobrpy_installed": True,
                "python_version": platform.python_version(),
            }
        except importlib_metadata.PackageNotFoundError:
            version_info = {
                "iobrpy_version": None,
                "iobrpy_installed": False,
                "python_version": platform.python_version(),
            }

    detected = _detect_input_summary(Path(path)) if path else None
    pipeline_map = analyze_pipeline_path(Path(path)) if path else None
    recommendation = _recommend_for_summary(detected, task) if detected else None
    payload = {
        "success": True,
        "entrypoint": "iobrpy-cli",
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "cwd": os.getcwd(),
        "console_scripts": _console_script_names(),
        "commands_on_path": {
            name: _command_on_path(name)
            for name in ("iobrpy", "iobrpy-cli", "iobrpy-cli-mcp")
        },
        "installed_distributions": {
            "iobrpy": root_version,
            "iobrpy-cli": standalone_harness_version,
        },
        "package_conflicts": (
            [
                "A standalone iobrpy-cli distribution is also installed. "
                "If its code is stale, it can shadow the harness bundled with iobrpy."
            ]
            if standalone_harness_version
            else []
        ),
        "version_info": version_info,
        "detected_input": detected,
        "pipeline_map": pipeline_map,
        "recommendation": recommendation,
        "agent_manifest": str(Path.cwd() / "agent-manifest.json"),
        "guidance": [
            "Prefer iobrpy-cli <native_command> or iobrpy-cli-mcp before scanning source files.",
            "Do not substitute the R IOBR package when Python iobrpy is available.",
            "Run iobrpy-cli map --path <path> --json first for path-driven requests so you know whether to continue downstream or rerun.",
            "When recommend returns confirm_parameters, confirm those before execution.",
        ],
    }
    _print_structured(payload)


@cli.command()
@click.option('--path', type=click.Path(), help='Input path to inspect')
@click.option('--task', default='', help='Natural-language task description')
def recommend(path: str, task: str):
    """Classify an input path and recommend the correct native iobrpy workflow using the required-params JSON."""
    if not path and not task:
        _print_structured(
            error_response("Provide at least one of --path or --task so the harness can recommend a workflow.")
        )
        raise click.exceptions.Exit(1)

    detected = _detect_input_summary(Path(path)) if path else {"kind": "unspecified", "path": None}
    pipeline_map = analyze_pipeline_path(Path(path)) if path else None
    recommendation = _recommend_for_summary(detected, task)
    payload = {
        "success": True,
        "entrypoint": "iobrpy-cli",
        "task": task,
        "detected_input": detected,
        "pipeline_map": pipeline_map,
        "recommendation": recommendation,
        "next_steps": [
            "Use iobrpy-cli map --path <path> --json first when you need to decide between continuing an existing analysis and rerunning it.",
            "Use iobrpy-cli commands --json to inspect the available native commands.",
            "Review recommendation.confirm_parameters before running a command so the JSON-defined confirmations are not skipped.",
            "Prefer iobrpy-cli <native_command> ... rather than analyze/quantify wrappers when possible.",
            "Use iobrpy-cli-mcp when your agent supports MCP tool registration.",
        ],
    }
    _print_structured(payload)


# ============================================================================
# REPL Mode
# ============================================================================

@cli.command('repl')
@click.option('--root', '-r', default='.', type=click.Path(), help='Root directory for projects')
@click.option('--sessions', type=click.Path(), help='Session directory')
def repl(root: str, sessions: str):
    """Start interactive REPL mode."""
    import time

    from prompt_toolkit import PromptSession
    from prompt_toolkit.history import FileHistory
    from prompt_toolkit.auto_suggest import AutoSuggestFromHistory

    click.echo(click.style("IOBRpy CLI Harness - REPL Mode", fg='blue', bold=True))
    click.echo("Type 'help' for available commands, 'exit' to quit.")
    click.echo()

    # Initialize managers
    proj_manager = ProjectManager(Path(root))
    session_dir = Path(sessions) if sessions else Path.home() / '.iobrpy' / 'sessions'
    session_manager = SessionManager(session_dir)
    current_session = session_manager.create_session()

    # Set up prompt session
    prompt_session = PromptSession(
        history=FileHistory(session_dir / 'history'),
        auto_suggest=AutoSuggestFromHistory(),
    )

    active_project = None

    while True:
        try:
            # Build prompt
            project_part = f"[{active_project}] " if active_project else ""
            prompt_text = f"{project_part}iobrpy> "

            user_input = prompt_session.prompt(prompt_text).strip()

            if not user_input:
                continue

            if user_input in ('exit', 'quit', 'q'):
                click.echo("Goodbye!")
                break

            if user_input == 'help':
                print_help()
                continue

            if user_input == 'status':
                if active_project:
                    proj = proj_manager.get_project(active_project)
                    click.echo(f"Project: {active_project}")
                    click.echo(f"  Path: {proj.project_dir}")
                    click.echo(f"  Input type: {proj.config.input_type or 'Not set'}")
                else:
                    click.echo("No active project. Use 'use <name>' to set one.")
                continue

            if user_input.startswith('use '):
                active_project = user_input[4:].strip()
                if active_project:
                    try:
                        proj_manager.get_project(active_project)
                        click.echo(f"Now using project: {active_project}")
                    except FileNotFoundError:
                        click.echo(f"Project '{active_project}' not found")
                        active_project = None
                continue

            if user_input == 'projects':
                projects = proj_manager.list_projects()
                for p in projects:
                    marker = "* " if p['name'] == active_project else "  "
                    click.echo(f"{marker}{p['name']}")
                continue

            if user_input == 'history':
                for i, entry in enumerate(current_session.get_history(limit=10), 1):
                    click.echo(f"  [{i}] {entry.command}")
                continue

            # Record and execute command
            start_time = time.time()
            exit_code = 0

            try:
                # Execute as shell command
                import subprocess
                result = subprocess.run(user_input, shell=True, cwd=Path.cwd())
                exit_code = result.returncode
            except Exception as e:
                click.echo(f"Error: {e}")
                exit_code = 1

            duration = int((time.time() - start_time) * 1000)
            current_session.record_command(user_input, str(Path.cwd()), exit_code, None, duration)

        except KeyboardInterrupt:
            click.echo("\nUse 'exit' to quit.")
        except EOFError:
            click.echo("\nGoodbye!")
            break


def print_help():
    """Print REPL help."""
    help_text = """
Available REPL commands:
  use <name>      Set active project
  projects         List all projects
  status           Show current status
  history          Show command history
  help             Show this help
  exit, quit, q    Exit REPL

Any IOBRpy command can be entered directly (e.g., 'analyze tme-profile input.csv').
"""
    click.echo(help_text)


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main entry point."""
    raise SystemExit(_dispatch_cli(sys.argv[1:]))


if __name__ == '__main__':
    main()

