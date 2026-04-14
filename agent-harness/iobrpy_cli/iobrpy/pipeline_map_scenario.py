"""
Scenario and recommendation helpers for pipeline_map.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set


def _next_nodes(next_stage_ids: List[str], stage_to_node: Dict[str, str]) -> List[str]:
    nodes: List[str] = []
    for stage_id in next_stage_ids:
        node = stage_to_node.get(stage_id, stage_id)
        if node not in nodes:
            nodes.append(node)
    return nodes


def _classify_scenario(
    completed: Set[str],
    detected: Set[str],
    current_stage: Optional[str],
) -> str:
    downstream = {"signature_scoring", "deconvolution", "ligand_receptor", "clustering"}
    immune = {"trust4", "spechla", "hla_typing"}
    upstream = {"fastq_qc", "salmon_quant", "salmon_merge", "prepare_salmon", "star_quant", "star_merge", "count2tpm"}

    has_raw = "raw_fastq" in detected
    has_upstream = bool(completed.intersection(upstream))
    has_tpm = "tpm_matrix" in completed
    has_downstream = bool(completed.intersection(downstream))
    has_immune = bool(completed.intersection(immune))

    if has_raw and (has_tpm or has_downstream or has_upstream):
        if has_tpm or has_downstream:
            return "mixed_raw_and_processed"
    if has_downstream:
        if {"deconvolution", "ligand_receptor"}.issubset(completed):
            return "downstream_complete"
        return "downstream_partial"
    if has_tpm:
        return "tpm_ready"
    if has_upstream:
        return "upstream_partial"
    if has_raw:
        return "raw_fastq_only"
    if has_immune or current_stage in immune:
        return "immune_only"
    return "unknown"


def _roadmap_progress_text(
    completed: Set[str],
    current_stage: Optional[str],
    next_stage_ids: List[str],
    *,
    ordered_stage_ids: callable,
    stage_to_node: Dict[str, str],
) -> Dict[str, str]:
    completed_nodes = []
    for stage_id in ordered_stage_ids(completed):
        node = stage_to_node.get(stage_id, stage_id)
        if node not in completed_nodes:
            completed_nodes.append(node)

    current_node = stage_to_node.get(current_stage, "Unclassified")
    next_nodes = _next_nodes(next_stage_ids, stage_to_node)

    completed_text = " -> ".join(completed_nodes) if completed_nodes else "No confirmed completed roadmap nodes yet"
    next_text = " | ".join(next_nodes) if next_nodes else "No clear next node identified yet"
    sentence = f"You are currently at {current_node}. Next possible nodes: {next_text}."
    return {
        "completed_text": completed_text,
        "next_text": next_text,
        "sentence": sentence,
    }


def _choice_details(
    choices: List[str],
    current_stage: Optional[str],
    suggested_commands: List[str],
    *,
    choice_labels: Dict[str, str],
    next_commands: Dict[str, str],
) -> List[Dict[str, Any]]:
    details: List[Dict[str, Any]] = []
    rerun_command: List[str] = []
    if current_stage and current_stage in next_commands:
        rerun_command = [next_commands[current_stage]]

    for choice in choices:
        item: Dict[str, Any] = {
            "id": choice,
            "label": choice_labels.get(choice, choice.replace("_", " ").title()),
        }
        if choice == "continue_downstream":
            item["when_to_choose"] = "Choose this when the directory already contains reusable intermediate or TPM-ready outputs."
            item["suggested_commands"] = suggested_commands
        elif choice == "rerun_current_stage":
            item["when_to_choose"] = "Choose this when the current results look incomplete, low quality, or were generated with the wrong parameters."
            item["suggested_commands"] = rerun_command or suggested_commands
        elif choice == "rerun_full_pipeline":
            item["when_to_choose"] = "Choose this when you want a clean end-to-end rerun, ideally in a new output directory."
            item["suggested_commands"] = [
                "iobrpy-cli runall --fastq <path_fastq_dir> --outdir <path_output_dir> --mode <salmon|star> --index <path_salmon_index|path_star_index> --threads <threads> --batch_size <batch_size> --project <project_name>"
            ]
        elif choice == "start_pipeline":
            item["when_to_choose"] = "Choose this when no reusable intermediate outputs were detected."
            item["suggested_commands"] = [
                "iobrpy-cli recommend --path <path> --task '<task>' --json"
            ]
        elif choice == "inspect_commands":
            item["when_to_choose"] = "Choose this when you need to inspect the command surface manually."
            item["suggested_commands"] = [
                "iobrpy-cli commands --json"
            ]
        details.append(item)
    return details


def _recommended_action(
    scenario: Dict[str, Any],
    choice_details: List[Dict[str, Any]],
    *,
    choice_labels: Dict[str, str],
) -> Dict[str, Any]:
    recommended_id = scenario.get("recommended_choice")
    for detail in choice_details:
        if detail.get("id") == recommended_id:
            return {
                "id": recommended_id,
                "label": detail.get("label"),
                "why": detail.get("when_to_choose"),
                "suggested_commands": detail.get("suggested_commands", []),
                "suggested_command_details": detail.get("suggested_command_details", []),
            }
    return {
        "id": recommended_id,
        "label": choice_labels.get(recommended_id, recommended_id or "unknown"),
        "why": "",
        "suggested_commands": [],
    }


def _scenario_payload(
    scenario_id: str,
    *,
    completed_labels: List[str],
    next_labels: List[str],
    choices: List[str],
    current_stage: Optional[str],
    suggested_commands: List[str],
    scenario_labels: Dict[str, str],
    choice_labels: Dict[str, str],
    next_commands: Dict[str, str],
) -> Dict[str, Any]:
    if scenario_id == "mixed_raw_and_processed":
        summary = "The directory contains both raw inputs and processed outputs. Reuse the processed outputs when possible, and use a fresh output directory if you decide to rerun end-to-end."
        communication_goals = [
            "Explain that raw FASTQ and processed outputs coexist in the same directory",
            "Summarize the current roadmap position and make it clear that existing results can often be reused",
            "Present the difference between continuing downstream, rerunning the current stage, and doing a clean full rerun",
        ]
        required_user_decisions = [
            "Whether to continue from existing results or rerun part or all of the workflow",
            "If rerunning end-to-end, which fresh output directory should be used",
            "If rerunning upstream quantification, which index path and mode should be used",
        ]
        recommended_choice = "continue_downstream"
    elif scenario_id == "tpm_ready":
        summary = "The directory is ready for downstream TME analysis from a TPM-like matrix."
        communication_goals = [
            "Explain that upstream quantification is already complete enough for downstream TME analysis",
            "Tell the user they are currently at the TPM Matrix node in the roadmap",
            "Offer a choice between continuing downstream analysis and rerunning upstream quantification",
        ]
        required_user_decisions = [
            "Whether to continue with downstream TME analysis or rerun upstream quantification",
            "If continuing downstream, whether to run all downstream analyses or only selected ones",
            "If rerunning upstream, which output directory and reference index should be used",
        ]
        recommended_choice = "continue_downstream"
    elif scenario_id == "downstream_partial":
        summary = "Some downstream TME analyses already exist, but the roadmap is not complete."
        communication_goals = [
            "List which downstream outputs already exist and which ones appear to be missing",
            "Make it clear that missing analyses can usually be filled in without rerunning everything",
            "Only suggest reruns after explaining what can be reused",
        ]
        required_user_decisions = [
            "Whether to fill in missing downstream analyses or rerun a downstream stage",
            "If rerunning, which stage should be rerun and where outputs should go",
        ]
        recommended_choice = "continue_downstream"
    elif scenario_id == "downstream_complete":
        summary = "Most downstream TME outputs already exist."
        communication_goals = [
            "Explain that most downstream roadmap nodes already have outputs",
            "Prioritize interpretation, summarization, export, clustering, or reporting before suggesting reruns",
            "Only discuss reruns if the user explicitly wants new parameters or new upstream processing",
        ]
        required_user_decisions = [
            "Whether to summarize existing results, add clustering or reporting, or rerun part of the workflow",
        ]
        recommended_choice = "continue_downstream"
    elif scenario_id == "upstream_partial":
        summary = "Intermediate upstream outputs were detected, so the directory looks resumable."
        communication_goals = [
            "Summarize which upstream stage appears to be complete and which one comes next",
            "Avoid defaulting to runall when reusable quantification outputs already exist",
            "Explain the difference between continuing from the current stage and rerunning it",
        ]
        required_user_decisions = [
            "Whether to continue from the current stage or rerun it",
            "If rerunning, which output directory should be used",
            "If continuing quantification, which mode and index path should be used if still missing",
        ]
        recommended_choice = "continue_downstream"
    elif scenario_id == "raw_fastq_only":
        summary = "Only raw FASTQ-like inputs were detected."
        communication_goals = [
            "Explain that only raw FASTQ input was detected and no reusable intermediate results were found",
            "Present runall as the default full-pipeline entrypoint",
            "Mention focused immune-side workflows only if they match the user task",
        ]
        required_user_decisions = [
            "Whether to start the full runall pipeline or a narrower workflow",
            "Which quantification mode, output directory, and index path should be used",
        ]
        recommended_choice = "start_pipeline"
    elif scenario_id == "immune_only":
        summary = "Immune-focused outputs were detected without enough evidence for the full TME roadmap."
        communication_goals = [
            "Clarify that immune-focused outputs were found",
            "Explain that the full TME roadmap may still require a TPM-like expression matrix",
            "Offer to continue the immune workflow or switch to a TPM/TME workflow if expression data exists",
        ]
        required_user_decisions = [
            "Whether to continue the immune-focused workflow or switch to a TPM/TME workflow",
            "If switching to TPM/TME, which matrix or expression directory should be used",
        ]
        recommended_choice = "continue_downstream"
    else:
        summary = "The directory could not be mapped confidently to the IOBRpy roadmap."
        communication_goals = [
            "Be explicit about uncertainty instead of guessing the stage",
            "Show the roadmap link and explain that a more specific entrypoint path would help",
            "Suggest inspecting commands or pointing to the key FASTQ or matrix file",
        ]
        required_user_decisions = [
            "Whether to inspect commands manually or provide a more specific matrix or FASTQ path",
        ]
        recommended_choice = "inspect_commands"

    return {
        "id": scenario_id,
        "label": scenario_labels.get(scenario_id, scenario_id.replace("_", " ").title()),
        "summary": summary,
        "completed_stage_labels": completed_labels,
        "next_stage_labels": next_labels,
        "communication_goals": communication_goals,
        "required_user_decisions": required_user_decisions,
        "recommended_choice": recommended_choice,
        "choice_details": _choice_details(
            choices,
            current_stage,
            suggested_commands,
            choice_labels=choice_labels,
            next_commands=next_commands,
        ),
    }
