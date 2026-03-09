#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IOBRpy BYOK assistant server.

Refactor boundaries:
1) Planner (LLM): understands user intent and proposes structured next action (JSON plan only).
2) Guardrails (Python): validates planner output, enforces schema/rules, updates state safely, builds/executes command.
3) UI (ai.py): only renders responses and handles terminal input/output.
"""

import json
import os
import re
import subprocess
import sys
import traceback
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import requests

REQUIRED_FILE = os.getenv(
    "IOBRPY_REQUIRED_PARAMS_FILE",
    os.path.join(os.path.dirname(__file__), "iobrpy_required_params.json"),
)
DEFAULTS_FILE = os.getenv(
    "IOBRPY_DEFAULTS_FILE",
    os.path.join(os.path.dirname(__file__), "iobrpy_defaults.json"),
)
CONDA_EXE = os.getenv("CONDA_EXE", "conda")
IOBRPY_CONDA_ENV = os.getenv("IOBRPY_CONDA_ENV", "iobrpy")
IOBRPY_SRC_DIR = os.getenv("IOBRPY_SRC_DIR", "")

SUPPORTED_PROTOCOLS = ["2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"]
SERVER_INFO = {
    "name": "iobrpy-ai-agentic",
    "title": "iobrpy AI Agentic (BYOK cloud LLM planner)",
    "version": "0.7.0",
    "description": "Planner-first natural language assistant for iobrpy.",
}
CAPABILITIES = {"tools": {"listChanged": False}}

VALID_ACTIONS = {
    "reply_only",
    "select_command",
    "clarify_intent",
    "update_params",
    "remove_params",
    "switch_command",
    "confirm_switch",
    "ask_missing_params",
    "ready_to_run",
    "execute",
    "undo",
    "restart_session",
}


@dataclass
class AIConfig:
    llm_alias: str
    api_key: str
    model: str
    base_url: str
    logdir: str


_AI_CONFIG: Optional[AIConfig] = None


def configure_runtime(config: AIConfig) -> None:
    global _AI_CONFIG
    _AI_CONFIG = config


def _require_ai_config() -> AIConfig:
    if _AI_CONFIG is None:
        raise RuntimeError("AI runtime config is not initialized. Call configure_runtime() first.")
    return _AI_CONFIG


def _runtime_identity() -> Dict[str, str]:
    cfg = _AI_CONFIG
    if cfg is None:
        return {"llm_alias": "unknown", "model": "unknown", "base_url": "unknown"}
    return {"llm_alias": cfg.llm_alias, "model": cfg.model, "base_url": cfg.base_url}


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def _debug_enabled() -> bool:
    return os.getenv("IOBRPY_AI_DEBUG", "0").strip().lower() in {"1", "true", "yes", "on"}


def dlog(*a):
    if _debug_enabled():
        log("[ai-debug]", *a)


def send(obj: Dict[str, Any]):
    sys.stdout.write(json.dumps(obj, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def jsonrpc_error(_id, code, message, data=None):
    err = {"code": code, "message": message}
    if data is not None:
        err["data"] = data
    send({"jsonrpc": "2.0", "id": _id, "error": err})


def _read_json(path: str, default):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return default


def _write_json(path: str, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def _normalize_text(s: str) -> str:
    if not s:
        return ""
    trans = {
        "；": ";", "，": ",", "。": ".", "：": ":", "＝": "=", "（": "(", "）": ")",
        "【": "[", "】": "]", "“": '"', "”": '"', "‘": "'", "’": "'",
    }
    for k, v in trans.items():
        s = s.replace(k, v)
    s = s.replace("\n", " ").replace("\r", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _contains_cjk(text: str) -> bool:
    return re.search(r"[\u4e00-\u9fff]", text or "") is not None


def _extract_json_block(text: str) -> Dict[str, Any]:
    if not text:
        return {}
    text = text.strip()
    try:
        return json.loads(text)
    except Exception:
        pass
    m = re.search(r"\{.*\}", text, flags=re.S)
    if not m:
        return {}
    try:
        return json.loads(m.group(0))
    except Exception:
        return {}


def _post_openai_compatible(base_url: str, api_key: str, model: str, prompt: str) -> str:
    url = base_url.rstrip("/") + "/chat/completions"
    payload = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0,
        "response_format": {"type": "json_object"},
    }
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    r = requests.post(url, headers=headers, json=payload, timeout=120)
    r.raise_for_status()
    data = r.json()
    return (((data.get("choices") or [{}])[0].get("message") or {}).get("content") or "").strip()


def _post_anthropic(base_url: str, api_key: str, model: str, prompt: str) -> str:
    url = base_url.rstrip("/") + "/messages"
    payload = {
        "model": model,
        "max_tokens": 1200,
        "temperature": 0,
        "messages": [{"role": "user", "content": prompt}],
    }
    headers = {
        "x-api-key": api_key,
        "anthropic-version": "2023-06-01",
        "content-type": "application/json",
    }
    r = requests.post(url, headers=headers, json=payload, timeout=120)
    r.raise_for_status()
    data = r.json()
    out = ""
    for p in data.get("content") or []:
        if isinstance(p, dict) and p.get("type") == "text":
            out += p.get("text") or ""
    return out.strip()


def _post_gemini(base_url: str, api_key: str, model: str, prompt: str) -> str:
    url = base_url.rstrip("/") + f"/models/{model}:generateContent?key={api_key}"
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": 0, "responseMimeType": "application/json"},
    }
    r = requests.post(url, headers={"Content-Type": "application/json"}, json=payload, timeout=120)
    r.raise_for_status()
    data = r.json()
    cands = data.get("candidates") or []
    if not cands:
        return ""
    parts = (((cands[0] or {}).get("content") or {}).get("parts") or [])
    return "\n".join([p.get("text", "") for p in parts if isinstance(p, dict)]).strip()


def _llm_generate_json(prompt: str) -> Dict[str, Any]:
    cfg = _require_ai_config()
    if cfg.llm_alias == "claude":
        txt = _post_anthropic(cfg.base_url, cfg.api_key, cfg.model, prompt)
    elif cfg.llm_alias == "gemini":
        txt = _post_gemini(cfg.base_url, cfg.api_key, cfg.model, prompt)
    else:
        txt = _post_openai_compatible(cfg.base_url, cfg.api_key, cfg.model, prompt)
    return _extract_json_block(txt)


def _translate_to_english(text: str) -> str:
    if not text or not _contains_cjk(text):
        return text or ""
    prompt = f"""
Translate user text to English faithfully.
Keep command names, parameter names, flags, file paths, and numbers unchanged.
Return JSON only: {{"english": "..."}}
Text: {text}
""".strip()
    try:
        j = _llm_generate_json(prompt)
        out = j.get("english") if isinstance(j, dict) else None
        if isinstance(out, str) and out.strip():
            return out.strip()
    except Exception:
        pass
    return text


def _translate_to_chinese(text: str) -> str:
    if not text:
        return ""
    prompt = f"""
Translate text to Chinese faithfully.
Keep command names, parameter names, flags, file paths, and numbers unchanged.
Return JSON only: {{"chinese": "..."}}
Text: {text}
""".strip()
    try:
        j = _llm_generate_json(prompt)
        out = j.get("chinese") if isinstance(j, dict) else None
        if isinstance(out, str) and out.strip():
            return out.strip()
    except Exception:
        pass
    return text


def enforce_output_language(text: str, prefer_chinese: bool) -> str:
    """Best-effort output language guardrail for assistant-facing text."""
    safe_en = "Please describe your task."
    if not text:
        return ""
    if prefer_chinese:
        return text
    if not _contains_cjk(text):
        return text

    translated = _translate_to_english(text)
    if translated and not _contains_cjk(translated):
        return translated
    return safe_en


DEFAULT_INTENT_KEYWORDS: Dict[str, List[str]] = {
    "runall": ["fastq", "workflow", "pipeline", "bulk", "rna-seq"],
    "trust4": ["tcr", "bcr", "vdj", "repertoire", "clonotype"],
    "spechla": ["hla", "typing"],
    "extract_hla_read": ["hla", "extract hla read"],
    "hla_typing": ["hla", "typing"],
}
DEFAULT_RULES = {
    "runall": {
        "required": ["fastq", "outdir", "mode", "index", "threads", "batch_size", "project"],
        "confirm": ["mode", "threads", "batch_size"],
        "choices": {"mode": ["salmon", "star"]},
    }
}


def load_rules() -> Dict[str, Any]:
    rules = _read_json(REQUIRED_FILE, DEFAULT_RULES)
    if not isinstance(rules, dict) or not rules:
        return DEFAULT_RULES
    out: Dict[str, Any] = {}
    for k, v in rules.items():
        if not isinstance(v, dict):
            continue
        req = v.get("required", []) if isinstance(v.get("required", []), list) else []
        conf = v.get("confirm", []) if isinstance(v.get("confirm", []), list) else []
        opt = v.get("optional", []) if isinstance(v.get("optional", []), list) else []
        choices = v.get("choices", {}) if isinstance(v.get("choices", {}), dict) else {}
        req_one = v.get("required_one_of", []) if isinstance(v.get("required_one_of", []), list) else []
        notes = v.get("notes", {}) if isinstance(v.get("notes", {}), dict) else {}
        optional_defaults = v.get("optional_defaults", {}) if isinstance(v.get("optional_defaults", {}), dict) else {}
        cli_flags = v.get("cli_flags", {}) if isinstance(v.get("cli_flags", {}), dict) else {}
        param_aliases = v.get("param_aliases", {}) if isinstance(v.get("param_aliases", {}), dict) else {}
        param_types = v.get("param_types", {}) if isinstance(v.get("param_types", {}), dict) else {}
        intent_keywords = v.get("intent_keywords", [])
        param_hints = v.get("param_hints", {}) if isinstance(v.get("param_hints", {}), dict) else {}
        function_summary = v.get("function_summary", "")
        if not isinstance(intent_keywords, (list, dict)):
            intent_keywords = []
        if (not intent_keywords) and k in DEFAULT_INTENT_KEYWORDS:
            intent_keywords = list(DEFAULT_INTENT_KEYWORDS[k])
        out[k] = {
            "required": req,
            "confirm": conf,
            "optional": opt,
            "choices": choices,
            "required_one_of": req_one,
            "notes": notes,
            "optional_defaults": optional_defaults,
            "cli_flags": cli_flags,
            "param_aliases": param_aliases,
            "param_types": param_types,
            "intent_keywords": intent_keywords,
            "param_hints": param_hints,
            "function_summary": function_summary,
        }
    return out or DEFAULT_RULES


def allowed_keys_for(subcommand: str, rules: Dict[str, Any]) -> List[str]:
    r = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    keys: List[str] = []
    for k in r.get("required", []):
        if k not in keys:
            keys.append(k)
    for k in r.get("optional", []):
        if k not in keys:
            keys.append(k)
    if "confirm" not in keys:
        keys.append("confirm")
    return keys


BUILTIN_OPTIONAL_DEFAULTS = {"threads": 8, "batch_size": 2, "project": "iobrpy_ai", "mode": "salmon"}


def load_defaults() -> Dict[str, Any]:
    d = _read_json(DEFAULTS_FILE, {})
    return d if isinstance(d, dict) else {}


def set_defaults(new_defaults: Dict[str, Any], rules: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(new_defaults, dict):
        return {"ok": False, "error": "defaults must be object"}
    allowed = set()
    for cmd in rules:
        for k in allowed_keys_for(cmd, rules):
            if k != "confirm":
                allowed.add(k)
    cur = load_defaults()
    for k, v in new_defaults.items():
        if k in allowed:
            cur[k] = v
    _write_json(DEFAULTS_FILE, cur)
    return {"ok": True, "defaults": cur}


def _shell_quote(s: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./=-]+", s):
        return s
    return "'" + s.replace("'", "'\"'\"'") + "'"


FLAG_MAP = {
    "fastq": "--fastq", "outdir": "--outdir", "mode": "--mode", "index": "--index", "threads": "--threads",
    "batch_size": "--batch_size", "project": "--project", "input": "--input", "output": "--output", "bam": "-b",
    "r1": "-1", "r2": "-2", "ru": "-u", "fqdir": "--fqdir", "f": "-f", "ref": "--ref", "o": "-o",
    "od": "--od", "t": "-t", "k": "-k", "bam_dir": "--bam-dir", "use_exon": "--use-exon",
    "sample": "--sample", "name": "--name", "read1": "--read1", "read2": "--read2",
}

PARAM_FLAGS_BY_COMMAND: Dict[str, Dict[str, str]] = {}

BOOL_FLAGS = {
    "repseq": "--repseq", "skipMateExtension": "--skipMateExtension", "abnormalUnmapFlag": "--abnormalUnmapFlag",
    "assembleWithRef": "--assembleWithRef", "noExtraction": "--noExtraction", "outputReadAssignment": "--outputReadAssignment",
    "no_auto_install": "--no-auto-install",
}

BOOL_FLAGS_BY_COMMAND: Dict[str, Dict[str, str]] = {}



CLI_FLAGS_BY_COMMAND: Dict[str, set] = {
    "hla_typing": {"-b", "--bam-dir", "-r", "--ref", "-o", "--outdir", "-j", "--threads", "-u", "--use-exon"},
    "spechla": {"-n", "--name", "-1", "--read1", "-2", "--read2", "-o", "--outdir", "-j", "--threads", "-u", "--use-exon"},
    "extract_hla_read": {"-s", "--sample", "-b", "--bam", "-r", "--ref", "-o", "--outdir", "--no-auto-install"},
    "runall": {"--mode", "--outdir", "--fastq", "--resume", "--dry_run", "--index", "--threads", "--batch_size", "--project"},
    "tme_profile": {"-i", "--input", "-o", "--output", "--threads"},
    "trust4": {
        "-b", "-1", "-2", "-u", "--fqdir", "-o", "--od", "-t", "-f", "-k",
        "--repseq", "--skipMateExtension", "--abnormalUnmapFlag", "--assembleWithRef", "--noExtraction", "--outputReadAssignment",
    },
}


def _default_cli_flag_for_key(key: str) -> Optional[str]:
    if not isinstance(key, str) or not key:
        return None
    if re.fullmatch(r"[A-Za-z0-9_]+", key) is None:
        return None
    return "--" + key


def _resolve_cli_flag(subcommand: str, key: str, rules: Dict[str, Any], is_bool: bool = False) -> Optional[str]:
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    cli_flags = rule.get("cli_flags", {}) if isinstance(rule.get("cli_flags", {}), dict) else {}
    if key in cli_flags and isinstance(cli_flags.get(key), str) and cli_flags.get(key):
        return str(cli_flags[key])

    per_command = (BOOL_FLAGS_BY_COMMAND if is_bool else PARAM_FLAGS_BY_COMMAND).get(subcommand, {})
    if key in per_command:
        return per_command[key]
    if is_bool:
        return BOOL_FLAGS.get(key) or _default_cli_flag_for_key(key)
    return FLAG_MAP.get(key) or _default_cli_flag_for_key(key)


def _unsupported_serialized_flags(subcommand: str, argv: List[str]) -> List[str]:
    allowed = CLI_FLAGS_BY_COMMAND.get(subcommand)
    if not allowed:
        return []
    used = [a for a in argv[1:] if isinstance(a, str) and a.startswith("-")]
    return sorted({f for f in used if f not in allowed})


def audit_parameter_link_consistency(rules: Dict[str, Any]) -> Dict[str, Any]:
    report: Dict[str, Any] = {}
    for cmd, rule in rules.items():
        if not isinstance(rule, dict):
            continue
        required = rule.get("required", []) if isinstance(rule.get("required", []), list) else []
        optional = rule.get("optional", []) if isinstance(rule.get("optional", []), list) else []
        required_unserializable: List[str] = []
        optional_unserializable: List[str] = []
        for k in required:
            if _resolve_cli_flag(cmd, k, rules, False) is None:
                required_unserializable.append(k)
        for k in optional:
            if _resolve_cli_flag(cmd, k, rules, False) is None and _resolve_cli_flag(cmd, k, rules, True) is None:
                optional_unserializable.append(k)
        report[cmd] = {
            "required_unserializable": required_unserializable,
            "optional_unserializable": optional_unserializable,
        }
    return report


def _required_without_serializable_flag(subcommand: str, rules: Dict[str, Any]) -> List[str]:
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    required = rule.get("required", []) if isinstance(rule.get("required", []), list) else []
    return [k for k in required if _resolve_cli_flag(subcommand, k, rules, False) is None and _resolve_cli_flag(subcommand, k, rules, True) is None]


def build_command(subcommand: str, params: Dict[str, Any], rules: Dict[str, Any]) -> Tuple[str, List[str]]:
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    argv = [subcommand]
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    unserializable_required = _required_without_serializable_flag(subcommand, rules)
    if unserializable_required:
        raise ValueError(
            f"Required params missing CLI flag mapping for '{subcommand}': {', '.join(unserializable_required)}"
        )
    ordered = list(rule.get("required", [])) + [k for k in rule.get("optional", []) if k not in rule.get("required", [])]
    for k in ordered:
        if k not in allowed_set or k not in params:
            continue
        v = params[k]
        if v in (None, "", []):
            continue
        if isinstance(v, bool):
            bf = _resolve_cli_flag(subcommand, k, rules, is_bool=True)
            if bf and v:
                argv.append(bf)
            continue
        flag = _resolve_cli_flag(subcommand, k, rules, is_bool=False)
        if not flag:
            raise ValueError(f"Parameter '{k}' for '{subcommand}' cannot be serialized to a CLI flag")
        argv.extend([flag, str(v)])

    unsupported_flags = _unsupported_serialized_flags(subcommand, argv)
    if unsupported_flags:
        raise ValueError(
            f"Serialized flags not supported by CLI for '{subcommand}': {', '.join(unsupported_flags)}"
        )
    return "iobrpy " + " ".join([_shell_quote(a) for a in argv]), argv


def merge_defaults(params: Dict[str, Any], subcommand: str, rules: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, str]]:
    d = load_defaults()
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    rule_defaults = rule.get("optional_defaults", {}) if isinstance(rule.get("optional_defaults", {}), dict) else {}
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    out: Dict[str, Any] = {}
    src: Dict[str, str] = {}

    for k in allowed_set:
        if k in params and params[k] in (None, "", []):
            continue
        if k in rule_defaults and rule_defaults[k] not in (None, "", []):
            out[k] = rule_defaults[k]
            src[k] = "default"
        elif k in BUILTIN_OPTIONAL_DEFAULTS:
            out[k] = BUILTIN_OPTIONAL_DEFAULTS[k]
            src[k] = "default"
        if k in d and d[k] not in (None, "", []):
            out[k] = d[k]
            src[k] = "default"

    for k, v in params.items():
        if k in allowed_set and v not in (None, "", []):
            out[k] = v
            src[k] = "user"
        elif k in allowed_set and v in (None, "", []) and k in out:
            out.pop(k, None)
            src.pop(k, None)

    return out, src


def compute_needs(subcommand: str, rules: Dict[str, Any], params: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    required = rule.get("required", []) if isinstance(rule.get("required", []), list) else []
    confirm = rule.get("confirm", []) if isinstance(rule.get("confirm", []), list) else []
    required_one_of = rule.get("required_one_of", []) if isinstance(rule.get("required_one_of", []), list) else []

    missing = [k for k in required if k not in params or params[k] in (None, "", [])]

    if required_one_of:
        def _group_ok(group):
            return isinstance(group, list) and group and all((k in params and params[k] not in (None, "", [])) for k in group)
        if not any(_group_ok(g) for g in required_one_of):
            missing.append("__one_of_group")

    # keep confirm as soft recommendation, not hard gate
    need_confirm = [k for k in confirm if k in params and params.get(k) not in (None, "", [])]
    return missing, need_confirm


def build_command_catalog(rules: Dict[str, Any]) -> Dict[str, Any]:
    tools: List[Dict[str, Any]] = []
    for cmd, r in rules.items():
        if not isinstance(r, dict):
            continue
        summary = r.get("function_summary", "")
        if isinstance(summary, dict):
            summary_txt = summary.get("zh") or summary.get("en") or ""
        else:
            summary_txt = summary or ""
        tools.append(
            {
                "name": cmd,
                "summary": summary_txt,
                "required": r.get("required", []),
                "optional": r.get("optional", []),
                "choices": r.get("choices", {}),
                "required_one_of": r.get("required_one_of", []),
                "optional_defaults": r.get("optional_defaults", {}),
                "intent_keywords": r.get("intent_keywords", []),
                "param_hints": r.get("param_hints", {}),
                "notes": r.get("notes", {}),
            }
        )
    return {"commands": tools}


PLANNER_SCHEMA = {
    "action": "one of: reply_only/select_command/clarify_intent/update_params/remove_params/switch_command/confirm_switch/ask_missing_params/ready_to_run/execute/undo/restart_session",
    "intent_type": "one of: side_question/execution_request/parameter_update/switch_request/clarification/control/other",
    "is_side_question": "boolean",
    "keep_current_command": "boolean",
    "switch_intent_strength": "0..1",
    "needs_confirmation": "boolean",
    "message": "natural-language assistant reply",
    "subcommand": "command name or null",
    "param_updates": "object of key->value",
    "param_removals": "array of parameter keys",
    "switch_to": "command name or null",
    "ask_for": "array of missing parameter keys",
    "confirm": "null or object with type/details",
    "confidence": "0..1",
    "reason": "short internal reason",
}


def build_planner_prompt(state: "SessionState", user_text: str, catalog: Dict[str, Any], merged_params: Dict[str, Any], missing: List[str], draft_command: Optional[str]) -> str:
    language_policy = (
        "Respond in Chinese."
        if state.prefer_chinese
        else "Respond in English only. Do not output Chinese unless the user has already used Chinese in this session."
    )
    return f"""
You are the iobrpy planner assistant (not a generic chatbot).
Your job: be a normal helpful LLM for conversation, while keeping iobrpy execution planning safe and catalog-bound.

Language policy (highest priority):
- {language_policy}

Decision principles:
- Understand user meaning semantically from context and conversation history.
- Keep current selected_command by default when user is asking explanation/comparison/API/method/parameter/workflow questions.
- Do not switch command merely because another tool/method/command is mentioned.
- Use switch_command or confirm_switch only when replacement intent is explicit or strongly implied.
- For side questions, prefer action=reply_only with keep_current_command=true and a direct helpful answer.
- If user asks identity/model/provider/runtime questions, answer truthfully using Runtime identity context.
- Distinguish role (IOBRpy assistant) from runtime model/provider.

Execution guardrails:
- Only use command names and parameter keys that exist in catalog.
- Never invent command or parameter.
- If intent is ambiguous: use action=clarify_intent.
- If user modifies/removes parameters for selected command: use update_params/remove_params.
- If required inputs are missing: do NOT execute.
- If user says run/execute but missing inputs exist: action=ask_missing_params.
- If everything required is satisfied and user asks to run: action=execute.
- Output JSON only, exactly following schema fields.

Planner JSON schema:
{json.dumps(PLANNER_SCHEMA, ensure_ascii=False)}

Current session state:
{json.dumps(state.to_planner_state(), ensure_ascii=False)}

Runtime identity:
{json.dumps(_runtime_identity(), ensure_ascii=False)}

Catalog:
{json.dumps(catalog, ensure_ascii=False)}

Current merged params:
{json.dumps(merged_params, ensure_ascii=False)}
Current missing params:
{json.dumps(missing, ensure_ascii=False)}
Current draft command:
{json.dumps(draft_command, ensure_ascii=False)}

User message:
{json.dumps(user_text, ensure_ascii=False)}

Return JSON only.
""".strip()




def build_selector_catalog(rules: Dict[str, Any]) -> Dict[str, Any]:
    commands: List[Dict[str, Any]] = []
    for cmd, r in rules.items():
        if not isinstance(r, dict):
            continue
        summary = r.get("function_summary", "")
        if isinstance(summary, dict):
            summary_txt = summary.get("zh") or summary.get("en") or ""
        else:
            summary_txt = summary or ""
        commands.append(
            {
                "name": cmd,
                "summary": summary_txt,
                "required": r.get("required", []),
                "optional": r.get("optional", []),
                "cli_flags": r.get("cli_flags", {}),
            }
        )
    return {"commands": commands}


def llm_select_command_from_catalog(
    user_text: str,
    rules: Dict[str, Any],
    current_command: Optional[str] = None,
    prefer_chinese: bool = False,
) -> Dict[str, Any]:
    catalog = build_selector_catalog(rules)
    prompt = f"""
You are an iobrpy command selector, not a chatbot.
Task: choose the single best matching command from catalog.commands for the user request.
If uncertain, return subcommand as null.
Never invent command names.
Only output JSON in this exact shape:
{{"subcommand": null, "confidence": 0.0, "reason": "..."}}

current_command:
{json.dumps(current_command, ensure_ascii=False)}

catalog:
{json.dumps(catalog, ensure_ascii=False)}

user_text:
{json.dumps(user_text, ensure_ascii=False)}
""".strip()

    try:
        out = _llm_generate_json(prompt)
    except Exception:
        return {"subcommand": None, "confidence": 0.0, "reason": "selector_call_failed"}

    if not isinstance(out, dict):
        return {"subcommand": None, "confidence": 0.0, "reason": "selector_invalid_output"}

    sub = out.get("subcommand")
    conf_raw = out.get("confidence", 0.0)
    try:
        conf = max(0.0, min(1.0, float(conf_raw)))
    except Exception:
        conf = 0.0

    if sub not in rules:
        sub = None
        conf = 0.0

    reason = str(out.get("reason") or "")
    return {"subcommand": sub, "confidence": conf, "reason": reason}

def _safe_fallback_plan(state: "SessionState", user_text: str) -> Dict[str, Any]:
    default_msg = "Please continue." if not state.prefer_chinese else "请继续。"
    return {
        "action": "reply_only",
        "intent_type": "other",
        "is_side_question": True,
        "keep_current_command": True,
        "switch_intent_strength": 0.0,
        "needs_confirmation": False,
        "message": default_msg,
        "subcommand": None,
        "switch_to": None,
        "param_updates": {},
        "param_removals": [],
        "ask_for": [],
        "confirm": None,
        "confidence": 0.0,
        "reason": "planner_fallback",
    }


def llm_plan_next_action(state: "SessionState", user_text: str, rules: Dict[str, Any]) -> Dict[str, Any]:
    catalog = build_command_catalog(rules)
    sub = state.selected_command or ""
    merged, _ = merge_defaults(state.params, sub, rules) if sub else ({}, {})
    missing, _ = compute_needs(sub, rules, merged) if sub else ([], [])
    try:
        draft = build_command(sub, merged, rules)[0] if sub else None
    except ValueError:
        draft = None
    prompt = build_planner_prompt(state, user_text, catalog, merged, missing, draft)

    try:
        plan = _llm_generate_json(prompt)
        if not isinstance(plan, dict) or not plan:
            plan = _safe_fallback_plan(state, user_text)
    except Exception as e:
        dlog("planner_call_failed", repr(e))
        plan = _safe_fallback_plan(state, user_text)

    dlog("planner_raw_plan", plan)
    return plan


def _normalize_identifier(text: str) -> str:
    t = (text or "").strip().lower()
    t = re.sub(r"[^a-z0-9_-]", "", t)
    return t


def _param_type_for_key(subcommand: str, key: str, rules: Dict[str, Any]) -> str:
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    param_types = rule.get("param_types", {}) if isinstance(rule.get("param_types", {}), dict) else {}
    if isinstance(param_types.get(key), str) and param_types.get(key):
        return str(param_types[key]).strip().lower()
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}
    if key in choices:
        return "choice"
    if key in {"threads", "batch_size", "t", "k"}:
        return "int"
    if key.startswith("no_"):
        return "bool"
    if key in {"index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od", "bam_dir", "read1", "read2"}:
        return "path"
    return "string"


def _is_path_like_value(val: str) -> bool:
    if not isinstance(val, str):
        return False
    v = val.strip()
    return bool(v) and (
        v.startswith("/") or
        v.startswith("./") or
        v.startswith("../") or
        bool(re.search(r"[A-Za-z0-9._-]+/[A-Za-z0-9._/-]+", v))
    )


def _truthy_bool_text(val: str) -> Optional[bool]:
    if not isinstance(val, str):
        return None
    t = val.strip().lower()
    if t in {"1", "true", "yes", "y", "on"}:
        return True
    if t in {"0", "false", "no", "n", "off"}:
        return False
    return None


def _build_alias_map(subcommand: str, rules: Dict[str, Any]) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    allowed = [k for k in allowed_keys_for(subcommand, rules) if k != "confirm"]
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    cfg = rule.get("param_aliases", {}) if isinstance(rule.get("param_aliases", {}), dict) else {}
    for key in allowed:
        aliases[_normalize_identifier(key)] = key
        aliases[_normalize_identifier(key.replace("_", "-"))] = key
        for a in cfg.get(key, []) if isinstance(cfg.get(key, []), list) else []:
            if isinstance(a, str) and a.strip():
                aliases[_normalize_identifier(a)] = key
    return aliases


def _split_structured_segments(text: str) -> List[str]:
    s = _normalize_text(text or "")
    if not s:
        return []
    return [x.strip() for x in re.split(r"[,;，；、]+", s) if x.strip()]


def _infer_unkeyed_param(candidates: List[str], value: str, subcommand: str, rules: Dict[str, Any], missing: List[str]) -> Optional[str]:
    if not candidates:
        return None
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}
    v = value.strip()

    choice_keys = [k for k in candidates if isinstance(choices.get(k), list) and v in choices.get(k, [])]
    if len(choice_keys) == 1:
        return choice_keys[0]

    if re.fullmatch(r"\d+", v):
        int_keys = [k for k in candidates if _param_type_for_key(subcommand, k, rules) == "int"]
        if len(int_keys) == 1:
            return int_keys[0]

    if _is_path_like_value(v):
        path_keys = [k for k in candidates if _param_type_for_key(subcommand, k, rules) == "path"]
        if len(path_keys) == 1:
            return path_keys[0]

    if len(missing) == 1 and missing[0] in candidates:
        return missing[0]
    return None


def extract_param_updates_from_text(
    user_text: str,
    subcommand: Optional[str],
    rules: Dict[str, Any],
    current_params: Dict[str, Any],
    missing_params: List[str],
) -> Dict[str, Any]:
    if not subcommand or subcommand not in rules:
        return {}

    allowed = [k for k in allowed_keys_for(subcommand, rules) if k != "confirm"]
    if not allowed:
        return {}

    alias_map = _build_alias_map(subcommand, rules)
    segments = _split_structured_segments(user_text)
    if not segments:
        return {}

    updates: Dict[str, Any] = {}
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}

    for seg in segments:
        m = re.match(r"^\s*([A-Za-z0-9_\-]+)\s*(=|:)\s*(.+?)\s*$", seg)
        key = None
        raw_val = None
        if m:
            key = alias_map.get(_normalize_identifier(m.group(1)))
            raw_val = m.group(3).strip()
        else:
            m2 = re.match(r"^\s*([A-Za-z0-9_\-]+)\s+(.+?)\s*$", seg)
            if m2:
                key = alias_map.get(_normalize_identifier(m2.group(1)))
                raw_val = m2.group(2).strip()

        if key and raw_val is not None:
            val = raw_val.strip().strip('"').strip("'")
            try:
                updates[key] = _convert_value_for_key(key, val, rules, subcommand)
            except Exception:
                pass
            continue

        # unkeyed inference with strong schema constraints
        tokens = [t.strip("\"'") for t in re.findall(r"[A-Za-z0-9_./\-]+", seg)]
        unresolved = [k for k in allowed if k not in updates]
        used_tokens = set()

        alias_hits = [alias_map.get(_normalize_identifier(tok)) for tok in tokens]
        alias_keys = [k for k in alias_hits if k]
        alias_keys = list(dict.fromkeys(alias_keys))
        if len(alias_keys) == 1:
            alias_key = alias_keys[0]
            expected_type = _param_type_for_key(subcommand, alias_key, rules)
            value_tokens = [t for t in tokens if alias_map.get(_normalize_identifier(t)) != alias_key]
            if not value_tokens:
                value_tokens = [t for t in tokens if not alias_map.get(_normalize_identifier(t))]
            for candidate in value_tokens:
                if expected_type == "path" and not _is_path_like_value(candidate):
                    continue
                if expected_type == "int" and not re.fullmatch(r"\d+", candidate):
                    continue
                if expected_type == "choice":
                    opts = choices.get(alias_key, []) if isinstance(choices.get(alias_key, []), list) else []
                    if candidate not in opts:
                        continue
                try:
                    updates[alias_key] = _convert_value_for_key(alias_key, candidate, rules, subcommand)
                    used_tokens.add(candidate)
                    unresolved = [k for k in unresolved if k != alias_key]
                    break
                except Exception:
                    continue

        ambiguous_choice_keys = set()
        for k in unresolved:
            opts = choices.get(k, []) if isinstance(choices.get(k, []), list) else []
            if not opts:
                continue
            matched = {tok for tok in tokens if tok in opts}
            if len(matched) > 1:
                ambiguous_choice_keys.add(k)

        for tok in tokens:
            if not tok or tok in used_tokens:
                continue
            inferred_key = _infer_unkeyed_param(
                [k for k in unresolved if k not in ambiguous_choice_keys],
                tok,
                subcommand,
                rules,
                missing_params,
            )
            if not inferred_key:
                continue
            try:
                updates[inferred_key] = _convert_value_for_key(inferred_key, tok, rules, subcommand)
                used_tokens.add(tok)
                unresolved = [k for k in unresolved if k != inferred_key]
            except Exception:
                continue

    return updates


def _text_has_parameter_shape(user_text: str, subcommand: str, rules: Dict[str, Any]) -> bool:
    text = _normalize_text(user_text or "")
    if not text:
        return False

    if re.search(r"\s[\w\-]+\s*[=:]\s*\S+", text):
        return True
    if re.search(r"(?:^|\s)/[^\s,;]+", text):
        return True
    if re.search(r"\b\d+\b", text):
        return True

    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}
    for vals in choices.values():
        if not isinstance(vals, list):
            continue
        for v in vals:
            if isinstance(v, str) and re.search(rf"\b{re.escape(v)}\b", text):
                return True
    return False


def _can_commit_command_from_plan(validated: Dict[str, Any], rules: Dict[str, Any]) -> bool:
    action = str(validated.get("action") or "")
    sub = validated.get("subcommand")
    switch_to = validated.get("switch_to")

    if action == "select_command" and sub in rules:
        return True
    if action in {"switch_command", "confirm_switch"} and (switch_to in rules or sub in rules):
        return True
    return False


def normalize_command_commit(validated: Dict[str, Any], state: "SessionState", rules: Dict[str, Any]) -> Dict[str, Any]:
    if state.selected_command is not None:
        return validated
    if _can_commit_command_from_plan(validated, rules):
        return validated

    sub = validated.get("subcommand")
    if sub not in rules:
        return validated

    action = str(validated.get("action") or "")
    intent_type = str(validated.get("intent_type") or "")
    is_side_question = bool(validated.get("is_side_question", False))
    keep_current = bool(validated.get("keep_current_command", False))

    if action in {"execute", "undo", "restart_session"}:
        return validated
    if is_side_question:
        return validated

    if intent_type in {"execution_request", "switch_request"} or (not keep_current and intent_type != "side_question"):
        out = dict(validated)
        out["action"] = "select_command"
        out["switch_to"] = None
        out["is_side_question"] = False
        out["keep_current_command"] = False
        out["reason"] = "normalized_select_command_from_planner_subcommand"
        return out

    return validated


def _should_try_command_selector(user_text: str, validated: Dict[str, Any], rules: Dict[str, Any]) -> bool:
    if not _normalize_text(user_text):
        return False

    action = str(validated.get("action") or "")
    intent_type = str(validated.get("intent_type") or "")
    if action in {"restart_session", "undo", "execute"}:
        return False
    if _can_commit_command_from_plan(validated, rules):
        return False
    if bool(validated.get("is_side_question", False)):
        return False
    if intent_type in {"side_question", "control"}:
        return False
    return action in {"reply_only", "clarify_intent", "ask_missing_params", "ready_to_run", "update_params"}


def llm_extract_param_updates(
    user_text: str,
    subcommand: str,
    rules: Dict[str, Any],
    current_params: Dict[str, Any],
    missing_params: List[str],
) -> Dict[str, Any]:
    if not subcommand or subcommand not in rules:
        return {}

    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    allowed = [k for k in allowed_keys_for(subcommand, rules) if k != "confirm"]
    if not allowed:
        return {}

    constrained_rule = {
        "required": rule.get("required", []),
        "optional": rule.get("optional", []),
        "choices": rule.get("choices", {}),
        "param_types": rule.get("param_types", {}),
        "param_hints": rule.get("param_hints", {}),
        "param_aliases": rule.get("param_aliases", {}),
    }

    prompt = f"""
You are an iobrpy parameter extractor, not a chatbot.
Current command: {subcommand}
You MUST only extract keys from allowed_keys.
Use only the provided command schema.
Do not invent parameters.
Do not include explanations in param_updates values.
If unsure, return empty param_updates.

Return JSON only in this exact shape:
{{"param_updates": {{}}, "confidence": 0.0, "reason": "..."}}

allowed_keys:
{json.dumps(allowed, ensure_ascii=False)}

command_schema:
{json.dumps(constrained_rule, ensure_ascii=False)}

current_params:
{json.dumps(current_params, ensure_ascii=False)}

missing_params:
{json.dumps(missing_params, ensure_ascii=False)}

user_text:
{json.dumps(user_text, ensure_ascii=False)}
""".strip()

    try:
        out = _llm_generate_json(prompt)
    except Exception:
        return {}
    if not isinstance(out, dict):
        return {}

    raw_updates = out.get("param_updates") if isinstance(out.get("param_updates"), dict) else {}
    if not raw_updates:
        return {}

    safe: Dict[str, Any] = {}
    allowed_set = set(allowed)
    for k, v in raw_updates.items():
        if k not in allowed_set:
            continue
        try:
            safe[k] = _convert_value_for_key(k, v, rules, subcommand)
        except Exception:
            continue
    return safe


def _convert_value_for_key(key: str, value: Any, rules: Dict[str, Any], subcommand: str) -> Any:
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    choices = rule.get("choices", {}) if isinstance(rule.get("choices", {}), dict) else {}
    if key in {"threads", "batch_size", "t", "k"}:
        if isinstance(value, int):
            return value
        if isinstance(value, str) and value.strip().isdigit():
            return int(value.strip())
    if key in {"mode"}:
        if isinstance(value, str):
            vl = value.strip().lower()
            if vl in {"salmon", "star"}:
                return vl
    if key in {"index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od"}:
        if isinstance(value, str):
            v = value.strip()
            if v.startswith("/"):
                return v
    if key in choices and isinstance(choices.get(key), list) and choices[key]:
        if value not in choices[key]:
            raise ValueError(f"invalid choice for {key}: {value}, allowed={choices[key]}")
    return value


def validate_planner_output(plan: Dict[str, Any], rules: Dict[str, Any], state: "SessionState") -> Dict[str, Any]:
    if not isinstance(plan, dict):
        return {"action": "reply_only", "message": "I could not parse planner output.", "confidence": 0.0}

    action = str(plan.get("action") or "reply_only")
    if action not in VALID_ACTIONS:
        dlog("invalid_action", action)
        action = "reply_only"

    out = {
        "action": action,
        "intent_type": str(plan.get("intent_type") or "other"),
        "is_side_question": bool(plan.get("is_side_question", False)),
        "keep_current_command": bool(plan.get("keep_current_command", False)),
        "switch_intent_strength": max(0.0, min(1.0, float(plan.get("switch_intent_strength") or 0.0))),
        "needs_confirmation": bool(plan.get("needs_confirmation", False)),
        "message": str(plan.get("message") or ""),
        "subcommand": plan.get("subcommand"),
        "param_updates": plan.get("param_updates") if isinstance(plan.get("param_updates"), dict) else {},
        "param_removals": plan.get("param_removals") if isinstance(plan.get("param_removals"), list) else [],
        "switch_to": plan.get("switch_to"),
        "ask_for": plan.get("ask_for") if isinstance(plan.get("ask_for"), list) else [],
        "confirm": plan.get("confirm"),
        "confidence": float(plan.get("confidence") or 0.0),
        "reason": str(plan.get("reason") or ""),
    }

    out["message"] = enforce_output_language(out["message"], state.prefer_chinese)

    # command guardrails
    if out["subcommand"] is not None and out["subcommand"] not in rules:
        dlog("invalid_subcommand", out["subcommand"])
        out["subcommand"] = None
        if out["action"] in {"select_command", "switch_command", "confirm_switch"}:
            out["action"] = "clarify_intent"

    if out["switch_to"] is not None and out["switch_to"] not in rules:
        dlog("invalid_switch_to", out["switch_to"])
        out["switch_to"] = None
        if out["action"] in {"switch_command", "confirm_switch"}:
            out["action"] = "clarify_intent"

    # parameter guardrails
    target_cmd = state.selected_command or out.get("subcommand") or out.get("switch_to")
    if out["action"] in {"update_params", "remove_params"} and not target_cmd:
        out["action"] = "clarify_intent"
        out["message"] = out["message"] or ("请先告诉我要运行哪个功能。" if state.prefer_chinese else "Please tell me which command you want first.")
        out["param_updates"] = {}
        out["param_removals"] = []

    if target_cmd and target_cmd in rules:
        allowed = set([k for k in allowed_keys_for(target_cmd, rules) if k != "confirm"])
        safe_updates: Dict[str, Any] = {}
        invalid_keys: List[str] = []
        invalid_choices: List[str] = []
        for k, v in out["param_updates"].items():
            if k not in allowed:
                invalid_keys.append(k)
                continue
            try:
                safe_updates[k] = _convert_value_for_key(k, v, rules, target_cmd)
            except Exception:
                invalid_choices.append(k)
        out["param_updates"] = safe_updates
        out["param_removals"] = [k for k in out["param_removals"] if k in allowed]

        if invalid_keys:
            out["message"] = out["message"] or (
                ("这些参数不存在，将忽略: " + ", ".join(invalid_keys)) if state.prefer_chinese else ("Ignored unknown params: " + ", ".join(invalid_keys))
            )
        if invalid_choices:
            out["message"] = out["message"] or (
                ("这些参数取值不合法，将忽略: " + ", ".join(invalid_choices)) if state.prefer_chinese else ("Ignored invalid choices for: " + ", ".join(invalid_choices))
            )

    return out


@dataclass
class SessionState:
    session_id: str
    selected_command: Optional[str] = None
    params: Dict[str, Any] = field(default_factory=dict)
    param_sources: Dict[str, str] = field(default_factory=dict)
    pending_confirmation: Optional[Dict[str, Any]] = None
    pending_switch: Optional[Dict[str, Any]] = None
    history: List[Dict[str, Any]] = field(default_factory=list)
    last_plan: Optional[Dict[str, Any]] = None
    last_message: Optional[str] = None
    prefer_chinese: bool = False
    phase: str = "idle"  # idle / clarifying / collecting / ready / executing

    def snapshot(self):
        self.history.append(
            {
                "selected_command": self.selected_command,
                "params": deepcopy(self.params),
                "param_sources": deepcopy(self.param_sources),
                "pending_confirmation": deepcopy(self.pending_confirmation),
                "pending_switch": deepcopy(self.pending_switch),
                "last_plan": deepcopy(self.last_plan),
                "last_message": self.last_message,
                "prefer_chinese": self.prefer_chinese,
                "phase": self.phase,
            }
        )
        if len(self.history) > 100:
            self.history = self.history[-100:]

    def undo(self) -> bool:
        if not self.history:
            return False
        snap = self.history.pop()
        self.selected_command = snap.get("selected_command")
        self.params = deepcopy(snap.get("params", {}))
        self.param_sources = deepcopy(snap.get("param_sources", {}))
        self.pending_confirmation = deepcopy(snap.get("pending_confirmation"))
        self.pending_switch = deepcopy(snap.get("pending_switch"))
        self.last_plan = deepcopy(snap.get("last_plan"))
        self.last_message = snap.get("last_message")
        self.prefer_chinese = bool(snap.get("prefer_chinese", self.prefer_chinese))
        self.phase = snap.get("phase", self.phase)
        return True

    def to_planner_state(self) -> Dict[str, Any]:
        return {
            "selected_command": self.selected_command,
            "params": self.params,
            "param_sources": self.param_sources,
            "pending_confirmation": self.pending_confirmation,
            "pending_switch": self.pending_switch,
            "last_message": self.last_message,
            "prefer_chinese": self.prefer_chinese,
            "phase": self.phase,
        }


SESSIONS: Dict[str, SessionState] = {}


def get_session(session_id: str) -> SessionState:
    if session_id not in SESSIONS:
        SESSIONS[session_id] = SessionState(session_id=session_id)
    return SESSIONS[session_id]


def _resolve_iobrpy_pythonpath() -> Optional[str]:
    cands: List[str] = []
    if IOBRPY_SRC_DIR:
        cands.append(IOBRPY_SRC_DIR)
    cands.extend(["/mnt/tmp-ext4/iobrpy/src", "/mnt/tmp-ext4/iobrpy"])
    for p in cands:
        try:
            if p and os.path.isdir(os.path.join(p, "iobrpy")):
                return p
        except Exception:
            pass
    return None


def _state_summary_text(state: SessionState, merged: Dict[str, Any], missing: List[str], sources: Dict[str, str], prefer_chinese: bool) -> str:
    items = []
    for k in sorted(merged.keys()):
        v = merged.get(k)
        if v in (None, "", []):
            continue
        src = sources.get(k, "user")
        tag = " (默认)" if src == "default" else ""
        items.append(f"{k}={v}{tag}")
    recognized = ", ".join(items) if items else ("无" if prefer_chinese else "none")
    miss = [m for m in missing if m != "__one_of_group"]
    missing_txt = ", ".join(miss) if miss else ("无" if prefer_chinese else "none")
    if prefer_chinese:
        return f"当前命令: iobrpy {state.selected_command}\n已识别: {recognized}\n缺失: {missing_txt}"
    return f"Current command: iobrpy {state.selected_command}\nRecognized: {recognized}\nMissing: {missing_txt}"


def _current_command_context(state: SessionState, rules: Dict[str, Any], message: str, needs: Optional[List[str]] = None) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "status": "need_info",
        "question": enforce_output_language(message, state.prefer_chinese),
        "needs": needs or [],
    }
    if not state.selected_command:
        return out

    merged, sources = merge_defaults(state.params, state.selected_command, rules)
    missing, _ = compute_needs(state.selected_command, rules, merged)
    try:
        draft_cmd, _ = build_command(state.selected_command, merged, rules)
    except ValueError:
        draft_cmd = None

    out.update(
        {
            "subcommand": state.selected_command,
            "draft_command": draft_cmd,
            "params": merged,
            "needs": missing if needs is None else needs,
            "state_summary": _state_summary_text(state, merged, missing, sources, state.prefer_chinese),
        }
    )
    return out


def apply_planner_output(state: SessionState, plan: Dict[str, Any], rules: Dict[str, Any], user_text: str) -> Dict[str, Any]:
    action = plan.get("action", "reply_only")
    msg = enforce_output_language(plan.get("message") or "", state.prefer_chinese)
    is_side_question = bool(plan.get("is_side_question", False))
    keep_current = bool(plan.get("keep_current_command", False))
    switch_strength = max(0.0, min(1.0, float(plan.get("switch_intent_strength") or 0.0)))
    needs_confirmation = bool(plan.get("needs_confirmation", False))

    if action == "restart_session":
        state.snapshot()
        state.selected_command = None
        state.params = {}
        state.param_sources = {}
        state.pending_confirmation = None
        state.pending_switch = None
        state.phase = "idle"
        return {"status": "need_info", "question": enforce_output_language(msg or ("请描述你的任务。" if state.prefer_chinese else "Please describe your task."), state.prefer_chinese), "needs": ["task"]}

    if action == "undo":
        ok = state.undo()
        if not ok:
            return {"status": "need_info", "question": enforce_output_language(msg or ("没有可撤销的历史。" if state.prefer_chinese else "No history to undo."), state.prefer_chinese), "needs": []}

    if is_side_question or (action == "reply_only" and keep_current):
        if state.selected_command:
            state.phase = "collecting" if state.phase in {"idle", "clarifying"} else state.phase
        else:
            state.phase = "idle"
        return _current_command_context(
            state,
            rules,
            msg or ("请继续告诉我你的任务。" if state.prefer_chinese else "Please continue with your task details."),
        )

    if state.pending_switch and action == "confirm_switch":
        target = plan.get("switch_to") or plan.get("subcommand") or state.pending_switch.get("target")
        if target and target in rules:
            state.snapshot()
            old_params = deepcopy(state.params)
            old_sources = deepcopy(state.param_sources)
            allowed = set([k for k in allowed_keys_for(target, rules) if k != "confirm"])
            state.selected_command = target
            state.params = {k: v for k, v in old_params.items() if k in allowed}
            state.param_sources = {k: old_sources.get(k, "user") for k in state.params.keys()}
            state.pending_switch = None
            state.phase = "collecting"

    if action in {"switch_command", "select_command"}:
        target = plan.get("switch_to") or plan.get("subcommand")
        if state.selected_command and target and target != state.selected_command:
            if keep_current:
                action = "reply_only"
            elif needs_confirmation or switch_strength < 0.75:
                state.pending_switch = {"from": state.selected_command, "target": target}
                state.phase = "clarifying"
                return _current_command_context(
                    state,
                    rules,
                    msg or (
                        f"检测到你可能想从 {state.selected_command} 切换到 {target}，是否切换？"
                        if state.prefer_chinese
                        else f"I detected you may want to switch from {state.selected_command} to {target}. Switch now?"
                    ),
                    needs=["confirm_switch"],
                )
            elif target in rules:
                state.snapshot()
                old_params = deepcopy(state.params)
                old_sources = deepcopy(state.param_sources)
                allowed = set([k for k in allowed_keys_for(target, rules) if k != "confirm"])
                state.selected_command = target
                state.params = {k: v for k, v in old_params.items() if k in allowed}
                state.param_sources = {k: old_sources.get(k, "user") for k in state.params.keys()}
                state.pending_switch = None
                state.phase = "collecting"
        elif action == "select_command":
            sub = plan.get("subcommand")
            if sub and sub in rules:
                state.snapshot()
                state.selected_command = sub
                state.params = {}
                state.param_sources = {}
                state.phase = "collecting"
            else:
                state.phase = "clarifying"
                return {"status": "need_info", "question": enforce_output_language(msg or ("我还不确定你要哪个功能，请再描述。" if state.prefer_chinese else "I am not sure which command you need yet. Please clarify."), state.prefer_chinese), "needs": ["clarify_intent"]}

    if action == "clarify_intent":
        state.phase = "clarifying"
        return _current_command_context(
            state,
            rules,
            msg or ("请再具体描述你的分析目标。" if state.prefer_chinese else "Please describe your analysis goal more specifically."),
            needs=["clarify_intent"],
        )

    if action == "reply_only":
        if state.selected_command:
            state.phase = "collecting" if state.phase in {"idle", "clarifying"} else state.phase
        else:
            state.phase = "idle"
        return _current_command_context(
            state,
            rules,
            msg or ("请继续告诉我你的任务。" if state.prefer_chinese else "Please continue with your task details."),
        )

    if action == "update_params":
        if not state.selected_command:
            return {"status": "need_info", "question": enforce_output_language(msg or ("请先告诉我要运行哪个功能。" if state.prefer_chinese else "Please tell me which command you want first."), state.prefer_chinese), "needs": ["task"]}
        if plan.get("param_updates"):
            state.snapshot()
            for k, v in plan.get("param_updates", {}).items():
                state.params[k] = v
                state.param_sources[k] = "user"
            state.phase = "collecting"

    if action == "remove_params":
        if not state.selected_command:
            return {"status": "need_info", "question": enforce_output_language(msg or ("请先告诉我要运行哪个功能。" if state.prefer_chinese else "Please tell me which command you want first."), state.prefer_chinese), "needs": ["task"]}
        if plan.get("param_removals"):
            state.snapshot()
            for k in plan.get("param_removals", []):
                state.params[k] = None
                state.param_sources[k] = "user"
            state.phase = "collecting"

    if action == "ask_missing_params":
        state.phase = "collecting"

    if action == "ready_to_run":
        state.phase = "ready"

    if action == "execute":
        state.phase = "executing"

    # compute final guarded status
    if not state.selected_command:
        state.phase = "idle"
        return {"status": "need_info", "question": enforce_output_language(msg or ("请先告诉我你想运行的分析功能。" if state.prefer_chinese else "Please first tell me which analysis command you want."), state.prefer_chinese), "needs": ["task"]}

    merged, sources = merge_defaults(state.params, state.selected_command, rules)
    missing, need_confirm = compute_needs(state.selected_command, rules, merged)
    mapping_errors = _required_without_serializable_flag(state.selected_command, rules)
    try:
        draft_cmd, _ = build_command(state.selected_command, merged, rules)
    except ValueError:
        draft_cmd = None

    if mapping_errors:
        state.phase = "collecting"
        q = msg or (
            "无法执行：命令参数映射配置不完整，请补齐以下必填参数的CLI flag映射："
            if state.prefer_chinese
            else "Cannot execute: command mapping is incomplete. Missing CLI flag mappings for required parameters:"
        )
        q = f"{q} {', '.join(mapping_errors)}"
        return {
            "status": "need_info",
            "subcommand": state.selected_command,
            "question": enforce_output_language(q, state.prefer_chinese),
            "needs": mapping_errors,
            "draft_command": draft_cmd,
            "params": merged,
            "state_summary": _state_summary_text(state, merged, missing, sources, state.prefer_chinese),
        }

    consistency = audit_parameter_link_consistency(rules).get(state.selected_command, {})
    optional_unserializable = consistency.get("optional_unserializable", []) if isinstance(consistency, dict) else []

    if action == "execute" and (missing or need_confirm):
        action = "ask_missing_params"

    if missing or action == "ask_missing_params":
        state.phase = "collecting"
        q = enforce_output_language(
            msg or ("还不能执行，参数尚未完整。" if state.prefer_chinese else "Cannot execute yet, parameters are incomplete."),
            state.prefer_chinese,
        )
        return {
            "status": "need_info",
            "subcommand": state.selected_command,
            "question": q,
            "needs": missing,
            "draft_command": draft_cmd,
            "params": merged,
            "warnings": optional_unserializable,
            "state_summary": _state_summary_text(state, merged, missing, sources, state.prefer_chinese),
        }

    # ready gate
    state.phase = "ready"
    if action == "execute":
        return {
            "status": "ready",
            "subcommand": state.selected_command,
            "question": enforce_output_language(msg or ("已准备好执行。" if state.prefer_chinese else "Ready to execute."), state.prefer_chinese),
            "draft_command": draft_cmd,
            "params": merged,
            "needs": missing,
            "warnings": optional_unserializable,
            "state_summary": _state_summary_text(state, merged, missing, sources, state.prefer_chinese),
        }

    return {
        "status": "ready",
        "subcommand": state.selected_command,
        "question": enforce_output_language(msg or ("参数已满足，可以执行。" if state.prefer_chinese else "All required parameters are satisfied. Ready to run."), state.prefer_chinese),
        "draft_command": draft_cmd,
        "params": merged,
        "needs": missing,
        "warnings": optional_unserializable,
        "state_summary": _state_summary_text(state, merged, missing, sources, state.prefer_chinese),
    }


def tool_iobrpy_assistant(session_id: str, task: Optional[str] = None, answer_text: Optional[str] = None, run: bool = False) -> Dict[str, Any]:
    rules = load_rules()
    state = get_session(session_id)

    text = answer_text if answer_text is not None else task
    user_text = _normalize_text(text or "")

    # Language policy:
    # - default English
    # - once user sends Chinese text, switch to Chinese mode
    # - do not auto-switch back to English in later turns
    if user_text and _contains_cjk(user_text):
        state.prefer_chinese = True

    if user_text in {":restart", "restart", "restart session"}:
        plan = {"action": "restart_session", "message": ""}
    else:
        plan = llm_plan_next_action(state, user_text, rules)

    validated = validate_planner_output(plan, rules, state)
    validated = normalize_command_commit(validated, state, rules)

    selector_fallback_used = False
    selector_out: Dict[str, Any] = {"subcommand": None, "confidence": 0.0, "reason": ""}
    selector_threshold = float(os.getenv("IOBRPY_SELECTOR_CONFIDENCE", "0.7"))
    if state.selected_command is None and _should_try_command_selector(user_text, validated, rules):
        selector_out = llm_select_command_from_catalog(
            user_text=user_text,
            rules=rules,
            current_command=state.selected_command,
            prefer_chinese=state.prefer_chinese,
        )
        selected = selector_out.get("subcommand")
        conf = float(selector_out.get("confidence") or 0.0)
        if selected in rules and conf >= selector_threshold:
            selector_fallback_used = True
            validated = validate_planner_output(
                {
                    "action": "select_command",
                    "intent_type": "execution_request",
                    "is_side_question": False,
                    "keep_current_command": False,
                    "switch_intent_strength": 1.0,
                    "needs_confirmation": False,
                    "message": "",
                    "subcommand": selected,
                    "param_updates": {},
                    "param_removals": [],
                    "switch_to": None,
                    "ask_for": [],
                    "confirm": None,
                    "confidence": conf,
                    "reason": "catalog_selector_fallback",
                },
                rules,
                state,
            )
            validated = normalize_command_commit(validated, state, rules)

    fallback_triggered = False
    schema_updates: Dict[str, Any] = {}
    fallback_updates: Dict[str, Any] = {}
    if state.selected_command:
        merged_now, _ = merge_defaults(state.params, state.selected_command, rules)
        missing_now, _ = compute_needs(state.selected_command, rules, merged_now)
        planner_updates = validated.get("param_updates") if isinstance(validated.get("param_updates"), dict) else {}
        looks_like_backfill = (
            validated.get("intent_type") != "side_question"
            and _text_has_parameter_shape(user_text, state.selected_command, rules)
        )

        if looks_like_backfill:
            planner_needs_help = (
                validated.get("action") != "update_params"
                or not planner_updates
            )
            if planner_needs_help:
                schema_updates = llm_extract_param_updates(
                    user_text=user_text,
                    subcommand=state.selected_command,
                    rules=rules,
                    current_params=merged_now,
                    missing_params=missing_now,
                )

            fallback_updates = extract_param_updates_from_text(
                user_text=user_text,
                subcommand=state.selected_command,
                rules=rules,
                current_params=merged_now,
                missing_params=missing_now,
            )

            merged_updates = dict(planner_updates)
            merged_updates.update(schema_updates)
            merged_updates.update(fallback_updates)
            if merged_updates:
                fallback_triggered = bool(schema_updates or fallback_updates)
                validated["param_updates"] = merged_updates
                validated["action"] = "update_params"
                validated["is_side_question"] = False
                validated["keep_current_command"] = True
                validated = validate_planner_output(validated, rules, state)

    dlog(
        "phase_before", state.phase,
        "planner_action", plan.get("action"),
        "planner_param_updates", plan.get("param_updates"),
        "validated_action", validated.get("action"),
        "validated_param_updates", validated.get("param_updates"),
        "selector_fallback_used", selector_fallback_used,
        "selector_out", selector_out,
        "fallback_triggered", fallback_triggered,
        "schema_updates", schema_updates,
        "fallback_updates", fallback_updates,
    )

    state.last_plan = deepcopy(validated)
    state.last_message = user_text

    out = apply_planner_output(state, validated, rules, user_text)

    if state.selected_command and out.get("status") not in {"done", "error"}:
        ctx = _current_command_context(
            state,
            rules,
            out.get("question") or "",
            needs=out.get("needs") if isinstance(out.get("needs"), list) else None,
        )
        for k in ["subcommand", "params", "draft_command", "state_summary", "needs"]:
            if k in ctx and (out.get(k) is None or k not in out):
                out[k] = ctx[k]

    out["prefer_chinese"] = state.prefer_chinese
    out["phase"] = state.phase
    out["intent_category"] = validated.get("action")

    # runtime execute path is guarded and optional (ai.py usually calls with run=False)
    if out.get("status") == "ready" and run:
        sub = out.get("subcommand")
        params = out.get("params") or {}
        try:
            draft_cmd, argv = build_command(str(sub), dict(params), rules)
        except ValueError as e:
            return {
                "status": "error",
                "error": str(e),
                "subcommand": sub,
                "params": params,
            }

        run_dir = os.getenv("IOBRPY_RUN_LOG_DIR", os.path.join(os.path.dirname(__file__), "mcp_runs"))
        os.makedirs(run_dir, exist_ok=True)
        log_path = os.path.join(run_dir, f"{session_id}_{sub}.log")
        env = os.environ.copy()
        pp = _resolve_iobrpy_pythonpath()
        if pp:
            env["PYTHONPATH"] = pp + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")

        cmd = [CONDA_EXE, "run", "-n", IOBRPY_CONDA_ENV, "python", "-m", "iobrpy.main"] + argv
        try:
            with open(log_path, "w", encoding="utf-8") as f:
                p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, env=env)
                rc = p.wait()
        except Exception as e:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write("\n\nEXCEPTION:\n" + repr(e) + "\n" + traceback.format_exc())
            rc = 1

        tail = ""
        try:
            with open(log_path, "r", encoding="utf-8") as f:
                tail_lines = 30 if rc == 0 else 80
                tail = "".join(f.readlines()[-tail_lines:])
        except Exception:
            pass

        out = {
            "status": "done" if rc == 0 else "error",
            "returncode": rc,
            "draft_command": draft_cmd,
            "log_path": log_path,
            "tail": tail,
            "prefer_chinese": state.prefer_chinese,
            "phase": state.phase,
            "intent_category": validated.get("action"),
        }

    dlog("phase_after", state.phase, "status", out.get("status"))
    return out


TOOLS = [
    {
        "name": "iobrpy_assistant",
        "title": "iobrpy assistant (planner-first BYOK)",
        "description": "Planner-first assistant: LLM plan + Python guardrails + local execution.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "task": {"type": "string"},
                "answer_text": {"type": "string"},
                "run": {"type": "boolean", "default": False},
            },
            "required": ["session_id"],
            "additionalProperties": False,
        },
    },
    {
        "name": "set_defaults",
        "title": "Set defaults (strict)",
        "description": "Persist defaults; only keys present in REQUIRED_PARAMS file are stored.",
        "inputSchema": {
            "type": "object",
            "properties": {"defaults": {"type": "object"}},
            "required": ["defaults"],
            "additionalProperties": False,
        },
    },
    {
        "name": "get_defaults",
        "title": "Get defaults",
        "description": "Read stored defaults.",
        "inputSchema": {"type": "object", "properties": {}, "additionalProperties": False},
    },
]


def handle_initialize(_id, params):
    requested = (params or {}).get("protocolVersion")
    chosen = requested if requested in SUPPORTED_PROTOCOLS else SUPPORTED_PROTOCOLS[-1]
    send({"jsonrpc": "2.0", "id": _id, "result": {"protocolVersion": chosen, "capabilities": CAPABILITIES, "serverInfo": SERVER_INFO}})


def handle_tools_list(_id, params):
    send({"jsonrpc": "2.0", "id": _id, "result": {"tools": TOOLS}})


def handle_tools_call(_id, params):
    name = (params or {}).get("name")
    args = (params or {}).get("arguments") or {}
    try:
        if name == "set_defaults":
            out = set_defaults(args.get("defaults") or {}, load_rules())
            send({"jsonrpc": "2.0", "id": _id, "result": {"content": [{"type": "text", "text": json.dumps(out, ensure_ascii=False)}], "isError": False}})
            return
        if name == "get_defaults":
            out = load_defaults()
            send({"jsonrpc": "2.0", "id": _id, "result": {"content": [{"type": "text", "text": json.dumps(out, ensure_ascii=False)}], "isError": False}})
            return
        if name == "iobrpy_assistant":
            sid = args.get("session_id")
            if not sid:
                raise ValueError("session_id is required")
            out = tool_iobrpy_assistant(sid, task=args.get("task"), answer_text=args.get("answer_text"), run=bool(args.get("run", False)))
            send({"jsonrpc": "2.0", "id": _id, "result": {"content": [{"type": "text", "text": json.dumps(out, ensure_ascii=False)}], "isError": False}})
            return
        jsonrpc_error(_id, -32601, f"Unknown tool: {name}")
    except Exception as e:
        dlog("tools/call error", repr(e))
        dlog(traceback.format_exc())
        jsonrpc_error(_id, -32603, "Internal error", data={"exception": repr(e)})


def main():
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        msg = {}
        try:
            msg = json.loads(line)
            method = msg.get("method")
            _id = msg.get("id", None)
            params = msg.get("params", None)
            if method == "initialize":
                handle_initialize(_id, params)
            elif method == "tools/list":
                handle_tools_list(_id, params)
            elif method == "tools/call":
                handle_tools_call(_id, params)
            elif method == "notifications/initialized":
                continue
            else:
                if _id is not None:
                    jsonrpc_error(_id, -32601, f"Method not found: {method}")
        except Exception:
            dlog(traceback.format_exc())
            if isinstance(msg, dict) and msg.get("id") is not None:
                jsonrpc_error(msg["id"], -32603, "Internal error")


if __name__ == "__main__":
    main()
