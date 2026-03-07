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
    "od": "--od", "t": "-t", "k": "-k",
}
BOOL_FLAGS = {
    "repseq": "--repseq", "skipMateExtension": "--skipMateExtension", "abnormalUnmapFlag": "--abnormalUnmapFlag",
    "assembleWithRef": "--assembleWithRef", "noExtraction": "--noExtraction", "outputReadAssignment": "--outputReadAssignment",
}


def build_command(subcommand: str, params: Dict[str, Any], rules: Dict[str, Any]) -> Tuple[str, List[str]]:
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    argv = [subcommand]
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    ordered = list(rule.get("required", [])) + [k for k in rule.get("optional", []) if k not in rule.get("required", [])]
    for k in ordered:
        if k not in allowed_set or k not in params:
            continue
        v = params[k]
        if v in (None, "", []):
            continue
        if isinstance(v, bool):
            bf = BOOL_FLAGS.get(k)
            if bf and v:
                argv.append(bf)
            continue
        flag = FLAG_MAP.get(k)
        if flag:
            argv.extend([flag, str(v)])
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
    draft = build_command(sub, merged, rules)[0] if sub else None
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
        state.phase = "idle" if not state.selected_command else state.phase
        return {"status": "need_info", "question": enforce_output_language(msg or ("请继续告诉我你的任务。" if state.prefer_chinese else "Please continue with your task details."), state.prefer_chinese), "needs": [] if state.selected_command else ["task"]}

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
                return {
                    "status": "need_info",
                    "question": enforce_output_language(msg or (
                        f"检测到你可能想从 {state.selected_command} 切换到 {target}，是否切换？"
                        if state.prefer_chinese
                        else f"I detected you may want to switch from {state.selected_command} to {target}. Switch now?"
                    ), state.prefer_chinese),
                    "needs": ["confirm_switch"],
                }
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
        return {"status": "need_info", "question": enforce_output_language(msg or ("请再具体描述你的分析目标。" if state.prefer_chinese else "Please describe your analysis goal more specifically."), state.prefer_chinese), "needs": ["clarify_intent"]}

    if action == "reply_only":
        state.phase = "idle" if not state.selected_command else state.phase
        return {"status": "need_info", "question": enforce_output_language(msg or ("请继续告诉我你的任务。" if state.prefer_chinese else "Please continue with your task details."), state.prefer_chinese), "needs": ["task"] if not state.selected_command else []}

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
    draft_cmd, _ = build_command(state.selected_command, merged, rules)

    if action == "execute" and (missing or need_confirm):
        action = "ask_missing_params"

    if missing or action == "ask_missing_params":
        state.phase = "collecting"
        q = msg or ("还不能执行，参数尚未完整。" if state.prefer_chinese else "Cannot execute yet, parameters are incomplete.")
        q += "\n\n" + _state_summary_text(state, merged, missing, sources, state.prefer_chinese)
        q = enforce_output_language(q, state.prefer_chinese)
        return {
            "status": "need_info",
            "subcommand": state.selected_command,
            "question": q,
            "needs": missing,
            "draft_command": draft_cmd,
            "params": merged,
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
        }

    return {
        "status": "ready",
        "subcommand": state.selected_command,
        "question": enforce_output_language(msg or ("参数已满足，可以执行。" if state.prefer_chinese else "All required parameters are satisfied. Ready to run."), state.prefer_chinese),
        "draft_command": draft_cmd,
        "params": merged,
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
    dlog("phase_before", state.phase, "validated_action", validated.get("action"))

    state.last_plan = deepcopy(validated)
    state.last_message = user_text

    out = apply_planner_output(state, validated, rules, user_text)

    out["prefer_chinese"] = state.prefer_chinese
    out["phase"] = state.phase
    out["intent_category"] = validated.get("action")

    # runtime execute path is guarded and optional (ai.py usually calls with run=False)
    if out.get("status") == "ready" and run:
        sub = out.get("subcommand")
        params = out.get("params") or {}
        draft_cmd, argv = build_command(str(sub), dict(params), rules)

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
                tail = "".join(f.readlines()[-80:])
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
