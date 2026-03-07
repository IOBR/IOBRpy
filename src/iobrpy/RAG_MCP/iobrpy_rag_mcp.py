#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BYOK iobrpy assistant server.

No local embedding / Chroma / Ollama dependency.
"""

import json
import os
import re
import subprocess
import sys
import traceback
from copy import deepcopy
from dataclasses import dataclass
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
    "title": "iobrpy AI Agentic (BYOK cloud LLM)",
    "version": "0.4.0",
    "description": "Natural language → choose iobrpy function → fill allowed params only → ask missing → run when ready.",
}
CAPABILITIES = {"tools": {"listChanged": False}}


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
    return os.getenv("IOBRPY_AI_DEBUG", "0").strip() in {"1", "true", "TRUE", "yes", "on"}


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
    if not text:
        return False
    return re.search(r"[\u4e00-\u9fff]", text) is not None


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
    parts = data.get("content") or []
    out = ""
    for p in parts:
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
    texts = [p.get("text", "") for p in parts if isinstance(p, dict)]
    return "\n".join(texts).strip()


def _llm_generate_json(prompt: str) -> Dict[str, Any]:
    cfg = _require_ai_config()
    alias = cfg.llm_alias
    if alias == "claude":
        txt = _post_anthropic(cfg.base_url, cfg.api_key, cfg.model, prompt)
    elif alias == "gemini":
        txt = _post_gemini(cfg.base_url, cfg.api_key, cfg.model, prompt)
    else:
        txt = _post_openai_compatible(cfg.base_url, cfg.api_key, cfg.model, prompt)
    return _extract_json_block(txt)


def _translate_to_english(text: str) -> str:
    if not text:
        return ""
    if not _contains_cjk(text):
        return text
    prompt = f"""
Translate the following user text to English faithfully.
Keep subcommand names, parameter names, CLI flags, file paths, numbers, and enum values unchanged.
Return JSON only: {{"english": "..."}}

Text:
{text}
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
Translate the following UI text to Chinese faithfully.
Keep subcommand names, parameter names, CLI flags, file paths, numbers, and enum values unchanged.
Return JSON only: {{"chinese": "..."}}

Text:
{text}
""".strip()
    try:
        j = _llm_generate_json(prompt)
        out = j.get("chinese") if isinstance(j, dict) else None
        if isinstance(out, str) and out.strip():
            return out.strip()
    except Exception:
        pass
    return text


DEFAULT_INTENT_KEYWORDS: Dict[str, List[str]] = {
    "runall": ["fastq", "workflow", "pipeline", "from fastq to tme", "bulk workflow"],
    "trust4": ["tcr", "bcr", "vdj", "repertoire", "clonotype", "immune receptor"],
    "spechla": ["hla", "typing"],
    "extract_hla_read": ["hla", "extract hla read"],
    "hla_typing": ["hla", "typing", "hla analysis"],
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
    out = {}
    for k, v in rules.items():
        if not isinstance(v, dict):
            continue
        req = v.get("required", [])
        conf = v.get("confirm", [])
        opt = v.get("optional", [])
        choices = v.get("choices", {})
        req_one = v.get("required_one_of", [])
        notes = v.get("notes", {})
        optional_defaults = v.get("optional_defaults", {})
        intent_keywords = v.get("intent_keywords", [])
        param_hints = v.get("param_hints", {})
        function_summary = v.get("function_summary", "")
        if not isinstance(req, list):
            req = []
        if not isinstance(conf, list):
            conf = []
        if not isinstance(opt, list):
            opt = []
        if not isinstance(choices, dict):
            choices = {}
        if not isinstance(req_one, list):
            req_one = []
        if not isinstance(notes, dict):
            notes = {}
        if not isinstance(optional_defaults, dict):
            optional_defaults = {}
        if not isinstance(intent_keywords, (list, dict)):
            intent_keywords = []
        if not isinstance(param_hints, dict):
            param_hints = {}
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
    r = rules.get(subcommand, {})
    if not isinstance(r, dict):
        return []
    keys = []
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


def _has_token(text_l: str, token: str) -> bool:
    return re.search(rf"\b{re.escape(token)}\b", text_l) is not None


def _sanitize_path(p: str) -> str:
    if not p:
        return p
    p = p.strip().strip('"').strip("'")
    p = re.sub(r"[,\.;:\)\]\}]+$", "", p)
    return p


POSITIVE_CONFIRM_WORDS = {
    "y", "yes", "confirm", "run", "execute", "ok", "go", "sure",
    "确认", "执行", "是", "好", "开始", "现在开始", "确认运行",
}
NEGATIVE_CONFIRM_WORDS = {
    "n", "no", "cancel", "stop", "nope", "取消", "否", "不", "不执行",
}
UNDO_WORDS = {
    "undo", "rollback", "back", "undo last", "撤销", "回退", "上一步",
}
RESET_COMMAND_WORDS = {
    "reset command", "rechoose command", "change command", "重新选命令", "重选命令", "重来",
}


def _is_positive_confirm_text(text: str) -> bool:
    t = _normalize_text(text).lower()
    if t in POSITIVE_CONFIRM_WORDS:
        return True
    return any(x in t for x in ["确认运行", "现在开始", "好 执行", "好,执行", "start now"])


def _is_negative_confirm_text(text: str) -> bool:
    t = _normalize_text(text).lower()
    if t in NEGATIVE_CONFIRM_WORDS:
        return True
    return any(x in t for x in ["不执行", "先不", "暂不"])


def _wants_undo(text: str) -> bool:
    t = _normalize_text(text).lower()
    return t in UNDO_WORDS or any(x in t for x in ["撤销", "回退", "上一步"])


def _wants_reset_command(text: str) -> bool:
    t = _normalize_text(text).lower()
    if t in RESET_COMMAND_WORDS:
        return True
    return any(x in t for x in ["重新选命令", "重选命令", "change command", "rechoose command"])


def _regex_extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    t = _normalize_text(text)
    tl = t.lower()
    out: Dict[str, Any] = {}

    if "mode" in allowed:
        if _has_token(tl, "salmon"):
            out["mode"] = "salmon"
        elif _has_token(tl, "star"):
            out["mode"] = "star"
        elif re.search(r"改成\s*star|改为\s*star", t, flags=re.I):
            out["mode"] = "star"
        elif re.search(r"改成\s*salmon|改为\s*salmon", t, flags=re.I):
            out["mode"] = "salmon"
        elif re.search(r"改回\s*star", t, flags=re.I):
            out["mode"] = "star"
        elif re.search(r"改回\s*salmon", t, flags=re.I):
            out["mode"] = "salmon"
    for key in ["index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od"]:
        if key in allowed:
            m = re.search(rf"\b{re.escape(key)}\b\s*(?:is|=|:)?\s*(/[^ \t;,，；。]+)", t, flags=re.I)
            if m:
                out[key] = _sanitize_path(m.group(1))
                continue
            m_cn_as = re.search(rf"用\s*(/[^ \t;,，；。]+)\s*作为\s*{re.escape(key)}", t, flags=re.I)
            if m_cn_as:
                out[key] = _sanitize_path(m_cn_as.group(1))

    for key in ["threads", "batch_size", "t", "k"]:
        if key in allowed:
            m = re.search(rf"\b{re.escape(key)}s?\b\s*(?:is|=|:)?\s*(\d+)\b", tl)
            if m:
                out[key] = int(m.group(1))
                continue
            m_cn = re.search(rf"(?:把\s*)?{re.escape(key)}\s*(?:改成|改为|设为|设置为|修改为)\s*(\d+)\b", t, flags=re.I)
            if m_cn:
                out[key] = int(m_cn.group(1))

    for key in ["project", "o", "f"]:
        if key in allowed:
            m = re.search(rf"\b{re.escape(key)}\b\s*(?:is|=|:)?\s*([A-Za-z0-9_.-]+)", t, flags=re.I)
            if m:
                out[key] = m.group(1)

    if "confirm" in allowed:
        confirms = re.findall(r"\bconfirm\s+([A-Za-z0-9_]+)\b", tl)
        if confirms:
            out["confirm"] = sorted(set(confirms))

    m_reset = re.findall(r"\b(?:reset|clear|unset|remove)\s+([A-Za-z0-9_]+)\b", tl)
    m_reset_cn = re.findall(r"(?:删除|去掉|清空|重置)\s*([A-Za-z0-9_]+)", t)
    to_reset = sorted(set([k for k in (m_reset + m_reset_cn) if k in allowed]))
    if to_reset:
        out["_reset"] = to_reset

    # generic key=value edits
    for key in allowed:
        if key in ("confirm",):
            continue
        mv = re.search(rf"\b{re.escape(key)}\b\s*=\s*([^,;\s]+)", t, flags=re.I)
        if mv:
            raw = mv.group(1)
            if key in ("threads", "batch_size", "t", "k"):
                if raw.isdigit():
                    out[key] = int(raw)
            elif key in ("mode",):
                if raw.lower() in ("salmon", "star"):
                    out[key] = raw.lower()
            elif key in ("index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od"):
                sv = _sanitize_path(raw)
                if sv.startswith("/"):
                    out[key] = sv
            else:
                out[key] = raw

    # Chinese/English set patterns with explicit target key
    for key in allowed:
        if key in ("confirm",):
            continue
        ms = re.search(rf"(?:set|change|modify|update)\s+{re.escape(key)}\s+(?:to|as)?\s*([^,;]+)", t, flags=re.I)
        ms_cn = re.search(rf"(?:把\s*)?{re.escape(key)}\s*(?:改成|改为|设置为|设为|修改为|用)\s*([^,;]+)", t, flags=re.I)
        g = (ms.group(1).strip() if ms else None) or (ms_cn.group(1).strip() if ms_cn else None)
        if not g:
            continue
        if key in ("threads", "batch_size", "t", "k"):
            mnum = re.search(r"\d+", g)
            if mnum:
                out[key] = int(mnum.group(0))
        elif key == "mode":
            gl = g.lower()
            if "star" in gl:
                out[key] = "star"
            elif "salmon" in gl:
                out[key] = "salmon"
        elif key in ("index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od"):
            mpath = re.search(r"/[^ \t;,，；。]+", g)
            if mpath:
                out[key] = _sanitize_path(mpath.group(0))
        else:
            out[key] = g

    if re.search(r"\b(?:undo|rollback|back)\b|撤销|回退|上一步", tl):
        out["_undo"] = True
    if re.search(r"重新选命令|重选命令|change\s+command|rechoose\s+command", tl):
        out["_reset_command"] = True

    return out


def _llm_extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    prompt = f"""
You are a strict information extraction system.
Extract ONLY fields explicitly present in the user's message.
Return JSON object only.

Allowed keys: {allowed}
Rules:
- Path-like fields must be absolute Linux paths (start with /).
- threads/batch_size/t/k must be integers.
- mode can only be salmon or star.
- If a field is not explicitly provided, omit it.

User message:
{text}
""".strip()
    j = _llm_generate_json(prompt)
    return j if isinstance(j, dict) else {}


def _evidence_filter(text: str, extracted: Dict[str, Any], allowed: List[str]) -> Dict[str, Any]:
    t = _normalize_text(text)
    out: Dict[str, Any] = {}
    for k, v in extracted.items():
        if k not in allowed and k not in ("confirm", "_reset"):
            continue
        if isinstance(v, str):
            if k in ("index", "input", "output", "fastq", "outdir", "bam", "fqdir", "ru", "r1", "r2", "ref", "od"):
                sv = _sanitize_path(v)
                if sv.startswith("/"):
                    out[k] = sv
            elif v and v.lower() in t.lower():
                out[k] = v
        elif isinstance(v, (int, float, bool, list)):
            out[k] = v
    return out


def extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    out = _regex_extract_slots(text, allowed)
    try:
        llm = _llm_extract_slots(text, allowed)
        out.update(_evidence_filter(text, llm, allowed))
    except Exception as e:
        log("LLM extraction failed:", e)
    return out


def extract_slots_multisource(raw_text: str, translated_text: Optional[str], allowed: List[str]) -> Dict[str, Any]:
    primary = extract_slots(raw_text, allowed)
    t = translated_text or ""
    if t and t != raw_text:
        secondary = extract_slots(t, allowed)
        primary.update(secondary)
    return primary


def _intent_tokenize(text: str) -> set:
    text = text.lower().replace("_", " ")
    toks = set(re.findall(r"[a-z][a-z0-9]+", text))
    toks.update([x for x in re.findall(r"[一-鿿]{2,}", text)])
    return toks


def _intent_keywords_list(rule: Dict[str, Any]) -> List[str]:
    kws = rule.get("intent_keywords", []) if isinstance(rule, dict) else []
    if isinstance(kws, dict):
        arr = kws.get("en", [])
        return arr if isinstance(arr, list) else []
    return kws if isinstance(kws, list) else []


def _intent_keyword_score(task_l: str, task_tokens: set, rule: Dict[str, Any]) -> int:
    score = 0
    for kw in _intent_keywords_list(rule):
        k = str(kw).strip().lower()
        if not k:
            continue
        if " " in k and _normalize_text(k) in _normalize_text(task_l):
            score += 2
        elif k in task_tokens or k in task_l:
            score += 1
    return score


def _is_function_discovery_query(text: str) -> bool:
    t = _normalize_text(text).lower()
    patterns = [
        r"有哪些\s*(function|functions|命令|功能)",
        r"which\s+(function|functions|command|commands)",
        r"what\s+can\s+i\s+use",
    ]
    return any(re.search(p, t) for p in patterns)


def classify_user_input(text: str, state: "SessionState") -> Dict[str, Any]:
    t = _normalize_text(text)
    tl = t.lower()
    if not t:
        return {"category": "unknown", "confidence": 0.2, "normalized_text": ""}

    if re.search(r"\b(restart|start over|new task|reset session)\b|重来|重新开始", tl):
        return {"category": "restart_session", "confidence": 0.99, "normalized_text": t}
    if _wants_undo(t):
        return {"category": "undo", "confidence": 0.99, "normalized_text": t}
    if _wants_reset_command(t):
        return {"category": "reset_command", "confidence": 0.99, "normalized_text": t}
    if _is_positive_confirm_text(t):
        return {"category": "confirm_yes", "confidence": 0.98, "normalized_text": t}
    if _is_negative_confirm_text(t):
        return {"category": "confirm_no", "confidence": 0.98, "normalized_text": t}

    greetings = ["你好", "您好", "hello", "hi", "hey", "早上好", "晚上好"]
    if any(g == tl for g in greetings) or any(tl.startswith(g + " ") for g in ["hello", "hi", "hey"]):
        return {"category": "greeting", "confidence": 0.95, "normalized_text": t}

    helps = ["你能做什么", "你可以做什么", "help", "what can you do", "有哪些功能", "有哪些命令"]
    if any(h in tl for h in helps):
        return {"category": "help", "confidence": 0.95, "normalized_text": t}

    chitchat = ["谢谢", "thanks", "thank you", "辛苦了"]
    if any(x in tl for x in chitchat):
        return {"category": "chit_chat", "confidence": 0.9, "normalized_text": t}

    switch_task_markers = ["runall", "trust4", "spechla", "hla typing", "hla_typing", "extract_hla_read"]
    if any(x in tl for x in switch_task_markers) and any(x in tl for x in ["改成", "改为", "换", "不做", "change", "switch"]):
        return {"category": "task_description", "confidence": 0.9, "normalized_text": t}

    execute_words = ["执行", "run", "execute", "开始", "现在开始", "确认运行", "现在运行", "确认执行"]
    if any(x == tl for x in execute_words):
        return {"category": "execute", "confidence": 0.95, "normalized_text": t}

    has_slot_tokens = bool(re.search(r"\b\w+\s*=\s*[^\s]+", t)) or bool(
        re.search(r"reset|clear|unset|remove|change\s+\w+|modify\s+\w+|update\s+\w+|set\s+\w+|删除\s*[A-Za-z_]|去掉\s*[A-Za-z_]|清空\s*[A-Za-z_]|重置\s*[A-Za-z_]|改成\s*(star|salmon)|改为\s*(star|salmon)", tl)
    )
    if has_slot_tokens:
        return {"category": "slot_update", "confidence": 0.85, "normalized_text": t}

    task_markers = [
        "分析", "analysis", "bulk", "rna", "fastq", "pipeline", "workflow", "tcr", "bcr", "hla",
        "runall", "trust4", "spechla", "hla_typing", "extract_hla_read",
    ]
    if any(x in tl for x in task_markers):
        return {"category": "task_description", "confidence": 0.8, "normalized_text": t}

    if state.subcommand:
        return {"category": "slot_update", "confidence": 0.55, "normalized_text": t}
    return {"category": "unknown", "confidence": 0.35, "normalized_text": t}


def _function_summary(rule: Dict[str, Any], prefer_chinese: bool) -> str:
    summary = rule.get("function_summary") if isinstance(rule, dict) else None
    if isinstance(summary, dict):
        txt = summary.get("zh" if prefer_chinese else "en")
        if isinstance(txt, str) and txt.strip():
            return txt.strip()
    if isinstance(summary, str) and summary.strip():
        return summary.strip()
    req = rule.get("required", []) if isinstance(rule.get("required", []), list) else []
    return (f"需要参数: {', '.join(req[:4]) if req else '无'}。" if prefer_chinese else f"Requires: {', '.join(req[:4]) if req else 'none'}.")


def _compose_function_suggestions(task: str, rules: Dict[str, Any], prefer_chinese: bool) -> str:
    ranked = _rank_subcommands(task, rules, top_n=5)
    if prefer_chinese:
        lines = ["你可以先从这些 function 里选择（按相关度排序）："]
    else:
        lines = ["You can choose from these functions first (ranked by relevance):"]
    for i, (cmd, _) in enumerate(ranked, 1):
        lines.append(f"{i}) {cmd}: {_function_summary(rules.get(cmd, {}), prefer_chinese)}")
    if prefer_chinese:
        lines.append("请回复 function 名称，我将继续补全参数。")
    else:
        lines.append("Reply with a function name, and I will continue parameter completion.")
    return "\n".join(lines)


def _rank_subcommands(task: str, rules: Dict[str, Any], top_n: Optional[int] = None):
    tl = _normalize_text(task).lower()
    bag = _intent_tokenize(tl)
    ranked: List[Tuple[Tuple[int, int, int], str]] = []
    for cmd, rule in rules.items():
        required = rule.get("required", []) if isinstance(rule, dict) else []
        req_tokens = set(re.findall(r"[a-z][a-z0-9]*", " ".join([str(x) for x in required]).lower()))
        name_tokens = set(re.findall(r"[a-z][a-z0-9]*", cmd.replace("_", " ").lower()))
        sc = (
            _intent_keyword_score(tl, bag, rule if isinstance(rule, dict) else {}),
            len(bag & req_tokens),
            len(bag & name_tokens),
        )
        ranked.append((sc, cmd))
    ranked.sort(key=lambda x: x[0], reverse=True)
    out = [(cmd, score) for score, cmd in ranked]
    return out[:top_n] if top_n else out


def choose_subcommand_with_confidence(task: str, rules: Dict[str, Any]) -> Dict[str, Any]:
    tl = _normalize_text(task).lower()
    trivial = ["你好", "您好", "hello", "hi", "hey", "thanks", "thank you", "谢谢"]
    if tl in trivial:
        return {
            "top1": None,
            "top1_score": 0.0,
            "top2": None,
            "top2_score": 0.0,
            "candidates": [],
            "is_clear": False,
        }

    # strong domain heuristics first
    heuristic_hits: List[Tuple[str, float]] = []
    heuristics = [
        ("runall", ["bulk rna", "rna-seq", "rnaseq", "fastq", "pipeline", "workflow", "bulk"], 0.86),
        ("trust4", ["tcr", "bcr", "vdj", "clonotype", "repertoire"], 0.9),
        ("spechla", ["spechla", "hla typing", "hla分型", "hla typing"], 0.88),
        ("hla_typing", ["hla_typing", "hla typing", "hla 分型"], 0.82),
    ]
    for cmd, keys, score in heuristics:
        if cmd not in rules:
            continue
        if any(k in tl for k in keys):
            heuristic_hits.append((cmd, score))
    if heuristic_hits:
        heuristic_hits = sorted(heuristic_hits, key=lambda x: x[1], reverse=True)
        top1 = heuristic_hits[0]
        top2 = heuristic_hits[1] if len(heuristic_hits) > 1 else (None, 0.0)
        clear = top1[1] >= 0.60 and (top1[1] - float(top2[1])) >= 0.15
        cands = [{"name": n, "score": float(s)} for n, s in heuristic_hits]
        return {
            "top1": top1[0] if clear else None,
            "top1_score": float(top1[1]),
            "top2": top2[0],
            "top2_score": float(top2[1]),
            "candidates": cands,
            "is_clear": clear,
        }

    ranked = _rank_subcommands(task, rules, top_n=5)
    if not ranked:
        return {
            "top1": None,
            "top1_score": 0.0,
            "top2": None,
            "top2_score": 0.0,
            "candidates": [],
            "is_clear": False,
        }

    def _score_norm(sc: Tuple[int, int, int]) -> float:
        raw = float((3 * sc[0]) + (2 * sc[1]) + sc[2])
        if raw <= 0:
            return 0.0
        return min(1.0, raw / (raw + 2.0))

    cands = []
    for cmd, sc in ranked:
        cands.append({"name": cmd, "score": _score_norm(sc)})

    top1 = cands[0] if cands else {"name": None, "score": 0.0}
    top2 = cands[1] if len(cands) > 1 else {"name": None, "score": 0.0}
    clear = bool(top1["name"]) and top1["score"] >= 0.60 and (top1["score"] - top2["score"]) >= 0.15

    return {
        "top1": top1["name"] if clear else None,
        "top1_score": float(top1["score"]),
        "top2": top2["name"],
        "top2_score": float(top2["score"]),
        "candidates": cands,
        "is_clear": clear,
    }


def choose_subcommand(task: str, rules: Dict[str, Any]) -> str:
    ranked = _rank_subcommands(task, rules, top_n=5)
    candidates = [c for c, _ in ranked]
    fallback = candidates[0] if candidates else list(rules.keys())[0]

    prompt = f"""
Select the best iobrpy subcommand from candidates: {candidates}
Return JSON only: {{"subcommand": "<candidate>"}}
User request: {task}
""".strip()
    try:
        j = _llm_generate_json(prompt)
    except Exception:
        j = {}
    sub = j.get("subcommand") if isinstance(j, dict) else None
    return sub if isinstance(sub, str) and sub in candidates else fallback


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


def _shell_quote(s: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./=-]+", s):
        return s
    return "'" + s.replace("'", "'\"'\"'") + "'"


def build_command(subcommand: str, params: Dict[str, Any], rules: Dict[str, Any]) -> Tuple[str, List[str]]:
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    argv = [subcommand]
    rule = rules.get(subcommand, {})
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


def compute_needs(subcommand: str, rules: Dict[str, Any], params: Dict[str, Any], confirmed: set) -> Tuple[List[str], List[str]]:
    rule = rules.get(subcommand, {})
    required = rule.get("required", [])
    confirm = rule.get("confirm", [])
    required_one_of = rule.get("required_one_of", [])
    missing = [k for k in required if k not in params or params[k] in (None, "", [])]
    if isinstance(required_one_of, list) and required_one_of:
        def _group_ok(group):
            return isinstance(group, list) and group and all((k in params and params[k] not in (None, "", [])) for k in group)
        if not any(_group_ok(g) for g in required_one_of):
            missing.append("__one_of_group")
    need_confirm = [k for k in confirm if k not in confirmed]
    return missing, need_confirm


def compose_questions(subcommand: str, missing: List[str], need_confirm: List[str], params: Dict[str, Any], rules: Dict[str, Any], prefer_chinese: bool = False) -> str:
    lines = [f"当前目标命令: iobrpy {subcommand}" if prefer_chinese else f"Current target command: iobrpy {subcommand}"]
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    notes = rule.get("notes", {}) if isinstance(rule.get("notes", {}), dict) else {}
    for m in missing:
        if m == "__one_of_group":
            default_msg = "请提供一种有效输入模式。" if prefer_chinese else "Please provide one valid input mode."
            msg = str(notes.get("required_one_of") or default_msg)
            if prefer_chinese and not _contains_cjk(msg):
                msg = _translate_to_chinese(msg)
            lines.append(msg)
        else:
            lines.append(f"- {m}（必填）" if prefer_chinese else f"- {m} (required)")
    for k in need_confirm:
        v = params.get(k)
        if v not in (None, "", []):
            if prefer_chinese:
                lines.append(f"请确认参数: {k}={v}（回复: confirm {k}）")
            else:
                lines.append(f"Please confirm parameter: {k}={v} (reply: confirm {k})")
    return "\n".join(lines)


def merge_defaults(params: Dict[str, Any], subcommand: str, rules: Dict[str, Any]) -> Dict[str, Any]:
    d = load_defaults()
    rule = rules.get(subcommand, {}) if isinstance(rules.get(subcommand, {}), dict) else {}
    rule_defaults = rule.get("optional_defaults", {}) if isinstance(rule.get("optional_defaults", {}), dict) else {}
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    out = {}
    for k in allowed_set:
        if k in params and params[k] in (None, "", []):
            continue
        if k in rule_defaults and rule_defaults[k] not in (None, "", []):
            out[k] = rule_defaults[k]
        elif k in BUILTIN_OPTIONAL_DEFAULTS:
            out[k] = BUILTIN_OPTIONAL_DEFAULTS[k]
        if k in d and d[k] not in (None, "", []):
            out[k] = d[k]
    for k, v in params.items():
        if k in allowed_set and v not in (None, "", []):
            out[k] = v
        elif k in allowed_set and v in (None, "", []) and k in out:
            out.pop(k, None)
    return out


class SessionState:
    def __init__(self):
        self.phase: str = "idle"
        self.task: str = ""
        self.subcommand: Optional[str] = None
        self.params: Dict[str, Any] = {}
        self.confirmed: set = set()
        self.prefer_chinese: bool = False
        self.history: List[Dict[str, Any]] = []
        self.pending_candidates: List[Dict[str, Any]] = []
        self.pending_switch_candidates: List[str] = []
        self.last_user_intent: Optional[str] = None
        self.subcommand_confidence: Optional[float] = None
        self.snapshots: List[Dict[str, Any]] = []
        self.pending_switch_subcommand: Optional[str] = None

    def push_snapshot(self):
        self.history.append({
            "phase": self.phase,
            "subcommand": self.subcommand,
            "params": deepcopy(self.params),
        })
        if len(self.history) > 80:
            self.history = self.history[-80:]
        self.snapshots.append({
            "phase": self.phase,
            "task": self.task,
            "subcommand": self.subcommand,
            "params": deepcopy(self.params),
            "confirmed": set(self.confirmed),
            "prefer_chinese": self.prefer_chinese,
            "pending_candidates": deepcopy(self.pending_candidates),
            "pending_switch_candidates": list(self.pending_switch_candidates),
            "last_user_intent": self.last_user_intent,
            "subcommand_confidence": self.subcommand_confidence,
            "pending_switch_subcommand": self.pending_switch_subcommand,
        })
        if len(self.snapshots) > 50:
            self.snapshots = self.snapshots[-50:]

    def rollback(self) -> bool:
        if not self.snapshots:
            return False
        snap = self.snapshots.pop()
        self.phase = snap.get("phase", "idle")
        self.task = snap.get("task", "")
        self.subcommand = snap.get("subcommand")
        self.params = deepcopy(snap.get("params", {}))
        self.confirmed = set(snap.get("confirmed", set()))
        self.prefer_chinese = bool(snap.get("prefer_chinese", self.prefer_chinese))
        self.pending_candidates = deepcopy(snap.get("pending_candidates", []))
        self.pending_switch_candidates = list(snap.get("pending_switch_candidates", []))
        self.last_user_intent = snap.get("last_user_intent")
        self.subcommand_confidence = snap.get("subcommand_confidence")
        self.pending_switch_subcommand = snap.get("pending_switch_subcommand")
        return True


SESSIONS: Dict[str, SessionState] = {}


def get_session(session_id: str) -> SessionState:
    if session_id not in SESSIONS:
        SESSIONS[session_id] = SessionState()
    return SESSIONS[session_id]


def apply_confirmations(state: SessionState, extracted: Dict[str, Any]):
    conf = extracted.get("confirm")
    if isinstance(conf, list):
        for k in conf:
            state.confirmed.add(str(k).strip())


def auto_confirm_filled_values(state: SessionState, extracted: Dict[str, Any], rules: Dict[str, Any]):
    if not state.subcommand:
        return
    rule = rules.get(state.subcommand, {}) if isinstance(rules.get(state.subcommand, {}), dict) else {}
    for k in rule.get("confirm", []) if isinstance(rule.get("confirm", []), list) else []:
        if k in extracted and extracted.get(k) not in (None, "", []):
            state.confirmed.add(str(k).strip())


def _detect_switch_candidate(text: str, current_subcommand: Optional[str], rules: Dict[str, Any]) -> Optional[str]:
    if not text or not current_subcommand:
        return None
    tl = _normalize_text(text).lower()
    # explicit phrasing: "switch/change to ..."
    mapped = [
        ("hla_typing", ["hla typing", "hla分型", "做 hla", "换 hla"]),
        ("spechla", ["spechla", "换 spechla", "做 spechla"]),
        ("trust4", ["trust4", "tcr", "bcr", "vdj"]),
        ("runall", ["runall", "bulk rna", "rna-seq", "workflow", "pipeline"]),
    ]
    for cmd, kws in mapped:
        if cmd == current_subcommand or cmd not in rules:
            continue
        if any(k in tl for k in kws):
            return cmd

    for cmd in rules.keys():
        if cmd == current_subcommand:
            continue
        if re.search(rf"\b{re.escape(cmd.lower())}\b", tl):
            return cmd
    # fallback: high-confidence re-routing signal
    sel = choose_subcommand_with_confidence(text, rules)
    top1 = sel.get("top1")
    if isinstance(top1, str) and top1 != current_subcommand:
        return top1
    return None


def _state_summary_text(state: SessionState, merged: Dict[str, Any], missing: List[str], prefer_chinese: bool) -> str:
    non_empty = []
    for k in sorted(merged.keys()):
        v = merged.get(k)
        if v in (None, "", []):
            continue
        tag = " (默认)" if k in BUILTIN_OPTIONAL_DEFAULTS and state.params.get(k) in (None, "", []) else ""
        non_empty.append(f"{k}={v}{tag}")
    recognized = ", ".join(non_empty) if non_empty else ("无" if prefer_chinese else "none")
    miss_show = ", ".join([m for m in missing if m != "__one_of_group"]) if missing else ("无" if prefer_chinese else "none")
    if prefer_chinese:
        return (
            f"当前命令: iobrpy {state.subcommand}\n"
            f"已识别: {recognized}\n"
            f"缺失: {miss_show}\n"
            "你可以继续补充参数，也可以说“把 threads 改成 16”“删除 project”“撤销上一条”“重新选命令”。"
        )
    return (
        f"Current command: iobrpy {state.subcommand}\n"
        f"Recognized: {recognized}\n"
        f"Missing: {miss_show}\n"
        "You can keep adding parameters, or say 'set threads to 16', 'remove project', 'undo', or 'change command'."
    )


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


def tool_iobrpy_assistant(session_id: str, task: Optional[str] = None, answer_text: Optional[str] = None, run: bool = False) -> Dict[str, Any]:
    rules = load_rules()
    state = get_session(session_id)
    phase_before = state.phase

    text = answer_text if answer_text is not None else task
    cls = classify_user_input(text or "", state)
    state.last_user_intent = cls.get("category")
    dlog("input_category=", cls.get("category"), "confidence=", cls.get("confidence"), "phase_before=", phase_before)

    if text:
        if _contains_cjk(text):
            state.prefer_chinese = True
        elif any(ch.isalpha() for ch in text):
            state.prefer_chinese = False

    if cls["category"] == "restart_session":
        state.push_snapshot()
        state.phase = "idle"
        state.task = ""
        state.subcommand = None
        state.params = {}
        state.confirmed = set()
        state.pending_candidates = []
        state.pending_switch_candidates = []
        state.pending_switch_subcommand = None
        state.subcommand_confidence = None
        q = "请用自然语言描述你的需求。" if state.prefer_chinese else "Please describe what you want to do (natural language)."
        return {"status": "need_info", "question": q, "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    if cls["category"] in ("greeting", "help", "chit_chat"):
        state.phase = "idle"
        state.subcommand = None
        state.pending_candidates = []
        if cls["category"] == "help":
            if state.prefer_chinese:
                q = (
                    "我可以帮助你运行 iobrpy 中的多个分析流程，例如 runall、trust4、spechla、hla_typing。"
                    "你可以直接告诉我你的任务，比如“我想做 bulk RNA-seq 分析”或“我想做 HLA typing”。"
                )
            else:
                q = (
                    "I can help run multiple iobrpy workflows, such as runall, trust4, spechla, and hla_typing. "
                    "You can directly describe your task, e.g. 'I want to run bulk RNA-seq analysis' or 'I want HLA typing'."
                )
        elif cls["category"] == "chit_chat":
            q = "不客气。你可以直接告诉我想做什么分析。" if state.prefer_chinese else "You're welcome. Tell me what analysis you want to run."
        else:
            q = "你好，我可以帮你选择合适的 iobrpy 功能，并一步步补全参数。你想做什么分析？" if state.prefer_chinese else "Hi, I can help choose the right iobrpy function and fill parameters step by step. What analysis do you want to run?"
        return {"status": "need_info", "question": q, "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    if cls["category"] == "undo":
        ok = state.rollback()
        msg = ("已撤销上一条。" if ok else "没有可撤销的历史。") if state.prefer_chinese else ("Undid last change." if ok else "No history to undo.")
        if not state.subcommand:
            state.phase = "idle"
            dlog("phase_after=", state.phase, "undo_ok=", ok)
            return {"status": "need_info", "question": msg, "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}
        merged_u = merge_defaults(state.params, state.subcommand, rules)
        missing_u, need_confirm_u = compute_needs(state.subcommand, rules, merged_u, state.confirmed)
        summary_u = _state_summary_text(state, merged_u, missing_u, state.prefer_chinese)
        state.phase = "parameter_filling" if (missing_u or need_confirm_u) else "ready_to_run"
        dlog("phase_after=", state.phase, "undo_ok=", ok)
        return {
            "status": "need_info" if state.phase == "parameter_filling" else "ready",
            "subcommand": state.subcommand,
            "question": msg + "\n\n" + summary_u,
            "draft_command": build_command(state.subcommand, merged_u, rules)[0],
            "params": merged_u,
            "prefer_chinese": state.prefer_chinese,
            "phase": state.phase,
            "intent_category": cls["category"],
        }

    if cls["category"] == "reset_command":
        state.push_snapshot()
        state.phase = "idle"
        state.subcommand = None
        state.subcommand_confidence = None
        state.pending_candidates = []
        state.pending_switch_candidates = []
        state.pending_switch_subcommand = None
        q = "已重置命令选择。请描述你要做什么。" if state.prefer_chinese else "Command selection reset. Please describe what you want to do."
        return {"status": "need_info", "question": q, "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    if state.pending_switch_subcommand and text:
        if cls["category"] == "confirm_yes":
            target = state.pending_switch_subcommand
            state.push_snapshot()
            old_params = dict(state.params)
            state.subcommand = target
            allowed = set([k for k in allowed_keys_for(target, rules) if k != "confirm"])
            state.params = {k: v for k, v in old_params.items() if k in allowed}
            state.confirmed = set([k for k in state.confirmed if k in allowed])
            state.pending_switch_subcommand = None
            state.pending_switch_candidates = []
            state.phase = "command_selected"
        elif cls["category"] == "confirm_no":
            state.pending_switch_subcommand = None
            state.pending_switch_candidates = []
            msg = "保持当前命令，继续补充参数。" if state.prefer_chinese else "Keep current command and continue filling parameters."
            return {"status": "need_info", "question": msg, "needs": [], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    # No selected command yet: only task_description can trigger routing
    if not state.subcommand:
        if cls["category"] == "task_description":
            sel = choose_subcommand_with_confidence(text or "", rules)
            dlog("subcommand_candidates=", sel.get("candidates"), "top1_score=", sel.get("top1_score"), "top2_score=", sel.get("top2_score"))
            state.pending_candidates = list(sel.get("candidates") or [])
            state.subcommand_confidence = sel.get("top1_score")
            if not sel.get("is_clear"):
                state.phase = "intent_clarification"
                cand_show = ", ".join([c.get("name") for c in state.pending_candidates[:4] if c.get("name")]) or "runall, trust4, spechla"
                if state.prefer_chinese:
                    q = (
                        "我需要再确认一下你的目标。你更接近下面哪一种？\n"
                        "1. bulk RNA-seq / 常规流程\n"
                        "2. HLA typing\n"
                        "3. TRUST4 相关分析\n"
                        f"4. 其他（可选候选: {cand_show}），请直接描述"
                    )
                else:
                    q = (
                        "I need one more clarification about your goal. Which is closest?\n"
                        "1. bulk RNA-seq / standard workflow\n"
                        "2. HLA typing\n"
                        "3. TRUST4-related analysis\n"
                        f"4. Other (candidate hints: {cand_show}), please describe"
                    )
                return {"status": "need_info", "question": q, "needs": ["clarify_intent"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "candidates": state.pending_candidates, "intent_category": cls["category"]}
            state.push_snapshot()
            state.task = text or ""
            state.subcommand = sel.get("top1")
            state.phase = "command_selected"
            state.params = {}
            state.confirmed = set()
            dlog("state transition -> command_selected", "subcommand=", state.subcommand)
        elif cls["category"] in ("slot_update", "execute", "unknown", "confirm_yes", "confirm_no"):
            state.phase = "idle"
            mkey = re.search(r"\b([A-Za-z_][A-Za-z0-9_]*)\b\s*=", text or "") if text else None
            key_hint = mkey.group(1) if mkey else None
            if state.prefer_chinese:
                q = (
                    f"请先告诉我你想运行哪个功能，我再帮你设置 {key_hint}。"
                    if key_hint
                    else "请先告诉我你要执行哪类分析任务。"
                )
            else:
                q = (
                    f"Please tell me which function you want to run first, then I can set {key_hint}."
                    if key_hint
                    else "Please first tell me what analysis task you want to run."
                )
            return {"status": "need_info", "question": q, "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    # lightweight switch detection only for task-like new description
    if answer_text and state.subcommand and cls["category"] == "task_description":
        switch_cand = _detect_switch_candidate(answer_text, state.subcommand, rules)
        if switch_cand:
            state.pending_switch_subcommand = switch_cand
            state.pending_switch_candidates = [switch_cand]
            dlog("switch_candidate=", switch_cand, "current=", state.subcommand)
            q = f"检测到你可能想切换到 {switch_cand}，是否切换？" if state.prefer_chinese else f"I detected you may want to switch to {switch_cand}. Switch now?"
            return {"status": "need_info", "question": q, "needs": ["confirm_switch"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    if not state.subcommand:
        state.phase = "idle"
        return {"status": "need_info", "question": ("请用自然语言描述你的需求。" if state.prefer_chinese else "Please describe what you want to do (natural language)."), "needs": ["task"], "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    state.phase = "parameter_filling"
    allowed = allowed_keys_for(state.subcommand, rules)
    parse_text = text or ""
    translated = _translate_to_english(parse_text) if parse_text else ""
    extracted = extract_slots_multisource(parse_text, translated, allowed) if parse_text else {}

    if extracted:
        state.push_snapshot()
    apply_confirmations(state, extracted)
    auto_confirm_filled_values(state, extracted, rules)

    resets = extracted.get("_reset")
    if isinstance(resets, list):
        for rk in resets:
            if rk in allowed:
                state.params[rk] = None
                state.confirmed.discard(rk)

    for k, v in extracted.items():
        if k in ("confirm", "_reset", "_undo", "_reset_command"):
            continue
        if k in allowed and v not in (None, "", []):
            state.params[k] = v
            state.confirmed.add(str(k).strip())

    merged = merge_defaults(state.params, state.subcommand, rules)
    missing, need_confirm = compute_needs(state.subcommand, rules, merged, state.confirmed)
    need_confirm = [k for k in need_confirm if state.params.get(k) not in (None, "", [])]

    draft_cmd, argv = build_command(state.subcommand, merged, rules)

    if cls["category"] == "execute" and (missing or need_confirm):
        q = "参数还不完整，暂时不能执行。" if state.prefer_chinese else "Parameters are still incomplete, cannot execute yet."
        q += "\n\n" + _state_summary_text(state, merged, missing, state.prefer_chinese)
        return {"status": "need_info", "subcommand": state.subcommand, "needs": missing + [f"confirm:{c}" for c in need_confirm], "question": q, "draft_command": draft_cmd, "params": merged, "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    if missing or need_confirm:
        q1 = compose_questions(state.subcommand, missing, need_confirm, merged, rules, prefer_chinese=state.prefer_chinese)
        q2 = _state_summary_text(state, merged, missing, state.prefer_chinese)
        return {"status": "need_info", "subcommand": state.subcommand, "needs": missing + [f"confirm:{c}" for c in need_confirm], "question": q1 + "\n\n" + q2, "draft_command": draft_cmd, "params": merged, "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    state.phase = "ready_to_run"
    dlog("phase_after=", state.phase, "missing=", missing, "need_confirm=", need_confirm)
    if not run:
        return {"status": "ready", "subcommand": state.subcommand, "draft_command": draft_cmd, "params": merged, "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}

    run_dir = os.getenv("IOBRPY_RUN_LOG_DIR", os.path.join(os.path.dirname(__file__), "mcp_runs"))
    os.makedirs(run_dir, exist_ok=True)
    log_path = os.path.join(run_dir, f"{session_id}_{state.subcommand}.log")
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
    return {"status": "done" if rc == 0 else "error", "returncode": rc, "draft_command": draft_cmd, "log_path": log_path, "tail": tail, "prefer_chinese": state.prefer_chinese, "phase": state.phase, "intent_category": cls["category"]}


TOOLS = [
    {
        "name": "iobrpy_assistant",
        "title": "iobrpy assistant (BYOK cloud LLM, strict params)",
        "description": "Natural language planner: choose iobrpy subcommand, fill allowed params only, ask missing, and run when ready.",
        "inputSchema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}, "task": {"type": "string"}, "answer_text": {"type": "string"}, "run": {"type": "boolean", "default": False}},
            "required": ["session_id"],
            "additionalProperties": False,
        },
    },
    {
        "name": "set_defaults",
        "title": "Set defaults (strict)",
        "description": "Persist defaults; only keys present in REQUIRED_PARAMS file are stored.",
        "inputSchema": {"type": "object", "properties": {"defaults": {"type": "object"}}, "required": ["defaults"], "additionalProperties": False},
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
        log("tools/call error:", repr(e))
        log(traceback.format_exc())
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
            log(traceback.format_exc())
            if isinstance(msg, dict) and msg.get("id") is not None:
                jsonrpc_error(msg["id"], -32603, "Internal error")


if __name__ == "__main__":
    main()
