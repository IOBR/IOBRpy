#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Agentic MCP stdio server (raw JSON-RPC) for iobrpy RAG + parameter completion.

V3 changes (per user request):
- STRICT: Only parameters that appear in REQUIRED_PARAMS json are ever extracted/asked/used.
- No guessing: values proposed by LLM are accepted only if they are explicitly evidenced in the user's text.
- Draft command only includes allowed params for the chosen subcommand (no cross-contamination like --input/--output for runall).
- Better path sanitization (avoid "/path.Please" etc).
"""

import os, sys, json, re, subprocess, traceback
from typing import Dict, Any, Optional, List, Tuple

import requests
import chromadb

CHROMA_DIR = os.getenv("CHROMA_DIR", "/mnt/tmp-ext4/RAG1/chroma")
CHROMA_COLLECTION = os.getenv("CHROMA_COLLECTION", "iobrpy_rag1")
OLLAMA_HOST = os.getenv("OLLAMA_HOST", "http://127.0.0.1:11434")
EMBED_MODEL = os.getenv("EMBED_MODEL", "qwen3-embedding:8b")
CHAT_MODEL = os.getenv("CHAT_MODEL", "qwen3:8b")

IOBRPY_CONDA_ENV = os.getenv("IOBRPY_CONDA_ENV", "iobrpy")
IOBRPY_SRC_DIR = os.getenv("IOBRPY_SRC_DIR", "/mnt/tmp-ext4/iobrpy/src")

REQUIRED_FILE = os.getenv("IOBRPY_REQUIRED_PARAMS_FILE", os.path.join(os.path.dirname(__file__), "iobrpy_required_params.json"))
DEFAULTS_FILE = os.getenv("IOBRPY_DEFAULTS_FILE", os.path.join(os.path.dirname(__file__), "iobrpy_defaults.json"))

CONDA_EXE = os.getenv("CONDA_EXE", "conda")

SUPPORTED_PROTOCOLS = ["2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"]

SERVER_INFO = {
    "name": "iobrpy-rag-agentic",
    "title": "iobrpy RAG Agentic (NL, strict params)",
    "version": "0.3.0",
    "description": "Natural language → choose iobrpy function → fill allowed params only → ask missing → run when ready.",
}

CAPABILITIES = {"tools": {"listChanged": False}}

def log(*a):
    print(*a, file=sys.stderr, flush=True)

def send(obj: Dict[str, Any]):
    sys.stdout.write(json.dumps(obj, ensure_ascii=False) + "\n")
    sys.stdout.flush()

def jsonrpc_error(_id, code, message, data=None):
    err = {"code": code, "message": message}
    if data is not None:
        err["data"] = data
    send({"jsonrpc": "2.0", "id": _id, "error": err})

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

def _chroma_collection():
    client = chromadb.PersistentClient(path=CHROMA_DIR)
    return client.get_or_create_collection(CHROMA_COLLECTION)

def _ollama_embed_one(text: str) -> List[float]:
    r = requests.post(f"{OLLAMA_HOST}/api/embed", json={"model": EMBED_MODEL, "input": text})
    if r.status_code == 200:
        j = r.json()
        if "embeddings" in j and isinstance(j["embeddings"], list) and len(j["embeddings"]) > 0:
            return j["embeddings"][0]
    r = requests.post(f"{OLLAMA_HOST}/api/embeddings", json={"model": EMBED_MODEL, "prompt": text})
    r.raise_for_status()
    return r.json()["embedding"]

def rag_search(query: str, top_k: int = 10) -> List[Dict[str, Any]]:
    col = _chroma_collection()
    qv = _ollama_embed_one(query)
    res = col.query(query_embeddings=[qv], n_results=top_k, include=["documents", "metadatas", "distances"])
    out = []
    docs = res.get("documents", [[]])[0]
    metas = res.get("metadatas", [[]])[0]
    dists = res.get("distances", [[]])[0]
    for doc, meta, dist in zip(docs, metas, dists):
        out.append({
            "distance": dist,
            "path": (meta or {}).get("path"),
            "start": (meta or {}).get("start"),
            "end": (meta or {}).get("end"),
            "text": doc,
        })
    return out

def _ollama_generate_json(prompt: str, model: str) -> Dict[str, Any]:
    payload = {"model": model, "prompt": prompt, "stream": False, "format": "json", "options": {"temperature": 0}}
    r = requests.post(f"{OLLAMA_HOST}/api/generate", json=payload)
    r.raise_for_status()
    txt = r.json().get("response", "").strip()
    try:
        return json.loads(txt) if txt else {}
    except Exception:
        m = re.search(r"\{.*\}", txt, flags=re.S)
        if m:
            try:
                return json.loads(m.group(0))
            except Exception:
                return {}
        return {}


def _contains_cjk(text: str) -> bool:
    if not text:
        return False
    return re.search(r"[㐀-䶿一-鿿豈-﫿]", text) is not None


def _translate_text(text: str, target_lang: str, key: str) -> str:
    if not text:
        return text
    prompt = f"""
Translate the following text to {target_lang}.
Keep file paths, numbers, parameter names, and flags exactly unchanged.
Return JSON only: {{{json.dumps(key)}: "..."}}

Text:
{text}
""".strip()
    try:
        j = _ollama_generate_json(prompt, CHAT_MODEL)
        out = j.get(key) if isinstance(j, dict) else None
        if isinstance(out, str) and out.strip():
            return out.strip()
    except Exception:
        pass
    return text


def _translate_to_english(text: str) -> str:
    """Best-effort translation to English for rule triggering; fallback to original text."""
    if not text or not _contains_cjk(text):
        return text
    return _translate_text(text, "English", "english")


def _translate_to_chinese(text: str) -> str:
    """Best-effort translation to Chinese for user-facing prompts; fallback to original text."""
    if not text:
        return text
    return _translate_text(text, "Chinese", "chinese")


def _resolve_iobrpy_pythonpath() -> Optional[str]:
    """Return a PYTHONPATH prefix that makes 'import iobrpy' work (best-effort)."""
    cand: List[str] = []
    if IOBRPY_SRC_DIR:
        cand.append(IOBRPY_SRC_DIR)
    # common fallback locations
    cand.extend(["/mnt/tmp-ext4/iobrpy/src", "/mnt/tmp-ext4/iobrpy"])
    for p in cand:
        try:
            if p and os.path.isdir(os.path.join(p, "iobrpy")):
                return p
        except Exception:
            pass
    return None

def _conda_run(args: List[str]) -> Tuple[int, str]:
    env = os.environ.copy()
    pp = _resolve_iobrpy_pythonpath()
    if pp:
        env["PYTHONPATH"] = pp + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    cmd = [CONDA_EXE, "run", "-n", IOBRPY_CONDA_ENV] + args
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, p.stdout

def iobrpy_help(subcommand: Optional[str] = None) -> str:
    base = ["python", "-m", "iobrpy.main"]
    if subcommand:
        base.append(subcommand)
    base.append("--help")
    _, out = _conda_run(base)
    return out[-30000:]

def _parse_help_options(help_text: str) -> set:
    return set(re.findall(r"--[A-Za-z0-9][A-Za-z0-9_-]*", help_text or ""))

def validate_command_options(subcommand: str, argv: List[str]) -> Dict[str, Any]:
    help_text = iobrpy_help(subcommand)
    # If help failed (e.g., iobrpy not importable in the target conda env), do not block with false 'unknown options'.
    if not help_text or ("usage" not in help_text.lower() and "options" not in help_text.lower()):
        return {"unknown_options": [], "help_unavailable": True, "help_excerpt": help_text[-1200:] if help_text else ""}
    opts = _parse_help_options(help_text)
    used = [a for a in argv if a.startswith("--")]
    unknown = [a for a in used if a not in opts]
    return {"unknown_options": unknown, "help_unavailable": False, "help_excerpt": help_text[-1200:]}

DEFAULT_RULES = {
  "runall": {
    "required": ["fastq", "outdir", "mode", "index", "threads", "batch_size", "project"],
    "confirm": ["mode", "threads", "batch_size"],
    "choices": {"mode": ["salmon", "star"]}
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
        if not isinstance(req, list): req = []
        if not isinstance(conf, list): conf = []
        if not isinstance(opt, list): opt = []
        if not isinstance(choices, dict): choices = {}
        if not isinstance(req_one, list): req_one = []
        if not isinstance(notes, dict): notes = {}
        out[k] = {"required": req, "confirm": conf, "optional": opt, "choices": choices, "required_one_of": req_one, "notes": notes}
    return out or DEFAULT_RULES

def allowed_keys_for(subcommand: str, rules: Dict[str, Any]) -> List[str]:
    rule = rules.get(subcommand, {})
    allowed = set(rule.get("required", []) + rule.get("optional", []))
    for k in (rule.get("choices", {}) or {}).keys():
        allowed.add(k)
    return sorted(list(allowed)) + ["confirm"]

def load_defaults() -> Dict[str, Any]:
    return _read_json(DEFAULTS_FILE, {})

def set_defaults(new_defaults: Dict[str, Any], rules: Dict[str, Any]) -> Dict[str, Any]:
    allowed_all = set()
    for sub in rules:
        for k in allowed_keys_for(sub, rules):
            if k != "confirm":
                allowed_all.add(k)
    cur = load_defaults()
    for k, v in (new_defaults or {}).items():
        if k in allowed_all:
            cur[k] = v
    _write_json(DEFAULTS_FILE, cur)
    return cur

PATH_KEYS = {"fastq", "outdir", "index", "input", "output", "bam", "r1", "r2", "ru", "fqdir", "o", "od"}
INT_KEYS = {"threads", "batch_size", "t", "k", "stage", "clean", "use_exon"}

def _sanitize_path(p: str) -> str:
    if not p:
        return p
    p = p.strip()
    p = re.sub(r"[,\.;:\)\]\}]+$", "", p)
    p = re.sub(r"\.?(?:Please|please)$", "", p).strip()
    p = p.strip('"').strip("'")
    return p

def _recover_full_path(short_path: str, text: str) -> str:
    """If model extracted a truncated absolute path like '/test/x', try to recover the full path from the user's text."""
    if not short_path or not short_path.startswith('/'):
        return short_path
    # collect absolute-path-like tokens from the original text
    cands = re.findall(r"(/[^ \t;,]+)", text)
    cands = [_sanitize_path(c) for c in cands]
    # prefer the longest candidate that endswith the short_path
    best = short_path
    for c in cands:
        if c == short_path:
            if len(c) > len(best):
                best = c
            continue
        if c.endswith(short_path) and len(c) > len(best):
            best = c
    return best


def _has_token(text_l: str, token: str) -> bool:
    return re.search(rf"\b{re.escape(token)}\b", text_l) is not None

def _regex_extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    t = _normalize_text(text)
    tl = t.lower()
    out: Dict[str, Any] = {}

    if "mode" in allowed:
        if _has_token(tl, "salmon"):
            out["mode"] = "salmon"
        elif _has_token(tl, "star"):
            out["mode"] = "star"

    if "index" in allowed:
        mi = re.search(r"\bindex\b\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if mi:
            out["index"] = _sanitize_path(mi.group(1))

    if "threads" in allowed:
        mt = re.search(r"\bthreads?\b\s*(?:is|=|:)?\s*(\d+)\b", tl)
        if mt:
            out["threads"] = int(mt.group(1))

    if "batch_size" in allowed:
        mb = re.search(r"\b(batch[_\s]?size|batchsize|parallel(?:\s*samples?)?)\b\s*(?:is|=|:)?\s*(\d+)\b", tl)
        if mb:
            out["batch_size"] = int(mb.group(2))

    if "project" in allowed:
        mp = re.search(r"\bproject\b\s*(?:is|=|:)?\s*([A-Za-z0-9_.-]+)\b", t, flags=re.I)
        if mp:
            val = mp.group(1).strip()
            if val.lower() != "is":
                out["project"] = val

    if "fastq" in allowed:
        # Accept: "FASTQ files located at /path", "fastq at /path", "fastq /path"
        mfq = re.search(
            r"\bfastq\b(?:\s+files?)?(?:\s+(?:located\s+at|located\s+in|at|in))?\s+(/[^ \t;,]+)",
            t,
            flags=re.I,
        )
        if mfq:
            out["fastq"] = _sanitize_path(mfq.group(1))

    if "outdir" in allowed:
        # Accept: "outdir /path", "outputs in /path", "put outputs in /path", "to /path"
        mout = re.search(
            r"\b(outdir|output\s*dir(?:ectory)?|outputs?\s*in|put\s*outputs?\s*in|to)\b\s*(?:is|=|:|to|in)?\s*(/[^ \t;,]+)",
            t,
            flags=re.I,
        )
        if mout:
            out["outdir"] = _sanitize_path(mout.group(2))

    if "bam" in allowed:
        mb = re.search(r"(?:\b-b\b|\bbam\b)\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if mb:
            out["bam"] = _sanitize_path(mb.group(1))

    if "fqdir" in allowed:
        mf = re.search(r"(?:\b--fqdir\b|\bfqdir\b|\bfastq\s*dir(?:ectory)?\b)\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if mf:
            out["fqdir"] = _sanitize_path(mf.group(1))

    if "ru" in allowed:
        mu = re.search(r"(?:\b-u\b|\bsingle(?:-end)?\b|\bru\b)\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if mu:
            out["ru"] = _sanitize_path(mu.group(1))

    if "r1" in allowed:
        m1 = re.search(r"(?:\b-1\b|\bread1\b|\br1\b)\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if m1:
            out["r1"] = _sanitize_path(m1.group(1))

    if "r2" in allowed:
        m2 = re.search(r"(?:\b-2\b|\bread2\b|\br2\b)\s*(?:is|=|:)?\s*(/[^ \t;,]+)", t, flags=re.I)
        if m2:
            out["r2"] = _sanitize_path(m2.group(1))

    if "confirm" in allowed:
        confirms = re.findall(r"\bconfirm\s+([A-Za-z0-9_]+)\b", tl)
        if confirms:
            out["confirm"] = sorted(set(confirms))

    # reset/unset keys: e.g. 'reset fastq', 'clear index', 'unset threads'
    m_reset = re.findall(r"\b(?:reset|clear|unset|remove)\s+([A-Za-z0-9_]+)\b", tl)
    if m_reset:
        out["_reset"] = sorted(set([k for k in m_reset if k in allowed]))

    return out

def _llm_extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    prompt = f"""
You are a strict information extraction system.
Extract ONLY fields that are explicitly present in the user's message. Do NOT guess.
Return JSON only.

Allowed keys: {allowed}

Conventions:
- mode must be exactly \"salmon\" or \"star\" if present.
- Path keys must be absolute Linux paths starting with \"/\".
- threads and batch_size must be integers if present.
- project must be an identifier like \"PRJNA320473\".
- confirm must be an array of keys explicitly confirmed (e.g. user says \"confirm threads\").

User message:
{text}
""".strip()
    j = _ollama_generate_json(prompt, CHAT_MODEL)
    return j if isinstance(j, dict) else {}

def _evidence_filter(text: str, extracted: Dict[str, Any], allowed: List[str]) -> Dict[str, Any]:
    t = _normalize_text(text)
    tl = t.lower()
    out: Dict[str, Any] = {}
    for k, v in extracted.items():
        if k == "confirm":
            if isinstance(v, list):
                out["confirm"] = [str(x) for x in v]
            continue
        if k == "_reset":
            if isinstance(v, list):
                out["_reset"] = [str(x) for x in v if str(x) in allowed]
            continue
        if k not in allowed:
            continue
        if k == "confirm":
            if isinstance(v, list):
                out["confirm"] = [str(x) for x in v]
            continue
        if k in PATH_KEYS:
            if isinstance(v, str) and v.startswith("/"):
                sv = _sanitize_path(v)
                if sv in t:
                    out[k] = sv
            continue
        if k == "mode":
            if isinstance(v, str) and v in ("salmon", "star") and _has_token(tl, v):
                out["mode"] = v
            continue
        if k in INT_KEYS:
            try:
                iv = int(v)
            except Exception:
                continue
            if k == "threads":
                if re.search(r"\bthreads?\b", tl) and re.search(rf"\b{iv}\b", tl):
                    out["threads"] = iv
            if k == "batch_size":
                if re.search(r"batch|parallel", tl) and re.search(rf"\b{iv}\b", tl):
                    out["batch_size"] = iv
            continue
        if k == "project":
            if isinstance(v, str) and re.fullmatch(r"[A-Za-z0-9_.-]+", v):
                if v in t and re.search(r"\bproject\b", tl):
                    out["project"] = v
            continue
    return out

def extract_slots(text: str, allowed: List[str]) -> Dict[str, Any]:
    r = _regex_extract_slots(text, allowed)
    missing = [k for k in allowed if k not in r and k != "confirm"]
    if missing:
        try:
            llm = _llm_extract_slots(text, allowed)
            llm2 = _evidence_filter(text, llm, allowed)
        except Exception as e:
            log("LLM extraction failed:", e)
            llm2 = {}
        for k in missing:
            if k in llm2 and llm2[k] not in (None, "", [], {}):
                r[k] = llm2[k]
        if "confirm" not in r and "confirm" in llm2:
            r["confirm"] = llm2["confirm"]
    return r

class SessionState:
    def __init__(self):
        self.task = ""
        self.subcommand = None
        self.params: Dict[str, Any] = {}
        self.confirmed = set()
        self.prefer_chinese = False

SESSIONS: Dict[str, SessionState] = {}

def get_session(session_id: str) -> SessionState:
    if session_id not in SESSIONS:
        SESSIONS[session_id] = SessionState()
    return SESSIONS[session_id]

def _intent_tokenize(text: str) -> set:
    """Tokenize user intent while avoiding noisy path fragments (e.g., '/.../star_2sample')."""
    if not text:
        return set()
    t = _normalize_text(text).lower()
    t = re.sub(r"/[A-Za-z0-9._/-]+", " ", t)
    tokens = re.findall(r"[a-z][a-z0-9]*", t.replace('_', ' '))
    return set(tokens)


def choose_subcommand(task: str, rules: Dict[str, Any]) -> str:
    tl = _normalize_text(task).lower()
    tools = list(rules.keys())

    # Build tool profile from rule metadata so the model can reason even when RAG is weak.
    tool_profiles = []
    for t in tools:
        r = rules.get(t, {}) if isinstance(rules.get(t, {}), dict) else {}
        tool_profiles.append({
            "name": t,
            "required": r.get("required", []),
            "required_one_of": r.get("required_one_of", []),
            "notes": r.get("notes", {}),
        })

    # Use both original query and English translation to improve retrieval/selection quality.
    q_en = _translate_to_english(task)
    query_for_rag = q_en if isinstance(q_en, str) and q_en.strip() else task

    try:
        ctx = rag_search(query_for_rag, top_k=10)
    except Exception:
        ctx = []

    ctx_text = "\n\n".join([
        f"[{i}] {c.get('path')}:{c.get('start')}-{c.get('end')}\n{c.get('text','')[:700]}"
        for i, c in enumerate(ctx)
    ])

    prompt = f"""
You are selecting the best iobrpy subcommand for the user intent.
Pick exactly ONE subcommand name from this list: {tools}

Tool profiles (from machine rules):
{json.dumps(tool_profiles, ensure_ascii=False)}

User request (original):
{task}

User request (english translation for intent matching):
{q_en}

Helpful retrieved docs/snippets:
{ctx_text}

Return JSON only: {{"subcommand": "<one of the list>", "reason": "<short>"}}
""".strip()

    try:
        j = _ollama_generate_json(prompt, CHAT_MODEL)
    except Exception:
        j = {}

    sub = j.get("subcommand") if isinstance(j, dict) else None
    if isinstance(sub, str) and sub in tools:
        return sub

    # Last-resort fallback: lexical overlap against tool names and notes (not hardcoded routing).
    bag = _intent_tokenize(tl)
    best = tools[0]
    best_score = (-1, -1, -1)
    has_abs_path = re.search(r"/[A-Za-z0-9._/-]+", _normalize_text(task)) is not None
    for t in tools:
        r = rules.get(t, {}) if isinstance(rules.get(t, {}), dict) else {}
        name_tokens = set(re.findall(r"[a-z][a-z0-9]*", t.lower().replace('_', ' ')))
        txt = " ".join([
            " ".join(r.get("required", []) if isinstance(r.get("required", []), list) else []),
            json.dumps(r.get("notes", {}), ensure_ascii=False),
        ]).lower()
        cand = set(re.findall(r"[a-z][a-z0-9]*", txt.replace('_', ' ')))
        required_keys = r.get("required", []) if isinstance(r.get("required", []), list) else []
        dir_bonus = 1 if has_abs_path and any(str(k).endswith("_dir") for k in required_keys) else 0
        score = (len(bag & name_tokens), len(bag & (cand | name_tokens)), dir_bonus)
        if score > best_score:
            best_score = score
            best = t

    return best
FLAG_MAP = {
    "fastq": "--fastq",
    "outdir": "--outdir",
    "mode": "--mode",
    "index": "--index",
    "threads": "--threads",
    "batch_size": "--batch_size",
    "project": "--project",
    "input": "--input",
    "output": "--output",
    "bam": "-b",
    "r1": "-1",
    "r2": "-2",
    "ru": "-u",
    "fqdir": "--fqdir",
    "f": "-f",
    "ref": "--ref",
    "o": "-o",
    "od": "--od",
    "t": "-t",
    "k": "-k",
}
BOOL_FLAGS = {
    "repseq": "--repseq",
    "skipMateExtension": "--skipMateExtension",
    "abnormalUnmapFlag": "--abnormalUnmapFlag",
    "assembleWithRef": "--assembleWithRef",
    "noExtraction": "--noExtraction",
    "outputReadAssignment": "--outputReadAssignment",
}

def _shell_quote(s: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./=-]+", s):
        return s
    # POSIX-safe single-quote escaping: close, escape, reopen
    return "'" + s.replace("'", "'\"'\"'") + "'"

def build_command(subcommand: str, params: Dict[str, Any], rules: Dict[str, Any]) -> Tuple[str, List[str]]:
    allowed = allowed_keys_for(subcommand, rules)
    allowed_set = set([k for k in allowed if k != "confirm"])
    argv = [subcommand]
    rule = rules.get(subcommand, {})
    ordered = list(rule.get("required", [])) + [k for k in rule.get("optional", []) if k not in rule.get("required", [])]
    for k in ordered:
        if k not in allowed_set:
            continue
        if k not in params:
            continue
        v = params[k]
        if v in (None, "", []):
            continue
        if isinstance(v, bool):
            bflag = BOOL_FLAGS.get(k)
            if bflag and v:
                argv.append(bflag)
            continue
        flag = FLAG_MAP.get(k)
        if flag:
            argv.extend([flag, str(v)])
    cmd = "iobrpy " + " ".join([_shell_quote(a) for a in argv])
    return cmd, argv

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
    lines = []
    rule = rules.get(subcommand, {})
    notes = rule.get("notes", {}) if isinstance(rule, dict) else {}

    for m in missing:
        if m == "__one_of_group":
            msg = notes.get("required_one_of") or notes.get("input_mode")
            lines.append((_translate_to_chinese(msg) if prefer_chinese else msg) if isinstance(msg, str) and msg.strip() else ("请提供一种有效输入模式。" if prefer_chinese else "Please provide one valid input mode."))
            continue
        msg = notes.get(m) if isinstance(notes, dict) else None
        if isinstance(msg, str) and msg.strip():
            lines.append(_translate_to_chinese(msg) if prefer_chinese else msg)
        else:
            lines.append((f"请提供 {m}。") if prefer_chinese else (f"Please provide {m}."))

    for k in need_confirm:
        v = params.get(k)
        if v not in (None, "", []):
            lines.append((f"请确认 {k}={v}。可回复：confirm {k}") if prefer_chinese else (f"Please confirm {k}={v}. Reply like: confirm {k}"))
    return "\n".join(lines) if lines else ("请补充缺失参数。" if prefer_chinese else "Please provide the missing parameters.")

def merge_defaults(params: Dict[str, Any], subcommand: str, rules: Dict[str, Any]) -> Dict[str, Any]:
    d = load_defaults()
    allowed_set = set([k for k in allowed_keys_for(subcommand, rules) if k != "confirm"])
    out = {}
    for k in allowed_set:
        if k in d and d[k] not in (None, "", []):
            out[k] = d[k]
    for k, v in params.items():
        if k in allowed_set and v not in (None, "", []):
            out[k] = v
    return out

def apply_confirmations(state: SessionState, extracted: Dict[str, Any]):
    conf = extracted.get("confirm")
    if isinstance(conf, list):
        for k in conf:
            state.confirmed.add(str(k).strip())

def tool_iobrpy_assistant(session_id: str, task: Optional[str] = None, answer_text: Optional[str] = None, run: bool = False) -> Dict[str, Any]:
    rules = load_rules()
    state = get_session(session_id)

    task_en = _translate_to_english(task) if task else None
    answer_text_en = _translate_to_english(answer_text) if answer_text else None

    if (task and _contains_cjk(task)) or (answer_text and _contains_cjk(answer_text)):
        state.prefer_chinese = True

    # Restart session on user request (match both original text and translated text)
    if answer_text:
        atl = _normalize_text(answer_text).lower()
        atl_en = _normalize_text(answer_text_en or "").lower()
        if re.search(r"\b(restart|start over|new task|reset session)\b", atl) or re.search(r"\b(restart|start over|new task|reset session)\b", atl_en):
            state.task = ""
            state.subcommand = None
            state.params = {}
            state.confirmed = set()
            state.prefer_chinese = False

    if task:
        parse_task = task_en or task
        state.task = parse_task
        state.subcommand = choose_subcommand(parse_task, rules)
        state.params = {}
        state.confirmed = set()

        allowed = allowed_keys_for(state.subcommand, rules)
        extracted = extract_slots(parse_task, allowed)
        apply_confirmations(state, extracted)
        # Apply reset/unset requests
        resets = extracted.get("_reset")
        if isinstance(resets, list):
            for rk in resets:
                state.params.pop(rk, None)
                state.confirmed.discard(rk)
        for k, v in extracted.items():
            if k != "confirm" and k in allowed and v not in (None, "", []):
                state.params[k] = v

    if answer_text:
        # Treat user reply as additional information for the current session.
        # If the user skipped the initial Task prompt, accept the first reply as the task.
        if not state.subcommand:
            parse_answer = answer_text_en or answer_text
            state.task = parse_answer
            state.subcommand = choose_subcommand(parse_answer, rules)
            state.params = {}
            state.confirmed = set()

        allowed = allowed_keys_for(state.subcommand, rules)
        parse_answer = answer_text_en or answer_text
        extracted = extract_slots(parse_answer, allowed)
        apply_confirmations(state, extracted)

        # Apply reset/unset requests (natural language): 'reset fastq', 'clear index', 'unset threads'
        resets = extracted.get("_reset")
        if isinstance(resets, list):
            for rk in resets:
                if rk in allowed:
                    state.params.pop(rk, None)
                    state.confirmed.discard(rk)

        # Merge extracted values (only allowed keys; no guessing)
        for k, v in extracted.items():
            if k in ("confirm", "_reset"):
                continue
            if k in allowed and v not in (None, "", []):
                state.params[k] = v


    if not state.subcommand:
        return {"status": "need_info", "question": ("请用自然语言描述你的需求。" if state.prefer_chinese else "Please describe what you want to do (natural language)."), "needs": ["task"]}

    merged = merge_defaults(state.params, state.subcommand, rules)
    missing, need_confirm = compute_needs(state.subcommand, rules, merged, state.confirmed)
    draft_cmd, argv = build_command(state.subcommand, merged, rules)

    # Optional validation against `iobrpy <cmd> --help` (DISABLED by default).
    # JSON rules are the source of truth; enable only if you really want strict CLI verification.
    if os.environ.get("IOBRPY_VALIDATE_HELP", "0") == "1" and not missing:
        try:
            val = validate_command_options(state.subcommand, argv)
        except Exception as e:
            val = {"unknown_options": [], "warning": f"help validation failed: {e}"}
        if val.get("unknown_options"):
            q = "I found options that do not appear in iobrpy --help for this subcommand: " + ", ".join(val["unknown_options"]) + ". Please revise your request."
            return {"status": "need_info", "question": q, "needs": ["revise_options"], "draft_command": draft_cmd, "validation": val}

    if missing or need_confirm:
        q = compose_questions(state.subcommand, missing, need_confirm, merged, rules, prefer_chinese=state.prefer_chinese)
        return {"status": "need_info", "subcommand": state.subcommand, "needs": missing + [f"confirm:{c}" for c in need_confirm],
                "question": q, "draft_command": draft_cmd, "params": merged}

    if not run:
        return {"status": "ready", "subcommand": state.subcommand, "draft_command": draft_cmd, "params": merged}

    run_dir = os.getenv("IOBRPY_RUN_LOG_DIR", os.path.join(os.path.dirname(__file__), "mcp_runs"))
    os.makedirs(run_dir, exist_ok=True)
    log_path = os.path.join(run_dir, f"{session_id}_{state.subcommand}.log")

    env = os.environ.copy()
    pp = _resolve_iobrpy_pythonpath()
    if pp:
        env["PYTHONPATH"] = pp + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")

    cmd1 = [CONDA_EXE, "run", "-n", IOBRPY_CONDA_ENV, "iobrpy"] + argv
    cmd2 = [CONDA_EXE, "run", "-n", IOBRPY_CONDA_ENV, "python", "-m", "iobrpy.main"] + argv

    def _run_cmd(cmd, mode="w"):
        with open(log_path, mode, encoding="utf-8") as f:
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, env=env)
            return p.wait()

    rc = None
    try:
        rc = _run_cmd(cmd1, "w")
        if rc != 0:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write("\n\n--- entrypoint failed, trying python -m iobrpy.main ---\n\n")
            rc = _run_cmd(cmd2, "a")
    except Exception as e:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write("\n\nEXCEPTION:\n" + repr(e) + "\n" + traceback.format_exc())

    tail = ""
    try:
        with open(log_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
            tail = "".join(lines[-80:])
    except Exception:
        pass

    return {"status": "done" if rc == 0 else "error", "returncode": rc, "draft_command": draft_cmd, "log_path": log_path, "tail": tail}

TOOLS = [
    {
        "name": "iobrpy_assistant",
        "title": "iobrpy assistant (agentic NL, strict params)",
        "description": "Natural language planner: choose iobrpy subcommand, fill allowed params only, ask missing, and run when ready.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "task": {"type": "string"},
                "answer_text": {"type": "string"},
                "run": {"type": "boolean", "default": False}
            },
            "required": ["session_id"],
            "additionalProperties": False
        }
    },
    {
        "name": "rag_search",
        "title": "RAG search",
        "description": "Search iobrpy code/docs chunks using embeddings + Chroma.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "query": {"type": "string"},
                "top_k": {"type": "integer", "default": 10, "minimum": 1, "maximum": 50}
            },
            "required": ["query"],
            "additionalProperties": False
        }
    },
    {
        "name": "set_defaults",
        "title": "Set defaults (strict)",
        "description": "Persist defaults; only keys present in REQUIRED_PARAMS file are stored.",
        "inputSchema": {
            "type": "object",
            "properties": {"defaults": {"type": "object"}},
            "required": ["defaults"],
            "additionalProperties": False
        }
    },
    {
        "name": "get_defaults",
        "title": "Get defaults",
        "description": "Read stored defaults.",
        "inputSchema": {"type": "object", "properties": {}, "additionalProperties": False}
    }
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
        if name == "rag_search":
            q = args["query"]
            top_k = int(args.get("top_k", 10))
            data = rag_search(q, top_k=top_k)
            send({"jsonrpc":"2.0","id":_id,"result":{"content":[{"type":"text","text":json.dumps(data,ensure_ascii=False)}],"isError":False}})
            return
        if name == "set_defaults":
            rules = load_rules()
            d = args.get("defaults") or {}
            out = set_defaults(d, rules)
            send({"jsonrpc":"2.0","id":_id,"result":{"content":[{"type":"text","text":json.dumps(out,ensure_ascii=False)}],"isError":False}})
            return
        if name == "get_defaults":
            out = load_defaults()
            send({"jsonrpc":"2.0","id":_id,"result":{"content":[{"type":"text","text":json.dumps(out,ensure_ascii=False)}],"isError":False}})
            return
        if name == "iobrpy_assistant":
            sid = args.get("session_id")
            if not sid:
                raise ValueError("session_id is required")
            out = tool_iobrpy_assistant(sid, task=args.get("task"), answer_text=args.get("answer_text"), run=bool(args.get("run", False)))
            send({"jsonrpc":"2.0","id":_id,"result":{"content":[{"type":"text","text":json.dumps(out,ensure_ascii=False)}],"isError":False}})
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
        except Exception as e:
            log("Server error:", repr(e))
            log(traceback.format_exc())
            if isinstance(msg, dict) and msg.get("id") is not None:
                jsonrpc_error(msg["id"], -32603, "Internal error")

if __name__ == "__main__":
    main()
