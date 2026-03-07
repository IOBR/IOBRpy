#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""iobrpy ai: BYOK cloud LLM planner + local iobrpy execution."""

from __future__ import annotations

import argparse
import getpass
import os
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from prompt_toolkit import PromptSession  # type: ignore
    from prompt_toolkit.history import InMemoryHistory  # type: ignore
except Exception:
    PromptSession = None
    InMemoryHistory = None

try:
    import readline  # noqa: F401
except Exception:
    readline = None  # type: ignore


@dataclass
class ProviderProfile:
    env_key: str
    default_model: str
    base_url: str


PROVIDERS: Dict[str, ProviderProfile] = {
    "qwen": ProviderProfile("DASHSCOPE_API_KEY", "qwen-plus", "https://dashscope.aliyuncs.com/compatible-mode/v1"),
    "kimi": ProviderProfile("MOONSHOT_API_KEY", "moonshot-v1-8k", "https://api.moonshot.cn/v1"),
    "deepseek": ProviderProfile("DEEPSEEK_API_KEY", "deepseek-chat", "https://api.deepseek.com/v1"),
    "glm": ProviderProfile("ZHIPUAI_API_KEY", "glm-4-plus", "https://open.bigmodel.cn/api/paas/v4"),
    "openai": ProviderProfile("OPENAI_API_KEY", "gpt-4o-mini", "https://api.openai.com/v1"),
    "claude": ProviderProfile("ANTHROPIC_API_KEY", "claude-3-5-sonnet-20241022", "https://api.anthropic.com/v1"),
    "gemini": ProviderProfile("GEMINI_API_KEY", "gemini-1.5-pro", "https://generativelanguage.googleapis.com/v1beta"),
}


@dataclass
class AIConfig:
    llm_alias: str
    api_key: str
    model: str
    base_url: str
    logdir: str


def _pkg_dir() -> Path:
    return Path(__file__).resolve().parent


def _required_params_file() -> Path:
    return _pkg_dir() / "iobrpy_required_params.json"


def _new_session_id() -> str:
    return f"s{int(time.time())}_{os.getpid()}_{uuid.uuid4().hex[:6]}"


def _contains_cjk(text: str) -> bool:
    return any("\u4e00" <= ch <= "\u9fff" for ch in (text or ""))


def _ui_text(key: str, prefer_chinese: bool, **kwargs: Any) -> str:
    zh = {
        "welcome": "IOBRpy AI（BYOK 云端 LLM）",
        "type_request": "请输入自然语言需求。",
        "commands": "命令: :exit  :quit  :restart",
        "bye": "再见",
        "ready": "就绪。草稿命令: {draft}",
        "confirm_run": "现在执行这个命令吗？[y/N] ",
        "cancelled": "已取消执行。你可以继续修改参数，或输入“执行”来运行。",
        "done": "[完成] rc={rc}",
        "error": "[失败] rc={rc}",
    }
    en = {
        "welcome": "IOBRpy AI (BYOK cloud LLM)",
        "type_request": "Type your request in natural language.",
        "commands": "Commands: :exit  :quit  :restart",
        "bye": "bye",
        "ready": "Ready. Draft: {draft}",
        "confirm_run": "Execute this command now? [y/N] ",
        "cancelled": "Execution cancelled. You can keep editing parameters, or reply 'run' to execute.",
        "done": "[done] rc={rc}",
        "error": "[error] rc={rc}",
    }
    t = (zh if prefer_chinese else en).get(key, key)
    return t.format(**kwargs)


def _resolve_api_key(alias: str, api_key_arg: Optional[str]) -> str:
    if api_key_arg:
        return api_key_arg
    profile = PROVIDERS[alias]
    env_val = os.getenv(profile.env_key, "").strip()
    if env_val:
        return env_val
    key = getpass.getpass(f"Enter API key for {alias} ({profile.env_key}): ").strip()
    if not key:
        raise ValueError(f"API key is required for --llm {alias}.")
    return key


_PROMPT_SESSION: Optional[Any] = None


def _init_prompt_session() -> Optional[Any]:
    global _PROMPT_SESSION
    if _PROMPT_SESSION is not None:
        return _PROMPT_SESSION
    if PromptSession is None or InMemoryHistory is None:
        return None
    try:
        _PROMPT_SESSION = PromptSession(history=InMemoryHistory(), multiline=False)
    except Exception:
        _PROMPT_SESSION = None
    return _PROMPT_SESSION


def read_user_input(prompt_text: str) -> str:
    session = _init_prompt_session()
    if session is not None:
        try:
            return session.prompt(prompt_text)
        except KeyboardInterrupt:
            raise
        except EOFError:
            raise
        except Exception:
            pass
    try:
        return input(prompt_text)
    except KeyboardInterrupt:
        raise
    except EOFError:
        raise


def _is_positive_confirm(text: str) -> bool:
    t = (text or "").strip().lower()
    return t in {"y", "yes", "confirm", "run", "execute", "ok", "确认", "执行", "是"}


def _is_negative_confirm(text: str) -> bool:
    t = (text or "").strip().lower()
    return t in {"n", "no", "cancel", "stop", "取消", "不执行", "否"}




def should_show_draft(assistant_result: Dict[str, Any]) -> bool:
    phase = str(assistant_result.get("phase") or "")
    category = str(assistant_result.get("intent_category") or "")
    if category in {"greeting", "help", "chit_chat"}:
        return False
    if phase in {"idle", "intent_clarification"}:
        return False
    if phase == "command_selected":
        return False
    if phase in {"parameter_filling", "ready_to_run"}:
        return bool(assistant_result.get("subcommand") and assistant_result.get("draft_command"))
    return False


def _should_show_state_summary(assistant_result: Dict[str, Any]) -> bool:
    phase = str(assistant_result.get("phase") or "")
    return phase in {"parameter_filling", "ready_to_run"}

def _print_state(obj: Dict[str, Any], prefer_chinese: bool) -> None:
    status = obj.get("status")
    phase = str(obj.get("phase") or "")

    if status == "need_info":
        q = obj.get("question") or ""
        if q:
            print(q)
        if phase == "command_selected" and obj.get("subcommand"):
            if prefer_chinese:
                print(f"已识别命令: {obj.get('subcommand')}")
            else:
                print(f"Recognized command: {obj.get('subcommand')}")
        if should_show_draft(obj):
            print(f"\nDraft: {obj.get('draft_command')}")
        return

    if status == "ready":
        print("\n" + _ui_text("ready", prefer_chinese, draft=obj.get("draft_command")))
        return

    if status in ("done", "error"):
        print("\n" + _ui_text("done" if status == "done" else "error", prefer_chinese, rc=obj.get("returncode")))
        if obj.get("draft_command"):
            print(f"Cmd : {obj.get('draft_command')}")
        if obj.get("log_path"):
            print(f"Log : {obj.get('log_path')}")
        tail = obj.get("tail")
        if tail:
            print("\n--- log tail ---")
            print(tail, end="" if str(tail).endswith("\n") else "\n")
        return

    print(obj)


def _run_iobrpy_current_env(*, session_id: str, subcommand: str, params: Dict[str, Any], server: Any, logdir: Path) -> Dict[str, Any]:
    rules = server.load_rules()
    draft_cmd, argv = server.build_command(subcommand, params, rules)
    log_path = logdir / f"{session_id}_{subcommand}.log"
    cmd = [sys.executable, "-m", "iobrpy.main"] + argv

    try:
        with log_path.open("w", encoding="utf-8") as f:
            p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, text=True, env=os.environ.copy())
            rc = p.wait()
    except Exception as e:
        with log_path.open("a", encoding="utf-8") as f:
            f.write("\n\nEXCEPTION:\n" + repr(e) + "\n")
        rc = 1

    tail = ""
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines(True)
        tail = "".join(lines[-80:])
    except Exception:
        pass
    return {"status": "done" if rc == 0 else "error", "returncode": rc, "draft_command": draft_cmd, "log_path": str(log_path), "tail": tail}


def run_interactive(logdir: str, *, llm: str, api_key: Optional[str] = None, model: Optional[str] = None, base_url: Optional[str] = None) -> None:
    logdir_p = Path(logdir).expanduser().resolve()
    logdir_p.mkdir(parents=True, exist_ok=True)

    req_file = _required_params_file()
    if not req_file.exists():
        raise FileNotFoundError(f"Required params JSON not found: {req_file}")

    profile = PROVIDERS[llm]
    resolved_key = _resolve_api_key(llm, api_key)
    resolved_model = model or profile.default_model
    resolved_base_url = (base_url or profile.base_url).strip()
    cfg = AIConfig(llm_alias=llm, api_key=resolved_key, model=resolved_model, base_url=resolved_base_url, logdir=str(logdir_p))

    os.environ["IOBRPY_REQUIRED_PARAMS_FILE"] = str(req_file)
    os.environ["IOBRPY_RUN_LOG_DIR"] = str(logdir_p)
    os.environ["IOBRPY_DEFAULTS_FILE"] = str(logdir_p / "iobrpy_defaults.json")

    from . import iobrpy_rag_mcp as server

    server.configure_runtime(server.AIConfig(**cfg.__dict__))
    session_id = _new_session_id()
    prefer_chinese = False
    pending_ready_plan: Optional[Dict[str, Any]] = None

    print(_ui_text("welcome", prefer_chinese))
    print(f"logdir : {logdir_p}")
    print(f"session: {session_id}")
    print(f"llm    : {llm}")
    print(f"model  : {resolved_model}")
    print(_ui_text("type_request", prefer_chinese))
    print(_ui_text("commands", prefer_chinese) + "\n")

    def call(answer: Optional[str] = None) -> Dict[str, Any]:
        return server.tool_iobrpy_assistant(session_id, task=None, answer_text=answer, run=False)

    last = call()
    prefer_chinese = bool(last.get("prefer_chinese", prefer_chinese))
    _print_state(last, prefer_chinese)

    while True:
        try:
            user_in = read_user_input("AI> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n" + _ui_text("bye", prefer_chinese))
            return

        if not user_in:
            continue

        if _contains_cjk(user_in):
            prefer_chinese = True
        elif not _contains_cjk(user_in) and any(c.isalpha() for c in user_in):
            prefer_chinese = False

        low = user_in.lower()
        if low in (":exit", ":quit", "exit", "quit"):
            print(_ui_text("bye", prefer_chinese))
            return
        if low in (":restart", "restart"):
            pending_ready_plan = None
            last = call("restart session")
            prefer_chinese = bool(last.get("prefer_chinese", prefer_chinese))
            _print_state(last, prefer_chinese)
            continue

        if pending_ready_plan is not None and _is_positive_confirm(user_in):
            out = _run_iobrpy_current_env(
                session_id=session_id,
                subcommand=str(pending_ready_plan.get("subcommand")),
                params=dict(pending_ready_plan.get("params") or {}),
                server=server,
                logdir=logdir_p,
            )
            _print_state(out, prefer_chinese)
            pending_ready_plan = None
            continue
        if pending_ready_plan is not None and _is_negative_confirm(user_in):
            print(_ui_text("cancelled", prefer_chinese))
            continue

        last = call(user_in)
        prefer_chinese = bool(last.get("prefer_chinese", prefer_chinese))
        _print_state(last, prefer_chinese)
        if last.get("status") == "ready":
            pending_ready_plan = {
                "subcommand": last.get("subcommand"),
                "params": dict(last.get("params") or {}),
                "draft_command": last.get("draft_command"),
            }
            confirm_in = read_user_input(_ui_text("confirm_run", prefer_chinese)).strip()
            if _is_positive_confirm(confirm_in):
                out = _run_iobrpy_current_env(
                    session_id=session_id,
                    subcommand=str(pending_ready_plan.get("subcommand")),
                    params=dict(pending_ready_plan.get("params") or {}),
                    server=server,
                    logdir=logdir_p,
                )
                _print_state(out, prefer_chinese)
                pending_ready_plan = None
            else:
                print(_ui_text("cancelled", prefer_chinese))


def main(argv: Optional[list[str]] = None) -> None:
    p = argparse.ArgumentParser(prog="iobrpy ai", add_help=True)
    p.add_argument("--logdir", required=True, help="Directory to store AI logs and defaults")
    p.add_argument("--llm", required=True, choices=list(PROVIDERS.keys()), help="LLM provider alias")
    p.add_argument("--api-key", default=None, help="BYOK API key; if omitted, reads env then hidden prompt")
    p.add_argument("--model", default=None, help="Override model name")
    p.add_argument("--base-url", default=None, help="Override provider base URL")
    args = p.parse_args(argv)
    run_interactive(args.logdir, llm=args.llm, api_key=args.api_key, model=args.model, base_url=args.base_url)


if __name__ == "__main__":
    main()
