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


def _print_state(obj: Dict[str, Any]) -> None:
    status = obj.get("status")
    if status == "need_info":
        if obj.get("draft_command"):
            print(f"\nDraft: {obj.get('draft_command')}")
        q = obj.get("question") or ""
        if q:
            print(q)
        return
    if status == "ready":
        print(f"\nReady. Draft: {obj.get('draft_command')}")
        return
    if status in ("done", "error"):
        print(f"\n[{status}] rc={obj.get('returncode')}")
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


def run_interactive(logdir: str, *, llm: str, api_key: Optional[str] = None, model: Optional[str] = None) -> None:
    logdir_p = Path(logdir).expanduser().resolve()
    logdir_p.mkdir(parents=True, exist_ok=True)

    req_file = _required_params_file()
    if not req_file.exists():
        raise FileNotFoundError(f"Required params JSON not found: {req_file}")

    profile = PROVIDERS[llm]
    resolved_key = _resolve_api_key(llm, api_key)
    resolved_model = model or profile.default_model
    cfg = AIConfig(llm_alias=llm, api_key=resolved_key, model=resolved_model, base_url=profile.base_url, logdir=str(logdir_p))

    os.environ["IOBRPY_REQUIRED_PARAMS_FILE"] = str(req_file)
    os.environ["IOBRPY_RUN_LOG_DIR"] = str(logdir_p)
    os.environ["IOBRPY_DEFAULTS_FILE"] = str(logdir_p / "iobrpy_defaults.json")

    from . import iobrpy_rag_mcp as server

    server.configure_runtime(server.AIConfig(**cfg.__dict__))
    session_id = _new_session_id()

    print("IOBRpy AI (BYOK cloud LLM)")
    print(f"logdir : {logdir_p}")
    print(f"session: {session_id}")
    print(f"llm    : {llm}")
    print(f"model  : {resolved_model}")
    print("\nType your request in natural language.")
    print("Commands: :exit  :quit  :restart\n")

    def call(answer: Optional[str] = None) -> Dict[str, Any]:
        return server.tool_iobrpy_assistant(session_id, task=None, answer_text=answer, run=False)

    last = call()
    _print_state(last)

    while True:
        try:
            user_in = input("AI> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nbye")
            return
        if not user_in:
            continue
        low = user_in.lower()
        if low in (":exit", ":quit", "exit", "quit"):
            print("bye")
            return
        if low in (":restart", "restart"):
            last = call("restart session")
            _print_state(last)
            continue

        last = call(user_in)
        _print_state(last)
        if last.get("status") == "ready":
            out = _run_iobrpy_current_env(session_id=session_id, subcommand=str(last.get("subcommand")), params=dict(last.get("params") or {}), server=server, logdir=logdir_p)
            _print_state(out)
            last = call("restart session")
            _print_state(last)


def main(argv: Optional[list[str]] = None) -> None:
    p = argparse.ArgumentParser(prog="iobrpy ai", add_help=True)
    p.add_argument("--logdir", required=True, help="Directory to store AI logs and defaults")
    p.add_argument("--llm", required=True, choices=list(PROVIDERS.keys()), help="LLM provider alias")
    p.add_argument("--api-key", default=None, help="BYOK API key; if omitted, reads env then hidden prompt")
    p.add_argument("--model", default=None, help="Override model name")
    args = p.parse_args(argv)
    run_interactive(args.logdir, llm=args.llm, api_key=args.api_key, model=args.model)


if __name__ == "__main__":
    main()
