#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""iobrpy ai (embedded RAG-MCP interactive)

Placement:
  src/iobrpy/RAG_MCP/ai.py

Run:
  iobrpy ai --logdir /path/to/logdir
Optional:
  iobrpy ai --logdir /path/to/logdir --ollama-host http://127.0.0.1:11434 --chat-model qwen3:8b --embed-model qwen3-embedding:8b

Design goals (per your requirements):
  - `iobrpy ai --logdir ...` enters interactive immediately
  - Natural language happens in the interactive prompt (no --task)
  - Chroma is embedded (auto-resolved), no chroma-dir flag
  - Runs iobrpy in CURRENT environment (no conda-related args)
  - When params are satisfied (and confirmations done), auto-run
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import uuid
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional


def _pkg_dir() -> Path:
    return Path(__file__).resolve().parent


def _embedded_chroma_dir() -> Path:
    return _pkg_dir() / "chroma"


def _required_params_file() -> Path:
    return _pkg_dir() / "iobrpy_required_params.json"


def _new_session_id() -> str:
    return f"s{int(time.time())}_{os.getpid()}_{uuid.uuid4().hex[:6]}"


def _set_env_for_embedded_assets(
    logdir: Path,
    *,
    ollama_host: Optional[str] = None,
    chat_model: Optional[str] = None,
    embed_model: Optional[str] = None,
) -> None:
    """Set env vars BEFORE importing the server module (it reads env vars at import-time)."""
    chroma_dir = _embedded_chroma_dir()
    req_file = _required_params_file()

    if not chroma_dir.exists():
        raise FileNotFoundError(f"Embedded chroma dir not found: {chroma_dir}")
    if not req_file.exists():
        raise FileNotFoundError(f"Required params JSON not found: {req_file}")

    # Embedded assets
    os.environ["CHROMA_DIR"] = str(chroma_dir)
    os.environ["IOBRPY_REQUIRED_PARAMS_FILE"] = str(req_file)

    # User-writable state/logs (avoid writing into site-packages)
    os.environ["IOBRPY_RUN_LOG_DIR"] = str(logdir)
    os.environ["IOBRPY_DEFAULTS_FILE"] = str(logdir / "iobrpy_defaults.json")

    # Optional: override Ollama endpoints/models
    if ollama_host:
        os.environ["OLLAMA_HOST"] = ollama_host
    if chat_model:
        os.environ["CHAT_MODEL"] = chat_model
    if embed_model:
        os.environ["EMBED_MODEL"] = embed_model

    # Optional help validation (keep off by default)
    os.environ.setdefault("IOBRPY_VALIDATE_HELP", "0")


def _import_server() -> Any:
    # Import AFTER env vars are set.
    from . import iobrpy_rag_mcp as server  # type: ignore

    return server


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
        if obj.get("draft_command"):
            print(f"\nReady. Draft: {obj.get('draft_command')}")
        else:
            print("\nReady.")
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


def _run_iobrpy_current_env(
    *,
    session_id: str,
    subcommand: str,
    params: Dict[str, Any],
    server: Any,
    logdir: Path,
) -> Dict[str, Any]:
    """Run planned iobrpy command using current python environment."""
    rules = server.load_rules()
    draft_cmd, argv = server.build_command(subcommand, params, rules)

    log_path = logdir / f"{session_id}_{subcommand}.log"

    # Always use current env python
    cmd = [sys.executable, "-m", "iobrpy.main"] + argv

    rc: Optional[int] = None
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

    return {
        "status": "done" if rc == 0 else "error",
        "returncode": rc,
        "draft_command": draft_cmd,
        "log_path": str(log_path),
        "tail": tail,
    }


def run_interactive(
    logdir: str,
    *,
    ollama_host: Optional[str] = None,
    chat_model: Optional[str] = None,
    embed_model: Optional[str] = None,
) -> None:
    """Called by iobrpy.main -> `iobrpy ai --logdir ...`."""
    logdir_p = Path(logdir).expanduser().resolve()
    logdir_p.mkdir(parents=True, exist_ok=True)

    _set_env_for_embedded_assets(
        logdir_p,
        ollama_host=ollama_host,
        chat_model=chat_model,
        embed_model=embed_model,
    )
    server = _import_server()
    server_file = Path(getattr(server, "__file__", "<unknown>")).resolve()
    rules_file = Path(os.environ.get("IOBRPY_REQUIRED_PARAMS_FILE", "<unknown>"))

    session_id = _new_session_id()

    print("IOBRpy AI (embedded RAG-MCP)")
    print(f"logdir : {logdir_p}")
    print(f"session: {session_id}")
    print(f"server : {server_file}")
    print(f"rules  : {rules_file}")
    if ollama_host:
        print(f"ollama : {ollama_host}")
    if chat_model:
        print(f"chat   : {chat_model}")
    if embed_model:
        print(f"embed  : {embed_model}")

    print("\nType your request in natural language.")
    print("Commands: :exit  :quit  :restart\n")

    def call(answer: Optional[str] = None) -> Dict[str, Any]:
        # run=False here; we run ourselves in current env.
        return server.tool_iobrpy_assistant(session_id, task=None, answer_text=answer, run=False)

    # initial prompt
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

        # feed user message
        last = call(user_in)
        _print_state(last)

        # Auto-run when ready
        if last.get("status") == "ready":
            out = _run_iobrpy_current_env(
                session_id=session_id,
                subcommand=str(last.get("subcommand")),
                params=dict(last.get("params") or {}),
                server=server,
                logdir=logdir_p,
            )
            _print_state(out)

            # reset for next task
            last = call("restart session")
            _print_state(last)


def main(argv: Optional[list[str]] = None) -> None:
    p = argparse.ArgumentParser(prog="iobrpy ai", add_help=True)
    p.add_argument("--logdir", required=True, help="Directory to store AI logs and defaults")
    p.add_argument("--ollama-host", default='http://127.0.0.1:11434', help="Override Ollama host, e.g. http://127.0.0.1:11434")
    p.add_argument("--chat-model", default='qwen3:8b', help="Override chat model, e.g. qwen3:8b")
    p.add_argument("--embed-model", default='qwen3-embedding:8b', help="Override embedding model, e.g. qwen3-embedding:8b")
    args = p.parse_args(argv)

    run_interactive(
        args.logdir,
        ollama_host=args.ollama_host,
        chat_model=args.chat_model,
        embed_model=args.embed_model,
    )


if __name__ == "__main__":
    main()
