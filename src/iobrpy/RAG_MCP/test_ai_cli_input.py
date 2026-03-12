import io
import sys
import types
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from iobrpy.RAG_MCP import ai


class _FakeSession:
    def __init__(self, out: str):
        self._out = out
        self.last_prompt = None
        self.calls = 0

    def prompt(self, *args, **kwargs):
        self.calls += 1
        if args:
            self.last_prompt = args[0]
        return self._out


class _FailingSession:
    def prompt(self, *args, **kwargs):
        raise RuntimeError("boom")


class AICliInputTests(unittest.TestCase):
    def test_main_prompt_prefix_is_iobrpy(self):
        self.assertEqual(ai.build_main_prompt_prefix(), "IOBRpy> ")

    def test_read_main_user_input_prefers_high_level_session(self):
        sess = _FakeSession("hello")
        out = ai.read_main_user_input(session=sess)
        self.assertEqual(out, "hello")
        self.assertEqual(sess.calls, 1)
        self.assertEqual(sess.last_prompt, "IOBRpy> ")

    def test_fallback_to_single_input_when_session_fails(self):
        prompts = []

        def _fake_input(prompt):
            prompts.append(prompt)
            return "fallback line"

        old_input = ai.__dict__.get("input", input)
        ai.input = _fake_input
        try:
            out = ai.read_main_user_input(session=_FailingSession())
        finally:
            ai.input = old_input

        self.assertEqual(out, "fallback line")
        self.assertEqual(prompts, ["IOBRpy> "])

    def test_normalize_multiline_payload(self):
        self.assertEqual(ai.normalize_main_input_text("你好\n你是什么模型\n"), "你好 你是什么模型")

    def test_help_text_does_not_overclaim_paste_guarantee(self):
        text = ai.build_main_input_help_text()
        self.assertIn("Single-line input mode", text)
        self.assertIn("depends on terminal behavior", text)
        self.assertNotIn("will be sent as one request", text)

    def test_confirmation_prompt_builder_exists_and_is_localized(self):
        self.assertIn("[y/n]", ai.build_confirmation_prompt(False))
        self.assertIn("[y/n]", ai.build_confirmation_prompt(True))

    def test_read_confirmation_input_prefers_confirmation_session(self):
        sess = _FakeSession("yes")
        out = ai.read_confirmation_input("Execute this command now? [y/n] ", session=sess)
        self.assertEqual(out, "yes")
        self.assertEqual(sess.calls, 1)


class AIRunInteractiveOutputTests(unittest.TestCase):
    def test_run_interactive_does_not_print_main_input_help_text(self):
        fake_server = types.SimpleNamespace()

        class _Cfg:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

        fake_server.AIConfig = _Cfg
        fake_server.configure_runtime = lambda _cfg: None
        fake_server.tool_iobrpy_assistant = lambda *args, **kwargs: {
            "status": "need_info",
            "question": "ok",
            "prefer_chinese": False,
        }

        old_module = sys.modules.get("iobrpy.RAG_MCP.iobrpy_rag_mcp")
        sys.modules["iobrpy.RAG_MCP.iobrpy_rag_mcp"] = fake_server

        old_read = ai.read_main_user_input
        ai.read_main_user_input = lambda *args, **kwargs: (_ for _ in ()).throw(KeyboardInterrupt())

        old_stdout = sys.stdout
        sys.stdout = io.StringIO()
        try:
            ai.run_interactive(logdir="./.tmp_ai_logs", llm="openai", api_key="k")
            out = sys.stdout.getvalue()
        finally:
            sys.stdout = old_stdout
            ai.read_main_user_input = old_read
            if old_module is None:
                sys.modules.pop("iobrpy.RAG_MCP.iobrpy_rag_mcp", None)
            else:
                sys.modules["iobrpy.RAG_MCP.iobrpy_rag_mcp"] = old_module

        self.assertNotIn("Single-line input mode", out)


class AIConfirmationFlowTests(unittest.TestCase):
    def test_confirmation_other_text_routes_back_to_main_flow(self):
        calls = []

        def _call_fn(x):
            calls.append(x)
            return {"status": "need_info", "question": "ok"}

        pending = {"subcommand": "runall", "params": {}}
        pending2, out, handled = ai._handle_confirmation_response(
            "free form update",
            pending_ready_plan=pending,
            session_id="s1",
            server=None,
            logdir_p=Path("."),
            prefer_chinese=False,
            call_fn=_call_fn,
        )
        self.assertFalse(handled)
        self.assertIsNotNone(pending2)
        self.assertEqual(calls, ["free form update"])
        self.assertEqual(out.get("status"), "need_info")

    def test_confirmation_yes_executes(self):
        old = ai._run_iobrpy_current_env
        ai._run_iobrpy_current_env = lambda **kwargs: {"status": "done", "returncode": 0}
        try:
            pending2, out, handled = ai._handle_confirmation_response(
                "yes",
                pending_ready_plan={"subcommand": "runall", "params": {}},
                session_id="s1",
                server=None,
                logdir_p=Path("."),
                prefer_chinese=False,
                call_fn=lambda _: {"status": "need_info"},
            )
        finally:
            ai._run_iobrpy_current_env = old

        self.assertTrue(handled)
        self.assertIsNone(pending2)
        self.assertEqual(out.get("status"), "done")

    def test_confirmation_no_cancels(self):
        pending2, out, handled = ai._handle_confirmation_response(
            "no",
            pending_ready_plan={"subcommand": "runall", "params": {}},
            session_id="s1",
            server=None,
            logdir_p=Path("."),
            prefer_chinese=False,
            call_fn=lambda _: {"status": "need_info"},
        )
        self.assertTrue(handled)
        self.assertIsNone(pending2)
        self.assertIsNone(out)


if __name__ == "__main__":
    unittest.main()
