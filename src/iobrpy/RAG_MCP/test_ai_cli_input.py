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

    def test_normalize_multiline_payload(self):
        payload = " line1\nline2 \n\n line3\r\n"
        self.assertEqual(ai.normalize_main_input_text(payload), "line1 line2 line3")

    def test_single_payload_with_newlines_is_normalized_once(self):
        payload = "你好\n你是什么模型\n"
        sess = _FakeSession(payload)
        out = ai.read_main_user_input(session=sess, drain_fn=lambda: "")
        self.assertEqual(out, "你好 你是什么模型")
        self.assertEqual(sess.calls, 1)

    def test_split_arrival_is_coalesced_to_single_request(self):
        sess = _FakeSession("你好")
        out = ai.read_main_user_input(session=sess, drain_fn=lambda: "\n你是什么模型\n")
        self.assertEqual(out, "你好 你是什么模型")
        self.assertEqual(sess.calls, 1)

    def test_normal_two_turn_input_is_not_merged(self):
        sess = _FakeSession("你好")
        out = ai.read_main_user_input(session=sess, drain_fn=lambda: "")
        self.assertEqual(out, "你好")

    def test_drain_reader_aggregates_chunks_until_empty(self):
        chunks = ["\n第二行", "\n第三行", ""]

        def _reader(_timeout):
            return chunks.pop(0)

        drained = ai._drain_immediately_available_paste_lines(read_chunk_fn=_reader)
        self.assertEqual(drained, "\n第二行\n第三行")

    def test_fallback_is_single_read_and_no_continuation_prompt(self):
        prompts = []

        def _fake_input(prompt):
            prompts.append(prompt)
            return "line1\nline2"

        old_input = ai.__dict__.get("input", input)
        ai.input = _fake_input
        try:
            out = ai.read_main_user_input(session=_FailingSession(), drain_fn=lambda: "")
        finally:
            ai.input = old_input

        self.assertEqual(out, "line1 line2")
        self.assertEqual(prompts, ["IOBRpy> "])
        self.assertNotIn("... ", prompts)


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

        self.assertNotIn(
            "Single-line input mode. Press Enter to submit. Multi-line paste is supported and will be sent as one request.",
            out,
        )


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
