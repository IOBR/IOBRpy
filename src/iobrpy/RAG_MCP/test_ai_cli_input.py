import sys
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


class AICliInputTests(unittest.TestCase):
    def test_main_prompt_prefix_is_iobrpy(self):
        self.assertEqual(ai.build_main_prompt_prefix(), "IOBRpy> ")

    def test_main_help_text_reflects_single_line_mode(self):
        help_text = ai.build_main_input_help_text()
        self.assertIn("Single-line input mode", help_text)
        self.assertIn("Enter to submit", help_text)
        self.assertIn("Multi-line paste is supported", help_text)
        self.assertNotIn("Shift+Enter", help_text)
        self.assertNotIn("Ctrl+J", help_text)
        self.assertNotIn("Esc+Enter", help_text)

    def test_normalize_single_line_input(self):
        self.assertEqual(ai.normalize_main_input_text("  hello world  "), "hello world")

    def test_normalize_multiline_paste(self):
        payload = " line1\nline2 \n\n line3\r\n"
        self.assertEqual(ai.normalize_main_input_text(payload), "line1 line2 line3")

    def test_read_main_user_input_keeps_multiline_paste_as_one_request(self):
        payload = "fastq路径在/a\n输出在/b\nmode用salmon"
        sess = _FakeSession(payload)
        out = ai.read_main_user_input(session=sess)
        self.assertEqual(out, "fastq路径在/a 输出在/b mode用salmon")
        self.assertEqual(sess.calls, 1)
        self.assertEqual(sess.last_prompt, "IOBRpy> ")


class AIConfirmationFlowTests(unittest.TestCase):
    def test_confirmation_other_text_routes_back_to_main_flow(self):
        calls = []

        def _call_fn(x):
            calls.append(x)
            return {"status": "need_info", "question": "ok"}

        pending = {"subcommand": "runall", "params": {}}
        pending2, out, handled = ai._handle_confirmation_response(  # pylint: disable=protected-access
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
        def _call_fn(_):
            return {"status": "need_info"}

        def _runner(**kwargs):
            return {"status": "done", "returncode": 0}

        old = ai._run_iobrpy_current_env
        ai._run_iobrpy_current_env = _runner
        try:
            pending2, out, handled = ai._handle_confirmation_response(  # pylint: disable=protected-access
                "yes",
                pending_ready_plan={"subcommand": "runall", "params": {}},
                session_id="s1",
                server=None,
                logdir_p=Path("."),
                prefer_chinese=False,
                call_fn=_call_fn,
            )
        finally:
            ai._run_iobrpy_current_env = old

        self.assertTrue(handled)
        self.assertIsNone(pending2)
        self.assertEqual(out.get("status"), "done")

    def test_confirmation_no_cancels(self):
        def _call_fn(_):
            return {"status": "need_info"}

        pending2, out, handled = ai._handle_confirmation_response(  # pylint: disable=protected-access
            "no",
            pending_ready_plan={"subcommand": "runall", "params": {}},
            session_id="s1",
            server=None,
            logdir_p=Path("."),
            prefer_chinese=False,
            call_fn=_call_fn,
        )
        self.assertTrue(handled)
        self.assertIsNone(pending2)
        self.assertIsNone(out)


if __name__ == "__main__":
    unittest.main()
