import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from iobrpy.RAG_MCP import ai


class _FakeSession:
    def __init__(self, out: str):
        self._out = out

    def prompt(self, *args, **kwargs):
        return self._out


class AICliInputTests(unittest.TestCase):
    def test_main_input_supports_multiline_payload_as_single_submission(self):
        payload = "line1\nline2\nline3"
        out = ai.read_main_user_input("AI> ", session=_FakeSession(payload))
        self.assertEqual(out, payload)

    def test_main_shortcut_submit_classification(self):
        self.assertEqual(ai.classify_main_shortcut("escape+enter"), "submit")
        self.assertEqual(ai.classify_main_shortcut("c-s"), "submit")

    def test_main_shortcut_newline_inserts_newline_not_submit(self):
        self.assertEqual(ai.classify_main_shortcut("c-j"), "newline")
        text, pos = ai.apply_newline_to_text("abc", 1)
        self.assertEqual(text, "a\nbc")
        self.assertEqual(pos, 2)

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
