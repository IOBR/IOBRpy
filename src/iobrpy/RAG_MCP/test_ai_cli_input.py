import sys
import threading
import unittest
from pathlib import Path

try:
    from prompt_toolkit.input import create_pipe_input
    from prompt_toolkit.output import DummyOutput
except Exception:
    create_pipe_input = None
    DummyOutput = None

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from iobrpy.RAG_MCP import ai


@unittest.skipIf(create_pipe_input is None or DummyOutput is None, "prompt_toolkit is not installed")
class AIComposerBehaviorTests(unittest.TestCase):
    def _run_composer_with_keys(self, keys: str) -> str:
        with create_pipe_input() as pipe_input:
            composer = ai.build_main_composer(input_obj=pipe_input, output_obj=DummyOutput())
            self.assertIsNotNone(composer)

            result = {}

            def _runner():
                result["text"] = composer.run()

            t = threading.Thread(target=_runner)
            t.start()
            pipe_input.send_text(keys)
            t.join(timeout=3)
            self.assertFalse(t.is_alive(), "composer.run did not finish")
            return result["text"]

    def test_main_prompt_prefix_is_iobrpy(self):
        self.assertEqual(ai.build_main_prompt_prefix(), "IOBRpy> ")

    def test_main_composer_renders_iobrpy_prefix(self):
        with create_pipe_input() as pipe_input:
            composer = ai.build_main_composer(input_obj=pipe_input, output_obj=DummyOutput())
            self.assertIsNotNone(composer)
            root_children = composer.app.layout.container.children
            self.assertEqual(root_children[0].content.text, "IOBRpy> ")

    def test_main_input_help_text_matches_runtime_shift_enter_capability(self):
        help_text = ai.build_main_input_help_text()
        if ai.main_input_supports_shift_enter():
            self.assertIn("Shift+Enter=newline", help_text)
        else:
            self.assertNotIn("Shift+Enter=newline", help_text)
            self.assertIn("cannot distinguish Shift+Enter", help_text)

    def test_enter_submits_main_composer(self):
        out = self._run_composer_with_keys("hello world\r")
        self.assertEqual(out, "hello world")

    def test_multiline_buffer_can_merge_lines_with_backspace_at_line_start(self):
        buf = ai.create_main_buffer()
        self.assertIsNotNone(buf)
        buf.text = "a\nb"
        buf.cursor_position = 2
        buf.delete_before_cursor(count=1)
        self.assertEqual(buf.text, "ab")
        self.assertEqual(buf.cursor_position, 1)

    def test_paste_multiline_keeps_single_buffer_submission(self):
        payload = "line1\nline2\nline3"
        out = self._run_composer_with_keys(payload + "\r")
        self.assertEqual(out, payload)


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
