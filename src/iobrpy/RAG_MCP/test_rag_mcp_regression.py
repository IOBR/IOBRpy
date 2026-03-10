import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from iobrpy.RAG_MCP import iobrpy_rag_mcp as mcp


def make_planner_output(
    *,
    action: str,
    intent_type: str = "other",
    is_side_question: bool = False,
    keep_current_command: bool = False,
    subcommand=None,
    param_updates=None,
    param_removals=None,
    switch_to=None,
    message: str = "",
    confidence: float = 0.9,
    reason: str = "test",
):
    return {
        "action": action,
        "intent_type": intent_type,
        "is_side_question": is_side_question,
        "keep_current_command": keep_current_command,
        "switch_intent_strength": 0.0,
        "needs_confirmation": False,
        "message": message,
        "subcommand": subcommand,
        "param_updates": param_updates or {},
        "param_removals": param_removals or [],
        "switch_to": switch_to,
        "ask_for": [],
        "confirm": None,
        "confidence": confidence,
        "reason": reason,
    }


def planner_sequence(*responses):
    it = iter(responses)

    def _fake_planner(state, user_text, rules):  # pylint: disable=unused-argument
        return next(it)

    return _fake_planner


class RagMcpRegressionTests(unittest.TestCase):
    def setUp(self):
        mcp.SESSIONS.clear()

    def test_param_updates_merge_into_existing_runall_state(self):
        plan1 = make_planner_output(
            action="select_command",
            intent_type="execution_request",
            subcommand="runall",
            reason="select_runall",
        )
        # non-update action to verify orchestration falls back to extractor merge path
        plan2 = make_planner_output(
            action="reply_only",
            intent_type="parameter_update",
            keep_current_command=True,
            reason="allow_backfill",
        )

        extractor_updates = {
            "fastq": "/mnt/cephfs/s2z5/iobrpy/01-fq",
            "outdir": "/mnt/cephfs/s2z5/iobrpy/03-result/runall",
            "mode": "salmon",
            "index": "/mnt/cephfs/s2z1/db/humandb/salmon/gencodev44/salmon_index",
            "threads": 32,
            "batch_size": 4,
            "project": "runall",
        }

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(plan1, plan2)), patch.object(
            mcp, "llm_extract_param_updates", return_value=extractor_updates
        ):
            mcp.tool_iobrpy_assistant("s1", answer_text="turn_1", run=False)
            out = mcp.tool_iobrpy_assistant("s1", answer_text="/tmp/p 32 salmon", run=False)

        for key, value in extractor_updates.items():
            self.assertEqual(out["params"].get(key), value)
        self.assertEqual(out.get("needs"), [])
        self.assertIn("iobrpy runall", out.get("draft_command") or "")

    def test_planner_subcommand_reply_only_is_normalized_and_committed(self):
        planner_out = make_planner_output(
            action="reply_only",
            intent_type="execution_request",
            is_side_question=False,
            keep_current_command=False,
            subcommand="runall",
            message="planner_message",
            reason="planner_has_subcommand_but_non_commit_action",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(planner_out)), patch.object(
            mcp,
            "llm_select_command_from_catalog",
            side_effect=AssertionError("selector should not run when normalization can commit"),
        ):
            out = mcp.tool_iobrpy_assistant("s_norm_commit", answer_text="turn_1", run=False)

        self.assertEqual(mcp.get_session("s_norm_commit").selected_command, "runall")
        self.assertEqual(out.get("subcommand"), "runall")
        self.assertTrue(out.get("draft_command"))
        self.assertTrue(out.get("state_summary"))

    def test_side_question_with_subcommand_does_not_commit_when_no_selected_command(self):
        planner_out = make_planner_output(
            action="reply_only",
            intent_type="side_question",
            is_side_question=True,
            keep_current_command=True,
            subcommand="runall",
            reason="side_question_with_subcommand",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(planner_out)), patch.object(
            mcp,
            "llm_select_command_from_catalog",
            return_value={"subcommand": None, "confidence": 0.0, "reason": "not_applicable"},
        ):
            out = mcp.tool_iobrpy_assistant("s_side_no_commit", answer_text="turn_1", run=False)

        self.assertIsNone(mcp.get_session("s_side_no_commit").selected_command)
        self.assertNotEqual(out.get("subcommand"), "runall")

    def test_first_turn_selector_fallback_sets_runall_when_planner_has_no_command(self):
        planner_no_command = make_planner_output(
            action="reply_only",
            intent_type="execution_request",
            keep_current_command=False,
            subcommand=None,
            reason="planner_miss",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(planner_no_command)), patch.object(
            mcp,
            "llm_select_command_from_catalog",
            return_value={"subcommand": "runall", "confidence": 0.9, "reason": "selector_match"},
        ):
            out = mcp.tool_iobrpy_assistant("s_selector_ok", answer_text="turn_1", run=False)

        self.assertEqual(out.get("subcommand"), "runall")
        self.assertTrue(out.get("draft_command"))
        self.assertEqual(mcp.get_session("s_selector_ok").selected_command, "runall")

    def test_selector_low_confidence_does_not_force_command_selection(self):
        planner_no_command = make_planner_output(
            action="clarify_intent",
            intent_type="clarification",
            subcommand=None,
            reason="planner_uncertain",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(planner_no_command)), patch.object(
            mcp,
            "llm_select_command_from_catalog",
            return_value={"subcommand": None, "confidence": 0.3, "reason": "low_confidence"},
        ):
            out = mcp.tool_iobrpy_assistant("s_selector_low", answer_text="turn_1", run=False)

        self.assertNotEqual(out.get("subcommand"), "runall")
        self.assertIsNone(mcp.get_session("s_selector_low").selected_command)

    def test_selected_command_context_always_returns_draft(self):
        plan = make_planner_output(
            action="select_command",
            intent_type="execution_request",
            subcommand="runall",
            reason="select_runall",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(plan)):
            out = mcp.tool_iobrpy_assistant("s2", answer_text="turn_1", run=False)

        self.assertEqual(out.get("subcommand"), "runall")
        self.assertTrue(out.get("draft_command"))
        self.assertTrue(out.get("state_summary"))

    def test_need_info_response_keeps_command_context(self):
        plan1 = make_planner_output(
            action="select_command",
            intent_type="execution_request",
            subcommand="runall",
            reason="select_runall",
        )
        plan2 = make_planner_output(
            action="reply_only",
            intent_type="parameter_update",
            keep_current_command=True,
            reason="partial_backfill",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(plan1, plan2)), patch.object(
            mcp,
            "llm_extract_param_updates",
            return_value={
                "fastq": "/mnt/cephfs/s2z5/iobrpy/01-fq",
                "outdir": "/mnt/cephfs/s2z5/iobrpy/03-result/runall",
            },
        ):
            mcp.tool_iobrpy_assistant("s3", answer_text="turn_1", run=False)
            out = mcp.tool_iobrpy_assistant("s3", answer_text="/tmp/p 8", run=False)

        self.assertEqual(out.get("status"), "need_info")
        self.assertEqual(out.get("subcommand"), "runall")
        self.assertTrue(out.get("draft_command"))
        self.assertTrue(out.get("state_summary"))

    def test_side_question_preserves_selected_command(self):
        plan1 = make_planner_output(
            action="select_command",
            intent_type="execution_request",
            subcommand="runall",
            reason="select_runall",
        )
        plan2 = make_planner_output(
            action="reply_only",
            intent_type="side_question",
            is_side_question=True,
            keep_current_command=True,
            message="context-only",
            reason="side_question",
        )

        with patch.object(mcp, "llm_plan_next_action", side_effect=planner_sequence(plan1, plan2)):
            mcp.tool_iobrpy_assistant("s4", answer_text="turn_1", run=False)
            out = mcp.tool_iobrpy_assistant("s4", answer_text="turn_2", run=False)

        self.assertEqual(mcp.get_session("s4").selected_command, "runall")
        self.assertEqual(out.get("subcommand"), "runall")
        self.assertTrue(out.get("draft_command"))
        self.assertTrue(out.get("state_summary"))


if __name__ == "__main__":
    unittest.main()
