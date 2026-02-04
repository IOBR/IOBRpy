from iobrpy.ai.plan import PlanStep, PlanResult, plan_schema, validate_plan_dict


def test_plan_schema():
    step = PlanStep(
        id=1,
        tool_name="count2tpm",
        mode="kwargs",
        args={"count_mat": "./counts.tsv", "output_path": "./tpm.tsv"},
        rationale="Convert counts to TPM.",
        expected_outputs=["./tpm.tsv"],
        depends_on=[],
    )
    result = PlanResult(steps=[step])
    plan_dict = result.to_dict()
    schema = plan_schema()

    assert isinstance(schema, dict)
    validation = validate_plan_dict(plan_dict)
    assert validation["ok"], validation["errors"]
