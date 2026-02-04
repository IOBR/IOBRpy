from iobrpy.ai.planner import SmartPlanner
from iobrpy.ai.registry import ToolRegistry


def test_smart_planner_runall_salmon_extraction():
    registry = ToolRegistry.from_main()
    planner = SmartPlanner()
    prompt = (
        "I have FASTQ files at /mnt/data/01-fq. "
        "Run the complete pipeline with salmon and output to /mnt/out/ai, threads 16"
    )
    result = planner.plan(prompt, registry)
    assert not result.need_clarification
    assert result.steps
    step = result.steps[0]
    assert step.tool_name == "runall"
    assert step.args.get("mode") == "salmon"
    assert step.args.get("fastq") == "/mnt/data/01-fq"
    assert step.args.get("outdir") == "/mnt/out/ai"
    assert "--threads" in step.argv
    assert "16" in step.argv
