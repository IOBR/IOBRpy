from iobrpy.ai.registry import ToolRegistry


def test_registry_finds_commands():
    registry = ToolRegistry.from_main()
    tools = {tool.name for tool in registry.list_tools()}
    assert len(tools) >= 25
    assert "spechla" in tools
    assert "extract_hla_read" in tools
    assert "bayesprism" in tools
