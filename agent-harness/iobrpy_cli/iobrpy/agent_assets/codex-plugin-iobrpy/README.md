# iobrpy Codex Plugin

This local plugin packages the iobrpy MCP server and an `iobrpy` command namespace for Codex.

## What it adds

- A local `iobrpy` MCP server via `iobrpy-cli-mcp`
- A plugin command namespace with `/iobrpy:iobrpy`

## Recommended usage

The `iobrpy-cli agent install --client codex` bootstrap also installs a top-level `~/.codex/skills/iobrpy` skill so users can type `/iobrpy` directly.

Use the plugin namespace when you specifically want the plugin-scoped command surface:

```text
/iobrpy:iobrpy parse /path/to/project
```
