---
title: "Agent Guide"
page-layout: article
toc: true
---

## What This Page Helps You Do

This guide is for users who want to connect **IOBRpy** to **Codex** or **Claude Code**. It focuses on three practical goals:

1. How to install the `/iobrpy` plugin and entrypoints.
2. How to verify that the installation worked.
3. What kinds of prompts help `/iobrpy` reach the right result faster.

---

## Tutorial 1: How to Install the Plugin

### Step 0: Activate the IOBRpy Environment

This page assumes that IOBRpy is already installed. If not, please read installation guide first. Here, the only environment preparation you need is:

```bash
conda activate iobrpy
```

---

### Step 1: Run the Installation Command

Choose the one command that matches your actual client:

```bash
# Claude Code only
iobrpy-cli agent install --client claude-code

# Codex only
iobrpy-cli agent install --client codex

# Both Codex and Claude Code
iobrpy-cli agent install --client all
```

What this installs:

- `codex`: the Codex skill, Codex plugin, and MCP integration.
- `claude-code`: the Claude Code memory, `/iobrpy` command entrypoint, and MCP integration.
- `all`: both of the above.

**Reference standard:**

- The command finishes successfully instead of exiting with an error.
- If a component was already present, an `already_installed`-style result is still acceptable.
- You should not need to manually copy plugin, skill, or MCP files.

---

### Step 2: Verify the Installation

Use the status command that matches what you installed:

```bash
# Check Codex only
iobrpy-cli agent status --client codex

# Check Claude Code only
iobrpy-cli agent status --client claude-code

# Check everything
iobrpy-cli agent status
```

**Reference standard:**

- Ideally, the output should include `Agent status: healthy`.
- You should no longer see states such as `not_installed`, `not_configured`, or `invalid_config`.
- After installation, the relevant client should be able to use `/iobrpy`.

---

### Step 3: Call the Plugin Inside the Agent

After installation, open your agent client and call the IOBRpy entrypoint directly in the chat input.

For both Codex and Claude Code, the normal entrypoint is:

```text
/iobrpy
```

If you are using Codex and specifically want the plugin-scoped namespace, you can also use:

```text
/iobrpy:iobrpy
```

**Reference standard:**

- Typing `/iobrpy` in the agent should show or trigger the IOBRpy entrypoint.
- After you send a message with `/iobrpy`, the agent should respond as an IOBRpy workflow assistant rather than a generic chat assistant.
- In Codex, `/iobrpy:iobrpy` should also work when you want the plugin namespace explicitly.
- The best first test is to give a real path or a concrete analysis goal, not just `/iobrpy` alone.

---

## Common Problems and Fixes

| Symptom | Common cause | Suggested fix |
| --- | --- | --- |
| `iobrpy-cli` is not recognized | The current shell is not in the right environment, or IOBRpy is not installed yet | Activate the environment with `conda activate iobrpy`, and if needed follow the installation guide first |
| `agent status` shows `needs attention` | One client is missing part of the skill, plugin, or MCP setup | Check `iobrpy-cli --json agent status` first, then rerun the relevant `agent install --client ...` command |
| Codex does not show `/iobrpy` | The skill or plugin did not take effect | Rerun `iobrpy-cli agent install --client codex` |
| Claude Code installation fails | The `claude` CLI is missing, or the config path is not writable | Install or fix the Claude Code CLI setup, then retry `iobrpy-cli agent install --client claude-code` |
| You suspect the setup is outdated | Older local agent files may still be present | Rerun the matching `iobrpy-cli agent install --client ...` command and check status again |

---

## Tutorial 2: What Kinds of Prompts Work Best for `/iobrpy`

### 1. Path Scanning Prompt

Use this first when you are not yet sure what the directory contains or which workflow command is the right one.

```text
/iobrpy Please scan path/to/your/target/directory
Show me the full checklist first, then tell me the current stage and the best next step.
```

### 2. Choose the Command Suggested by the Agent

After the scan, the agent will usually tell you which workflow best matches your directory and your current stage.

Typical examples:

- If the directory contains raw FASTQ data, the agent will usually suggest `runall`.
- If the directory already contains a TPM matrix, the agent will usually suggest `tme_profile`.
- If the directory contains BAM files for HLA analysis, the agent will usually suggest `hla_typing`.
- If the goal is TCR/BCR repertoire analysis, the agent will usually suggest `trust4`.

This means you do not need to decide the command first. Start with the scan, then choose the command that matches the agent's recommendation and your real analysis goal.

**Reference standard:**

- Your first prompt includes a real path.
- The agent returns a checklist and a current-stage summary.
- You choose the next command based on the scan result instead of guessing before the scan.
