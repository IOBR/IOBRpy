# Skill: IOBRpy development

Use this skill when you need to make changes to the IOBRpy codebase, update docs, or prepare releases.

## Scope
- Repository: `/workspace/IOBRpy`
- Focus: CLI workflows, documentation, packaging metadata, and pipelines.

## Quick orientation
- CLI entry points: `src/iobrpy/workflow/`
- Core resources: `src/iobrpy/resources/`
- Documentation: `docs/` and `README.md`
- Packaging config: `pyproject.toml`, `requirements.txt`, `environment.yml`

## Common tasks
1. Find relevant code paths with `rg` (avoid recursive `ls -R`).
2. Keep command-line interfaces consistent with existing flags and output formats.
3. Update README/docs when behavior or output changes.

## Suggested checks (optional)
- `python -m pip install -e .`
- `python -m iobrpy --help`
- Run any relevant CLI subcommand for a smoke test.
