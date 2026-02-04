from __future__ import annotations

from abc import ABC, abstractmethod
import re
from typing import Dict, List, Optional, Tuple

from iobrpy.ai.plan import PlanResult, PlanStep
from iobrpy.ai.registry import ToolRegistry, ToolSpec


class Planner(ABC):
    @abstractmethod
    def plan(
        self,
        prompt: str,
        registry: ToolRegistry,
        context: Optional[Dict[str, object]] = None,
    ) -> PlanResult:
        raise NotImplementedError


class SimplePlanner(Planner):
    """A lightweight planner that maps prompts to a single tool call."""

    def plan(
        self,
        prompt: str,
        registry: ToolRegistry,
        context: Optional[Dict[str, object]] = None,
    ) -> PlanResult:
        prompt = (prompt or "").strip()
        if not prompt:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=["Please describe the analysis task you want to run."],
                required_fields=["tool_name", "inputs", "outputs"],
                example_payload={
                    "tool_name": "tme_profile",
                    "inputs": {"input": "./data/tpm.tsv"},
                    "outputs": {"output": "./results/tme"},
                },
            )

        tool = self._find_tool(prompt, registry)
        if tool is None:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=[
                    "Which iobrpy command should be used?",
                    "Provide the required input/output parameters for that command.",
                ],
                required_fields=["tool_name", "required_args"],
                example_payload={
                    "tool_name": "bayesprism",
                    "required_args": {
                        "input": "./data/bulk.tsv",
                        "output": "./results/bayesprism",
                    },
                },
            )

        required_args = [
            name
            for name, spec in tool.parameters.items()
            if spec.get("required") and spec.get("action") != "store_true"
        ]

        step = PlanStep(
            id=1,
            tool_name=tool.name,
            mode=tool.mode,
            args={},
            argv=[],
            rationale=f"Matched '{tool.name}' from prompt.",
            expected_outputs=[],
            depends_on=[],
        )

        if tool.mode == "kwargs" and required_args:
            return PlanResult(
                steps=[step],
                need_clarification=True,
                questions=[
                    f"Provide values for required arguments: {', '.join(required_args)}",
                ],
                required_fields=required_args,
                example_payload={
                    "tool_name": tool.name,
                    "args": {arg: f"<{arg}>" for arg in required_args},
                },
            )

        if tool.mode == "argv":
            return PlanResult(
                steps=[step],
                need_clarification=True,
                questions=[
                    "Provide the raw CLI arguments list for this tool.",
                ],
                required_fields=["argv"],
                example_payload={
                    "tool_name": tool.name,
                    "argv": ["--input", "./data/input.bam", "--out", "./results"],
                },
            )

        return PlanResult(steps=[step])

    @staticmethod
    def _find_tool(prompt: str, registry: ToolRegistry) -> Optional[object]:
        prompt_lower = prompt.lower()
        for tool in registry.list_tools():
            if tool.name.lower() in prompt_lower:
                return tool
        return None


class SmartPlanner(Planner):
    """Rule-based planner that selects tools and extracts parameters without LLMs."""

    _THREADS_RE = re.compile(r"(?:threads?|线程|thread)\s*[:=]?\s*(\d+)", re.IGNORECASE)
    _OUTDIR_RE = re.compile(
        r"(?:output to|outdir|output|save to|输出到|结果放到|结果输出到)\s+(\S+)",
        re.IGNORECASE,
    )
    _FASTQ_HINT_RE = re.compile(r"(fastq|fq)", re.IGNORECASE)
    _PIPELINE_RE = re.compile(r"(pipeline|run all|complete|全套|完整|runall)", re.IGNORECASE)
    _MODE_RE = re.compile(r"\b(salmon|star)\b", re.IGNORECASE)
    _TME_RE = re.compile(r"\btme\b|肿瘤微环境", re.IGNORECASE)
    _CLUSTER_RE = re.compile(r"cluster|聚类", re.IGNORECASE)
    _HLA_RE = re.compile(r"\bhla\b", re.IGNORECASE)
    _SPECHLA_RE = re.compile(r"spechla", re.IGNORECASE)
    _BAYESPRISM_RE = re.compile(r"bayesprism|\bbp\b", re.IGNORECASE)
    _QC_RE = re.compile(r"qc|质控|fastp", re.IGNORECASE)

    def plan(
        self,
        prompt: str,
        registry: ToolRegistry,
        context: Optional[Dict[str, object]] = None,
    ) -> PlanResult:
        prompt = (prompt or "").strip()
        context = context or {}

        if not prompt and not context:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=["Please describe the analysis task you want to run."],
                required_fields=["tool_name", "inputs", "outputs"],
                example_payload={
                    "tool_name": "tme_profile",
                    "inputs": {"input": "./data/tpm.tsv"},
                    "outputs": {"output": "./results/tme"},
                },
            )

        tool = self._select_tool(prompt, registry)
        if tool is None:
            return PlanResult(
                steps=[],
                need_clarification=True,
                questions=[
                    "无法识别要运行的子命令。请说明要使用的 iobrpy 命令名称。",
                ],
                required_fields=["tool_name"],
                example_payload={
                    "tool_name": "runall",
                    "args": {"mode": "salmon", "fastq": "./data/fastq", "outdir": "./out"},
                },
            )

        args = {}
        argv: List[str] = []
        extracted = self._extract_args(prompt, tool)
        extra_argv = extracted.pop("_unknown_args", [])
        args.update(extracted)
        if extra_argv:
            argv.extend([str(item) for item in extra_argv])

        if context.get("args") and isinstance(context["args"], dict):
            args.update(context["args"])
        if context.get("argv") and isinstance(context["argv"], list):
            argv.extend([str(item) for item in context["argv"]])

        missing = self._missing_required(tool, args)
        questions = []
        if missing:
            questions = self._build_questions(missing, tool, prompt)

        step = PlanStep(
            id=1,
            tool_name=tool.name,
            mode=tool.mode,
            args=args,
            argv=argv,
            rationale=f"Rule-based selection for '{tool.name}'.",
            expected_outputs=self._expected_outputs(tool, args),
            depends_on=[],
        )

        return PlanResult(
            steps=[step],
            need_clarification=bool(missing),
            questions=questions,
            required_fields=missing,
            example_payload=self._example_payload(tool, missing),
        )

    def _select_tool(self, prompt: str, registry: ToolRegistry) -> Optional[ToolSpec]:
        prompt_lower = prompt.lower()
        for tool in registry.list_tools():
            if tool.name == "ai":
                continue
            pattern = r"\\b" + re.escape(tool.name.lower()) + r"\\b"
            if re.search(pattern, prompt_lower):
                return tool

        if self._BAYESPRISM_RE.search(prompt):
            return registry.get("bayesprism")
        if self._HLA_RE.search(prompt):
            if self._SPECHLA_RE.search(prompt):
                return registry.get("spechla")
            if re.search(r"typing", prompt, re.IGNORECASE):
                return registry.get("hla_typing")
            return registry.get("hla_typing")
        if self._TME_RE.search(prompt):
            if self._CLUSTER_RE.search(prompt):
                return registry.get("tme_cluster")
            return registry.get("tme_profile")
        if self._QC_RE.search(prompt):
            return registry.get("fastq_qc")

        if self._PIPELINE_RE.search(prompt) and self._FASTQ_HINT_RE.search(prompt):
            return registry.get("runall")
        if self._FASTQ_HINT_RE.search(prompt):
            if re.search(r"salmon", prompt, re.IGNORECASE):
                return registry.get("batch_salmon")
            if re.search(r"star", prompt, re.IGNORECASE):
                return registry.get("batch_star_count")
        if re.search(r"merge\\s+salmon", prompt, re.IGNORECASE):
            return registry.get("merge_salmon")
        if re.search(r"merge\\s+star", prompt, re.IGNORECASE):
            return registry.get("merge_star_count")

        return None

    def _extract_args(self, prompt: str, tool: ToolSpec) -> Dict[str, object]:
        args: Dict[str, object] = {}
        outdir = self._extract_outdir(prompt)
        threads = self._extract_threads(prompt)
        mode = self._extract_mode(prompt)
        paths = self._extract_paths(prompt)

        if tool.name == "runall":
            if mode:
                args["mode"] = mode
            if outdir:
                args["outdir"] = outdir
            fastq = self._pick_path(paths, prompt, hint="fastq")
            if fastq:
                args["fastq"] = fastq
            if threads:
                args["_unknown_args"] = ["--threads", str(threads)]
        elif tool.name == "fastq_qc":
            if outdir:
                args["path2_fastp"] = outdir
            fastq = self._pick_path(paths, prompt, hint="fastq")
            if fastq:
                args["path1_fastq"] = fastq
            if threads:
                args["num_threads"] = threads
        elif tool.name == "batch_salmon":
            if outdir:
                args["path_out"] = outdir
            fastq = self._pick_path(paths, prompt, hint="fastq")
            if fastq:
                args["path_fq"] = fastq
            if threads:
                args["num_threads"] = threads
        elif tool.name == "merge_salmon":
            if outdir:
                args["project"] = outdir
            if paths:
                args["path_salmon"] = paths[0]
            if threads:
                args["num_processes"] = threads
        elif tool.name == "batch_star_count":
            if outdir:
                args["path_out"] = outdir
            fastq = self._pick_path(paths, prompt, hint="fastq")
            if fastq:
                args["path_fq"] = fastq
            if threads:
                args["num_threads"] = threads
        elif tool.name in {"tme_profile", "tme_cluster"}:
            if outdir:
                args["output"] = outdir
            if paths:
                args["input"] = paths[0]
            if threads and tool.name == "tme_profile":
                args["threads"] = threads
        elif tool.name == "bayesprism":
            if outdir:
                args["output"] = outdir
            if paths:
                args["input"] = paths[0]
            if threads:
                args["threads"] = threads

        return args

    def _missing_required(self, tool: ToolSpec, args: Dict[str, object]) -> List[str]:
        missing = []
        for name, spec in tool.parameters.items():
            if spec.get("required") and spec.get("action") != "store_true":
                if args.get(name) in {None, ""}:
                    missing.append(name)
        return missing

    def _build_questions(self, missing: List[str], tool: ToolSpec, prompt: str) -> List[str]:
        questions = []
        for name in missing:
            spec = tool.parameters.get(name, {})
            help_text = spec.get("help") or ""
            choices = spec.get("choices")
            default = spec.get("default")
            flags = ", ".join(spec.get("flags") or [])
            parts = [f"需要参数 '{name}' ({flags})。"]
            if help_text:
                parts.append(f"说明：{help_text}")
            if choices:
                parts.append(f"可选值：{', '.join(map(str, choices))}")
            if default not in (None, "None"):
                parts.append(f"默认值：{default}")
            questions.append(" ".join(parts))

        if "`" in prompt:
            questions.append(
                "注意：在 bash 中反引号会触发命令替换，请不要用 `...` 包路径；"
                "可改用普通引号 '...' 或直接写路径。"
            )

        return questions

    def _expected_outputs(self, tool: ToolSpec, args: Dict[str, object]) -> List[str]:
        if tool.name in {"runall", "tme_profile", "bayesprism", "batch_salmon"}:
            output = args.get("outdir") or args.get("output") or args.get("path_out")
            if output:
                return [str(output)]
        return []

    def _example_payload(self, tool: ToolSpec, missing: List[str]) -> Optional[Dict[str, object]]:
        if not missing:
            return None
        return {
            "tool_name": tool.name,
            "args": {name: f"<{name}>" for name in missing},
        }

    @classmethod
    def _extract_threads(cls, prompt: str) -> Optional[int]:
        match = cls._THREADS_RE.search(prompt)
        if not match:
            return None
        try:
            return int(match.group(1))
        except ValueError:
            return None

    @classmethod
    def _extract_outdir(cls, prompt: str) -> Optional[str]:
        match = cls._OUTDIR_RE.search(prompt)
        if match:
            return cls._clean_path(match.group(1))
        return None

    @classmethod
    def _extract_mode(cls, prompt: str) -> Optional[str]:
        match = cls._MODE_RE.search(prompt)
        if match:
            return match.group(1).lower()
        return None

    @staticmethod
    def _clean_path(value: str) -> str:
        value = value.strip().strip(",.;")
        if (value.startswith("'") and value.endswith("'")) or (
            value.startswith('"') and value.endswith('"')
        ):
            return value[1:-1]
        return value

    @classmethod
    def _extract_paths(cls, prompt: str) -> List[str]:
        candidates = re.findall(r"(?:(?:/|\./|\~/)[^\s]+|[A-Za-z]:\\[^\s]+)", prompt)
        return [cls._clean_path(item) for item in candidates]

    @classmethod
    def _pick_path(cls, paths: List[str], prompt: str, hint: str) -> Optional[str]:
        if not paths:
            return None
        if hint and hint in prompt.lower():
            return paths[0]
        return paths[0]
