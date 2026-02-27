# IOBRpy AI 架构说明（纯“用户 ↔ iobrpy ai”对话版）

> **用户**：我完全不会 iobrpy 命令，`iobrpy ai` 能做什么？  
> **iobrpy ai**：你只要描述目标，我负责把自然语言需求转成可执行命令，并在确认后执行。

---

> **用户**：我有 FASTQ，想做完整 TME 分析。  
> **iobrpy ai**：建议使用 `runall`。请先提供 FASTQ 目录。  
> **用户**：`/data/fastq`。  
> **iobrpy ai**：请提供输出目录。  
> **用户**：`/data/out`。  
> **iobrpy ai**：请选择模式：`salmon` 或 `star`。  
> **用户**：`salmon`。  
> **iobrpy ai**：请继续提供 `index`、`threads`、`batch_size`、`project`。  
> **用户**：`/ref/salmon_index`、`8`、`2`、`demo_proj`。  
> **iobrpy ai**：参数已满足，请确认草拟命令：  
> `python -m iobrpy.main runall --fastq /data/fastq --outdir /data/out --mode salmon --index /ref/salmon_index --threads 8 --batch_size 2 --project demo_proj`  
> **用户**：确认执行。  
> **iobrpy ai**：开始执行，日志路径：`<logdir>/<session>_runall.log`。

---

> **用户**：如果我只说“帮我跑 runall”，你会直接执行吗？  
> **iobrpy ai**：不会。必填参数不完整时，我会继续追问，直到满足执行条件。  
> **用户**：如果我给了这个子命令不支持的参数呢？  
> **iobrpy ai**：我会过滤该参数，不会写入草拟命令。  
> **用户**：如果执行失败呢？  
> **iobrpy ai**：我会返回失败状态，并建议你查看 `<logdir>/<session>_<subcommand>.log`。

---

> **用户**：你内部是怎么分工的？  
> **iobrpy ai**：我会根据 `src/iobrpy/RAG_MCP/iobrpy_required_params.json` 的参数契约做参数收敛，并通过 `src/iobrpy/RAG_MCP/iobrpy_rag_mcp.py` 完成意图识别、参数抽取、命令组装与状态推进，再由 `src/iobrpy/RAG_MCP/ai.py` 负责交互循环与真实执行。  
> **用户**：状态推进是什么意思？  
> **iobrpy ai**：典型路径是 `need_info -> ready -> done/error`。  
> **用户**：执行闭环怎么理解？  
> **iobrpy ai**：意图识别 -> 参数收敛 -> 草拟命令确认 -> 执行与日志回传。

---

> **用户**：部署时我主要要关注哪些配置？  
> **iobrpy ai**：重点包括 `CHROMA_DIR`、`IOBRPY_REQUIRED_PARAMS_FILE`、`IOBRPY_RUN_LOG_DIR`、`IOBRPY_DEFAULTS_FILE`、`OLLAMA_HOST`、`CHAT_MODEL`、`EMBED_MODEL`。  
> **用户**：这些配置的底线要求是什么？  
> **iobrpy ai**：模型版本建议固定，规则与代码版本建议联动发布，日志目录必须可写。

---

> **用户**：你怎么保证不乱执行？  
> **iobrpy ai**：我依赖三层控制：参数白名单、非猜测策略、执行前确认。  
> **用户**：非猜测策略是什么？  
> **iobrpy ai**：用户未明确提供的关键值，我不会主观补全。  
> **用户**：也就是说你是“先问清楚再跑”？  
> **iobrpy ai**：是的，这是默认执行原则。

---

> **用户**：下一步 API 化怎么做？  
> **iobrpy ai**：会先建设统一 API 接入层，再固化请求/响应 JSON Schema 与错误码规范。  
> **用户**：安全侧呢？  
> **iobrpy ai**：会增加 API Key/Token 鉴权与基础限流。  
> **用户**：可观测性呢？  
> **iobrpy ai**：会提供任务状态与日志查询 API，并把会话状态机映射到 API 层。  
> **用户**：最后给我一句总结。  
> **iobrpy ai**：你描述目标，我补齐参数、给出草拟命令、确认后执行，并保证全过程可审查、可追溯。
