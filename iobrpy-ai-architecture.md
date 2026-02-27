# IOBRpy AI 架构说明（全函数、参数按 site 文档）

> **用户**：这份文档怎么保证不乱写参数？  
> **iobrpy ai**：每个 function 的对话都只使用 site 文档里出现的参数名与用法；不编造额外参数。

---

> **用户**：统一入口和状态是什么？  
> **iobrpy ai**：入口：`iobrpy ai --logdir <dir>`；状态：`need_info -> ready -> done/error`。

---

## 1) prepare_salmon
> **用户**：我要清理 salmon 合并结果。  
> **iobrpy ai**：请给 `-i/--input` 和 `-o/--output`。  
> **用户**：`-i /data/merge_salmon.tsv -o /data/prepare_salmon.tsv`。  
> **iobrpy ai**：是否设置 `-r/--return_feature {ENST|ENSG|symbol}` 或 `--remove_version`？  
> **用户**：`-r symbol`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_prepare_salmon.log`。

## 2) count2tpm
> **用户**：我要把 count 转 TPM。  
> **iobrpy ai**：请给 `-i/--input` 与 `-o/--output`。  
> **用户**：`-i /data/count.tsv -o /data/tpm.tsv`。  
> **iobrpy ai**：可选：`--idtype`、`--org`、`--effLength_csv`、`--remove_version`。  
> **用户**：`--idtype ensembl --org hsa`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_count2tpm.log`。

## 3) anno_eset
> **用户**：我要注释表达矩阵。  
> **iobrpy ai**：请给 `-i/--input`、`-o/--output`，以及 `--annotation`（或 `--annotation-file`）。  
> **用户**：`-i /data/eset.tsv -o /data/eset_anno.tsv --annotation anno_rnaseq`。  
> **iobrpy ai**：是否指定 `--symbol`、`--probe`、`--method`？  
> **用户**：`--method mean`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_anno_eset.log`。

## 4) calculate_sig_score
> **用户**：我要做 signature 评分。  
> **iobrpy ai**：请给输入输出，并指定 `--method {integration|pca|zscore|ssgsea}`。  
> **用户**：`-i /data/tpm.tsv -o /data/sig.tsv --method ssgsea`。  
> **iobrpy ai**：是否指定 `--signature`？  
> **用户**：`--signature hallmark`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_calculate_sig_score.log`。

## 5) cibersort
> **用户**：我要跑 CIBERSORT。  
> **iobrpy ai**：请给 `-i/--input` 与 `-o/--output`。  
> **用户**：`-i /data/tpm.tsv -o /data/cibersort.tsv`。  
> **iobrpy ai**：是否设置 `--perm` 与 `--QN`？  
> **用户**：`--perm 100 --QN true`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_cibersort.log`。

## 6) IPS
> **用户**：我要算 IPS。  
> **iobrpy ai**：请给 `-i/--input` 和 `-o/--output`。  
> **用户**：`-i /data/tpm.tsv -o /data/ips.tsv`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_IPS.log`。

## 7) estimate
> **用户**：我要跑 ESTIMATE。  
> **iobrpy ai**：请给 `-i/--input`、`-o/--output`，可选 `-p/--platform`。  
> **用户**：`-i /data/tpm.tsv -o /data/estimate.tsv -p affymetrix`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_estimate.log`。

## 8) mcpcounter
> **用户**：我要跑 MCPcounter。  
> **iobrpy ai**：请给 `-i/--input`、`-f/--features`、`-o/--output`。  
> **用户**：`-i /data/tpm.tsv -f HUGO_symbols -o /data/mcpcounter.tsv`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_mcpcounter.log`。

## 9) quantiseq
> **用户**：我要跑 quanTIseq。  
> **iobrpy ai**：请给 `-i/--input` 和 `-o/--output`。  
> **用户**：`-i /data/tpm.tsv -o /data/quantiseq.tsv`。  
> **iobrpy ai**：是否设置 `--method`、`--arrays`、`--tumor`、`--scale_mrna`？  
> **用户**：`--method lsei`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_quantiseq.log`。

## 10) epic
> **用户**：我要跑 EPIC。  
> **iobrpy ai**：请给 `-i/--input`、`-o/--output`，可选 `--reference {TRef|BRef|both}`。  
> **用户**：`-i /data/tpm.tsv -o /data/epic.tsv --reference TRef`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_epic.log`。

## 11) deside
> **用户**：我要跑 DeSide。  
> **iobrpy ai**：必填是 `-m/--model_dir`、`-i/--input`、`-o/--output`。  
> **用户**：`-m /ref/deside_model -i /data/tpm.tsv -o /data/deside.csv`。  
> **iobrpy ai**：可选 `--exp_type`、`--gmt`、`--transpose`。  
> **用户**：`--exp_type TPM`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_deside.log`。

## 12) tme_cluster
> **用户**：我要做 tme_cluster。  
> **iobrpy ai**：请给 `-i/--input` 与 `-o/--output`。  
> **用户**：`-i /data/tme.tsv -o /data/tme_cluster.tsv`。  
> **iobrpy ai**：可选 `--features/--pattern`、`--min_nc`、`--max_nc`。  
> **用户**：`--features 2:20`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_tme_cluster.log`。

## 13) LR_cal
> **用户**：我要跑 LR_cal。  
> **iobrpy ai**：请给 `-i/--input` 和 `-o/--output`。  
> **用户**：`-i /data/tpm.tsv -o /data/lr.tsv`。  
> **iobrpy ai**：可选 `--data_type`、`--id_type`、`--cancer_type`。  
> **用户**：`--data_type tpm --id_type ensembl --cancer_type pancan`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_LR_cal.log`。

## 14) nmf
> **用户**：我要做 NMF。  
> **iobrpy ai**：请给 `-i/--input` 与 `-o/--output`。  
> **用户**：`-i /data/tme.tsv -o /data/nmf_out`。  
> **iobrpy ai**：可选 `--kmin`、`--kmax`、`--features`。  
> **用户**：`--kmin 2 --kmax 8`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_nmf.log`。

## 15) mouse2human_eset
> **用户**：我要做 mouse2human_eset。  
> **iobrpy ai**：请给 `-i/--input`、`-o/--output`。若表格模式还需 `--column_of_symbol`。  
> **用户**：`-i /data/mouse.tsv -o /data/human.tsv --is_matrix`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_mouse2human_eset.log`。

## 16) batch_salmon
> **用户**：我要批量跑 salmon。  
> **iobrpy ai**：请给 `--index`、`--path_fq`、`--path_out`。  
> **用户**：`--index /ref/salmon_idx --path_fq /data/fastq --path_out /data/salmon_out`。  
> **iobrpy ai**：可选 `--num_threads`、`--batch_size`、`--suffix1`。  
> **用户**：`--num_threads 8 --batch_size 2`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_batch_salmon.log`。

## 17) merge_salmon
> **用户**：我要合并 salmon quant 结果。  
> **iobrpy ai**：请给 `--path_salmon` 与 `--project`。  
> **用户**：`--path_salmon /data/salmon_out --project Proj1`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_merge_salmon.log`。

## 18) merge_star_count
> **用户**：我要合并 STAR count。  
> **iobrpy ai**：请给 `--path` 与 `--project`。  
> **用户**：`--path /data/star_out --project Proj1`。确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_merge_star_count.log`。

## 19) batch_star_count
> **用户**：我要批量跑 STAR count。  
> **iobrpy ai**：请给 `--index`、`--path_fq`、`--path_out`。  
> **用户**：`--index /ref/star_idx --path_fq /data/fastq --path_out /data/star_out`。  
> **iobrpy ai**：可选 `--num_threads`、`--batch_size`。  
> **用户**：`--num_threads 16 --batch_size 2`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_batch_star_count.log`。

## 20) fastq_qc
> **用户**：我只做 fastq 质控。  
> **iobrpy ai**：请给 `--path1_fastq` 与 `--path2_fastp`。  
> **用户**：`--path1_fastq /data/fastq --path2_fastp /data/qc_out`。  
> **iobrpy ai**：可选 `--num_threads`、`--batch_size`、`--suffix1`。  
> **用户**：`--num_threads 8 --batch_size 5`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_fastq_qc.log`。

## 21) log2_eset
> **用户**：我要做 log2_eset。  
> **iobrpy ai**：请给 `-i/--input` 和 `-o/--output`。  
> **用户**：`-i /data/eset.tsv -o /data/eset_log2.tsv`，确认执行。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_log2_eset.log`。

## 22) trust4
> **用户**：我要跑 trust4。  
> **iobrpy ai**：按 site 文档，可用 BAM 输入 `-b`，或 FASTQ 输入 `-1/-2`，或批量 `--fqdir`；并给 `-o`、`-t`。  
> **用户**：`-b /data/sample.bam -o /data/trust4_out -t 8`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_trust4.log`。

## 23) spechla
> **用户**：我要跑 spechla。  
> **iobrpy ai**：请给 `-n`、`-1`、`-2`、`-o`，可选 `-j`。  
> **用户**：`-n SAMPLE1 -1 /data/SAMPLE1_1.fq.gz -2 /data/SAMPLE1_2.fq.gz -o /data/spechla_out -j 8`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_spechla.log`。

## 24) extract_hla_read
> **用户**：我要先抽取 HLA reads。  
> **iobrpy ai**：请给 `-s`、`-b`、`-r`、`-o`。  
> **用户**：`-s SAMPLE1 -b /data/bam/SAMPLE1.bam -r hg38 -o /data/ExtractHLAread/SAMPLE1`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_extract_hla_read.log`。

## 25) hla_typing
> **用户**：我要批量做 HLA typing。  
> **iobrpy ai**：请给 `-b`、`-r`、`-o`，可选 `-j`。  
> **用户**：`-b /data/bam -r hg38 -o /data/hla_typing_results -j 8`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_hla_typing.log`。

## 26) tme_profile
> **用户**：我有 TPM，想跑 tme_profile。  
> **iobrpy ai**：请给 `-i/--input` 与 `-o/--output`，可选 `--threads`。  
> **用户**：`-i /data/tpm.tsv -o /data/tme_profile_out --threads 1`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_tme_profile.log`。

## 27) runall
> **用户**：我要从 FASTQ 一条龙跑完。  
> **iobrpy ai**：请给 `--mode`、`--outdir`、`--fastq`、`--threads`、`--batch_size`、`--index`、`--project`。  
> **用户**：`--mode salmon --outdir /data/runall_out --fastq /data/fastq --threads 8 --batch_size 1 --index /ref/salmon_idx --project MyProj`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_runall.log`。

## 28) bayesprism
> **用户**：我要跑 bayesprism。  
> **iobrpy ai**：请给 `-i/--input`、`-o/--output`；如需单细胞先验，再给 `--sc_dat`、`--cell_state_labels`、`--cell_type_labels`（可选 `--key`）。  
> **用户**：`-i /data/bulk.tsv -o /data/bp_out --sc_dat /data/sc.csv --cell_state_labels /data/state.csv --cell_type_labels /data/type.csv --key Tumor`。  
> **iobrpy ai**：请确认执行。  
> **用户**：确认。  
> **iobrpy ai**：已执行，日志：`<logdir>/<session>_bayesprism.log`。

## 29) ai
> **用户**：`ai` function 本身怎么用？  
> **iobrpy ai**：`iobrpy ai --logdir <dir>` 启动交互。  
> **用户**：会话控制呢？  
> **iobrpy ai**：可用 `:exit` / `:quit` 退出，`:restart` 重置会话。

---

> **用户**：最后总结下你引用的内部文件和配置。  
> **iobrpy ai**：参数契约来自 `src/iobrpy/RAG_MCP/iobrpy_required_params.json`；路由/组装在 `src/iobrpy/RAG_MCP/iobrpy_rag_mcp.py`；执行入口在 `src/iobrpy/RAG_MCP/ai.py`。  
> **用户**：部署配置？  
> **iobrpy ai**：`CHROMA_DIR`、`IOBRPY_REQUIRED_PARAMS_FILE`、`IOBRPY_RUN_LOG_DIR`、`IOBRPY_DEFAULTS_FILE`、`OLLAMA_HOST`、`CHAT_MODEL`、`EMBED_MODEL`。  
> **用户**：最终原则？  
> **iobrpy ai**：只按规则问参，不编参数；先确认再执行；全程日志可追溯。
