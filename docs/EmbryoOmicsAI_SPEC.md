# EmbryoOmics AI — 技术规格（SPEC v0.1）

日期：2026-09-27。对应：`EmbryoOmicsAI_PRD.md`。本文件规定 v0.1 可开发和可测试的契约，并约定后续质谱蛋白组接入标准；后续更改需记录原因。

## 1. 架构及边界

Python 分析包 `embryoomics` 独立于本地 Streamlit 界面。命令行与网页调用同一服务入口 `run_analysis(config)`，返回产物路径及运行摘要。v0.1 不依赖 GPU，不调用远程模型也可完整运行。

```text
embryoomics/
  io.py              # 读取和校验
  qc.py              # 指标和过滤预览
  pipeline.py        # Scanpy 主流程
  evidence.py        # marker 与人工注释记录
  reporting.py       # HTML/JSON 导出
  config.py          # 配置校验
  ai.py              # 可选，受约束的解释接口
app.py               # Streamlit
tests/               # 小矩阵契约测试与一次端到端冒烟测试
docs/                # 需求与使用说明
```

依赖：Python 3.11；固定兼容版本的 Scanpy、AnnData、NumPy、SciPy、Pandas、Matplotlib、Leiden 实现、Streamlit。安装前用目标环境核对版本并记录锁定文件。v0.2 才引入 MuData/scvi-tools。无需为 v0.1 安装 CellTypist 或强行套用人类免疫细胞模型给小鼠胚胎注释。

## 2. 输入契约

### 数据

- 10x 三件套：`.mtx[.gz]`、`features.tsv[.gz]`、`barcodes.tsv[.gz]`；解析后 `cells × genes`，稀疏非负整数矩阵。
- `.h5ad`：配置明确 `raw_count_source: X | layers/<name> | raw/X`；检测浮点/对数化数值时拒绝当作原始计数。
- 可选 `metadata.csv`：`cell_id,sample_id,stage,embryo_id,batch`；前三列必备，可由单样本配置补充 `sample_id` 与 `stage`，但不能猜出 `embryo_id`。细胞 ID 与矩阵条码精确匹配，重复和未匹配一律报告。
- 配置明确 `species: mus_musculus`、`genome_build`（若未知则 `unknown` 并在报告中披露）、`source_url`、`dataset_accession`、`sample_id`、`stage`、`seed`、QC 参数、Leiden resolution。

### 示例配置

```yaml
dataset_accession: GSE278981
sample_id: selected_e7_sample
stage: E7.0
species: mus_musculus
genome_build: unknown
source_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278981
raw_count_source: X
seed: 42
qc:
  min_genes: null   # null = 显示分布，由用户在数据审查后填写
  max_genes: null
  max_pct_mt: null # 线粒体基因匹配不可靠时禁用该过滤器
analysis:
  n_top_genes: 2000
  n_neighbors: 15
  leiden_resolution: 0.5
```

## 3. 处理流程与产物

1. **validate**：格式、ID 唯一性、计数整数性、维度、元数据匹配、非空样本；返回错误代码和说明。
2. **QC**：计算总计数和检测基因数；按小鼠基因符号或可核实的 gene ID 映射识别线粒体基因，记录识别方法和数量。展示按样本分布，产出过滤预览；用户提交阈值后过滤，不暗中使用固定值。
3. **RNA**：在 `layers['counts']` 保留稀疏原始计数；对工作对象 normalize_total、log1p，记录具体 HVG 方法与参数，再 PCA、neighbors、Leiden、UMAP；必要时对 HVG 子集建模但保留完整基因 marker 分析视图。
4. **证据**：每 cluster 输出 marker 统计（基因 ID、符号、组内比例、参考组比例、效应量和检验方法）；同时输出阶段、样本、胚胎组成，不因一张 UMAP 自动生成谱系结论。
5. **报告**：HTML 页面包含数据清单、QC、分析流程、图、候选注释、局限性。除图像外输出 `run_manifest.json`、`qc_metrics.csv`、`markers.csv`、`annotations.csv`、`processed.h5ad`。

## 4. 注释与可选 AI 规则

- v0.1 人工注释字段：`cluster_id`, `candidate_cell_type`, `support_genes`, `conflict_genes`, `reference_url`, `review_status`, `reviewer_note`。
- 如启用 LLM，只传递聚合的 marker 与图表说明，默认不上传逐细胞矩阵；显示发送内容预览；API 密钥只从环境变量读取。
- AI 输出必须标记来源 `observed / literature / inferred / unknown`。无法用提供的 marker 和引用支持的结论返回 `unknown`；不生成数字化可信度，除非另有独立校准数据。
- AI 出错时报告仍能生成；人工确认的注释优先，AI 文本保留原始建议和审核状态。

## 5. 执行入口及目录

```text
embryoomics validate --input DATA --config config.yaml
embryoomics run --input DATA --config config.yaml --out results/run_001
streamlit run app.py
```

`results/run_001` 含上述产物和 `figures/`。重复运行默认拒绝覆盖同名目录，允许新目录/显式 `--overwrite`。`run_manifest.json` 至少包含输入文件名与 SHA256、数据来源、输入维度、各 QC 阶段细胞数、参数、软件版本、随机种子、警告与实际生成的文件清单。

## 6. 验收用例

| 情况 | 预期结果 |
| --- | --- |
| 小型合成稀疏整数矩阵与匹配元数据 | 验证通过，QC 统计正确；可运行到报告生成或在样本太少时给清晰提示 |
| 重复条码或元数据缺失细胞 | 拒绝分析，指出重复或缺失数与示例 ID |
| 归一化浮点矩阵当原始计数 | 拒绝分析，并解释如何指定正确 counts 层 |
| 线粒体基因零匹配 | 报告检测失败，跳过线粒体百分比阈值，不输出伪零值 |
| 没有 API 密钥 | 基础分析与报告正常完成，AI 面板说明未启用 |
| 正式 GSE278981 单个样本 | 本地端到端运行，输出四页内容和所有可下载产物，记录真实耗时、峰值内存和原始文件版本 |

## 7. v0.2 质谱蛋白组接入闸门

- 数据候选为小鼠真实胚胎 E7.5、E8.0、E8.5、E9.5 蛋白组，论文 https://doi.org/10.1016/j.stem.2024.04.017，PRIDE PXD041309。先在仓库中核实是否有可直接使用的蛋白强度表和完整样本注释；不得把同项目的体外 gastruloid 样本混作胚胎。E7.5 蛋白组为整胚组织，多个胚胎合并成一个生物学重复。
- 接口输入为 `protein_abundance.csv`（行=蛋白 ID，列=质谱样本 ID）与 `proteomics_samples.csv`（质谱样本 ID、阶段、物种、整胚/组织区域、体内/体外、实验批次、生物学重复、合并胚胎数、定量方法）。蛋白组基因名映射要保留原蛋白组 ID、多重映射和未映射项；绝不把缺失值当成零表达。
- 分开质控样本总强度、检出蛋白数、缺失率、批次与重复相似性；归一化、对数化和缺失值处理遵从论文的定量方法并记录。磷酸化位点单设表，不与全蛋白丰度直接混算。
- 与单细胞 RNA 比较时，先按**取材区域、阶段、实际胚胎/样本**聚合 RNA，再标记对齐等级 `same_biological_sample / same_embryo_different_cells / same_stage_region / same_stage_only / incomparable`。不同研究的相同阶段仅代表探索性对照，不能计算假装配对的逐样本相关或差异统计；E7.5 后部 RNA 与 E7.5 整胚蛋白组最多为 `same_stage_only`。
- 导出 `cross_modality_manifest.csv`，逐对列出研究、物种、阶段、组织、胚胎 ID（若有）、样本 ID、重复、对齐等级、可执行分析、被禁止的推断，并在报告中将 RNA/蛋白不一致列为候选问题。

## 8. v0.3 RNA＋ATAC 设计闸门

- 先用仓库 https://github.com/czbiohub-sf/zebrahub-multiome-analysis 的处理数据核对细胞条码、模态匹配、阶段、胚胎和基因组构建；建立 `pairing_report.json`。
- 同细胞 RNA/ATAC 可做联合嵌入；样本级配对只做样本层级对照；仅同阶段的数据不得按行拼接。多批次 ATAC 要统一峰集合。ATAC 的同细胞配对不自动扩展到另一个研究的质谱蛋白组。
- 多组学潜在空间模型可按官方教程 https://docs.scvi-tools.org/en/latest/tutorials/notebooks/multimodal/MultiVI_tutorial.html 实施；模型预测或补齐数据必须与实测数据分开标注。
- 峰—基因关联作为候选假设，记录配对方式、距离/统计方法、负例或置换对照和重复结构，避免称为因果关系。

## 9. 明确的开发取舍

首版做本地 Streamlit＋CLI，不搭数据库和服务器；优先保留稀疏矩阵，避免把整个 RNA 矩阵转密集；导入大文件前提示预计资源；公开数据放在用户本地 `data/`，`.gitignore` 排除数据、密钥与输出。代码和报告附研究与数据集来源及许可信息。
