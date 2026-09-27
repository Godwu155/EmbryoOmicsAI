"""Presentation translations. Analysis fields and exported values remain stable."""

LANGUAGES = {"English": "en", "简体中文": "zh-CN"}

TEXT = {
    "en": {
        "subtitle": "Local mouse embryo scRNA-seq review. Cell types need human confirmation; UMAP is not a trajectory.",
        "data": "Data", "qc": "QC", "results": "Results", "report": "Report",
        "input_path": "10x folder or .h5ad path", "config_path": "Config YAML path",
        "metadata_path": "Optional metadata CSV path", "out_path": "New results folder",
        "configuration": "Configuration", "field": "Field", "value": "Value",
        "validate": "Validate input", "mt_detection": "Mitochondrial detection",
        "mt_method": "Detection method", "mt_matched": "Matched genes",
        "mt_symbol_method": "Mouse gene symbols, case-insensitive ^mt-",
        "mt_no_match": "No matching mouse gene symbols",
        "qc_intro": "Review QC distributions and choose thresholds. A run requires explicit confirmation.",
        "qc_thresholds": "Configured thresholds", "qc_edit": "Edit the YAML, then validate again to update the preview.",
        "qc_confirm": "I reviewed the QC preview and confirm these thresholds", "run": "Run analysis",
        "completed": "Completed", "validate_first": "Validate input on the Data tab first.",
        "annotation_hint": "Candidate annotations default to unknown. Export and review the table before reuse.",
        "annotation_download": "Download reviewed annotations CSV",
        "annotation_rerun": "Pass the reviewed file to the CLI with --annotations for a new reproducible run.",
        "run_first": "Run an analysis to see results.", "report_location": "HTML report",
        "ai_disabled": "AI interpretation is disabled in v0.1; all analysis and reporting run locally without an API key.",
        "report_first": "No report in the selected results folder yet.",
        "chart_count": "Cell count",
        "unknown": "unknown", "unset": "Not set",
        "report_title": "EmbryoOmics AI report", "report_results": "Results and marker evidence",
        "report_review": "Review and limitations", "downloads": "Downloads",
        "dataset": "Dataset", "sample": "Sample", "stage": "Stage", "source": "Source",
        "geo_record": "GEO record", "input_dimensions": "Input dimensions", "cells": "cells", "genes": "genes",
        "genome": "Genome", "license": "Data reuse license", "seed": "Analysis seed",
        "ai_status": "AI interpretation: disabled.", "qc_figure": "QC distributions before filtering",
        "umap_figure": "UMAP by cluster, stage and sample",
        "qc_axes": "Left to right: total counts, detected genes, mitochondrial percentage.",
        "umap_axes": "Left to right: cluster, stage, sample.",
        "umap_limit": "UMAP shows similarity in this analysis; it is not a developmental trajectory.",
        "review_limit": "Cell type labels require human review. Cluster markers are exploratory and cells are not biological replicates. No lineage, causal, or cross-modality claim is made.",
    },
    "zh-CN": {
        "subtitle": "本地小鼠胚胎单细胞 RNA 分析。细胞类型需人工确认；UMAP 不代表发育轨迹。",
        "data": "数据", "qc": "质控", "results": "结果", "report": "报告",
        "input_path": "10x 文件夹或 .h5ad 路径", "config_path": "配置 YAML 路径",
        "metadata_path": "可选元数据 CSV 路径", "out_path": "新结果文件夹",
        "configuration": "配置", "field": "字段", "value": "值",
        "validate": "验证输入", "mt_detection": "线粒体基因识别",
        "mt_method": "识别方法", "mt_matched": "匹配基因数",
        "mt_symbol_method": "小鼠基因符号，不区分大小写地匹配 ^mt-",
        "mt_no_match": "未匹配到小鼠线粒体基因符号",
        "qc_intro": "查看质控分布并选择阈值。运行分析前须明确确认。",
        "qc_thresholds": "当前阈值", "qc_edit": "修改 YAML 后重新验证，以更新预览。",
        "qc_confirm": "我已查看质控预览并确认这些阈值", "run": "运行分析",
        "completed": "已完成", "validate_first": "请先在“数据”页验证输入。",
        "annotation_hint": "候选注释默认是 unknown（未知）。复用前请导出并审核表格。",
        "annotation_download": "下载已审核注释 CSV",
        "annotation_rerun": "将已审核文件通过 CLI 的 --annotations 参数传入新运行，以保留可复现记录。",
        "run_first": "运行分析后可查看结果。", "report_location": "HTML 报告",
        "ai_disabled": "v0.1 未启用 AI 解释；无需 API 密钥即可在本地完成分析和报告。",
        "report_first": "所选结果文件夹中尚无报告。",
        "chart_count": "细胞数",
        "unknown": "未知", "unset": "未设置",
        "report_title": "EmbryoOmics AI 报告", "report_results": "结果与 marker 证据",
        "report_review": "人工审核与局限", "downloads": "下载文件",
        "dataset": "数据集", "sample": "样本", "stage": "阶段", "source": "来源",
        "geo_record": "GEO 记录", "input_dimensions": "输入维度", "cells": "个细胞", "genes": "个基因",
        "genome": "基因组版本", "license": "数据再利用许可", "seed": "分析随机种子",
        "ai_status": "AI 解释：未启用。", "qc_figure": "过滤前的质控分布",
        "umap_figure": "按聚类、阶段及样本显示的 UMAP",
        "qc_axes": "从左到右：总计数、检出基因数、线粒体计数比例。",
        "umap_axes": "从左到右：聚类、阶段、样本。",
        "umap_limit": "UMAP 显示本次分析中的相似性，不代表发育轨迹。",
        "review_limit": "细胞类型标签须人工审核。聚类 marker 仅供探索，细胞不能视为生物学重复。不据此作谱系、因果或跨模态结论。",
    },
}

COLUMNS = {
    "dataset_accession": "数据集编号", "species": "物种", "genome_build": "基因组版本",
    "source_url": "来源链接", "data_license": "数据再利用许可", "raw_count_source": "原始计数来源",
    "min_genes": "最少检出基因数", "max_genes": "最多检出基因数", "max_pct_mt": "最高线粒体计数比例 (%)",
    "sample_id": "样本 ID", "stage": "阶段", "embryo_id": "胚胎 ID", "cluster": "聚类",
    "input_cells": "输入细胞数", "retained_cells": "保留细胞数", "removed_cells": "去除细胞数",
    "total_counts": "总计数", "n_genes_by_counts": "检出基因数", "pct_counts_mt": "线粒体计数比例 (%)",
    "cluster_id": "聚类 ID", "gene_id": "基因 ID", "gene_symbol": "基因符号",
    "pct_in": "组内检出比例", "pct_reference": "其余细胞检出比例",
    "logfoldchange": "对数倍数变化", "pval_adj": "校正后 P 值", "test_method": "检验方法",
    "candidate_cell_type": "候选细胞类型", "support_genes": "支持基因", "conflict_genes": "冲突基因",
    "reference_url": "参考链接", "review_status": "审核状态", "reviewer_note": "审核备注",
}

WARNINGS = {
    "No mitochondrial gene symbols matched; pct_counts_mt unavailable and mitochondrial filtering disabled.": "未匹配到线粒体基因符号；无法计算线粒体计数比例，线粒体过滤已禁用。",
    "Single-sample clustering and marker tests are exploratory; cell counts are not biological replicate counts.": "单样本聚类与 marker 检验仅供探索；细胞数不能当作生物学重复数。",
    "Genome build is unknown; gene ID mapping is not verified.": "基因组版本未知；基因 ID 映射尚未核实。",
}


def t(language, key):
    return TEXT[language][key]


def columns_for(language):
    return COLUMNS if language == "zh-CN" else {}


def warning_for(language, message):
    return WARNINGS.get(message, message) if language == "zh-CN" else message
