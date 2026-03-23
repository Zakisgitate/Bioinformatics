# Bioinformatics
GSE204844 KDM6A scRNA-seq Pipeline

English

Overview

This repository contains an R-based single-cell RNA-seq analysis pipeline for GEO dataset GSE204844, focusing on KDM6A, inflammasome-related signaling (Nlrp3, Pycard, Casp1, Il1b), and NETosis / neutrophil-associated inflammatory features in a mouse chronic pancreatitis context.

The workflow starts from raw count matrices downloaded from GEO, converts mouse Ensembl IDs to gene symbols, performs standard Seurat v5 preprocessing and clustering, annotates major cell populations, and then carries out focused analyses on macrophage/myeloid cells and neutrophils. The pipeline also generates publication-oriented figures, including UMAP plots, feature plots, module score plots, heatmaps, and custom bubble plots.

Main Functions
	•	Import raw GEO single-cell count matrices
	•	Convert Ensembl gene IDs to mouse gene symbols
	•	Perform quality control and filtering
	•	Run normalization, highly variable gene selection, PCA, UMAP, and clustering
	•	Identify cluster markers and perform rough cell-type annotation
	•	Explore Kdm6a expression together with inflammasome and NETosis-related genes
	•	Focus on Macrophage_Myeloid and Neutrophil subsets for downstream visualization
	•	Export figures in both PNG and SVG formats
	•	Save intermediate Seurat objects and summary matrices for reproducibility

Dataset
	•	GEO accession: GSE204844
	•	Species: Mus musculus
	•	Input files used in this pipeline:
	•	GSM6198711_Case.counts.tsv
	•	GSM6198712_Control.counts.tsv

Analysis Workflow
	1.	Read raw count matrices
	2.	Convert Ensembl IDs to gene symbols
	3.	Build Seurat objects for Case and Control samples
	4.	Perform QC and filter low-quality cells
	5.	Merge filtered objects
	6.	Run standard single-cell workflow:
	•	NormalizeData()
	•	FindVariableFeatures()
	•	ScaleData()
	•	RunPCA()
	•	FindNeighbors()
	•	FindClusters()
	•	RunUMAP()
	7.	Explore target genes and calculate module scores
	8.	Run marker analysis with Seurat v5-compatible JoinLayers()
	9.	Perform rough annotation of major cell populations
	10.	Extract and analyze Macrophage_Myeloid and Neutrophil subsets
	11.	Generate heatmaps and custom bubble plots

Key Biological Focus

This pipeline is designed to address the following exploratory questions:
	•	In which cell populations is Kdm6a mainly expressed?
	•	Is Kdm6a associated with the NLRP3/IL-1β inflammasome axis?
	•	Are NETosis-like inflammatory signatures enriched in neutrophils?
	•	How do macrophage/myeloid cells and neutrophils differ in their inflammatory programs in chronic pancreatitis?

Software Environment

Recommended environment:
	•	R >= 4.3
	•	RStudio
	•	Seurat v5 compatible

Main R packages:
	•	Seurat
	•	SeuratObject
	•	data.table
	•	dplyr
	•	ggplot2
	•	org.Mm.eg.db
	•	AnnotationDbi
	•	pheatmap
	•	tidyr
	•	svglite
	•	grid

Repository Structure

GSE204844_analysis_output/
├── 01_QC/
├── 02_DimReduction/
├── 03_TargetGenes/
├── 04_Annotation/
├── 05_FocusSubset/
├── 06_Heatmap/
├── 07_BubblePlot/
├── obj_case_filtered.rds
├── obj_ctrl_filtered.rds
├── obj_merge_filtered.rds
├── obj_after_umap_cluster.rds
├── obj_after_target_module_score.rds
├── obj_marker_ready.rds
├── obj_after_rough_annotation.rds
├── obj_focus_macrophage_neutrophil.rds
└── obj_focus_joined_for_plot.rds

Notes
	•	This pipeline is written for Seurat v5 layered objects.
	•	JoinLayers() is used before marker analysis and selected downstream expression summarization steps.
	•	The current analysis is exploratory, because it is based on one Case sample and one Control sample.
	•	Conclusions should therefore be interpreted cautiously and ideally validated with additional datasets or experiments.

Citation

If you use this repository, please cite:
	•	The original GEO dataset GSE204844
	•	The Seurat framework
	•	Any additional software packages or downstream methods used in your work

Contact

For questions, suggestions, or collaboration, please open an issue in this repository.

⸻

中文

项目简介

本仓库提供了一套基于 R 的单细胞 RNA 测序分析流程，用于分析 GEO 数据集 GSE204844，重点关注 KDM6A、炎症小体相关信号通路（Nlrp3、Pycard、Casp1、Il1b）以及 NETosis / 中性粒细胞相关炎症特征 在小鼠慢性胰腺炎中的表达模式。

本流程从 GEO 下载的原始 count 矩阵出发，先将小鼠 Ensembl ID 转换为基因 symbol，随后使用 Seurat v5 完成标准单细胞分析，包括质量控制、降维、聚类和细胞类型粗注释。在此基础上，进一步聚焦于 Macrophage/Myeloid 和 Neutrophil 两类关键炎症细胞群，开展目标基因表达分析，并生成适合展示和论文整理的图形结果，如 UMAP 图、FeaturePlot、module score 图、热图以及定制气泡图。

主要功能
	•	读取 GEO 原始单细胞 count 矩阵
	•	将 Ensembl 基因 ID 转换为小鼠基因 symbol
	•	进行质量控制和细胞过滤
	•	执行标准单细胞分析流程：标准化、高变基因筛选、PCA、UMAP 和聚类
	•	提取 cluster marker 并进行粗略细胞类型注释
	•	分析 Kdm6a 与炎症小体及 NETosis 相关基因的关系
	•	聚焦 Macrophage_Myeloid 和 Neutrophil 子集进行下游可视化
	•	同时导出 PNG 和 SVG 图像
	•	保存中间 Seurat 对象和汇总表达矩阵，便于复现

数据来源
	•	GEO 编号： GSE204844
	•	物种： Mus musculus
	•	本流程使用的输入文件：
	•	GSM6198711_Case.counts.tsv
	•	GSM6198712_Control.counts.tsv

分析流程
	1.	读取原始 count 矩阵
	2.	将 Ensembl ID 转换为 gene symbol
	3.	分别构建 Case 和 Control 的 Seurat 对象
	4.	进行质量控制并过滤低质量细胞
	5.	合并过滤后的对象
	6.	运行标准单细胞流程：
	•	NormalizeData()
	•	FindVariableFeatures()
	•	ScaleData()
	•	RunPCA()
	•	FindNeighbors()
	•	FindClusters()
	•	RunUMAP()
	7.	检查目标基因并计算 module score
	8.	使用 Seurat v5 兼容的 JoinLayers() 进行 marker 分析
	9.	对主要细胞群进行粗注释
	10.	提取 Macrophage_Myeloid 和 Neutrophil 子集进行聚焦分析
	11.	输出热图和定制气泡图

主要生物学问题

本流程主要用于探索以下问题：
	•	Kdm6a 主要表达在哪些细胞群中？
	•	Kdm6a 是否与 NLRP3/IL-1β 炎症小体轴相关？
	•	NETosis 样炎症特征是否主要富集于中性粒细胞？
	•	慢性胰腺炎背景下，巨噬/髓系细胞与中性粒细胞的炎症程序有何不同？

软件环境

推荐环境：
	•	R >= 4.3
	•	RStudio
	•	Seurat v5 兼容版本

主要 R 包：
	•	Seurat
	•	SeuratObject
	•	data.table
	•	dplyr
	•	ggplot2
	•	org.Mm.eg.db
	•	AnnotationDbi
	•	pheatmap
	•	tidyr
	•	svglite
	•	grid

输出目录结构

GSE204844_analysis_output/
├── 01_QC/
├── 02_DimReduction/
├── 03_TargetGenes/
├── 04_Annotation/
├── 05_FocusSubset/
├── 06_Heatmap/
├── 07_BubblePlot/
├── obj_case_filtered.rds
├── obj_ctrl_filtered.rds
├── obj_merge_filtered.rds
├── obj_after_umap_cluster.rds
├── obj_after_target_module_score.rds
├── obj_marker_ready.rds
├── obj_after_rough_annotation.rds
├── obj_focus_macrophage_neutrophil.rds
└── obj_focus_joined_for_plot.rds

说明
	•	本流程基于 Seurat v5 layered object 编写。
	•	在 marker 分析及部分下游表达汇总步骤前，使用了 JoinLayers() 以兼容 Seurat v5 的多 layer 数据结构。
	•	当前分析属于探索性分析，因为仅基于 1 个 Case 样本和 1 个 Control 样本。
	•	因此结论应谨慎解释，并建议结合更多样本或实验进一步验证。

引用建议

如果你使用了本仓库，请同时引用：
	•	原始 GEO 数据集 GSE204844
	•	Seurat 方法框架
	•	以及你在研究中实际使用的其他软件包或分析方法

联系方式

如有问题、建议或合作需求，欢迎在本仓库中提交 issue。
