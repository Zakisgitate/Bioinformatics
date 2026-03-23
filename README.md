# Bioinformatics
# GSE204844_KDM6A_scRNA_pipeline

Single-cell RNA-seq analysis pipeline for exploring **KDM6A**, **inflammasome-related signaling**, and **NETosis-associated neutrophil programs** in a mouse chronic pancreatitis dataset from **GEO GSE204844**.

---

## Overview

This repository contains an R/Seurat-based workflow for analyzing the GEO dataset **GSE204844**, with a focus on:

- identifying major cell populations in chronic pancreatitis tissue
- locating **Kdm6a** expression across cell types
- assessing the relationship between **Kdm6a** and **inflammasome-related genes**
  - `Nlrp3`
  - `Pycard`
  - `Casp1`
  - `Il1b`
- evaluating **NETosis / neutrophil-associated inflammatory programs**
  - `Cybb`
  - `S100a8`
  - `S100a9`
  - `Camp`
  - `Ly6g`

The analysis is based on **single-cell RNA-seq count matrices** from:

- `GSM6198711_Case.counts.tsv`
- `GSM6198712_Control.counts.tsv`

---

## Project goals

The main goals of this project are:

1. build a reproducible single-cell analysis workflow for **GSE204844**
2. identify major cell populations in the dataset
3. explore whether **Kdm6a** is preferentially expressed in specific inflammatory cell subsets
4. characterize the relative distribution of:
   - **inflammasome-like macrophage/myeloid signatures**
   - **NETosis-like neutrophil signatures**
5. generate publication-friendly figures and summary matrices

---

## Main workflow

The pipeline includes the following analysis modules:

1. **Raw data import**
   - read GEO `counts.tsv` matrices
   - inspect matrix dimensions and structure

2. **Gene annotation conversion**
   - convert mouse Ensembl IDs to gene symbols
   - collapse duplicated symbols after mapping

3. **Seurat object construction and QC**
   - create Case and Control Seurat objects
   - calculate mitochondrial percentage
   - visualize QC metrics
   - filter low-quality cells

4. **Data integration at the object level**
   - merge filtered Case and Control objects

5. **Standard scRNA-seq workflow**
   - normalization
   - highly variable gene selection
   - scaling
   - PCA
   - neighborhood graph construction
   - clustering
   - UMAP visualization

6. **Target-gene exploration**
   - visualize `Kdm6a`, `Nlrp3`, `Il1b`, and related genes
   - compute inflammasome and NETosis module scores

7. **Marker analysis and rough annotation**
   - identify cluster markers
   - annotate major cell populations

8. **Focused myeloid/neutrophil analysis**
   - subset **Macrophage_Myeloid** and **Neutrophil**
   - compare key inflammatory signatures across these subsets

9. **Summary visualization**
   - heatmaps
   - custom bubble plot
   - exported summary matrices

---

## Biological focus

This analysis is designed to answer the following questions:

- Is **Kdm6a** enriched in macrophage/myeloid cells or neutrophils?
- Are **NLRP3 inflammasome-related genes** preferentially associated with macrophage-like populations?
- Are **NETosis-associated genes** mainly enriched in case-specific neutrophils?
- Does chronic pancreatitis exhibit marked immune cell remodeling at the single-cell level?

---

## Repository structure

A recommended repository structure is:

```text
.
├── README.md
├── GSE204844_KDM6A_scRNA_pipeline.R
├── data/
│   ├── GSM6198711_Case.counts.tsv
│   └── GSM6198712_Control.counts.tsv
├── results/
│   └── GSE204844_analysis_output/
│       ├── 01_QC/
│       ├── 02_DimReduction/
│       ├── 03_TargetGenes/
│       ├── 04_Annotation/
│       ├── 05_FocusSubset/
│       ├── 06_Heatmap/
│       ├── 07_BubblePlot/
│       ├── *.rds
│       └── *.csv
└── figures/
