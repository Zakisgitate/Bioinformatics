# ============================================================
# Script Name : GSE204844_KDM6A_scRNA_pipeline.R
# Project     : Single-cell RNA-seq analysis of chronic pancreatitis
# Dataset     : GEO GSE204844
# Species     : Mus musculus
# Author      : Jiangzhenmin
# Platform    : R / RStudio
#
# Purpose:
#   1. Read raw scRNA-seq count matrices from Case and Control samples
#   2. Convert mouse Ensembl gene IDs to gene symbols
#   3. Perform quality control, filtering, normalization, dimensional reduction,
#      clustering, and rough cell-type annotation
#   4. Focus on Kdm6a and its relationship with inflammasome-related genes
#      (Nlrp3, Pycard, Casp1, Il1b) and NETosis/neutrophil-related genes
#      (e.g., Cybb, S100a8, S100a9, Camp, Ly6g)
#   5. Generate publication-oriented plots and summary matrices
#
# Main Analysis Modules:
#   - Raw data import and matrix checking
#   - Ensembl-to-symbol conversion
#   - Seurat object construction and QC
#   - PCA / UMAP / clustering
#   - Marker analysis and rough annotation
#   - Focused analysis of macrophage/myeloid cells and neutrophils
#   - Heatmap and bubble plot visualization
#
# Environment:
#   - R version      : Recommended >= 4.3
#   - Seurat         : v5 compatible
#   - Required pkgs  : Seurat, SeuratObject, data.table, dplyr, ggplot2,
#                      org.Mm.eg.db, AnnotationDbi, pheatmap, tidyr,
#                      svglite, grid
#
# Input Files:
#   - GSM6198711_Case.counts.tsv
#   - GSM6198712_Control.counts.tsv
#
# Output:
#   - Filtered Seurat objects (.rds)
#   - Marker tables and summary matrices (.csv)
#   - Figures in PNG and SVG format, organized by module
#
# Notes:
#   1. This script is designed for Seurat v5 layered objects
#   2. JoinLayers() is used before marker analysis and selected downstream plots
#   3. Current analysis is exploratory because only one Case and one Control
#      sample are included
# ============================================================


# ============================================================
# GSE204844 | Chronic Pancreatitis scRNA-seq Analysis Pipeline
# Author: Jiangzhenmin
# Platform: R / RStudio
# Notes:
#   1) Compatible with Seurat v5
#   2) Focus on Kdm6a, inflammasome, and NETosis-related genes
#   3) Output figures are saved as both PNG and SVG
# ============================================================


# =========================
# 0. Load packages
# =========================
library(Seurat)
library(SeuratObject)
library(data.table)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(tidyr)
library(svglite)
library(grid)


# =========================
# 1. Global settings
# =========================

# ---------- 1.1 File paths ----------
case_file <- "/Users/jiangzhenmin/Desktop/bioinformatics/GSE204844_RAW/GSM6198711_Case.counts.tsv"
ctrl_file <- "/Users/jiangzhenmin/Desktop/bioinformatics/GSE204844_RAW/GSM6198712_Control.counts.tsv"

# ---------- 1.2 Output directory ----------
outdir <- "GSE204844_analysis_output"
if (!dir.exists(outdir)) dir.create(outdir)

plot_dirs <- c(
  "01_QC",
  "02_DimReduction",
  "03_TargetGenes",
  "04_Annotation",
  "05_FocusSubset",
  "06_Heatmap",
  "07_BubblePlot"
)

for (d in plot_dirs) {
  dir.create(file.path(outdir, d), showWarnings = FALSE, recursive = TRUE)
}

# ---------- 1.3 QC thresholds ----------
MIN_FEATURES <- 200
MAX_FEATURES <- 2000
MAX_MT <- 10

# ---------- 1.4 Dimension reduction parameters ----------
N_HVG <- 2000
USE_DIMS <- 1:20
CLUSTER_RES <- 0.5

# ---------- 1.5 Gene panels ----------
genes_interest <- c(
  "Kdm6a", "Nlrp3", "Il1b", "Pycard", "Casp1",
  "Padi4", "Cybb", "Tlr4", "Cxcl1", "Cxcl2"
)

inflam_genes <- c("Nlrp3", "Pycard", "Casp1", "Il1b")
netosis_genes <- c("Padi4", "Cybb", "Tlr4", "Cxcl1", "Cxcl2")

marker_panel <- c(
  # Endocrine
  "Ins1", "Ins2", "Iapp", "Gcg", "Sst", "Chga", "Ppy", "Pyy",
  # Acinar
  "Prss1", "Prss2", "Cela3b", "Cpa1", "Ctrb1",
  # Ductal
  "Krt19", "Krt8", "Krt18", "Spp1", "Sox9", "Mmp7", "Epcam", "Cldn4", "Pigr",
  # Fibroblast / Stellate
  "Col1a1", "Col1a2", "Col3a1", "Dcn", "Sparc", "Clec3b", "Pdgfra", "Dpt",
  # Endothelial
  "Kdr", "Pecam1", "Cdh5", "Emcn", "Plvap", "Flt1", "Ptprb",
  # Macrophage / Myeloid
  "Lyz2", "Adgre1", "Cd68", "Ctss", "Tyrobp", "Fcgr3", "C1qa", "C1qb", "C1qc",
  "Csf1r", "Aif1", "Mpeg1", "Clec7a", "Lgals3",
  # Neutrophil / NETosis
  "S100a8", "S100a9", "Ly6g", "Cxcr2", "Camp", "Cybb", "Cxcl1", "Cxcl2", "Tlr4", "Padi4",
  # T / NK
  "Cd3d", "Cd3e", "Trbc2", "Nkg7", "Ccl5", "Gzmb", "Gzma", "Il2rb",
  # B
  "Cd79a", "Cd79b", "Ms4a1", "Ighm", "Igkc", "Cd74", "H2-Ab1", "H2-Aa"
)

genes_focus <- c(
  "Kdm6a", "Nlrp3", "Il1b", "Pycard", "Casp1",
  "Padi4", "Cybb", "Tlr4", "Cxcl1", "Cxcl2",
  "S100a8", "S100a9", "Camp", "Ly6g",
  "Clec7a", "Aif1", "C1qb", "C1qc", "Lyz2", "Ctss"
)


# =========================
# 2. Utility functions
# =========================

# ---------- 2.1 Read GEO count matrix ----------
read_counts_tsv <- function(file_path) {
  df <- fread(file_path, data.table = FALSE)
  gene_ids <- df[, 1]
  expr_mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(expr_mat) <- gene_ids
  storage.mode(expr_mat) <- "numeric"
  expr_mat
}

# ---------- 2.2 Convert Ensembl to gene symbol ----------
convert_ensembl_to_symbol <- function(mat) {
  ens_ids <- rownames(mat)
  ens_ids_clean <- sub("\\..*$", "", ens_ids)
  
  symbol_map <- AnnotationDbi::mapIds(
    x = org.Mm.eg.db,
    keys = ens_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  symbol_map <- as.character(symbol_map)
  keep <- !is.na(symbol_map) & symbol_map != ""
  
  mat2 <- mat[keep, , drop = FALSE]
  rownames(mat2) <- symbol_map[keep]
  mat2 <- rowsum(mat2, group = rownames(mat2))
  as.matrix(mat2)
}

# ---------- 2.3 Print matrix summary ----------
print_matrix_summary <- function(mat, name = "Matrix") {
  cat("\n============================\n")
  cat(name, "\n")
  cat("============================\n")
  cat("Dim:", dim(mat)[1], "genes x", dim(mat)[2], "cells\n")
  cat("First 6 genes:\n")
  print(head(rownames(mat)))
  cat("First 6 barcodes:\n")
  print(head(colnames(mat)))
}

# ---------- 2.4 Create Seurat object ----------
create_sc_object <- function(mat, project_name, group_name) {
  obj <- CreateSeuratObject(
    counts = mat,
    project = project_name,
    min.cells = 3,
    min.features = MIN_FEATURES
  )
  obj$group <- group_name
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj
}

# ---------- 2.5 QC summary ----------
print_qc_summary <- function(obj, name = "Object") {
  cat("\n============================\n")
  cat(name, "QC Summary\n")
  cat("============================\n")
  print(summary(obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]))
}

# ---------- 2.6 QC filtering ----------
filter_sc_object <- function(obj) {
  subset(
    obj,
    subset = nFeature_RNA > MIN_FEATURES &
      nFeature_RNA < MAX_FEATURES &
      percent.mt < MAX_MT
  )
}

# ---------- 2.7 Safe gene subset ----------
keep_genes <- function(gene_vec, obj) {
  gene_vec[gene_vec %in% rownames(obj)]
}

# ---------- 2.8 Save ggplot as PNG + SVG ----------
save_plot_dual <- function(plot_obj,
                           module_dir,
                           filename_stub,
                           width = 8,
                           height = 6,
                           dpi = 300,
                           bg = "white") {
  ggsave(
    filename = file.path(outdir, module_dir, paste0(filename_stub, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )
  
  ggsave(
    filename = file.path(outdir, module_dir, paste0(filename_stub, ".svg")),
    plot = plot_obj,
    width = width,
    height = height,
    bg = bg
  )
}

# ---------- 2.9 Save pheatmap as PNG + SVG ----------
save_pheatmap_dual <- function(pheatmap_obj,
                               module_dir,
                               filename_stub,
                               width = 8,
                               height = 6,
                               dpi = 300,
                               bg = "white") {
  png(
    filename = file.path(outdir, module_dir, paste0(filename_stub, ".png")),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = bg
  )
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  dev.off()
  
  svglite::svglite(
    file = file.path(outdir, module_dir, paste0(filename_stub, ".svg")),
    width = width,
    height = height,
    bg = bg
  )
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  dev.off()
}

# ---------- 2.10 Average expression helper ----------
get_avg_expr <- function(obj, genes, group_by_col) {
  avg_res <- AverageExpression(
    obj,
    features = genes,
    group.by = group_by_col,
    assays = "RNA",
    layer = "data"
  )
  avg_res$RNA
}


# =========================
# 3. Read raw matrices
# =========================
mat_case <- read_counts_tsv(case_file)
mat_ctrl <- read_counts_tsv(ctrl_file)

print_matrix_summary(mat_case, "Case raw matrix")
print_matrix_summary(mat_ctrl, "Control raw matrix")


# =========================
# 4. Convert Ensembl ID to Symbol
# =========================
mat_case_sym <- convert_ensembl_to_symbol(mat_case)
mat_ctrl_sym <- convert_ensembl_to_symbol(mat_ctrl)

print_matrix_summary(mat_case_sym, "Case symbol matrix")
print_matrix_summary(mat_ctrl_sym, "Control symbol matrix")

cat("\nCase target genes present:\n")
print(genes_interest[genes_interest %in% rownames(mat_case_sym)])

cat("\nControl target genes present:\n")
print(genes_interest[genes_interest %in% rownames(mat_ctrl_sym)])


# =========================
# 5. Build Seurat objects and QC
# =========================
obj_case <- create_sc_object(mat_case_sym, "CP_Case", "Case")
obj_ctrl <- create_sc_object(mat_ctrl_sym, "CP_Control", "Control")

cat("\nCase object dim:\n")
print(dim(obj_case))
cat("\nControl object dim:\n")
print(dim(obj_ctrl))

print_qc_summary(obj_case, "Case")
print_qc_summary(obj_ctrl, "Control")

p_case_qc <- VlnPlot(
  obj_case,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + ggtitle("Case QC")

p_ctrl_qc <- VlnPlot(
  obj_ctrl,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + ggtitle("Control QC")

print(p_case_qc)
print(p_ctrl_qc)

save_plot_dual(p_case_qc, "01_QC", "01_case_qc", 10, 4)
save_plot_dual(p_ctrl_qc, "01_QC", "02_control_qc", 10, 4)


# =========================
# 6. QC filtering and merge
# =========================
cat("\nBefore filtering:\n")
cat("Case:", ncol(obj_case), "\n")
cat("Control:", ncol(obj_ctrl), "\n")

obj_case_filt <- filter_sc_object(obj_case)
obj_ctrl_filt <- filter_sc_object(obj_ctrl)

cat("\nAfter filtering:\n")
cat("Case:", ncol(obj_case_filt), "\n")
cat("Control:", ncol(obj_ctrl_filt), "\n")

print_qc_summary(obj_case_filt, "Filtered Case")
print_qc_summary(obj_ctrl_filt, "Filtered Control")

obj_merge <- merge(
  x = obj_ctrl_filt,
  y = obj_case_filt,
  add.cell.ids = c("Control", "Case"),
  project = "CP_scRNA"
)

cat("\nMerged object dim:\n")
print(dim(obj_merge))
cat("\nMerged group table:\n")
print(table(obj_merge$group))

saveRDS(obj_case_filt, file.path(outdir, "obj_case_filtered.rds"))
saveRDS(obj_ctrl_filt, file.path(outdir, "obj_ctrl_filtered.rds"))
saveRDS(obj_merge, file.path(outdir, "obj_merge_filtered.rds"))


# =========================
# 7. Standard scRNA-seq workflow
# =========================
obj <- NormalizeData(obj_merge)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = N_HVG
)

cat("\nTop 10 HVGs:\n")
print(head(VariableFeatures(obj), 10))

p_hvg <- VariableFeaturePlot(obj)
p_hvg_labeled <- LabelPoints(
  plot = p_hvg,
  points = head(VariableFeatures(obj), 10),
  repel = TRUE
)
print(p_hvg_labeled)
save_plot_dual(p_hvg_labeled, "02_DimReduction", "03_hvg_plot", 8, 6)

obj <- ScaleData(obj)
obj <- RunPCA(obj)

cat("\nPCA embeddings head:\n")
print(head(Embeddings(obj, "pca")))

cat("\nPCA loading summary:\n")
print(Loadings(obj[["pca"]])[1:20, 1:5])

p_elbow <- ElbowPlot(obj)
print(p_elbow)
save_plot_dual(p_elbow, "02_DimReduction", "04_elbow_plot", 6, 5)

obj <- FindNeighbors(obj, dims = USE_DIMS)
obj <- FindClusters(obj, resolution = CLUSTER_RES)
obj <- RunUMAP(obj, dims = USE_DIMS)

p_umap_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("UMAP by cluster")

p_umap_group <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "group"
) + ggtitle("UMAP by group")

print(p_umap_cluster)
print(p_umap_group)

save_plot_dual(p_umap_cluster, "02_DimReduction", "05_umap_cluster", 8, 6)
save_plot_dual(p_umap_group, "02_DimReduction", "06_umap_group", 8, 6)

cat("\nCluster distribution:\n")
print(table(obj$seurat_clusters))

saveRDS(obj, file.path(outdir, "obj_after_umap_cluster.rds"))


# =========================
# 8. Target gene exploration + module scores
# =========================
genes_present <- keep_genes(genes_interest, obj)

cat("\nTarget genes present in merged object:\n")
print(genes_present)

p_feature_core <- FeaturePlot(
  obj,
  features = keep_genes(c("Kdm6a", "Nlrp3", "Il1b"), obj),
  ncol = 3
)
print(p_feature_core)
save_plot_dual(p_feature_core, "03_TargetGenes", "07_feature_core", 12, 4)

p_feature_inflam <- FeaturePlot(
  obj,
  features = keep_genes(inflam_genes, obj),
  ncol = 2
)
print(p_feature_inflam)
save_plot_dual(p_feature_inflam, "03_TargetGenes", "08_feature_inflam", 8, 6)

p_feature_netosis <- FeaturePlot(
  obj,
  features = keep_genes(netosis_genes, obj),
  ncol = 3
)
print(p_feature_netosis)
save_plot_dual(p_feature_netosis, "03_TargetGenes", "09_feature_netosis", 12, 6)

p_dot_targets <- DotPlot(obj, features = genes_present) +
  RotatedAxis() +
  ggtitle("Target genes across clusters")
print(p_dot_targets)
save_plot_dual(p_dot_targets, "03_TargetGenes", "10_dot_targets_clusters", 10, 6)

p_vln_group_core <- VlnPlot(
  obj,
  features = keep_genes(c("Kdm6a", "Nlrp3", "Il1b"), obj),
  group.by = "group",
  pt.size = 0
)
print(p_vln_group_core)
save_plot_dual(p_vln_group_core, "03_TargetGenes", "11_vln_group_core", 10, 4)

obj <- AddModuleScore(
  obj,
  features = list(keep_genes(inflam_genes, obj)),
  name = "InflamScore"
)

obj <- AddModuleScore(
  obj,
  features = list(keep_genes(netosis_genes, obj)),
  name = "NetosisScore"
)

p_module_feature <- FeaturePlot(
  obj,
  features = c("InflamScore1", "NetosisScore1"),
  ncol = 2
)
print(p_module_feature)
save_plot_dual(p_module_feature, "03_TargetGenes", "12_module_feature", 10, 4)

p_module_vln_group <- VlnPlot(
  obj,
  features = c("InflamScore1", "NetosisScore1"),
  group.by = "group",
  pt.size = 0
)
print(p_module_vln_group)
save_plot_dual(p_module_vln_group, "03_TargetGenes", "13_module_vln_group", 10, 4)

p_module_vln_cluster <- VlnPlot(
  obj,
  features = c("InflamScore1", "NetosisScore1"),
  group.by = "seurat_clusters",
  pt.size = 0
)
print(p_module_vln_cluster)
save_plot_dual(p_module_vln_cluster, "03_TargetGenes", "14_module_vln_cluster", 10, 5)

saveRDS(obj, file.path(outdir, "obj_after_target_module_score.rds"))


# =========================
# 9. Marker analysis (Seurat v5 compatible)
# =========================
Idents(obj) <- "seurat_clusters"
obj_marker <- JoinLayers(obj)

markers_all <- FindAllMarkers(
  obj_marker,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cat("\nMarkers dim:\n")
print(dim(markers_all))
cat("\nMarkers columns:\n")
print(colnames(markers_all))
cat("\nMarkers head:\n")
print(head(markers_all))

write.csv(markers_all, file.path(outdir, "all_cluster_markers.csv"), row.names = FALSE)
saveRDS(obj_marker, file.path(outdir, "obj_marker_ready.rds"))

top10_markers <- markers_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 10) %>%
  dplyr::ungroup()

write.csv(top10_markers, file.path(outdir, "top10_markers_each_cluster.csv"), row.names = FALSE)

split_top10 <- split(top10_markers, top10_markers$cluster)

cat("\nCluster 1 top markers:\n")
print(split_top10[["1"]])
cat("\nCluster 2 top markers:\n")
print(split_top10[["2"]])
cat("\nCluster 6 top markers:\n")
print(split_top10[["6"]])
cat("\nCluster 10 top markers:\n")
print(split_top10[["10"]])
cat("\nCluster 11 top markers:\n")
print(split_top10[["11"]])


# =========================
# 10. Rough annotation
# =========================
obj$cluster_id <- as.character(obj$seurat_clusters)

cluster_annot <- c(
  "0"  = "Endothelial",
  "1"  = "Ductal",
  "2"  = "Fibroblast_Stellate",
  "6"  = "Macrophage_Myeloid",
  "10" = "Acinar",
  "11" = "Neutrophil"
)

all_clusters <- levels(obj$seurat_clusters)

annot_vector <- sapply(all_clusters, function(x) {
  if (x %in% names(cluster_annot)) cluster_annot[[x]] else paste0("Cluster_", x)
})

annot_vector <- unname(annot_vector)
names(annot_vector) <- all_clusters

obj$celltype_rough <- plyr::mapvalues(
  x = as.character(obj$seurat_clusters),
  from = names(annot_vector),
  to = annot_vector
)

obj$celltype_rough <- factor(obj$celltype_rough, levels = annot_vector)

cat("\nRough cell type distribution:\n")
print(table(obj$celltype_rough))

cat("\nRough cell type by group:\n")
print(table(obj$celltype_rough, obj$group))

p_umap_celltype <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "celltype_rough",
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP by rough cell type")

print(p_umap_celltype)
save_plot_dual(p_umap_celltype, "04_Annotation", "15_umap_rough_celltype", 10, 7)

saveRDS(obj, file.path(outdir, "obj_after_rough_annotation.rds"))


# =========================
# 11. Focus analysis: macrophage + neutrophil
# =========================
obj_focus <- subset(
  obj,
  subset = celltype_rough %in% c("Macrophage_Myeloid", "Neutrophil")
)

cat("\nFocus subset dim:\n")
print(dim(obj_focus))
cat("\nFocus subset by cell type:\n")
print(table(obj_focus$celltype_rough))
cat("\nFocus subset by cell type and group:\n")
print(table(obj_focus$celltype_rough, obj_focus$group))

p_focus_group <- DimPlot(
  obj_focus,
  reduction = "umap",
  group.by = "group"
) + ggtitle("Focus subset by group")

p_focus_celltype <- DimPlot(
  obj_focus,
  reduction = "umap",
  group.by = "celltype_rough",
  label = TRUE,
  repel = TRUE
) + ggtitle("Focus subset by cell type")

print(p_focus_group)
print(p_focus_celltype)

save_plot_dual(p_focus_group, "05_FocusSubset", "16_focus_umap_group", 8, 6)
save_plot_dual(p_focus_celltype, "05_FocusSubset", "17_focus_umap_celltype", 8, 6)

genes_focus_present <- keep_genes(genes_focus, obj_focus)

cat("\nFocus genes present:\n")
print(genes_focus_present)

p_dot_focus_celltype <- DotPlot(
  obj_focus,
  features = genes_focus_present,
  group.by = "celltype_rough"
) + RotatedAxis() + ggtitle("Focus genes in macrophage/neutrophil")
print(p_dot_focus_celltype)
save_plot_dual(p_dot_focus_celltype, "05_FocusSubset", "18_dot_focus_celltype", 12, 5)

p_dot_focus_group <- DotPlot(
  obj_focus,
  features = genes_focus_present,
  group.by = "group"
) + RotatedAxis() + ggtitle("Focus genes by group")
print(p_dot_focus_group)
save_plot_dual(p_dot_focus_group, "05_FocusSubset", "19_dot_focus_group", 12, 5)

core_genes <- keep_genes(
  c("Kdm6a", "Nlrp3", "Il1b", "Pycard", "Casp1", "Cybb", "Cxcl1", "Cxcl2"),
  obj_focus
)

p_feature_focus <- FeaturePlot(
  obj_focus,
  features = core_genes,
  ncol = 4
)
print(p_feature_focus)
save_plot_dual(p_feature_focus, "05_FocusSubset", "20_feature_focus", 14, 8)

p_vln_celltype <- VlnPlot(
  obj_focus,
  features = core_genes,
  group.by = "celltype_rough",
  pt.size = 0,
  ncol = 4
)
print(p_vln_celltype)
save_plot_dual(p_vln_celltype, "05_FocusSubset", "21_vln_focus_celltype", 14, 8)

p_vln_group <- VlnPlot(
  obj_focus,
  features = core_genes,
  group.by = "group",
  pt.size = 0,
  ncol = 4
)
print(p_vln_group)
save_plot_dual(p_vln_group, "05_FocusSubset", "22_vln_focus_group", 14, 8)

saveRDS(obj_focus, file.path(outdir, "obj_focus_macrophage_neutrophil.rds"))


# =========================
# 12. Heatmaps for focus subsets
# =========================
genes_heatmap <- keep_genes(genes_focus, obj_focus)

cat("\nGenes used in heatmaps:\n")
print(genes_heatmap)

obj_focus_heat <- JoinLayers(obj_focus)
obj_focus_heat <- NormalizeData(obj_focus_heat)

avg_celltype_mat <- get_avg_expr(
  obj_focus_heat,
  genes = genes_heatmap,
  group_by_col = "celltype_rough"
)

mat_celltype_scaled <- t(scale(t(avg_celltype_mat)))
mat_celltype_scaled[is.na(mat_celltype_scaled)] <- 0

ph_celltype <- pheatmap(
  mat_celltype_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Focus genes heatmap by cell type",
  fontsize_row = 10,
  fontsize_col = 11,
  border_color = NA
)

save_pheatmap_dual(
  ph_celltype,
  "06_Heatmap",
  "24_heatmap_by_celltype",
  width = 8,
  height = 6
)

write.csv(avg_celltype_mat, file.path(outdir, "focus_avg_expression_by_celltype.csv"))
write.csv(mat_celltype_scaled, file.path(outdir, "focus_heatmap_matrix_by_celltype_scaled.csv"))

obj_focus_heat$group_celltype <- paste(
  obj_focus_heat$group,
  obj_focus_heat$celltype_rough,
  sep = "-"
)

avg_group_celltype_mat <- get_avg_expr(
  obj_focus_heat,
  genes = genes_heatmap,
  group_by_col = "group_celltype"
)

mat_group_celltype_scaled <- t(scale(t(avg_group_celltype_mat)))
mat_group_celltype_scaled[is.na(mat_group_celltype_scaled)] <- 0

ph_group_celltype <- pheatmap(
  mat_group_celltype_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Focus genes heatmap by group and cell type",
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

save_pheatmap_dual(
  ph_group_celltype,
  "06_Heatmap",
  "25_heatmap_by_group_celltype",
  width = 8,
  height = 6
)

write.csv(avg_group_celltype_mat, file.path(outdir, "focus_avg_expression_by_group_celltype.csv"))
write.csv(mat_group_celltype_scaled, file.path(outdir, "focus_heatmap_matrix_by_group_celltype_scaled.csv"))


# =========================
# 13. Bubble plot for key myeloid subsets
# =========================
obj_focus2 <- JoinLayers(obj_focus)
obj_focus2 <- NormalizeData(obj_focus2)

genes_plot <- keep_genes(genes_focus, obj_focus2)

cat("genes_plot:\n")
print(genes_plot)
cat("genes_plot length =", length(genes_plot), "\n")

if (length(genes_plot) == 0) {
  stop("genes_plot 为空，请检查 obj_focus2 中是否保留了目标基因。")
}

obj_focus2$plot_group <- NA_character_

obj_focus2$plot_group[
  obj_focus2$celltype_rough == "Macrophage_Myeloid" & obj_focus2$group == "Case"
] <- "Case_Macrophage"

obj_focus2$plot_group[
  obj_focus2$celltype_rough == "Macrophage_Myeloid" & obj_focus2$group == "Control"
] <- "Control_Macrophage"

obj_focus2$plot_group[
  obj_focus2$celltype_rough == "Neutrophil" & obj_focus2$group == "Case"
] <- "Case_Neutrophil"

obj_focus2$plot_group <- factor(
  obj_focus2$plot_group,
  levels = c("Case_Macrophage", "Control_Macrophage", "Case_Neutrophil")
)

cat("\nplot_group 分布：\n")
print(table(obj_focus2$plot_group, useNA = "ifany"))

expr_mat <- GetAssayData(
  obj_focus2,
  assay = "RNA",
  layer = "data"
)[genes_plot, , drop = FALSE]

cat("\nexpr_mat dim:\n")
print(dim(expr_mat))

meta_df <- obj_focus2@meta.data[, c("plot_group"), drop = FALSE]
meta_df$cell_id <- rownames(meta_df)

plot_df <- data.frame()

for (g in genes_plot) {
  tmp <- data.frame(
    cell_id = colnames(expr_mat),
    gene = rep(g, ncol(expr_mat)),
    expr = as.numeric(expr_mat[g, ]),
    stringsAsFactors = FALSE
  )
  
  tmp <- merge(tmp, meta_df, by = "cell_id", all.x = TRUE)
  tmp <- tmp[!is.na(tmp$plot_group), , drop = FALSE]
  
  tmp_sum <- tmp %>%
    dplyr::group_by(plot_group, gene) %>%
    dplyr::summarise(
      avg_expr = mean(expr),
      pct_expr = mean(expr > 0) * 100,
      .groups = "drop"
    )
  
  plot_df <- rbind(plot_df, as.data.frame(tmp_sum))
}

cat("\nplot_df columns:\n")
print(colnames(plot_df))
cat("\nplot_df head:\n")
print(head(plot_df))

plot_df <- plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(avg_expr_scaled = as.numeric(scale(avg_expr))) %>%
  dplyr::ungroup()

plot_df$avg_expr_scaled[is.na(plot_df$avg_expr_scaled)] <- 0
plot_df$gene <- factor(plot_df$gene, levels = rev(genes_plot))
plot_df$plot_group <- factor(
  plot_df$plot_group,
  levels = c("Case_Macrophage", "Control_Macrophage", "Case_Neutrophil")
)

p_bubble <- ggplot(plot_df, aes(x = plot_group, y = gene)) +
  geom_point(aes(size = pct_expr, color = avg_expr_scaled)) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "#F3F4F6",
    high = "#EF4444",
    midpoint = 0
  ) +
  scale_size(range = c(1.5, 10)) +
  labs(
    title = "Inflammasome- and NETosis-related genes in key myeloid subsets",
    x = NULL,
    y = NULL,
    color = "Scaled\nexpression",
    size = "% expressed"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_bubble)

save_plot_dual(
  p_bubble,
  "07_BubblePlot",
  "26_bubble_plot_myeloid_neutrophil",
  width = 12,
  height = 8
)

write.csv(
  plot_df,
  file.path(outdir, "07_BubblePlot", "bubble_plot_data_myeloid_neutrophil.csv"),
  row.names = FALSE
)

saveRDS(
  obj_focus2,
  file.path(outdir, "obj_focus_joined_for_plot.rds")
)


# =========================
# 14. Output structure summary
# =========================
cat("\n============================\n")
cat("Analysis outputs saved to:\n")
cat(normalizePath(outdir), "\n")
cat("============================\n")

cat("\nSubfolders:\n")
print(list.dirs(outdir, recursive = FALSE))