library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)

setwd("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/submission/git")

# ---- Loading cell ranger objects ----
run_dir <- "/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/data/scShred/cellranger/tx_Shred_v2"
data_dir <- file.path(run_dir, "outs", "filtered_feature_bc_matrix")

# ---- Load counts ----
counts <- Read10X(data.dir = data_dir)

# ---- Create Seurat scShred_cellsect ----
scShred <- CreateSeuratObject(counts = counts, project = basename(run_dir), min.cells = 3, min.features = 200) 

# ---- Basic QC metrics and plots ----
scShred[["percent.mt"]] <- PercentageFeatureSet(scShred, pattern = "^MT-")

# ---- Simple filtering ----
# for transcripts
scShred_cells <- subset(
  scShred,
  subset = nFeature_RNA >= 400 & nFeature_RNA <= 3000 &
    nCount_RNA >= 500  & nCount_RNA <= 4000 &
    percent.mt <= 5
)

# ---- Standard Seurat workflow ----
scShred_cells <- NormalizeData(scShred_cells)
scShred_cells <- FindVariableFeatures(scShred_cells, selection.method = "vst", nfeatures = 2000)
scShred_cells <- ScaleData(scShred_cells, features = VariableFeatures(scShred_cells))
scShred_cells <- RunPCA(scShred_cells, features = VariableFeatures(scShred_cells))
scShred_cells <- CellCycleScoring(scShred_cells, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)

# Plot the PCs:
p_ncount <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_cells, "nCount_RNA", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
  p}), ncol = 3)

p_s     <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_cells, "S.Score", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
  p}),  ncol = 3)
p_g2m   <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_cells, "G2M.Score", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
  p}), ncol = 3)

pdf("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/data/scShred/plots/scShred_PC_correlations_v2.pdf", width = 11, height = 8.5)
print(p_ncount + plot_annotation(title = "nCount_RNA vs PCs"))
print(p_s      + plot_annotation(title = "PCs vs S.Score"))
print(p_g2m    + plot_annotation(title = "PCs vs G2M.Score"))
dev.off()

# Regress out PC1 (correlates with read number)
scShred_regressed <- ScaleData(scShred_cells, vars.to.regress=c("nCount_RNA"))
scShred_regressed <- RunPCA(scShred_regressed)

# Plot the PCs after regressing nCount:
p_ncount <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_regressed, "nCount_RNA", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
}), ncol = 3)

p_s <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_regressed, "S.Score", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
}), ncol = 3)

p_g2m <- wrap_plots(lapply(pcs, function(pc) {
  p <- FeatureScatter(scShred_regressed, "G2M.Score", pc)
  p$layers[[1]]$aes_params$size <- 0.1
  p + NoLegend()
}), ncol = 3)

pdf("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/data/scShred/plots/scShred_PC_correlations_regressed_v2.pdf", width = 11, height = 8.5)
print(p_ncount + plot_annotation(title = "nCount_RNA vs PCs"))
print(p_s      + plot_annotation(title = "PCs vs S.Score"))
print(p_g2m    + plot_annotation(title = "PCs vs G2M.Score"))
dev.off()

# Choose dims after looking at elbow plot
ElbowPlot(scShred_regressed)

# UMAP
dims_use <- 1:10
scShred_regressed <- FindNeighbors(scShred_regressed, dims = dims_use)
scShred_regressed <- FindClusters(scShred_regressed, resolution = 0.5)
scShred_regressed <- RunUMAP(scShred_regressed, dims = dims_use)

# ---- Save  ----
saveRDS(scShred_regressed, file = file.path(run_dir, "seurat_qc_umap.rds"))