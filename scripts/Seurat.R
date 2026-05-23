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

# --- Extract and save useful information ----
scShred_regressed_tibble <- scShred_regressed@meta.data %>% 
  rownames_to_column("barcode") %>%
  as_tibble() %>%
  mutate(barcode = str_remove(barcode, "-1"))

passing_cells <- scShred_regressed_tibble %>% pull(barcode)
write_tsv(tibble(passing_cells = passing_cells), "./scShred/passing_cells_v2.tsv")

# ---- Annotate with T7 IST metadata - generated in T7_IST_processing.R ----
passing_clones_v2 <- read_tsv("./scShred/passing_clones_v2.tsv.gz")
sc_deletions_filtered_v2 <- read_tsv("./scShred/scShred_v2_deletions_filtered.tsv.gz")

clone_metadata <- passing_clones_v2 %>% 
  mutate(cell = paste0(cell, "-1")) %>%
  column_to_rownames("cell")

deletion_metadata <- sc_deletions_filtered_v2 %>% 
  group_by(cell) %>%
  summarise(n_deletions = n_distinct(deletion_label), deletion = paste(deletion_label, collapse = "_"), n_umi = sum(n_umi), .groups = "drop") %>%
  dplyr::select(cell, deletion, n_umi, n_deletions) %>%
  mutate(cell = paste0(cell, "-1")) %>%
  column_to_rownames("cell")

# ---- Add to Seurat object ----
scShred_regressed <- AddMetaData(scShred_regressed, clone_metadata)
scShred_regressed <- AddMetaData(scShred_regressed, deletion_metadata)

# Extract a data set with the interesting genes and metadata
scShred_goi <- FetchData(scShred_regressed, vars = c("RBM3", "WDR13", "GFP", "Clone", "deletion", "nCount_RNA", "n_umi", "seurat_clusters")) %>%
  rownames_to_column("cell_barcode") %>%
  mutate(cell_barcode = str_remove(cell_barcode, "-1")) %>%
  as_tibble()

write_tsv(scShred_goi, "./scShred/scShred_goi.tsv.gz")

