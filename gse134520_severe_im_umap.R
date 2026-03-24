# ============================================================
# GSE134520 Severe Intestinal metaplasia scRNA-seq analysis
# Script 2/2: UMAP export and marker gene output
# Suitable for Windows + RStudio
# Recommended working directory: D:/Rdata/GSE134520
# ============================================================

# ------------------------------
# 0. User-adjustable parameters
# ------------------------------
project_dir <- "D:/Rdata/GSE134520"
results_dir <- file.path(project_dir, "results_severe_im")
seurat_rds_file <- file.path(results_dir, "severe_im_seurat_processed.rds")

# Marker-finding defaults.
# You can relax min.pct or logfc.threshold if too few markers are returned.
marker_only_positive <- TRUE
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25

# ------------------------------
# 1. Package installation / loading
# ------------------------------
cran_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "data.table"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in cran_packages) {
  install_if_missing(pkg)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------
# 2. Set working directory
# ------------------------------
if (dir.exists(project_dir)) {
  setwd(project_dir)
} else {
  warning("project_dir does not exist. Using current working directory instead: ", getwd())
  project_dir <- getwd()
  results_dir <- file.path(project_dir, "results_severe_im")
  seurat_rds_file <- file.path(results_dir, "severe_im_seurat_processed.rds")
}

if (!file.exists(seurat_rds_file)) {
  stop(
    "Processed Seurat object not found: ", seurat_rds_file,
    "\nPlease run gse134520_severe_im_rstudio.R first."
  )
}

message("Loading processed Seurat object from: ", seurat_rds_file)
seu <- readRDS(seurat_rds_file)

if (!"umap" %in% names(seu@reductions)) {
  stop("UMAP reduction was not found in the Seurat object. Please rerun script 1.")
}

# ------------------------------
# 3. UMAP plots
# ------------------------------
message("\n[Step 1] Saving UMAP plots ...")

p_umap_cluster <- DimPlot(
  object = seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 0.4
) +
  ggtitle("GSE134520 Severe Intestinal metaplasia - UMAP by cluster") +
  theme_bw()

p_umap_labeled <- DimPlot(
  object = seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle("GSE134520 Severe Intestinal metaplasia - UMAP with cluster labels") +
  theme_bw()

ggsave(
  filename = file.path(results_dir, "umap_by_cluster.png"),
  plot = p_umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(results_dir, "umap_by_cluster.pdf"),
  plot = p_umap_cluster,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(results_dir, "umap_by_cluster_labeled.png"),
  plot = p_umap_labeled,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(results_dir, "umap_by_cluster_labeled.pdf"),
  plot = p_umap_labeled,
  width = 8,
  height = 6
)

message("UMAP plots saved successfully.")

# ------------------------------
# 4. Find marker genes for each cluster
# ------------------------------
message("\n[Step 2] Finding marker genes for each cluster ...")

Idents(seu) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  object = seu,
  only.pos = marker_only_positive,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold,
  test.use = "wilcox"
)

if (nrow(all_markers) == 0) {
  warning(
    "No marker genes were found with the current settings.\n",
    "You may need to lower marker_min_pct or marker_logfc_threshold in this script."
  )
}

write.csv(
  all_markers,
  file = file.path(results_dir, "cluster_markers_all.csv"),
  row.names = FALSE
)

message("Full marker table saved: cluster_markers_all.csv")

# Top 10 markers per cluster, ordered by average log2 fold change.
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10_markers,
  file = file.path(results_dir, "cluster_markers_top10.csv"),
  row.names = FALSE
)

message("Top 10 marker table saved: cluster_markers_top10.csv")

# ------------------------------
# 5. Optional marker heatmap for top markers
# ------------------------------
message("\n[Step 3] Saving optional top-marker heatmap ...")

if (nrow(top10_markers) > 0) {
  top_features <- unique(top10_markers$gene)

  png(
    filename = file.path(results_dir, "cluster_top10_marker_heatmap.png"),
    width = 1800,
    height = 2400,
    res = 200
  )
  print(DoHeatmap(seu, features = top_features, size = 3) + NoLegend())
  dev.off()

  message("Top-marker heatmap saved: cluster_top10_marker_heatmap.png")
} else {
  message("Skipped heatmap because no top markers were available.")
}

# ------------------------------
# 6. Save session outputs and finish
# ------------------------------
cluster_sizes <- data.frame(
  cluster = names(table(Idents(seu))),
  n_cells = as.integer(table(Idents(seu))),
  row.names = NULL,
  stringsAsFactors = FALSE
)
write.csv(cluster_sizes, file = file.path(results_dir, "cluster_cell_counts.csv"), row.names = FALSE)

message("\nAnalysis finished.")
message("Outputs saved in: ", results_dir)
print(cluster_sizes)
