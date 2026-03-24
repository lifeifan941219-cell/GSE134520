# ============================================================
# GSE134520 Severe IM 后续分析脚本 1/2
# 功能：读取既有 Seurat 对象，进行 cluster 初步细胞类型注释，
#      重点识别上皮细胞群，并保存注释后的对象
# 适用环境：Windows + RStudio
# 默认工作目录：D:/Rdata/GSE134520
# ============================================================

# ------------------------------
# 0. 用户可修改参数
# ------------------------------
project_dir <- "D:/Rdata/GSE134520"

# 如果你已经明确知道 Seurat 对象文件名，可以直接在这里填写完整路径；
# 如果留空（NULL），脚本会自动搜索工作目录及常见子目录中的 .rds 文件。
input_seurat_file <- NULL

# 输出目录
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
tables_dir <- file.path(project_dir, "tables")
objects_dir <- file.path(project_dir, "objects")
annotation_figures_dir <- file.path(figures_dir, "annotation")
annotation_tables_dir <- file.path(tables_dir, "annotation")

# marker 查找参数
marker_only_positive <- TRUE
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25

# 如果你的对象中 cluster 列名不是 seurat_clusters，可以在这里手动指定；
# 留空则脚本自动判断。
cluster_column_override <- NULL

# 如果你想手动修改 cluster 注释名称，可在脚本运行第一次后，
# 查看输出的 cluster_annotation_summary.csv，再回到这里修改映射。
# 名称格式：c("0" = "Epithelial", "1" = "T cells")
manual_cluster_annotation <- c()

# ------------------------------
# 1. 安装并加载常用 R 包
# ------------------------------
cran_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork",
  "Matrix",
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
  library(patchwork)
  library(Matrix)
  library(data.table)
})

# ------------------------------
# 2. 设置工作目录与输出目录
# ------------------------------
if (dir.exists(project_dir)) {
  setwd(project_dir)
} else {
  warning("project_dir 不存在，改用当前工作目录：", getwd())
  project_dir <- getwd()
  results_dir <- file.path(project_dir, "results")
  figures_dir <- file.path(project_dir, "figures")
  tables_dir <- file.path(project_dir, "tables")
  objects_dir <- file.path(project_dir, "objects")
  annotation_figures_dir <- file.path(figures_dir, "annotation")
  annotation_tables_dir <- file.path(tables_dir, "annotation")
}

for (dir_i in c(results_dir, figures_dir, tables_dir, objects_dir,
                annotation_figures_dir, annotation_tables_dir)) {
  dir.create(dir_i, showWarnings = FALSE, recursive = TRUE)
}

message("当前工作目录：", getwd())
message("注释图输出目录：", annotation_figures_dir)
message("注释表输出目录：", annotation_tables_dir)

# ------------------------------
# 3. 工具函数：自动寻找 Seurat 对象 .rds 文件
# ------------------------------
find_candidate_rds_files <- function(root_dir) {
  search_dirs <- unique(c(
    root_dir,
    file.path(root_dir, "results"),
    file.path(root_dir, "results_severe_im"),
    file.path(root_dir, "objects"),
    file.path(root_dir, "output"),
    file.path(root_dir, "outputs")
  ))

  rds_files <- unlist(lapply(search_dirs, function(dd) {
    if (!dir.exists(dd)) return(character(0))
    list.files(dd, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  }), use.names = FALSE)

  unique(normalizePath(rds_files, winslash = "/", mustWork = FALSE))
}

choose_default_rds <- function(files) {
  if (length(files) == 0) return(NULL)

  score_pattern <- c(
    "annotated" = 100,
    "severe_im_seurat_processed" = 90,
    "seurat" = 60,
    "processed" = 50,
    "umap" = 20
  )

  scores <- rep(0, length(files))
  file_base <- tolower(basename(files))

  for (nm in names(score_pattern)) {
    scores <- scores + ifelse(grepl(nm, file_base, fixed = TRUE), score_pattern[[nm]], 0)
  }

  files[order(scores, decreasing = TRUE)][1]
}

candidate_rds_files <- find_candidate_rds_files(project_dir)

message("\n[Step 1] 检查工作目录中的 .rds 文件 ...")
if (length(candidate_rds_files) == 0) {
  stop(
    "未在工作目录及常见子目录中找到任何 .rds 文件。\n",
    "请先把之前保存的 Seurat 对象放到 D:/Rdata/GSE134520，或手动设置 input_seurat_file。"
  )
}

print(data.frame(index = seq_along(candidate_rds_files), file = candidate_rds_files))

if (is.null(input_seurat_file)) {
  input_seurat_file <- choose_default_rds(candidate_rds_files)
  message("自动选择的 Seurat 对象文件：", input_seurat_file)
} else {
  message("使用手动指定的 Seurat 对象文件：", input_seurat_file)
}

if (!file.exists(input_seurat_file)) {
  stop("指定的 Seurat 对象文件不存在：", input_seurat_file)
}

# ------------------------------
# 4. 读取对象并检查基本信息
# ------------------------------
message("\n[Step 2] 读取 Seurat 对象 ...")
seu <- readRDS(input_seurat_file)

if (!inherits(seu, "Seurat")) {
  stop("选中的 .rds 文件不是 Seurat 对象，请手动修改 input_seurat_file。")
}

message("Seurat 对象读取成功。")
message("对象维度（基因 x 细胞）：", nrow(seu), " x ", ncol(seu))
message("Assays：", paste(names(seu@assays), collapse = ", "))
message("Reductions：", paste(names(seu@reductions), collapse = ", "))
message("metadata 列名：")
print(colnames(seu@meta.data))

# ------------------------------
# 5. 检查 cluster 列名与 UMAP 是否存在
# ------------------------------
message("\n[Step 3] 检查 cluster 信息与已有 UMAP ...")

possible_cluster_columns <- c(
  cluster_column_override,
  "seurat_clusters",
  "cluster",
  "Cluster",
  "clusters",
  "RNA_snn_res.0.5"
)
possible_cluster_columns <- possible_cluster_columns[!is.na(possible_cluster_columns) & !is.null(possible_cluster_columns)]
possible_cluster_columns <- unique(possible_cluster_columns)

cluster_col <- possible_cluster_columns[possible_cluster_columns %in% colnames(seu@meta.data)][1]

# 如果 metadata 中没有显式 cluster 列，则尝试从 active identities 获取。
if (length(cluster_col) == 0 || is.na(cluster_col)) {
  current_ident <- Idents(seu)
  if (length(current_ident) == ncol(seu)) {
    seu$active_ident_cluster <- as.character(current_ident)
    cluster_col <- "active_ident_cluster"
    message("metadata 中未发现常见 cluster 列，已从 Idents(seu) 写入 active_ident_cluster。")
  } else {
    stop("未能识别 cluster 列。请检查对象是否已经完成聚类。")
  }
}

message("将使用以下列作为 cluster 列：", cluster_col)
message("cluster 唯一值：")
print(sort(unique(as.character(seu@meta.data[[cluster_col]]))))

# 为了后续统一，复制到 seurat_clusters_use
seu$seurat_clusters_use <- as.character(seu@meta.data[[cluster_col]])
Idents(seu) <- "seurat_clusters_use"

if (!"umap" %in% names(seu@reductions)) {
  warning("对象中未发现已有 UMAP。第一个脚本仍可继续注释，但注释前后 UMAP 图无法直接输出。")
} else {
  message("检测到已有 UMAP，可直接用于注释前后展示。")
}

# ------------------------------
# 6. 定义常见 marker panel
# ------------------------------
message("\n[Step 4] 准备 marker panel ...")

marker_panel <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"),
  Goblet_intestinal_like_epithelial = c("MUC2", "TFF3", "REG4", "GPA33", "SI", "FABP1", "SLC26A3"),
  Stem_like_progenitor_epithelial = c("OLFM4", "SOX9", "MKI67", "TOP2A"),
  T_cells = c("CD3D", "CD3E", "IL7R"),
  NK_cells = c("NKG7", "GNLY", "KLRD1"),
  B_cells = c("CD79A", "MS4A1", "CD74"),
  Plasma_cells = c("MZB1", "JCHAIN", "IGKC"),
  Myeloid_macrophage = c("LYZ", "FCN1", "C1QC", "APOE", "CD68"),
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM"),
  Endothelial = c("PECAM1", "EMCN", "VWF", "KDR"),
  Neuroendocrine = c("CHGA", "CHGB", "NEUROD1")
)

all_panel_genes <- unique(unlist(marker_panel))
genes_present <- intersect(all_panel_genes, rownames(seu))
genes_missing <- setdiff(all_panel_genes, rownames(seu))

message("marker panel 中存在于对象内的基因数：", length(genes_present))
message("marker panel 中缺失的基因数：", length(genes_missing))
if (length(genes_missing) > 0) {
  message("以下 marker 在对象中未找到（前 50 个）：")
  print(utils::head(genes_missing, 50))
}

write.csv(
  data.frame(gene_present = genes_present),
  file = file.path(annotation_tables_dir, "marker_panel_genes_present.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_missing = genes_missing),
  file = file.path(annotation_tables_dir, "marker_panel_genes_missing.csv"),
  row.names = FALSE
)

# ------------------------------
# 7. 计算 cluster marker
# ------------------------------
message("\n[Step 5] 计算 cluster marker genes ...")

all_markers <- FindAllMarkers(
  object = seu,
  only.pos = marker_only_positive,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold,
  test.use = "wilcox"
)

write.csv(
  all_markers,
  file = file.path(annotation_tables_dir, "cluster_markers_all.csv"),
  row.names = FALSE
)

message("已保存所有 cluster marker：cluster_markers_all.csv")

if (nrow(all_markers) > 0) {
  top10_markers <- all_markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
} else {
  top10_markers <- data.frame()
}

write.csv(
  top10_markers,
  file = file.path(annotation_tables_dir, "cluster_markers_top10.csv"),
  row.names = FALSE
)

message("已保存 top10 marker：cluster_markers_top10.csv")

# ------------------------------
# 8. 计算 marker panel 模块得分并给出自动注释建议
# ------------------------------
message("\n[Step 6] 基于 marker panel 进行自动注释建议 ...")

score_names <- names(marker_panel)
for (ct in names(marker_panel)) {
  features_use <- intersect(marker_panel[[ct]], rownames(seu))
  if (length(features_use) == 0) {
    warning("细胞类型面板没有一个基因出现在对象中：", ct)
    next
  }
  seu <- AddModuleScore(
    object = seu,
    features = list(features_use),
    name = paste0(ct, "_Score")
  )
}

# AddModuleScore 会生成 *_Score1，因此这里统一收集真正的得分列。
score_columns <- grep("_Score1$", colnames(seu@meta.data), value = TRUE)
message("实际生成的模块打分列：")
print(score_columns)

if (length(score_columns) == 0) {
  stop("未生成任何模块打分列，请检查 marker panel。")
}

# 为每个细胞找出得分最高的类型。
score_matrix <- as.matrix(seu@meta.data[, score_columns, drop = FALSE])
max_idx <- apply(score_matrix, 1, which.max)
celltype_auto <- score_columns[max_idx]
celltype_auto <- sub("_Score1$", "", celltype_auto)
celltype_auto <- gsub("_", " ", celltype_auto)
seu$celltype_auto <- celltype_auto

# cluster 层面：统计每个 cluster 各类模块得分平均值。
cluster_score_df <- seu@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, seurat_clusters_use, all_of(score_columns)) %>%
  group_by(seurat_clusters_use) %>%
  summarise(across(all_of(score_columns), mean), .groups = "drop")

# 根据 cluster 的最高平均得分给出 cluster 级别建议注释。
cluster_score_matrix <- as.matrix(cluster_score_df[, score_columns, drop = FALSE])
cluster_max_idx <- apply(cluster_score_matrix, 1, which.max)
cluster_pred <- score_columns[cluster_max_idx]
cluster_pred <- sub("_Score1$", "", cluster_pred)
cluster_pred <- gsub("_", " ", cluster_pred)
cluster_score_df$predicted_celltype <- cluster_pred
cluster_score_df$predicted_score_column <- score_columns[cluster_max_idx]

# 叠加 top markers，帮助人工判断。
cluster_top_marker_summary <- if (nrow(top10_markers) > 0) {
  top10_markers %>%
    group_by(cluster) %>%
    summarise(top10_genes = paste(gene, collapse = ", "), .groups = "drop")
} else {
  data.frame(cluster = character(0), top10_genes = character(0))
}

cluster_annotation_summary <- cluster_score_df %>%
  rename(cluster = seurat_clusters_use) %>%
  left_join(cluster_top_marker_summary, by = "cluster")

# 如果用户手动提供了 cluster 注释映射，则优先使用手动结果；否则使用自动建议。
cluster_annotation_summary$celltype_final <- cluster_annotation_summary$predicted_celltype
if (length(manual_cluster_annotation) > 0) {
  manual_map_df <- data.frame(
    cluster = names(manual_cluster_annotation),
    manual_celltype = as.character(manual_cluster_annotation),
    stringsAsFactors = FALSE
  )
  cluster_annotation_summary <- cluster_annotation_summary %>%
    left_join(manual_map_df, by = "cluster") %>%
    mutate(celltype_final = ifelse(!is.na(manual_celltype), manual_celltype, celltype_final))
}

write.csv(
  cluster_score_df,
  file = file.path(annotation_tables_dir, "cluster_marker_panel_scores.csv"),
  row.names = FALSE
)
write.csv(
  cluster_annotation_summary,
  file = file.path(annotation_tables_dir, "cluster_annotation_summary.csv"),
  row.names = FALSE
)

message("cluster 注释总结表已保存：cluster_annotation_summary.csv")
print(cluster_annotation_summary[, c("cluster", "predicted_celltype", "celltype_final")])

# 将 cluster 级别注释写回每个细胞 metadata。
cluster_to_celltype <- setNames(cluster_annotation_summary$celltype_final, cluster_annotation_summary$cluster)
seu$celltype_cluster <- unname(cluster_to_celltype[as.character(seu$seurat_clusters_use)])
seu$celltype <- seu$celltype_cluster

# 保存每个细胞的主导模块得分，便于后续排查。
seu$celltype_score <- apply(score_matrix, 1, max)

message("metadata 中已新增列：celltype_auto, celltype_cluster, celltype, celltype_score")
message("celltype 唯一值：")
print(sort(unique(as.character(seu$celltype))))

# ------------------------------
# 9. 输出注释前后的 UMAP 与 marker panel 图
# ------------------------------
message("\n[Step 7] 输出注释前后的 UMAP 与 marker panel 图 ...")

if ("umap" %in% names(seu@reductions)) {
  p_before <- DimPlot(
    object = seu,
    reduction = "umap",
    group.by = "seurat_clusters_use",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.4
  ) +
    ggtitle("注释前：按 cluster 显示") +
    theme_bw()

  p_after <- DimPlot(
    object = seu,
    reduction = "umap",
    group.by = "celltype",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.4
  ) +
    ggtitle("注释后：按 celltype 显示") +
    theme_bw()

  ggsave(
    filename = file.path(annotation_figures_dir, "umap_before_annotation.png"),
    plot = p_before,
    width = 8,
    height = 6,
    dpi = 300
  )
  ggsave(
    filename = file.path(annotation_figures_dir, "umap_after_annotation.png"),
    plot = p_after,
    width = 9,
    height = 6,
    dpi = 300
  )
  ggsave(
    filename = file.path(annotation_figures_dir, "umap_before_after_annotation_combined.png"),
    plot = p_before + p_after,
    width = 16,
    height = 6,
    dpi = 300
  )
}

# DotPlot：适合快速比较各 cluster 对 marker panel 的表达。
marker_panel_for_plot <- unique(unlist(marker_panel))
marker_panel_for_plot <- intersect(marker_panel_for_plot, rownames(seu))

if (length(marker_panel_for_plot) > 0) {
  p_dot <- DotPlot(
    object = seu,
    features = marker_panel_for_plot,
    group.by = "seurat_clusters_use"
  ) +
    RotatedAxis() +
    ggtitle("各 cluster 的常见 marker 表达 DotPlot")

  ggsave(
    filename = file.path(annotation_figures_dir, "cluster_marker_panel_dotplot.png"),
    plot = p_dot,
    width = 16,
    height = 8,
    dpi = 300
  )
}

# 单独投影上皮相关 marker，帮助识别 epithelial clusters。
epithelial_markers_to_plot <- intersect(
  c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "MUC2", "TFF3", "REG4", "GPA33", "OLFM4", "SOX9", "MKI67"),
  rownames(seu)
)

if ("umap" %in% names(seu@reductions) && length(epithelial_markers_to_plot) > 0) {
  p_feature <- FeaturePlot(
    object = seu,
    features = epithelial_markers_to_plot,
    reduction = "umap",
    ncol = 4,
    order = TRUE
  )

  ggsave(
    filename = file.path(annotation_figures_dir, "epithelial_related_featureplot.png"),
    plot = p_feature,
    width = 16,
    height = 12,
    dpi = 300
  )
}

# ------------------------------
# 10. 标记哪些 cluster 更可能是上皮细胞
# ------------------------------
message("\n[Step 8] 标记可能的上皮 cluster ...")

# 这里把包含以下关键词的注释视为上皮相关。
# 这样可以同时保留：
# - Epithelial
# - Goblet / intestinal-like epithelial
# - Stem-like / progenitor epithelial
cluster_annotation_summary$is_epithelial_related <- grepl(
  pattern = "epithelial|goblet|intestinal|enterocyte|progenitor|stem-like",
  x = tolower(cluster_annotation_summary$celltype_final)
)

write.csv(
  cluster_annotation_summary,
  file = file.path(annotation_tables_dir, "cluster_annotation_summary.csv"),
  row.names = FALSE
)

message("以下 cluster 被标记为上皮相关：")
print(cluster_annotation_summary[cluster_annotation_summary$is_epithelial_related, c("cluster", "celltype_final")])

# ------------------------------
# 11. 保存注释后的对象
# ------------------------------
message("\n[Step 9] 保存注释后的 Seurat 对象 ...")

annotated_rds_file <- file.path(objects_dir, "gse134520_severe_im_annotated.rds")
saveRDS(seu, file = annotated_rds_file)

message("注释后的 Seurat 对象已保存：", annotated_rds_file)
message("脚本 1 完成。接下来请运行：gse134520_epithelial_reg4_hnf4g_projection.R")
