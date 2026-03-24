# ============================================================
# GSE134520 Severe IM 后续分析脚本 2/2
# 功能：读取注释后的 Seurat 对象，提取上皮细胞，重新聚类，
#      并分析 REG4 / HNF4G 在上皮亚群中的投影与平均表达
# 适用环境：Windows + RStudio
# 默认工作目录：D:/Rdata/GSE134520
# ============================================================

# ------------------------------
# 0. 用户可修改参数
# ------------------------------
project_dir <- "D:/Rdata/GSE134520"

# 如果你知道脚本 1 输出对象的准确路径，可直接写在这里；
# 如果留空（NULL），脚本会自动优先寻找 objects/gse134520_severe_im_annotated.rds
input_annotated_rds <- NULL

results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")
tables_dir <- file.path(project_dir, "tables")
objects_dir <- file.path(project_dir, "objects")
epithelial_figures_dir <- file.path(figures_dir, "epithelial")
epithelial_tables_dir <- file.path(tables_dir, "epithelial")

# Seurat 再分析参数
variable_features_n <- 2000
pcs_to_compute <- 30
dims_use <- 1:20
cluster_resolution <- 0.4

# marker 查找参数
marker_only_positive <- TRUE
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25

# 如果 celltype 列名不是 celltype，可在这里手动指定
celltype_column_override <- NULL

# 用于识别“哪些注释应当被视为上皮细胞”的正则模式
# 如果你的注释名称不同，请优先修改这里。
epithelial_pattern <- "epithelial|goblet|intestinal|enterocyte|stem-like|progenitor"

# 目标 marker
projection_genes <- c(
  "REG4", "HNF4G", "EPCAM", "KRT20", "MUC2", "TFF3",
  "FABP1", "GPA33", "OLFM4", "SOX9", "MKI67"
)

# REG4 / HNF4G 阳性阈值：默认使用 log-normalized data > 0
# 如果你希望更严格，可改成 0.25、0.5 等
positive_cutoff <- 0

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


# 兼容 Seurat v4 / v5 的辅助函数：
# - Seurat v5 推荐使用 layer
# - 较早版本仍使用 slot
get_assay_layer_data <- function(obj, assay = NULL, layer = "data") {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }

  assay_obj <- obj[[assay]]

  if (inherits(assay_obj, "Assay5")) {
    return(LayerData(object = obj, assay = assay, layer = layer))
  }

  GetAssayData(object = obj, assay = assay, slot = layer)
}

average_expression_compatible <- function(obj, features, group_by, layer = "data") {
  avg_formals <- names(formals(AverageExpression))

  if ("layer" %in% avg_formals) {
    return(AverageExpression(
      object = obj,
      assays = DefaultAssay(obj),
      features = features,
      group.by = group_by,
      layer = layer,
      verbose = FALSE
    ))
  }

  AverageExpression(
    object = obj,
    assays = DefaultAssay(obj),
    features = features,
    group.by = group_by,
    slot = layer,
    verbose = FALSE
  )
}

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
  epithelial_figures_dir <- file.path(figures_dir, "epithelial")
  epithelial_tables_dir <- file.path(tables_dir, "epithelial")
}

for (dir_i in c(results_dir, figures_dir, tables_dir, objects_dir,
                epithelial_figures_dir, epithelial_tables_dir)) {
  dir.create(dir_i, showWarnings = FALSE, recursive = TRUE)
}

message("当前工作目录：", getwd())
message("上皮分析图输出目录：", epithelial_figures_dir)
message("上皮分析表输出目录：", epithelial_tables_dir)

# ------------------------------
# 3. 自动寻找注释后的对象
# ------------------------------
find_candidate_rds_files <- function(root_dir) {
  search_dirs <- unique(c(
    root_dir,
    file.path(root_dir, "objects"),
    file.path(root_dir, "results"),
    file.path(root_dir, "results_severe_im")
  ))

  rds_files <- unlist(lapply(search_dirs, function(dd) {
    if (!dir.exists(dd)) return(character(0))
    list.files(dd, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  }), use.names = FALSE)

  unique(normalizePath(rds_files, winslash = "/", mustWork = FALSE))
}

choose_default_annotated_rds <- function(files) {
  if (length(files) == 0) return(NULL)

  file_base <- tolower(basename(files))
  scores <- rep(0, length(files))
  scores <- scores + ifelse(grepl("annotated", file_base), 100, 0)
  scores <- scores + ifelse(grepl("severe_im", file_base), 50, 0)
  scores <- scores + ifelse(grepl("seurat", file_base), 30, 0)

  files[order(scores, decreasing = TRUE)][1]
}

candidate_rds_files <- find_candidate_rds_files(project_dir)
message("\n[Step 1] 检查可用的 .rds 文件 ...")
print(data.frame(index = seq_along(candidate_rds_files), file = candidate_rds_files))

if (is.null(input_annotated_rds)) {
  preferred_file <- file.path(objects_dir, "gse134520_severe_im_annotated.rds")
  if (file.exists(preferred_file)) {
    input_annotated_rds <- preferred_file
  } else {
    input_annotated_rds <- choose_default_annotated_rds(candidate_rds_files)
  }
  message("自动选择的输入对象：", input_annotated_rds)
} else {
  message("使用手动指定的输入对象：", input_annotated_rds)
}

if (is.null(input_annotated_rds) || !file.exists(input_annotated_rds)) {
  stop("未找到可用的注释后 Seurat 对象，请先运行 gse134520_epithelial_annotation.R。")
}

# ------------------------------
# 4. 读取注释后的对象并检查 celltype 列
# ------------------------------
message("\n[Step 2] 读取注释后的 Seurat 对象 ...")
seu <- readRDS(input_annotated_rds)

if (!inherits(seu, "Seurat")) {
  stop("读取到的对象不是 Seurat 对象，请检查 input_annotated_rds。")
}

message("对象读取成功。维度（基因 x 细胞）：", nrow(seu), " x ", ncol(seu))
message("metadata 列名：")
print(colnames(seu@meta.data))

possible_celltype_columns <- c(
  celltype_column_override,
  "celltype",
  "celltype_cluster",
  "celltype_auto",
  "CellType",
  "cell_type"
)
possible_celltype_columns <- possible_celltype_columns[!is.na(possible_celltype_columns) & !is.null(possible_celltype_columns)]
possible_celltype_columns <- unique(possible_celltype_columns)

celltype_col <- possible_celltype_columns[possible_celltype_columns %in% colnames(seu@meta.data)][1]
if (length(celltype_col) == 0 || is.na(celltype_col)) {
  stop(
    "未找到可用的 celltype 注释列。\n",
    "请先运行脚本 1，或手动指定 celltype_column_override。"
  )
}

message("将使用以下列作为 celltype 列：", celltype_col)
message("该列的唯一值：")
print(sort(unique(as.character(seu@meta.data[[celltype_col]]))))

# ------------------------------
# 5. 根据 celltype 提取上皮细胞
# ------------------------------
message("\n[Step 3] 提取上皮细胞 ...")

is_epithelial_cell <- grepl(
  pattern = epithelial_pattern,
  x = tolower(as.character(seu@meta.data[[celltype_col]]))
)

message("识别到的上皮细胞数：", sum(is_epithelial_cell), " / ", ncol(seu))

if (sum(is_epithelial_cell) == 0) {
  stop(
    "没有识别到任何上皮细胞。\n",
    "请检查 celltype 列内容，必要时修改 epithelial_pattern。"
  )
}

epi <- subset(seu, cells = colnames(seu)[is_epithelial_cell])
message("上皮子对象维度（基因 x 细胞）：", nrow(epi), " x ", ncol(epi))

# 保存原始上皮对象（提取后、再聚类前）
saveRDS(epi, file = file.path(objects_dir, "gse134520_epithelial_subset_raw.rds"))

# ------------------------------
# 6. 对上皮细胞重新运行 Seurat 标准流程
# ------------------------------
message("\n[Step 4] 对上皮细胞重新运行 Seurat 标准流程 ...")

DefaultAssay(epi) <- DefaultAssay(seu)

epi <- NormalizeData(epi, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = variable_features_n, verbose = FALSE)
epi <- ScaleData(epi, features = rownames(epi), verbose = FALSE)
epi <- RunPCA(epi, features = VariableFeatures(epi), npcs = pcs_to_compute, verbose = FALSE)
epi <- FindNeighbors(epi, dims = dims_use, verbose = FALSE)
epi <- FindClusters(epi, resolution = cluster_resolution, verbose = FALSE)
epi <- RunUMAP(epi, dims = dims_use, verbose = FALSE)

message("上皮细胞重新聚类完成。")
message("上皮 cluster 唯一值：")
print(sort(unique(as.character(epi$seurat_clusters))))

# ------------------------------
# 7. 输出上皮 UMAP 图
# ------------------------------
message("\n[Step 5] 输出上皮亚群 UMAP 图 ...")

p_epi_umap <- DimPlot(
  object = epi,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = FALSE,
  pt.size = 0.4
) +
  ggtitle("上皮细胞子群 UMAP") +
  theme_bw()

p_epi_umap_label <- DimPlot(
  object = epi,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle("上皮细胞子群 UMAP（带 cluster 编号）") +
  theme_bw()

ggsave(file.path(epithelial_figures_dir, "epithelial_umap_by_cluster.png"), p_epi_umap, width = 8, height = 6, dpi = 300)
ggsave(file.path(epithelial_figures_dir, "epithelial_umap_by_cluster_labeled.png"), p_epi_umap_label, width = 8, height = 6, dpi = 300)
ggsave(file.path(epithelial_figures_dir, "epithelial_umap_by_cluster.pdf"), p_epi_umap, width = 8, height = 6)
ggsave(file.path(epithelial_figures_dir, "epithelial_umap_by_cluster_labeled.pdf"), p_epi_umap_label, width = 8, height = 6)

# ------------------------------
# 8. 计算上皮亚群 marker
# ------------------------------
message("\n[Step 6] 计算上皮亚群 marker genes ...")

Idents(epi) <- "seurat_clusters"

epithelial_markers <- FindAllMarkers(
  object = epi,
  only.pos = marker_only_positive,
  min.pct = marker_min_pct,
  logfc.threshold = marker_logfc_threshold,
  test.use = "wilcox"
)

write.csv(
  epithelial_markers,
  file = file.path(epithelial_tables_dir, "epithelial_markers_all.csv"),
  row.names = FALSE
)

if (nrow(epithelial_markers) > 0) {
  epithelial_top10_markers <- epithelial_markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
} else {
  epithelial_top10_markers <- data.frame()
}

write.csv(
  epithelial_top10_markers,
  file = file.path(epithelial_tables_dir, "epithelial_markers_top10.csv"),
  row.names = FALSE
)

message("已保存上皮 marker 表和 top10 表。")

# ------------------------------
# 9. 输出上皮相关 marker 表达图
# ------------------------------
message("\n[Step 7] 输出上皮 marker 表达图 ...")

genes_present <- intersect(projection_genes, rownames(epi))
genes_missing <- setdiff(projection_genes, rownames(epi))

message("存在于上皮对象中的目标基因：")
print(genes_present)
if (length(genes_missing) > 0) {
  message("未在上皮对象中找到的目标基因：")
  print(genes_missing)
}

write.csv(data.frame(gene_present = genes_present), file = file.path(epithelial_tables_dir, "projection_genes_present.csv"), row.names = FALSE)
write.csv(data.frame(gene_missing = genes_missing), file = file.path(epithelial_tables_dir, "projection_genes_missing.csv"), row.names = FALSE)

if (length(genes_present) > 0) {
  p_dot_genes <- DotPlot(
    object = epi,
    features = genes_present,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("上皮亚群中的目标 marker DotPlot")

  ggsave(
    filename = file.path(epithelial_figures_dir, "epithelial_marker_dotplot.png"),
    plot = p_dot_genes,
    width = 14,
    height = 7,
    dpi = 300
  )

  p_vln_genes <- VlnPlot(
    object = epi,
    features = genes_present,
    group.by = "seurat_clusters",
    pt.size = 0,
    ncol = 3
  )

  ggsave(
    filename = file.path(epithelial_figures_dir, "epithelial_marker_violin.png"),
    plot = p_vln_genes,
    width = 14,
    height = 12,
    dpi = 300
  )

  p_feature_all <- FeaturePlot(
    object = epi,
    features = genes_present,
    reduction = "umap",
    order = TRUE,
    ncol = 3
  )

  ggsave(
    filename = file.path(epithelial_figures_dir, "epithelial_marker_featureplot_all.png"),
    plot = p_feature_all,
    width = 14,
    height = 12,
    dpi = 300
  )
}

# ------------------------------
# 10. REG4 / HNF4G 投影与比较分析
# ------------------------------
message("\n[Step 8] 进行 REG4 / HNF4G 投影分析 ...")

focus_genes <- intersect(c("REG4", "HNF4G"), rownames(epi))
if (length(focus_genes) == 0) {
  stop("上皮对象中既没有 REG4，也没有 HNF4G，无法继续投影分析。")
}

# 分别输出 FeaturePlot
for (gene_i in focus_genes) {
  p_feature_gene <- FeaturePlot(
    object = epi,
    features = gene_i,
    reduction = "umap",
    order = TRUE
  ) +
    ggtitle(paste0(gene_i, " 在上皮 UMAP 上的投影"))

  ggsave(
    filename = file.path(epithelial_figures_dir, paste0("featureplot_", gene_i, ".png")),
    plot = p_feature_gene,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# DotPlot：比较 REG4 / HNF4G 在各上皮亚群中的表达
if (length(focus_genes) > 0) {
  p_dot_focus <- DotPlot(
    object = epi,
    features = focus_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("REG4 / HNF4G 在各上皮亚群中的表达 DotPlot")

  ggsave(
    filename = file.path(epithelial_figures_dir, "reg4_hnf4g_dotplot.png"),
    plot = p_dot_focus,
    width = 7,
    height = 5,
    dpi = 300
  )

  p_vln_focus <- VlnPlot(
    object = epi,
    features = focus_genes,
    group.by = "seurat_clusters",
    pt.size = 0,
    ncol = 2
  )

  ggsave(
    filename = file.path(epithelial_figures_dir, "reg4_hnf4g_violin.png"),
    plot = p_vln_focus,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ------------------------------
# 11. 计算平均表达量并导出 csv
# ------------------------------
message("\n[Step 9] 计算上皮亚群平均表达量 ...")

avg_expr_all <- average_expression_compatible(
  obj = epi,
  features = genes_present,
  group_by = "seurat_clusters",
  layer = "data"
)

avg_expr_matrix <- avg_expr_all[[DefaultAssay(epi)]]
avg_expr_df <- as.data.frame(avg_expr_matrix)
avg_expr_df$gene <- rownames(avg_expr_df)
avg_expr_df <- avg_expr_df[, c("gene", setdiff(colnames(avg_expr_df), "gene")), drop = FALSE]

write.csv(
  avg_expr_df,
  file = file.path(epithelial_tables_dir, "epithelial_cluster_average_expression.csv"),
  row.names = FALSE
)

if (all(c("REG4", "HNF4G") %in% rownames(avg_expr_matrix))) {
  reg4_hnf4g_avg <- avg_expr_matrix[c("REG4", "HNF4G"), , drop = FALSE]
} else {
  reg4_hnf4g_avg <- avg_expr_matrix[intersect(c("REG4", "HNF4G"), rownames(avg_expr_matrix)), , drop = FALSE]
}

reg4_hnf4g_avg_df <- as.data.frame(reg4_hnf4g_avg)
reg4_hnf4g_avg_df$gene <- rownames(reg4_hnf4g_avg_df)
write.csv(
  reg4_hnf4g_avg_df,
  file = file.path(epithelial_tables_dir, "reg4_hnf4g_average_expression_by_cluster.csv"),
  row.names = FALSE
)

# 自动输出 REG4 / HNF4G 最高表达的上皮 cluster
if ("REG4" %in% rownames(avg_expr_matrix)) {
  reg4_values <- avg_expr_matrix["REG4", ]
  reg4_top_cluster <- names(which.max(reg4_values))
  message("REG4 最高表达的上皮 cluster：", reg4_top_cluster)
  message("REG4 在该 cluster 的平均表达量：", signif(max(reg4_values), 4))
}

if ("HNF4G" %in% rownames(avg_expr_matrix)) {
  hnf4g_values <- avg_expr_matrix["HNF4G", ]
  hnf4g_top_cluster <- names(which.max(hnf4g_values))
  message("HNF4G 最高表达的上皮 cluster：", hnf4g_top_cluster)
  message("HNF4G 在该 cluster 的平均表达量：", signif(max(hnf4g_values), 4))
}

# ------------------------------
# 12. 标记 REG4+ / HNF4G+ / Double positive 细胞并绘图
# ------------------------------
message("\n[Step 10] 标记 REG4 positive、HNF4G positive、double positive 细胞 ...")

expr_data <- get_assay_layer_data(epi, layer = "data")

reg4_expr <- if ("REG4" %in% rownames(expr_data)) as.numeric(expr_data["REG4", ]) else rep(0, ncol(epi))
hnf4g_expr <- if ("HNF4G" %in% rownames(expr_data)) as.numeric(expr_data["HNF4G", ]) else rep(0, ncol(epi))

names(reg4_expr) <- colnames(epi)
names(hnf4g_expr) <- colnames(epi)

epi$REG4_positive <- reg4_expr > positive_cutoff
epi$HNF4G_positive <- hnf4g_expr > positive_cutoff

epi$REG4_HNF4G_status <- dplyr::case_when(
  epi$REG4_positive & epi$HNF4G_positive ~ "REG4+/HNF4G+ double positive",
  epi$REG4_positive & !epi$HNF4G_positive ~ "REG4 positive only",
  !epi$REG4_positive & epi$HNF4G_positive ~ "HNF4G positive only",
  TRUE ~ "double negative"
)

epi$REG4_HNF4G_status <- factor(
  epi$REG4_HNF4G_status,
  levels = c(
    "double negative",
    "REG4 positive only",
    "HNF4G positive only",
    "REG4+/HNF4G+ double positive"
  )
)

status_count_df <- as.data.frame(table(epi$REG4_HNF4G_status))
colnames(status_count_df) <- c("group", "n_cells")
write.csv(
  status_count_df,
  file = file.path(epithelial_tables_dir, "reg4_hnf4g_positive_cell_counts.csv"),
  row.names = FALSE
)

if ("umap" %in% names(epi@reductions)) {
  p_status <- DimPlot(
    object = epi,
    reduction = "umap",
    group.by = "REG4_HNF4G_status",
    pt.size = 0.4
  ) +
    ggtitle("REG4 / HNF4G 阳性分组 UMAP") +
    theme_bw()

  ggsave(
    filename = file.path(epithelial_figures_dir, "reg4_hnf4g_status_umap.png"),
    plot = p_status,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# 进一步统计每个上皮 cluster 内不同阳性组别的细胞数
cluster_status_df <- epi@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  count(seurat_clusters, REG4_HNF4G_status, name = "n_cells")

write.csv(
  cluster_status_df,
  file = file.path(epithelial_tables_dir, "reg4_hnf4g_positive_cell_counts_by_cluster.csv"),
  row.names = FALSE
)

# ------------------------------
# 13. 保存对象与解释说明
# ------------------------------
message("\n[Step 11] 保存上皮分析对象 ...")

saveRDS(epi, file = file.path(objects_dir, "gse134520_epithelial_reclustered.rds"))

message("已保存对象：gse134520_epithelial_reclustered.rds")
message("\n解释提示：")
message("1) 如果 REG4 更偏 goblet-like / secretory epithelial 亚群，通常提示更偏分泌样、杯状细胞样或肠型分泌程序。")
message("2) 如果 HNF4G 更偏 intestinal differentiation / enterocyte-like 亚群，通常提示更偏成熟肠分化方向。")
message("3) 如果两者部分重叠，可能代表同一上皮谱系内部的连续过渡状态。")
message("4) 如果两者明显分离，可能代表不同的功能分支或分化路线。")
message("脚本 2 完成。所有结果已输出到 tables/epithelial、figures/epithelial、objects/。")
