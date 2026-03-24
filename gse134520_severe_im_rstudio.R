# ============================================================
# GSE134520 Severe Intestinal metaplasia scRNA-seq analysis
# Script 1/2: data import, Severe IM filtering, QC, Seurat workflow
# Suitable for Windows + RStudio
# Recommended working directory: D:/Rdata/GSE134520
# ============================================================

# ------------------------------
# 0. User-adjustable parameters
# ------------------------------
project_dir <- "D:/Rdata/GSE134520"
raw_dir <- file.path(project_dir, "GSE134520_RAW")
series_matrix_file <- file.path(project_dir, "GSE134520_series_matrix.txt.gz")
results_dir <- file.path(project_dir, "results_severe_im")

# Quality-control defaults.
# These are reasonable starting values and can be adjusted if needed.
min_features <- 400
max_features <- 7000
max_percent_mt <- 20

# Seurat analysis parameters.
variable_features_n <- 2000
pcs_to_compute <- 30
dims_use <- 1:20
cluster_resolution <- 0.5

# ------------------------------
# 1. Package installation / loading
# ------------------------------
cran_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "Matrix",
  "data.table",
  "stringr",
  "patchwork"
)

bioc_packages <- c("GEOquery")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in cran_packages) {
  install_if_missing(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(GEOquery)
})

# ------------------------------
# 2. Set working directory and create output folder
# ------------------------------
if (dir.exists(project_dir)) {
  setwd(project_dir)
} else {
  warning("project_dir does not exist. Using current working directory instead: ", getwd())
  project_dir <- getwd()
  raw_dir <- file.path(project_dir, "GSE134520_RAW")
  series_matrix_file <- file.path(project_dir, "GSE134520_series_matrix.txt.gz")
  results_dir <- file.path(project_dir, "results_severe_im")
}

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

message("Current working directory: ", getwd())
message("Results will be saved in: ", results_dir)

# ------------------------------
# 3. File existence checks
# ------------------------------
message("\n[Step 1] Checking files ...")
message("series_matrix_file: ", series_matrix_file)
message("raw_dir: ", raw_dir)

if (!file.exists(series_matrix_file)) {
  warning(
    "Series matrix file was not found: ", series_matrix_file,
    "\nYou can still proceed if you manually place the file here or edit series_matrix_file."
  )
}

if (!dir.exists(raw_dir)) {
  stop(
    "raw_dir does not exist: ", raw_dir,
    "\nPlease extract GSE134520_RAW.tar into this folder or modify raw_dir."
  )
}

all_expression_files <- list.files(
  raw_dir,
  pattern = "^GSM\\d+.*processed_.*\\.txt\\.gz$",
  full.names = TRUE
)

if (length(all_expression_files) == 0) {
  stop(
    "No processed expression files were found in: ", raw_dir,
    "\nPlease check whether GSE134520_RAW.tar has been extracted correctly."
  )
}

message("Number of processed sample files found: ", length(all_expression_files))
message("Detected files:")
print(basename(all_expression_files))

# ------------------------------
# 4. Build / read GEO metadata
# ------------------------------
message("\n[Step 2] Reading GEO metadata ...")

build_fallback_metadata <- function() {
  data.frame(
    geo_accession = c(
      "GSM3954946", "GSM3954947", "GSM3954948",
      "GSM3954949", "GSM3954950", "GSM3954951",
      "GSM3954952", "GSM3954953",
      "GSM3954954", "GSM3954955", "GSM3954956", "GSM3954957",
      "GSM3954958"
    ),
    title = c(
      "Non-atrophic gastritis 1",
      "Non-atrophic gastritis 2",
      "Non-atrophic gastritis 3",
      "Chronic atrophic gastritis 1",
      "Chronic atrophic gastritis 2",
      "Chronic atrophic gastritis 3",
      "Wild Intestinal metaplasia 1",
      "Wild Intestinal metaplasia 2",
      "Severe Intestinal metaplasia 1",
      "Severe Intestinal metaplasia 2",
      "Severe Intestinal metaplasia 3",
      "Severe Intestinal metaplasia 4",
      "Early gastric cancer"
    ),
    lesion_manual = c(
      "NAG", "NAG", "NAG",
      "CAG", "CAG", "CAG",
      "Wild IM", "Wild IM",
      "Severe IM", "Severe IM", "Severe IM", "Severe IM",
      "EGC"
    ),
    stringsAsFactors = FALSE
  )
}

sample_metadata <- NULL

if (file.exists(series_matrix_file)) {
  message("Found series matrix file. Attempting to parse with GEOquery ...")

  sample_metadata <- tryCatch({
    geo_object <- getGEO(filename = series_matrix_file, GSEMatrix = TRUE)
    if (is.list(geo_object)) {
      geo_object <- geo_object[[1]]
    }
    pData(geo_object) %>%
      as.data.frame(stringsAsFactors = FALSE)
  }, error = function(e) {
    warning("Failed to parse series matrix file with GEOquery: ", conditionMessage(e))
    NULL
  })
}

if (is.null(sample_metadata) || nrow(sample_metadata) == 0) {
  warning("Using fallback metadata table because series matrix metadata was unavailable.")
  sample_metadata <- build_fallback_metadata()
}

# Ensure there is a sample accession column.
if (!"geo_accession" %in% colnames(sample_metadata)) {
  if ("name" %in% colnames(sample_metadata)) {
    sample_metadata$geo_accession <- sample_metadata$name
  } else if ("title" %in% colnames(sample_metadata)) {
    # Try to infer accession from row names if present.
    if (!is.null(rownames(sample_metadata)) && all(grepl("^GSM", rownames(sample_metadata)))) {
      sample_metadata$geo_accession <- rownames(sample_metadata)
    }
  }
}

if (!"geo_accession" %in% colnames(sample_metadata)) {
  stop("Could not identify a geo_accession column in metadata. Please inspect your series matrix file.")
}

message("Metadata column names:")
print(colnames(sample_metadata))

candidate_text_columns <- colnames(sample_metadata)[
  str_detect(
    string = tolower(colnames(sample_metadata)),
    pattern = "title|source|characteristics|group|lesion|pathology|disease"
  )
]

if (length(candidate_text_columns) == 0) {
  message("No obvious metadata text columns were detected automatically.")
} else {
  message("\nPotential metadata columns for group checking:")
  print(candidate_text_columns)

  for (col_name in candidate_text_columns) {
    message("\nUnique values preview for column: ", col_name)
    unique_values <- unique(as.character(sample_metadata[[col_name]]))
    print(utils::head(unique_values, 20))
  }
}

# ------------------------------
# 5. Identify Severe IM samples as accurately as possible
# ------------------------------
message("\n[Step 3] Filtering Severe Intestinal metaplasia samples ...")

sample_metadata$geo_accession <- as.character(sample_metadata$geo_accession)

# Combine several text fields so we can search them jointly.
text_cols_present <- intersect(candidate_text_columns, colnames(sample_metadata))

if (length(text_cols_present) > 0) {
  metadata_text <- apply(
    sample_metadata[, text_cols_present, drop = FALSE],
    1,
    function(x) paste(na.omit(as.character(x)), collapse = " | ")
  )
} else {
  metadata_text <- rep("", nrow(sample_metadata))
}

metadata_text_lower <- tolower(metadata_text)

positive_pattern <- paste(
  c(
    "severe intestinal metaplasia",
    "intestinal metaplasia 4",   # optional fallback if title formatting is odd
    "gastric metaplasia \\(severe\\)",
    "metaplasia \\(severe\\)",
    "\\bsevere\\b"
  ),
  collapse = "|"
)

negative_pattern <- paste(
  c(
    "wild intestinal metaplasia",
    "non-atrophic gastritis",
    "chronic atrophic gastritis",
    "early gastric cancer",
    "\\bnag\\b",
    "\\bcag\\b",
    "\\begc\\b"
  ),
  collapse = "|"
)

is_severe_candidate <- str_detect(metadata_text_lower, positive_pattern) &
  !str_detect(metadata_text_lower, negative_pattern)

selected_metadata <- sample_metadata[is_severe_candidate, , drop = FALSE]

# If metadata matching fails, use the GEO-known Severe IM GSM IDs as a fallback.
known_severe_gsms <- c("GSM3954954", "GSM3954955", "GSM3954956", "GSM3954957")

if (nrow(selected_metadata) == 0) {
  warning(
    "No Severe IM samples were detected by metadata keyword matching. Falling back to known GEO Severe IM sample IDs."
  )
  selected_metadata <- sample_metadata[sample_metadata$geo_accession %in% known_severe_gsms, , drop = FALSE]
}

if (nrow(selected_metadata) == 0) {
  stop(
    "Severe Intestinal metaplasia samples could not be identified.\n",
    "Please inspect metadata column names and unique values printed above, then adjust the filtering logic."
  )
}

message("Selected Severe IM samples:")
print(selected_metadata[, intersect(c("geo_accession", "title", candidate_text_columns), colnames(selected_metadata)), drop = FALSE])

write.csv(
  selected_metadata,
  file = file.path(results_dir, "severe_im_metadata_selected_samples.csv"),
  row.names = FALSE
)

# ------------------------------
# 6. Match selected samples to local expression files
# ------------------------------
message("\n[Step 4] Matching selected metadata to local processed matrix files ...")

extract_gsm_from_path <- function(x) {
  sub("^(GSM\\d+).*$", "\\1", basename(x))
}

file_table <- data.frame(
  file_path = all_expression_files,
  file_name = basename(all_expression_files),
  geo_accession = extract_gsm_from_path(all_expression_files),
  stringsAsFactors = FALSE
)

selected_files <- file_table %>%
  filter(geo_accession %in% selected_metadata$geo_accession)

if (nrow(selected_files) == 0) {
  stop(
    "No local expression files matched the selected Severe IM samples.\n",
    "Please check your extracted file names in raw_dir."
  )
}

message("Matched Severe IM matrix files:")
print(selected_files)

# ------------------------------
# 7. Read each Severe IM expression matrix
# ------------------------------
message("\n[Step 5] Reading Severe IM expression matrices ...")

read_geo_processed_matrix <- function(file_path) {
  message("Reading file: ", basename(file_path))

  dt <- data.table::fread(file_path, data.table = FALSE)

  if (ncol(dt) < 2) {
    stop("The file appears to have fewer than 2 columns: ", file_path)
  }

  gene_col <- dt[[1]]
  expr_df <- dt[, -1, drop = FALSE]

  expr_matrix <- as.matrix(expr_df)
  mode(expr_matrix) <- "numeric"
  rownames(expr_matrix) <- make.unique(as.character(gene_col))

  Matrix::Matrix(expr_matrix, sparse = TRUE)
}

sample_objects <- list()
cell_metadata_list <- list()

for (i in seq_len(nrow(selected_files))) {
  sample_id <- selected_files$geo_accession[i]
  sample_file <- selected_files$file_path[i]
  sample_name <- selected_files$file_name[i]

  sample_title <- if ("title" %in% colnames(selected_metadata)) {
    selected_metadata$title[match(sample_id, selected_metadata$geo_accession)]
  } else {
    sample_id
  }

  counts_mat <- read_geo_processed_matrix(sample_file)

  original_barcodes <- colnames(counts_mat)
  new_barcodes <- paste(sample_id, original_barcodes, sep = "_")
  colnames(counts_mat) <- new_barcodes

  sample_objects[[sample_id]] <- counts_mat
  cell_metadata_list[[sample_id]] <- data.frame(
    cell_barcode = new_barcodes,
    original_barcode = original_barcodes,
    sample_id = sample_id,
    sample_title = sample_title,
    matrix_file = sample_name,
    stringsAsFactors = FALSE,
    row.names = new_barcodes
  )

  message("  - ", sample_id, ": ", nrow(counts_mat), " genes x ", ncol(counts_mat), " cells")
}

combined_counts <- do.call(cbind, sample_objects)
combined_cell_metadata <- bind_rows(cell_metadata_list)

message("Combined matrix dimension: ", nrow(combined_counts), " genes x ", ncol(combined_counts), " cells")

# ------------------------------
# 8. Construct Seurat object and basic QC
# ------------------------------
message("\n[Step 6] Creating Seurat object and calculating QC metrics ...")

seu <- CreateSeuratObject(
  counts = combined_counts,
  project = "GSE134520_Severe_IM",
  meta.data = combined_cell_metadata,
  min.cells = 3,
  min.features = 0
)

# Calculate mitochondrial percentage.
# Human mitochondrial genes usually start with "MT-".
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

qc_before <- seu@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")
write.csv(qc_before, file = file.path(results_dir, "severe_im_qc_metrics_before_filter.csv"), row.names = FALSE)

message("QC summary before filtering:")
print(summary(seu$nFeature_RNA))
print(summary(seu$nCount_RNA))
print(summary(seu$percent.mt))

# Save QC violin plots for initial inspection.
png(filename = file.path(results_dir, "qc_violin_before_filter.png"), width = 1800, height = 600, res = 150)
print(
  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
)
dev.off()

# Filter cells.
message(
  "Applying filters: nFeature_RNA >= ", min_features,
  "; nFeature_RNA <= ", max_features,
  "; percent.mt < ", max_percent_mt
)

seu <- subset(
  seu,
  subset = nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt < max_percent_mt
)

if (ncol(seu) == 0) {
  stop("No cells remain after filtering. Please relax the QC thresholds.")
}

qc_after <- seu@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")
write.csv(qc_after, file = file.path(results_dir, "severe_im_qc_metrics_after_filter.csv"), row.names = FALSE)

message("Cell number after filtering: ", ncol(seu))
message("Gene number after filtering: ", nrow(seu))

png(filename = file.path(results_dir, "qc_violin_after_filter.png"), width = 1800, height = 600, res = 150)
print(
  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
)
dev.off()

# ------------------------------
# 9. Standard Seurat workflow
# ------------------------------
message("\n[Step 7] Running standard Seurat workflow ...")

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = variable_features_n, verbose = FALSE)
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = pcs_to_compute, verbose = FALSE)
seu <- FindNeighbors(seu, dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)
seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)

# Save an elbow plot for reference.
png(filename = file.path(results_dir, "pca_elbow_plot.png"), width = 1200, height = 900, res = 150)
print(ElbowPlot(seu, ndims = pcs_to_compute))
dev.off()

# Save a quick UMAP after preprocessing so users can immediately inspect it.
png(filename = file.path(results_dir, "umap_quick_preview.png"), width = 1200, height = 900, res = 150)
print(DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Severe IM - quick UMAP preview"))
dev.off()

# ------------------------------
# 10. Save processed object for script 2
# ------------------------------
message("\n[Step 8] Saving processed objects ...")

saveRDS(seu, file = file.path(results_dir, "severe_im_seurat_processed.rds"))
saveRDS(selected_metadata, file = file.path(results_dir, "severe_im_selected_sample_metadata.rds"))

message("Done. Processed Seurat object saved to:")
message(file.path(results_dir, "severe_im_seurat_processed.rds"))
message("Now run the second script: gse134520_severe_im_umap.R")
