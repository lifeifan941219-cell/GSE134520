#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GEOquery)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

flatten_metadata <- function(x) {
  if (is.null(x)) {
    return("")
  }
  paste(as.character(x), collapse = " | ")
}

classify_lesion <- function(text) {
  lowered <- tolower(text)

  if (str_detect(lowered, "severe") && str_detect(lowered, "metaplasia")) {
    return("IMS")
  }
  if (str_detect(lowered, "\\bims\\b")) {
    return("IMS")
  }
  if (str_detect(lowered, "wild") && str_detect(lowered, "metaplasia")) {
    return("IMW")
  }
  if (str_detect(lowered, "\\bimw\\b")) {
    return("IMW")
  }
  if (str_detect(lowered, "early gastric cancer") || str_detect(lowered, "\\begc\\b")) {
    return("EGC")
  }
  if (str_detect(lowered, "atrophic gastritis") || str_detect(lowered, "\\bcag\\b")) {
    return("CAG")
  }
  if (str_detect(lowered, "non-atrophic") || str_detect(lowered, "\\bnag\\b")) {
    return("NAG")
  }
  if (str_detect(lowered, "intestinal metaplasia")) {
    return("IM_unspecified")
  }
  "unknown"
}

get_series_object <- function(geo_id, cache_dir) {
  message("Downloading GEO series metadata for ", geo_id, " ...")
  getGEO(filename = NULL, GEO = geo_id, destdir = cache_dir, GSEMatrix = FALSE)
}

download_supp_files <- function(geo_id, cache_dir) {
  supp_dir <- file.path(cache_dir, paste0(geo_id, "_supp"))
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  message("Downloading GEO supplementary files for ", geo_id, " ...")
  getGEOSuppFiles(GEO = geo_id, baseDir = supp_dir, fetch_files = TRUE)
  supp_dir
}

build_sample_metadata <- function(gse) {
  gsm_list <- GSMList(gse)
  metadata <- lapply(names(gsm_list), function(gsm_name) {
    gsm <- gsm_list[[gsm_name]]
    meta <- Meta(gsm)
    merged <- paste(
      c(gsm_name, unlist(lapply(names(meta), function(key) paste0(key, ": ", flatten_metadata(meta[[key]]))))),
      collapse = " || "
    )
    tibble(
      gsm = gsm_name,
      title = flatten_metadata(meta$title),
      source_name = flatten_metadata(meta$source_name_ch1),
      characteristics = flatten_metadata(meta$characteristics_ch1),
      lesion_class = classify_lesion(merged),
      metadata_blob = merged
    )
  })

  bind_rows(metadata) %>% arrange(lesion_class, gsm)
}

unpack_archives <- function(root_dir) {
  archives <- list.files(root_dir, pattern = "[.](tar|tar[.]gz|tgz)$", recursive = TRUE, full.names = TRUE)
  extracted_dirs <- character(0)

  for (archive in archives) {
    archive_name <- basename(archive)
    target_dir <- file.path(dirname(archive), str_replace(archive_name, "[.](tar|tar[.]gz|tgz)$", ""))
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

    existing_files <- list.files(target_dir, recursive = TRUE, full.names = TRUE)
    if (length(existing_files) == 0) {
      try(utils::untar(archive, exdir = target_dir, tar = "internal"), silent = TRUE)
      extracted_files <- list.files(target_dir, recursive = TRUE, full.names = TRUE)
      if (length(extracted_files) == 0) {
        try(utils::untar(archive, exdir = target_dir), silent = TRUE)
      }
    }

    extracted_dirs <- c(extracted_dirs, target_dir)
  }

  unique(c(root_dir, extracted_dirs))
}

write_directory_manifest <- function(search_roots, outdir) {
  manifest_path <- file.path(outdir, "input_file_manifest.txt")
  manifest_entries <- unlist(lapply(search_roots, function(path) {
    if (!dir.exists(path)) {
      return(character(0))
    }
    list.files(path, recursive = TRUE, full.names = TRUE)
  }))
  writeLines(sort(unique(manifest_entries)), manifest_path)
  manifest_path
}

find_10x_dirs <- function(root_dir) {
  all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  candidates <- vapply(all_dirs, function(dir_path) {
    files <- basename(list.files(dir_path, full.names = FALSE))
    has_matrix <- any(files %in% c("matrix.mtx", "matrix.mtx.gz"))
    has_barcodes <- any(files %in% c("barcodes.tsv", "barcodes.tsv.gz"))
    has_features <- any(files %in% c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))
    has_matrix && has_barcodes && has_features
  }, logical(1))
  unique(all_dirs[candidates])
}

find_10x_file_sets <- function(root_dir) {
  all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  tenx_files <- all_files[str_detect(basename(all_files), "(matrix[.]mtx|barcodes[.]tsv|genes[.]tsv|features[.]tsv)([.]gz)?$")]
  if (length(tenx_files) == 0) {
    return(list())
  }

  prefixes <- basename(tenx_files) %>%
    str_replace("_(matrix[.]mtx|barcodes[.]tsv|genes[.]tsv|features[.]tsv)([.]gz)?$", "")

  split_files <- split(tenx_files, prefixes)
  lapply(split_files, function(paths) {
    names(paths) <- basename(paths)
    list(
      matrix = unname(paths[str_detect(names(paths), "matrix[.]mtx([.]gz)?$")][1]),
      barcodes = unname(paths[str_detect(names(paths), "barcodes[.]tsv([.]gz)?$")][1]),
      features = unname(paths[str_detect(names(paths), "(features|genes)[.]tsv([.]gz)?$")][1])
    )
  })
}

choose_sample_input <- function(gsm, tenx_dirs, tenx_file_sets) {
  dir_matches <- tenx_dirs[str_detect(tenx_dirs, fixed(gsm))]
  if (length(dir_matches) > 0) {
    return(list(type = "dir", value = dir_matches[order(nchar(dir_matches))][1]))
  }

  prefix_matches <- names(tenx_file_sets)[str_detect(names(tenx_file_sets), fixed(gsm))]
  if (length(prefix_matches) > 0) {
    return(list(type = "files", value = tenx_file_sets[[prefix_matches[1]]]))
  }

  stop("Could not find a 10x directory or file set for sample ", gsm, call. = FALSE)
}

read_sample <- function(sample_input, gsm) {
  if (sample_input$type == "dir") {
    counts <- Read10X(data.dir = sample_input$value)
  } else {
    counts <- ReadMtx(
      mtx = sample_input$value$matrix,
      features = sample_input$value$features,
      cells = sample_input$value$barcodes,
      feature.column = 2
    )
  }
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  seu <- CreateSeuratObject(counts = counts, project = gsm, min.cells = opt$min_cells, min.features = 0)
  seu$gsm <- gsm
  seu
}

save_umap_plot <- function(seu, group_by, output_file, label = FALSE) {
  plot <- DimPlot(seu, reduction = "umap", group.by = group_by, label = label, repel = TRUE) +
    theme_bw(base_size = 12)
  ggsave(filename = output_file, plot = plot, width = 8, height = 6, dpi = 300)
}

option_list <- list(
  make_option("--geo-id", type = "character", default = "GSE134520", help = "GEO series accession [default %default]"),
  make_option("--input-dir", type = "character", default = NULL, help = "Optional directory with downloaded supplementary files."),
  make_option("--cache-dir", type = "character", default = "data/cache", help = "Directory for GEO downloads [default %default]"),
  make_option("--outdir", type = "character", default = "results/gse134520_severe_im", help = "Output directory [default %default]"),
  make_option("--min-features", type = "integer", default = 200, help = "Minimum detected genes per cell [default %default]"),
  make_option("--min-cells", type = "integer", default = 3, help = "Minimum detected cells per gene [default %default]"),
  make_option("--max-mt-pct", type = "double", default = 20, help = "Maximum mitochondrial percent [default %default]"),
  make_option("--nfeatures", type = "integer", default = 2000, help = "Number of variable genes [default %default]"),
  make_option("--dims", type = "integer", default = 30, help = "Number of PCA dimensions [default %default]"),
  make_option("--resolution", type = "double", default = 0.5, help = "Clustering resolution [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

gse <- get_series_object(opt$geo_id, opt$cache_dir)
sample_df <- build_sample_metadata(gse)
write_csv(sample_df, file.path(opt$outdir, "gse134520_sample_metadata_all.csv"))

severe_df <- sample_df %>% filter(lesion_class == "IMS")
if (nrow(severe_df) == 0) {
  stop(
    "No Severe Intestinal Metaplasia samples were identified from GEO metadata. Check gse134520_sample_metadata_all.csv and adjust classify_lesion() if necessary.",
    call. = FALSE
  )
}

input_dir <- opt$input_dir
if (is.null(input_dir) || identical(input_dir, "")) {
  input_dir <- download_supp_files(opt$geo_id, opt$cache_dir)
}

search_roots <- unpack_archives(input_dir)
manifest_path <- write_directory_manifest(search_roots, opt$outdir)
tenx_dirs <- unique(unlist(lapply(search_roots, find_10x_dirs)))
tenx_file_sets <- unlist(lapply(search_roots, find_10x_file_sets), recursive = FALSE)
if (length(tenx_dirs) == 0 && length(tenx_file_sets) == 0) {
  stop(
    paste0(
      "No 10x matrix directories or flat GEO matrix/barcode/feature files were found under input-dir. ",
      "Please verify the supplementary files. A manifest was written to: ", manifest_path
    ),
    call. = FALSE
  )
}

sample_objects <- lapply(severe_df$gsm, function(gsm) {
  sample_input <- choose_sample_input(gsm, tenx_dirs, tenx_file_sets)
  read_sample(sample_input, gsm)
})
names(sample_objects) <- severe_df$gsm

if (length(sample_objects) == 1) {
  seu <- sample_objects[[1]]
} else {
  seu <- merge(sample_objects[[1]], y = sample_objects[-1], add.cell.ids = names(sample_objects), project = opt$geo_id)
}
seu <- subset(seu, subset = nFeature_RNA >= opt$min_features)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = percent.mt <= opt$max_mt_pct)

sample_meta <- severe_df %>% select(gsm, title, source_name, characteristics, lesion_class)
seu@meta.data <- seu@meta.data %>% tibble::rownames_to_column("cell_id") %>%
  left_join(sample_meta, by = "gsm") %>%
  tibble::column_to_rownames("cell_id")

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = opt$nfeatures)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = opt$dims)
seu <- FindNeighbors(seu, dims = seq_len(opt$dims))
seu <- FindClusters(seu, resolution = opt$resolution)
seu <- RunUMAP(seu, dims = seq_len(opt$dims))

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv(markers, file.path(opt$outdir, "severe_im_cluster_markers.csv"))
write_csv(severe_df, file.path(opt$outdir, "severe_im_sample_metadata.csv"))
write_csv(seu@meta.data %>% tibble::rownames_to_column("cell_id"), file.path(opt$outdir, "severe_im_cell_qc_and_clusters.csv"))
saveRDS(seu, file.path(opt$outdir, "severe_im_seurat.rds"))

save_umap_plot(seu, "seurat_clusters", file.path(opt$outdir, "umap_clusters.png"), label = TRUE)
save_umap_plot(seu, "gsm", file.path(opt$outdir, "umap_sample.png"), label = FALSE)

message("Analysis completed successfully.")
message("Cells retained: ", ncol(seu))
message("Genes retained: ", nrow(seu))
message("Output directory: ", opt$outdir)
