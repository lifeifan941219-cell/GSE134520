# GSE134520 Severe Intestinal Metaplasia UMAP（R 版本）

这个仓库现在提供 **R 语言 / Seurat** 脚本，用于从 **GSE134520** 中筛选 **Severe Intestinal Metaplasia（IMS）** 细胞，完成质控、降维、聚类，并输出 UMAP 图。

## 文件说明

- `scripts/gse134520_severe_im_umap.R`：主分析脚本。
- `results/gse134520_severe_im/`：默认输出目录。

## 依赖安装

建议先在 R 中安装下面这些包：

```r
install.packages(c("optparse", "dplyr", "stringr", "readr", "ggplot2", "BiocManager"))
BiocManager::install(c("GEOquery", "Seurat"))
```

## 直接运行

```bash
Rscript scripts/gse134520_severe_im_umap.R
```

脚本默认会：

1. 自动下载 `GSE134520` 的 GEO 元数据；
2. 自动下载 GEO supplementary files；
3. 根据 GEO 样本注释识别 `Severe Intestinal Metaplasia` / `IMS` 样本；
4. 读取对应 10x 表达矩阵；
5. 用 Seurat 完成 QC、标准化、高变基因、PCA、邻近图、聚类和 UMAP；
6. 导出聚类结果、marker 基因和图片。

## 常用参数

```bash
Rscript scripts/gse134520_severe_im_umap.R \
  --cache-dir data/cache \
  --outdir results/gse134520_severe_im \
  --min-features 200 \
  --min-cells 3 \
  --max-mt-pct 20 \
  --nfeatures 2000 \
  --dims 30 \
  --resolution 0.5
```

如果你已经手动下载了 GEO 补充文件，可以指定本地目录：

```bash
Rscript scripts/gse134520_severe_im_umap.R --input-dir /path/to/GSE134520_supp
```

## 输出文件

- `gse134520_sample_metadata_all.csv`：所有 GEO 样本及自动识别分组。
- `severe_im_sample_metadata.csv`：筛选出的 IMS 样本信息。
- `severe_im_seurat.rds`：Seurat 对象。
- `severe_im_cluster_markers.csv`：每个 cluster 的 marker 基因。
- `severe_im_cell_qc_and_clusters.csv`：每个细胞的样本、QC 和 cluster 信息。
- `umap_clusters.png`：按 cluster 着色的 UMAP。
- `umap_sample.png`：按样本着色的 UMAP。

## 说明

根据公开资料，`GSE134520` 包含 NAG、CAG、IMW、IMS 和 EGC 等病理阶段。这个脚本优先从 GEO 元数据里自动识别 `Severe Intestinal Metaplasia` / `IMS` 标签，不把样本编号硬编码到脚本里，便于你后续复用和调整。

补充说明：GEO 的 `GSE134520_RAW.tar` 有时解压后不是标准的“每个样本一个文件夹”，而是所有 `GSM*_matrix/barcodes/genes` 文件平铺在同一目录。当前脚本已经同时支持这两种结构。


## 在 RStudio 中运行

如果你想直接在 **RStudio** 中点开脚本运行，推荐使用：

- `scripts/gse134520_severe_im_rstudio.R`：适合在 RStudio 中 `Source` 或逐段运行；
- `scripts/gse134520_severe_im_umap.R`：适合在终端里用 `Rscript` 命令运行。

### RStudio 运行步骤

1. 在 RStudio 中打开仓库目录；
2. 打开 `scripts/gse134520_severe_im_rstudio.R`；
3. 根据需要修改顶部的 `CONFIG` 参数；
4. 点击 **Source**，或者逐段运行；
5. 结果会输出到 `results/gse134520_severe_im/`。

### RStudio 脚本的特点

- 不依赖命令行参数解析；
- 顶部集中配置参数；
- 主流程封装在 `run_gse134520_severe_im(CONFIG)` 中；
- 运行完成后会返回一个 Seurat 对象 `seu`，方便你继续在 Console 中画图或查看 cluster。

例如运行完成后，你可以继续在 RStudio Console 里执行：

```r
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(seu, reduction = "umap", group.by = "gsm")
head(seu@meta.data)
```
