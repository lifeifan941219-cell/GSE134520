# GSE134520：Severe Intestinal metaplasia 结果基础上的上皮细胞注释与 REG4/HNF4G 投影分析

本说明对应以下两个脚本：

1. `gse134520_epithelial_annotation.R`
2. `gse134520_epithelial_reg4_hnf4g_projection.R`

这套代码默认你已经完成了 **Severe Intestinal metaplasia（Severe IM）** 样本的初始 Seurat 聚类和 UMAP，并且已经保存了一个可读取的 `.rds` Seurat 对象。新的脚本**不会重新从 GEO 原始矩阵开始**，而是尽量复用你之前的分析结果继续向下做：

- cluster 初步细胞类型注释；
- 识别并标注上皮细胞群；
- 提取上皮细胞；
- 上皮细胞重新聚类与 UMAP；
- REG4 / HNF4G 在上皮亚群中的投影、比较和平均表达统计。

---

## 一、工作目录

默认工作目录：

```r
D:/Rdata/GSE134520
```

### 在 RStudio 中设置工作目录

#### 方法 1：菜单设置
- `Session` → `Set Working Directory` → `Choose Directory...`
- 选择 `D:/Rdata/GSE134520`

#### 方法 2：代码设置

```r
setwd("D:/Rdata/GSE134520")
getwd()
```

正确时应显示：

```r
[1] "D:/Rdata/GSE134520"
```

---

## 二、运行顺序

请按以下顺序运行：

1. **先运行** `gse134520_epithelial_annotation.R`
2. **再运行** `gse134520_epithelial_reg4_hnf4g_projection.R`

### 推荐运行方式

#### 方式 A：整段运行

```r
source("gse134520_epithelial_annotation.R", echo = TRUE)
source("gse134520_epithelial_reg4_hnf4g_projection.R", echo = TRUE)
```

#### 方式 B：逐段运行
适合初学者：
- 打开脚本；
- 从上到下逐段运行；
- 每一步观察 Console 输出；
- 如果对象文件名、metadata 列名或注释列名与你的对象不一致，可按脚本前部参数提示修改。

---

## 三、你需要提前准备什么

你至少需要一个**已经完成 Severe IM 初始聚类与 UMAP 的 Seurat 对象 `.rds` 文件**。

例如它可能叫：

- `severe_im_seurat_processed.rds`
- `seu.rds`
- `seurat_object.rds`
- 或其他名称

脚本已经考虑到**文件名不确定**的情况：

- 会自动列出工作目录及常见结果目录中的 `.rds` 文件；
- 会尝试根据文件名自动优先选择更像 Seurat 对象的文件；
- 如果自动选择不符合你的情况，你可以直接在脚本顶部手动修改 `input_seurat_file`。

---

## 四、输出目录说明

脚本会自动在 `D:/Rdata/GSE134520` 下创建以下目录：

```text
results/
figures/
tables/
objects/
```

推荐理解如下：

- `figures/`：保存 UMAP、FeaturePlot、DotPlot、VlnPlot、Heatmap 等图
- `tables/`：保存 marker 表、top10 marker、平均表达矩阵、cluster 注释表等
- `objects/`：保存注释后的对象、上皮子对象等
- `results/`：保存一些综合性输出或日志型文件

---

## 五、第一个脚本做什么：cluster 注释

`gse134520_epithelial_annotation.R` 主要完成：

1. 读取你之前保存的 `.rds` Seurat 对象；
2. 自动检查：
   - metadata 列名；
   - 当前 cluster 信息；
   - 是否已有 UMAP；
3. 使用常见 marker panel 做初步细胞类型注释；
4. 输出：
   - 每个 cluster 的 marker genes；
   - 每个 cluster 的 top10 marker；
   - 注释前 UMAP；
   - 注释后 UMAP；
5. 在 metadata 中新增注释列，例如：
   - `celltype`
   - `celltype_cluster`
   - `celltype_score`
6. 保存注释后的 Seurat 对象。

### 脚本中提供的常见 marker panel

- **上皮细胞**：`EPCAM, KRT8, KRT18, KRT19, KRT20`
- **goblet / intestinal-like epithelial**：`MUC2, TFF3, REG4, GPA33, SI, FABP1, SLC26A3`
- **stem-like / progenitor-like epithelial**：`OLFM4, SOX9, MKI67, TOP2A`
- **T cells**：`CD3D, CD3E, IL7R`
- **NK cells**：`NKG7, GNLY, KLRD1`
- **B cells**：`CD79A, MS4A1, CD74`
- **plasma cells**：`MZB1, JCHAIN, IGKC`
- **myeloid / macrophage**：`LYZ, FCN1, C1QC, APOE, CD68`
- **fibroblasts**：`COL1A1, COL1A2, DCN, LUM`
- **endothelial cells**：`PECAM1, EMCN, VWF, KDR`
- **neuroendocrine cells**：`CHGA, CHGB, NEUROD1`

### 如何判断哪些 cluster 是上皮细胞

请优先结合以下几类信息：

1. **cluster marker genes**
   - 如果 cluster 高表达 `EPCAM/KRT8/KRT18/KRT19/KRT20`，通常提示上皮来源；
2. **marker panel 打分结果**
   - 脚本会按 marker panel 计算各类细胞得分；
3. **UMAP 上 marker 可视化**
   - 上皮 cluster 往往在 `FeaturePlot(EPCAM)`、`FeaturePlot(KRT19)` 上较明显；
4. **肠化生偏向 marker**
   - 若同时高表达 `MUC2/TFF3/REG4/GPA33/FABP1/SI`，更提示 intestinal-like / goblet-like epithelial。

> 注意：脚本会自动给出“建议注释”，但这属于**初步注释**。如果你结合 marker 结果后认为某个 cluster 应改名，可手动修改脚本里的 `manual_cluster_annotation` 映射表再运行一次。

---

## 六、第二个脚本做什么：上皮细胞再分析 + REG4/HNF4G 投影

`gse134520_epithelial_reg4_hnf4g_projection.R` 主要完成：

1. 读取第一个脚本输出的“注释后 Seurat 对象”；
2. 根据 metadata 中的 `celltype`（或候选注释列）提取上皮细胞；
3. 对上皮细胞重新运行 Seurat 标准流程：
   - `NormalizeData`
   - `FindVariableFeatures`
   - `ScaleData`
   - `RunPCA`
   - `FindNeighbors`
   - `FindClusters`
   - `RunUMAP`
4. 输出上皮亚群的：
   - UMAP 图；
   - 带 cluster 编号的 UMAP；
   - marker genes；
   - top10 marker；
5. 输出以下 marker 在上皮亚群中的表达：
   - `REG4`
   - `HNF4G`
   - `EPCAM`
   - `KRT20`
   - `MUC2`
   - `TFF3`
   - `FABP1`
   - `GPA33`
   - `OLFM4`
   - `SOX9`
   - `MKI67`
6. 完成 REG4/HNF4G 投影分析：
   - `FeaturePlot`
   - `DotPlot`
   - `VlnPlot`
   - 各上皮亚群平均表达量统计；
   - 自动输出 REG4 最高表达 cluster；
   - 自动输出 HNF4G 最高表达 cluster；
   - 标记 `REG4 positive`、`HNF4G positive`、`double positive` 细胞并绘图。

---

## 七、如何确认 Seurat 对象是否成功读取

运行第一个脚本后，请重点观察 Console：

1. 是否打印出找到的 `.rds` 文件列表；
2. 是否成功读入 Seurat 对象；
3. 是否显示对象维度，例如：
   - 基因数；
   - 细胞数；
4. 是否打印出 metadata 列名；
5. 是否提示已有 PCA / UMAP / cluster 信息。

如果脚本提示：

```r
The selected file is not a Seurat object
```

说明你选错了 `.rds` 文件，需要手动改脚本顶部的 `input_seurat_file`。

---

## 八、如何查看 cluster marker 和 top10 marker

运行第一个脚本后，到 `tables/annotation/` 目录查看：

- `cluster_markers_all.csv`
- `cluster_markers_top10.csv`
- `cluster_annotation_summary.csv`
- `cluster_marker_panel_scores.csv`

运行第二个脚本后，到 `tables/epithelial/` 目录查看：

- `epithelial_markers_all.csv`
- `epithelial_markers_top10.csv`
- `epithelial_cluster_average_expression.csv`
- `reg4_hnf4g_average_expression_by_cluster.csv`
- `reg4_hnf4g_positive_cell_counts.csv`

---

## 九、如何查看 REG4 和 HNF4G 分别在哪个上皮亚群高表达

请重点看以下输出：

### 1. DotPlot
文件通常在：

- `figures/epithelial/reg4_hnf4g_dotplot.png`

用途：
- 比较 `REG4` 和 `HNF4G` 在各上皮亚群的平均表达与阳性比例。

### 2. VlnPlot
文件通常在：

- `figures/epithelial/reg4_hnf4g_violin.png`

用途：
- 查看不同上皮亚群中 `REG4` / `HNF4G` 的表达分布。

### 3. 平均表达表
文件通常在：

- `tables/epithelial/reg4_hnf4g_average_expression_by_cluster.csv`

用途：
- 直接查看每个上皮 cluster 中 `REG4` 和 `HNF4G` 的平均表达量；
- 脚本还会在 Console 中自动打印：
  - `REG4 最高表达的上皮 cluster`
  - `HNF4G 最高表达的上皮 cluster`

### 4. FeaturePlot
文件通常在：

- `figures/epithelial/featureplot_REG4.png`
- `figures/epithelial/featureplot_HNF4G.png`

用途：
- 直观看这两个基因在 UMAP 上主要投影到哪些上皮亚群。

---

## 十、REG4 / HNF4G 的解释思路

### 1. 如果 REG4 更偏 goblet-like / secretory epithelial 亚群
这通常提示：
- 这些细胞更偏向肠型分化中的分泌性上皮程序；
- 常常可伴随 `MUC2`、`TFF3`、`GPA33` 等一起升高；
- 在肠化生背景下，可能代表更明显的肠型化、杯状细胞样或分泌型上皮特征。

### 2. 如果 HNF4G 更偏 intestinal differentiation / enterocyte-like 亚群
这通常提示：
- 这些细胞更偏向肠上皮吸收样/分化样程序；
- 常可与 `FABP1`、`SI`、`GPA33`、`KRT20` 等共同出现；
- 在肠化生中可理解为更成熟的肠分化方向。

### 3. 如果 REG4 和 HNF4G 部分重叠
这可能意味着：
- 同一类肠化生上皮细胞内部存在连续谱；
- 一部分细胞同时具有分泌样与肠分化样特征；
- 也可能代表从 stem/progenitor-like 向更成熟 intestinal-like 状态过渡。

### 4. 如果 REG4 和 HNF4G 明显分离
这可能意味着：
- 上皮细胞内部存在较清晰的功能分工或分化分支；
- 一支更偏 goblet/secretory-like；
- 另一支更偏 enterocyte-like / intestinal differentiation-like。

---

## 十一、常见报错与解决方法

### 报错 1：找不到 `.rds` 文件
可能提示：

```r
No .rds files found
```

解决方法：
- 确认工作目录是 `D:/Rdata/GSE134520`；
- 确认你之前保存的 Seurat 对象确实在该目录或其常见子目录中；
- 直接手动修改脚本中的：

```r
input_seurat_file <- "D:/Rdata/GSE134520/你的对象文件名.rds"
```

### 报错 2：读取的不是 Seurat 对象
可能提示：

```r
The selected file is not a Seurat object
```

解决方法：
- 说明该 `.rds` 可能只是一个表格或列表；
- 请换成真正保存 Seurat 对象的 `.rds` 文件。

### 报错 3：对象里没有 cluster 信息
解决方法：
- 检查对象是否已经完成 `FindClusters()`；
- 看 metadata 中是否存在 `seurat_clusters`；
- 如果没有，可先回到之前的 Severe IM 聚类脚本重新保存对象。

### 报错 4：对象里没有 UMAP
解决方法：
- 检查对象是否已经运行过 `RunUMAP()`；
- 若没有，脚本会尽量提示你重新运行上一步流程；
- 对于上皮子对象，第二个脚本会重新生成新的 UMAP。

### 报错 5：没有识别到上皮细胞
解决方法：
- 查看第一个脚本输出的：
  - marker 表；
  - cluster 注释总结表；
  - marker panel DotPlot；
- 检查 `celltype` 是否被标成：
  - `Epithelial`
  - `Goblet / intestinal-like epithelial`
  - `Stem-like / progenitor epithelial`
- 如命名不同，请修改第二个脚本中的 `epithelial_pattern`。

### 报错 6：某些 marker 不在对象中
解决方法：
- 单细胞对象可能因为过滤、基因名格式或低表达导致某些基因缺失；
- 脚本会自动打印哪些 marker 存在、哪些不存在；
- 若你的对象用的是别名或 Ensembl ID，需要先转换基因名。

### 报错 7：FindAllMarkers 太慢
解决方法：
- 可以先只运行部分 cluster；
- 或降低细胞数后测试；
- 也可以改脚本参数，放宽/收紧 `min.pct` 和 `logfc.threshold`。

### 报错 8：Seurat v5 提示 `GetAssayData()` 的 `slot` 参数已废弃
可能提示：

```r
The `slot` argument of `GetAssayData()` was deprecated in SeuratObject 5.0.0
```

解决方法：
- 当前版本脚本已经改成 **Seurat v4 / v5 兼容写法**；
- 如果你仍看到旧报错，请确认运行的是最新脚本，而不是之前保存的旧版本；
- 若你自己改代码读取表达矩阵，请优先使用 `layer = "data"`，而不是 `slot = "data"`。

---

## 十二、缺失 R 包如何安装

脚本已经包含自动安装逻辑，但如果你想手动安装，可在 RStudio 中运行：

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix", "data.table"))
```

如果某些依赖安装较慢，请耐心等待，或更换镜像：

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

---

## 十三、建议你第一次运行时重点关注什么

建议你第一次运行时，重点关注 Console 中以下输出：

1. 找到了哪些 `.rds` 文件；
2. 自动选择了哪个对象；
3. metadata 列名是什么；
4. cluster 列名是什么；
5. marker panel 是否能正常出图；
6. 哪些 cluster 被建议注释为上皮；
7. 是否成功提取到 epithelial cells；
8. REG4 和 HNF4G 最高表达分别落在哪个上皮 cluster。

如果这些步骤都正常，那么后续做更深入的：
- epithelial 再注释；
- goblet-like / enterocyte-like 细分；
- stem-like 上皮轨迹分析；
- REG4/HNF4G 双阳性群进一步 marker 分析；

都会比较顺利。
