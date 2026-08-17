# bioBulkR

`bioBulkR` 是一个用于常规 bulk RNA-seq 分析与可视化的个人 R 包。目前提供
GO/KEGG 富集分析以及条形图、气泡图、桑基图和环形图，后续可继续加入差异
表达分析、PPI 网络和其他下游分析函数。

## 安装

首先安装 `remotes`：

```r
install.packages("remotes")
```

然后从 GitHub 安装：

```r
remotes::install_github("xianyuyinghua/bioBulkR", dependencies = TRUE)
```

如果尚未配置 Bioconductor，请先运行：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler", "ComplexHeatmap", "org.Hs.eg.db", "org.Mm.eg.db"
))
```

`ggsankey` 由 GitHub 安装，`basicR` 会从作者的 GitHub 仓库安装。仅使用本地
KEGG 数据时还需要安装旧版 `KEGG.db`；通常建议使用默认的在线模式。

## 快速开始

```r
library(bioBulkR)

genes <- c("TP53", "EGFR", "MYC", "CDKN1A", "BAX")

result <- enrich_go_kegg(
  genes = genes,
  species = "9606",                  # 人：9606；小鼠：10090
  p_name = "p.adj",                  # p.value 或 p.adj
  filter_pvalue = 0.05,
  pvalueCutoff = 0.5,
  qvalueCutoff = 1,
  kegg_analysis_method = "online",
  plot_type = c("barplot", "dotplot"),
  go_kegg_topn = c(10, 10),
  description_wrap_width = 60
)

head(result$result)
result$plot$barplot$gokegg
result$plot$dotplot$go
```

## 参数说明

- `genes`：基因 SYMBOL 字符向量。函数会将其转换为 ENTREZID；无法映射的
  SYMBOL 不参与后续分析。
- `species`：NCBI taxonomy ID。当前支持人 `"9606"` 和小鼠 `"10090"`。
- `p_name`：最终筛选和作图使用的显著性指标。`"p.value"` 表示原始 P 值，
  `"p.adj"` 表示 BH 校正后的 P 值。
- `filter_pvalue`：最终保留通路的显著性阈值，要求所选指标严格小于该值。
- `pvalueCutoff`：传递给 clusterProfiler 的初始 P 值阈值。
- `qvalueCutoff`：传递给 clusterProfiler 的 q-value 阈值。
- `kegg_analysis_method`：`"online"` 使用在线 KEGG 注释；`"local"` 使用
  `KEGG.db` 的旧版内部注释。
- `plot_type`：`NULL` 表示不画图；也可选择 `"barplot"`、`"dotplot"`、
  `"sankeyplot"`、`"circularplot"`，并可一次传入多个类型。
- `scankeyplot_mode`：桑基图布局。`"sankey_buble"` 为基因—通路—气泡；
  `"buble_sankey"` 为气泡—通路—基因。参数保留了原脚本中的拼写。
- `fill_colors`：图形使用的颜色向量，默认由 `basicR::get_colors()` 生成。
- `go_kegg_topn`：长度为 1 时，GO 和 KEGG 使用同一展示数量；长度为 2 时，
  第一个数是每类 GO 展示数量，第二个数是 KEGG 展示数量。
- `description_wrap_width`：通路描述换行的近似字符宽度。

安装包后可以用以下命令查看更完整的参数、返回值和示例：

```r
?enrich_go_kegg
```

## 返回值

函数返回包含 `result` 和 `plot` 的列表。`result` 是合并后的 GO/KEGG 富集
结果表；`plot` 按图形类型组织，每种类型包含 `go`、`kegg` 和 `gokegg`。
环形图返回的是绘图函数，例如：

```r
result$plot$circularplot$gokegg()
```
