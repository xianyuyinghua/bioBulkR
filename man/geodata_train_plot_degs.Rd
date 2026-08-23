\name{geodata_train_plot_degs}
\alias{geodata_train_plot_degs}
\title{Run differential expression, Venn, enrichment, and PPI analyses for GEO data}
\usage{
geodata_train_plot_degs(
  workdir = "/project/project_analysis/26_26YLM199F/02_Result",
  Project_ID = "03_Expression_DEG",
  Data_ID = "GSE28829",
  exp_dir = "01_Expression_data",
  group_levels = NULL,
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05,
  enrich_figure_type = "barplot",
  enrich_pname = "p.value",
  go_kegg_topn = c(10, 10),
  description_wrap_width = 100,
  sample_color = basicR::get_colors(number = 10.1),
  confidence = 150,
  target_gene_file = NULL,
  target_file = NULL,
  target_gene_column = NULL,
  target_name = NULL,
  target_set = NULL
)
}
\arguments{
\item{workdir}{Parent analysis directory. The function changes the working
directory to this path.}
\item{Project_ID}{Project directory identifier used to find or create a
numbered project folder.}
\item{Data_ID}{GEO Series identifier naming the data subdirectory.}
\item{exp_dir}{Expression-data subdirectory within the GEO data directory.}
\item{group_levels}{Optional ordering of sample groups passed to
\code{\link[=group_order_set]{group_order_set()}}.}
\item{logFC}{Absolute log2 fold-change threshold for differential genes.}
\item{p_value}{Differential-expression significance measure:
\code{"p.value"} or \code{"p.adj"}.}
\item{p}{Differential-expression significance threshold.}
\item{enrich_figure_type}{Enrichment plot family requested from
\code{\link[=enrich_go_kegg]{enrich_go_kegg()}}, such as \code{"barplot"},
\code{"dotplot"}, \code{"sankeyplot"}, \code{"treeplot"}, or
\code{"circularplot"}.}
\item{enrich_pname}{Enrichment significance column, \code{"p.value"} or
\code{"p.adj"}.}
\item{go_kegg_topn}{Number of top GO and KEGG terms displayed in enrichment
plots. A length-two value controls GO and KEGG separately.}
\item{description_wrap_width}{Approximate character width for wrapping
enrichment descriptions.}
\item{sample_color}{Colors used for sample annotations in heatmaps.}
\item{confidence}{STRING required-score threshold used for the PPI query.}
\item{target_gene_file}{Response to the prompt asking whether a target-gene
CSV has been placed in \code{03_Venn}; \code{NULL} keeps the interactive
prompt.}
\item{target_file}{Path selected when multiple target-gene CSV files are
found; \code{NULL} keeps the interactive prompt.}
\item{target_gene_column}{Name of the target file's gene-symbol column when it
cannot be identified automatically; \code{NULL} prompts interactively.}
\item{target_name}{Name shown for the target set in Venn plots; \code{NULL}
prompts when needed.}
\item{target_set}{Response indicating whether
\code{Target_venn_result.csv} is available; \code{NULL} keeps the interactive
prompt.}
}
\value{
The differential-expression result table returned in \code{Res$result} after
volcano and heatmap processing. It includes differential statistics, the
classification column, gene symbols, and plot-label columns. The workflow also
saves its full intermediate \code{data} list to
\code{00_data_result/02_data.rds}.
}
\description{
Continue a processed GEO project through array or count-based differential
expression, volcano and heatmap plotting, target-gene intersections, GO and
KEGG enrichment, and STRING protein-interaction network visualization.
}
\details{
The function expects output from the GEO download and processing workflow,
including \code{01_data.rds}, an expression CSV, and
\code{02_Condition.csv} below the selected project and GEO directories. It
creates \code{02_DEGs}, \code{03_Venn}, \code{04_Enrich}, and \code{05_PPI}
directories as necessary.

Array data are analyzed with \code{\link[=Array_DEGs]{Array_DEGs()}}, whereas
high-throughput count data use
\code{\link[=High_Throughtput_Count_DEG]{High_Throughtput_Count_DEG()}}. Both
branches create volcano, rectangular heatmap, and circular heatmap outputs.
The target-gene file in \code{03_Venn} is intersected with all, upregulated,
and downregulated DEGs. Enrichment is then run separately for those three
intersections, followed by an online STRING query and network plot.

This is a file-writing and network-enabled workflow. It writes CSV, RDS, PDF,
and PNG outputs, can overwrite existing results, changes the process
working directory, and can request console input. It also loads numerous
analysis and plotting packages with \code{library()}. If genes in an existing
\code{Target_venn_result.csv} disagree with the newly calculated intersection,
the current implementation calls \code{quit(save = "no")}.
}
\examples{
\dontrun{
deg_result <- geodata_train_plot_degs(
  workdir = "/project/project_analysis/26_26YLM199F/02_Result",
  Project_ID = "03_Expression_DEG",
  Data_ID = "GSE28829",
  group_levels = c("Control", "Disease"),
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05,
  enrich_figure_type = "treeplot",
  enrich_pname = "p.value",
  go_kegg_topn = c(10, 10),
  confidence = 400,
  target_gene_column = "Symbol",
  target_name = "Target",
  target_set = "no"
)
}
}
\seealso{
\code{\link[=geodata_download_and_process]{geodata_download_and_process()}},
\code{\link[=Array_DEGs]{Array_DEGs()}},
\code{\link[=High_Throughtput_Count_DEG]{High_Throughtput_Count_DEG()}},
\code{\link[=Venn_plot]{Venn_plot()}},
\code{\link[=enrich_go_kegg]{enrich_go_kegg()}},
\code{\link[=string_ppi_get_data]{string_ppi_get_data()}},
\code{\link[=network_graph_plot]{network_graph_plot()}}
}
