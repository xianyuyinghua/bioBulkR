\name{gsea_analysis}
\alias{gsea_analysis}
\title{Run GSEA for multiple ranked feature groups and MSigDB collections}
\usage{
gsea_analysis(
  feature_value_df,
  feature_col,
  value_col,
  separate_col,
  background_gene_sets = c("GO", "KEGG", "Reactome", "Hallmark", "C2All", "Wikipathway"),
  colors = basicR::get_colors(number = 10.1),
  species = 9606,
  p_threshold_type = "pvalue",
  topn = 5,
  subplots = 1:3,
  base_size = 16,
  legend_text_size = 16,
  width_break = 60,
  gmt_MSigDB_dir = "/reference/functional/gene_sets/MSigDB"
)
}
\arguments{
\item{feature_value_df}{A data frame containing gene identifiers, numeric
ranking values, and a grouping column.}
\item{feature_col}{Name of the gene or feature identifier column.}
\item{value_col}{Name of the numeric ranking-statistic column.}
\item{separate_col}{Name of the column whose factor levels define separate
ranked lists and analyses.}
\item{background_gene_sets}{A character vector selecting collections by
case-insensitive patterns. Supported patterns include \code{"GO"},
\code{"BP"}, \code{"CC"}, \code{"MF"}, \code{"KEGG"},
\code{"Reactome"}, \code{"Hallmark"}, and \code{"Wikipathway"}.}
\item{colors}{A color vector for enrichment-score curves. Supply at least as
many colors as pathways shown together.}
\item{species}{NCBI taxonomy identifier: \code{9606} for human or
\code{10090} for mouse.}
\item{p_threshold_type}{Result column used for the final 0.05 filter, normally
\code{"pvalue"} or \code{"p.adjust"}.}
\item{topn}{Maximum number of significant pathways displayed per plot.}
\item{subplots}{An integer subset of \code{1:3} selecting the enrichment-score,
gene-hit, and ranked-metric panels.}
\item{base_size}{Base plotting text size.}
\item{legend_text_size}{Enrichment-curve legend text size.}
\item{width_break}{Approximate character width used to wrap pathway labels.}
\item{gmt_MSigDB_dir}{Root directory containing the expected human and mouse
MSigDB version 2025.1 GMT directory structure.}
}
\value{
A nested list whose first level is named by each level of
\code{separate_col}, and whose second level is named by each selected gene-set
collection. A successful element contains:
\describe{
\item{\code{result}}{The GSEA table filtered at 0.05 by
\code{p_threshold_type}, ordered by decreasing absolute NES, with internal
columns 8 through 10 removed by the current implementation.}
\item{\code{gsea}}{The combined enrichment plot for up to \code{topn}
pathways.}
\item{\code{topn}}{The requested maximum pathway count.}
}
A failed group-collection combination is returned as a \code{gsea_error}
object containing an \code{error} message rather than stopping all analyses.
}
\description{
Build a ranked gene vector for each level of a grouping column, load human or
mouse MSigDB GMT collections, run clusterProfiler GSEA for every
group-collection combination in parallel, and create enrichment plots for the
strongest significant pathways.
}
\details{
Within each separate group, ranking values are named by \code{feature_col} and
sorted decreasingly. \code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA()}}
is run with \code{pvalueCutoff = 1}; results are subsequently filtered at 0.05
using \code{p_threshold_type} and ordered by decreasing absolute normalized
enrichment score. Ranked vectors should contain finite numeric values and
appropriate gene symbols for the selected GMT files.

The function expects MSigDB 2025.1 symbol GMT files below
\code{human/v2025.1} or \code{mouse/v2025.1}. Selecting \code{"GO"} includes
GO BP, CC, and MF. Human collections use c5 GO, KEGG Medicus, Reactome,
Hallmark, C2 all, and WikiPathways files. Mouse collections use their m5, m2,
mh, and WikiPathways counterparts.

The internal all-C2 collection is currently named \code{"CM2All"}; therefore,
the default label \code{"C2All"} does not match that collection under the
current pattern selection. A pattern matching \code{"CM2All"}, such as
\code{"M2All"}, selects it.

Analyses use \code{future::multisession} with up to 15 workers and seeded
furrr execution. The function sets \code{future.globals.maxSize} to 10 GiB and
changes the future plan to sequential at the end; prior option and plan values
are not restored.
}
\examples{
\dontrun{
gsea_input <- data.frame(
  gene = rep(paste0("Gene", seq_len(1000)), 2),
  rank = rnorm(2000),
  target = rep(c("TargetA", "TargetB"), each = 1000)
)

gsea_result <- gsea_analysis(
  feature_value_df = gsea_input,
  feature_col = "gene",
  value_col = "rank",
  separate_col = "target",
  background_gene_sets = c("GO", "KEGG", "Hallmark"),
  species = 9606,
  p_threshold_type = "p.adjust",
  topn = 5,
  gmt_MSigDB_dir = "/reference/functional/gene_sets/MSigDB"
)

gsea_result$TargetA$KEGG$result
gsea_result$TargetA$KEGG$gsea
}
}
\seealso{
\code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA()}},
\code{\link[clusterProfiler:read-gmt]{clusterProfiler::read.gmt()}},
\code{\link[enrichplot:gseaplot2]{enrichplot::gseaplot2()}},
\code{\link[furrr:future_imap]{furrr::future_imap()}}
}
