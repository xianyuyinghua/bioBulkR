\name{Volcano_Heatmap_plot}
\alias{Volcano_Heatmap_plot}
\title{Create volcano, rectangular heatmap, and circular heatmap plots}
\usage{
Volcano_Heatmap_plot(
  Res = Res,
  expr = data$expr,
  condition = data$condition,
  sample_color = colors,
  logFC_name = "logFC",
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05,
  top_gene = 10,
  heatmap_scale_range = c(-2, 2),
  group_order = data$group_order,
  group_order_color = c("#FF7F00", "#1177b0"),
  change = c("Up", "Down", "Not"),
  change_colors = c("orange", "#1177b0", "#7b7c7d"),
  heatmap_value_color = viridis(10),
  prominent_genes = NULL,
  prominent_genes_color = c("red"),
  prominent_genes_size = 2,
  track_height = 0.3,
  dend_track_height = 0.2
)
}
\arguments{
\item{Res}{A differential-expression data frame. Row names must contain gene
symbols, and columns must include \code{logFC_name}, \code{p.value}, and
\code{p.adj}.}

\item{expr}{A numeric gene-by-sample expression matrix or data frame whose row
names match the gene symbols in \code{Res}.}

\item{condition}{A sample annotation data frame containing a \code{group} or
\code{groups} column, ignoring case. Sample names may be in row names or a
\code{sample} or \code{samples} column.}

\item{sample_color}{A character vector of sample colors. This argument is
retained by the interface but is not used in the current plotting code.}

\item{logFC_name}{The name of the fold-change column in \code{Res}.}

\item{logFC}{A non-negative absolute fold-change threshold on the log2 scale.}

\item{p_value}{The significance column used for classification. Use
\code{"p.adj"} for adjusted p-values; any other value selects
\code{"p.value"}.}

\item{p}{A numeric significance threshold.}

\item{top_gene}{The number of top upregulated and top downregulated genes used
for labels and heatmaps.}

\item{heatmap_scale_range}{A numeric vector of length two giving lower and
upper truncation limits for the scaled rectangular heatmap matrix.}

\item{group_order}{A character vector specifying the sample-group order and
the group names displayed in plot annotations.}

\item{group_order_color}{A color vector corresponding to
\code{group_order}.}

\item{change}{Three labels used, in order, for upregulated, downregulated, and
non-significant genes.}

\item{change_colors}{Three colors corresponding to \code{change}.}

\item{heatmap_value_color}{A color vector used by the rectangular and circular
heatmaps.}

\item{prominent_genes}{An optional character vector of genes to emphasize in
the volcano plot and include in the heatmaps.}

\item{prominent_genes_color}{Colors for the regulation classes represented by
\code{prominent_genes}.}

\item{prominent_genes_size}{Numeric point-size control for emphasized genes;
their text size is twice this value.}

\item{track_height}{Numeric height of the circular heatmap track.}

\item{dend_track_height}{Numeric height of the circular heatmap dendrogram
track.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{The updated differential-expression table, including
\code{sig}, \code{Symbols}, and gene-label columns.}
\item{\code{plot}}{A list containing \code{valcano}, the volcano
\code{ggplot}; \code{rHeatmap}, the rectangular ComplexHeatmap object; and
\code{cHeatmap}, a recorded circular heatmap plot. The historical name
\code{valcano} is retained by the function.}
}
}
\description{
Reclassify a differential-expression result using specified fold-change and
p-value thresholds, label top genes, and create a volcano plot plus
rectangular and circular heatmaps for the selected genes.
}
\details{
Genes are classified using strict p-value and absolute fold-change
thresholds. Within each direction, genes are ordered first by fold change and
then by the selected p-value. Up to \code{top_gene} genes from each direction
are labelled and included in the heatmaps.

The rectangular heatmap scales each selected gene across samples and clips
values to \code{heatmap_scale_range}. The circular heatmap clusters the
unscaled selected expression matrix, draws directly on the active graphics
device, and is captured with \code{\link[grDevices:recordplot]{grDevices::recordPlot()}}.

When \code{prominent_genes} is supplied, every requested gene must occur in
\code{Res}. The number of \code{prominent_genes_color} values must equal the
number of distinct regulation classes represented by those genes; otherwise,
the function prints a message and stops.
}
\examples{
\dontrun{
plots <- Volcano_Heatmap_plot(
  Res = deg_result,
  expr = expression_matrix,
  condition = sample_condition,
  group_order = c("Control", "Disease"),
  group_order_color = c("#FF7F00", "#1177B0"),
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05,
  top_gene = 10
)

plots$plot$valcano
ComplexHeatmap::draw(plots$plot$rHeatmap)
replayPlot(plots$plot$cHeatmap)
}
}
\seealso{
\code{\link[=Array_DEGs]{Array_DEGs()}},
\code{\link[ggplot2:ggplot]{ggplot2::ggplot()}},
\code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}},
\code{\link[circlize:circos.heatmap]{circlize::circos.heatmap()}}
}
