\name{Venn_plot}
\alias{Venn_plot}
\title{Draw a Venn diagram and return genes shared by all sets}
\usage{
Venn_plot(
  target_genes = list(data$DEGs, data$target_gene),
  target_genes_names = c("DEGs", data$target_gene_name),
  fill_colors = c("#F89FA8", "#F9E9A4"),
  base_size = 2,
  label_size = 2,
  sigdigs = 3,
  digits = NULL,
  print_mode = c("percent", "raw"),
  font = "Arial"
)
}
\arguments{
\item{target_genes}{A list of vectors containing the genes or other values in
each set. At least two sets are expected.}

\item{target_genes_names}{A character vector of display names corresponding to
the elements of \code{target_genes}.}

\item{fill_colors}{A character vector of fill colors. It must contain at least
one color for every element of \code{target_genes}.}

\item{base_size}{Numeric size passed to
\code{\link[VennDiagram:venn.diagram]{VennDiagram::venn.diagram()}} for
set-category labels.}

\item{label_size}{Numeric size used for region count or percentage labels.}

\item{sigdigs}{Numeric control retained by the function interface. The current
implementation uses two significant digits internally.}

\item{digits}{\code{NULL}, or a non-negative integer specifying the decimal
places used when recalculating percentage labels. This is applied only when
\code{print_mode} contains \code{"percent"}.}

\item{print_mode}{A character vector passed to
\code{\link[VennDiagram:venn.diagram]{VennDiagram::venn.diagram()}} to select
region-label content, commonly \code{"raw"}, \code{"percent"}, or
\code{c("percent", "raw")}.}

\item{font}{A font-family name used for region and set-category labels.}
}
\value{
A named list containing:
\describe{
\item{\code{common}}{The intersection of all vectors in
\code{target_genes}.}
\item{\code{plot}}{A grid grob containing the Venn diagram. It can be drawn
with \code{\link[grid:grid.draw]{grid::grid.draw()}}.}
}
}
\description{
Create a Venn diagram from a named collection of gene vectors, optionally
format region labels as percentages and counts, and calculate the
intersection shared by every supplied gene set.
}
\details{
Diagrams are drawn without proportional scaling. For two sets, category
labels are positioned at 180 degrees; diagrams with more sets use the
defaults from VennDiagram. When percentage output and \code{digits} are both
requested, the text grobs are rewritten to display the recalculated
percentage followed by the raw count on a new line. Label line height is set
to 0.7.

The function stops with \code{"colors less"} when fewer colors than gene sets
are supplied. It also sets the \code{VennDiagramLogger} threshold in
\code{futile.logger} to \code{ERROR} to suppress lower-priority VennDiagram
log messages.
}
\examples{
genes <- list(
  DEGs = c("TP53", "EGFR", "MYC"),
  Targets = c("EGFR", "MYC", "AKT1")
)

venn <- Venn_plot(
  target_genes = genes,
  target_genes_names = names(genes),
  fill_colors = c("#F89FA8", "#F9E9A4"),
  print_mode = c("percent", "raw"),
  digits = 1
)

venn$common
grid::grid.draw(venn$plot)
}
\seealso{
\code{\link[VennDiagram:venn.diagram]{VennDiagram::venn.diagram()}},
\code{\link[grid:grid.draw]{grid::grid.draw()}}
}
