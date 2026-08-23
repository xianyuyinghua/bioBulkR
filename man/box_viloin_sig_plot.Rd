\name{box_viloin_sig_plot}
\alias{box_viloin_sig_plot}
\title{Draw grouped box or violin plots with pairwise significance annotations}
\usage{
box_viloin_sig_plot(
  data_long,
  type_col = "Gene",
  group_col = "Group",
  value_col = "Value",
  test_method = "wilcox",
  type_level = NA,
  group_level = NA,
  x_text_angle = 0,
  add_point = FALSE,
  point_size = 1.8,
  point_alpha = 0.7,
  point_width = 0.15,
  group_color = basicR::get_colors(package = "ggsci", name = "jco", number = 5.1),
  plot_type = "box",
  diff_show_type = "signif",
  diff_show_name = "p.value",
  sig_bracket_vjust = 0.7,
  step_increase = 0,
  hide_ns = TRUE,
  bracket_size = 0.8,
  bracket_y = 0,
  axis_y_title = "Expression",
  axis_x_title = "Gene",
  legend_title = "Group",
  title_name = NULL
)
}
\arguments{
\item{data_long}{A long-format data frame containing feature, group, and
numeric value columns.}
\item{type_col}{Name of the feature column. Set to \code{NA} for data
containing only groups and values, without a feature dimension.}
\item{group_col}{Name of the group column.}
\item{value_col}{Name of the numeric measurement column.}
\item{test_method}{Pairwise test. A value beginning with \code{"wilcox"} uses
\code{\link[rstatix:wilcox_test]{rstatix::pairwise_wilcox_test()}};
\code{"t"} uses \code{\link[rstatix:t_test]{rstatix::pairwise_t_test()}}.}
\item{type_level}{Optional feature order. \code{NA} uses observed factor
order. Every supplied level must occur in the data.}
\item{group_level}{Optional group order. \code{NA} places \code{Normal} or
\code{Control} first when present, then keeps the remaining observed order.}
\item{x_text_angle}{Numeric angle for x-axis labels.}
\item{add_point}{Logical; overlay jittered observations.}
\item{point_size, point_alpha, point_width}{Point size, transparency, and
jitter width.}
\item{group_color}{A character vector of colors for groups. A manual scale is
added when more than one color is supplied.}
\item{plot_type}{Plot family. Values beginning with \code{"box"} select a box
plot; values beginning with \code{"vi"} select a violin plot.}
\item{diff_show_type}{Annotation form containing \code{"sig"} for significance
stars or \code{"value"} for numeric p-values.}
\item{diff_show_name}{Significance source containing \code{"value"} for raw
p-values or \code{"adj"} for Bonferroni-adjusted p-values.}
\item{sig_bracket_vjust}{Vertical justification of annotation text.}
\item{step_increase}{Incremental vertical spacing between comparison
brackets.}
\item{hide_ns}{Logical; hide non-significant comparisons.}
\item{bracket_size}{Bracket line width.}
\item{bracket_y}{Vertical bracket nudge.}
\item{axis_y_title, axis_x_title, legend_title}{Axis and color-legend titles.}
\item{title_name}{Optional plot title.}
}
\value{
A named list containing:
\describe{
\item{\code{plot}}{The box or violin ggplot with pairwise significance
annotations.}
\item{\code{stat}}{The pairwise-test table, including feature names, compared
groups, sample sizes, raw and Bonferroni-adjusted p-values, significance
labels, and x/y plotting coordinates.}
}
}
\description{
Standardize a long-format measurement table, perform all pairwise group
comparisons within each feature, and draw a box or violin plot with optional
sample points and raw or adjusted p-value annotations.
}
\details{
Input columns are renamed internally to \code{Type}, \code{Group}, and
\code{Value}. Existing columns with those standard names can be renamed with a
\code{_raw} suffix to avoid collisions. All pairwise comparisons among factor
levels of \code{Group} are tested separately for each \code{Type}, with
Bonferroni adjustment.

Features whose values are all zero are excluded from testing and receive
non-significant placeholder rows. For certain cases where at least two groups
are constant, the implementation adds 0.0001 to one value before testing. Raw
significance labels use \code{***}, \code{**}, and \code{*} for p-values below
0.001, 0.01, and 0.05. Four-star adjusted labels returned by rstatix are
reduced to three stars.

When \code{type_col = NA}, groups are drawn on the x-axis and the legend is
hidden. Otherwise, feature types are placed on the x-axis and groups control
color and dodging. \code{add_point = TRUE} suppresses boxplot outliers and
overlays individual observations.
}
\examples{
\dontrun{
plot_data <- data.frame(
  Gene = rep(c("GeneA", "GeneB"), each = 12),
  Group = rep(rep(c("Control", "Disease"), each = 6), 2),
  Value = c(
    rnorm(6, 5), rnorm(6, 7),
    rnorm(6, 4), rnorm(6, 4.5)
  )
)

result <- box_viloin_sig_plot(
  data_long = plot_data,
  group_level = c("Control", "Disease"),
  test_method = "wilcox",
  plot_type = "box",
  diff_show_type = "signif",
  diff_show_name = "p.adj",
  add_point = TRUE
)

result$plot
result$stat
}
}
\seealso{
\code{\link[ggpubr:ggboxplot]{ggpubr::ggboxplot()}},
\code{\link[ggpubr:ggviolin]{ggpubr::ggviolin()}},
\code{\link[ggpubr:stat_pvalue_manual]{ggpubr::stat_pvalue_manual()}},
\code{\link[rstatix:wilcox_test]{rstatix::pairwise_wilcox_test()}}
}
