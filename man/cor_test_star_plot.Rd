\name{cor_test_star_plot}
\alias{cor_test_star_plot}
\title{Calculate feature correlations and optionally create significance plots}
\usage{
cor_test_star_plot(
  cor_df,
  common_features,
  separate_features,
  separate_features_levles,
  cor_method = "spearman",
  plot = FALSE,
  plot_title = NULL,
  axis_x_angle = 45,
  facet_ncol = NA,
  scatter_position = NULL,
  heatmap_sig_height = 0.015
)
}
\arguments{
\item{cor_df}{A data frame containing observations in rows and all requested
numeric features in columns. Row names are used to align reshaped data.}

\item{common_features}{A character vector of features each compared with every
member of \code{separate_features}.}

\item{separate_features}{A character vector defining the other side of each
pairwise correlation.}

\item{separate_features_levles}{Factor levels controlling the order of separate
features in matrices, plot lists, and facets. The historical argument spelling
\code{levles} is retained.}

\item{cor_method}{The correlation method passed to
\code{\link[stats:cor.test]{stats::cor.test()}}: \code{"spearman"},
\code{"kendall"}, or \code{"pearson"}.}

\item{plot}{Logical; create plots in addition to result tables.}

\item{plot_title}{Optional title for lollipop plots. \code{NULL} uses
\code{"Correlation"}.}

\item{axis_x_angle}{Numeric angle for common-feature labels on lollipop
x-axes. Values of 0 and 90 receive specialized justification; other values use
right-aligned labels.}

\item{facet_ncol}{Number of columns in the faceted lollipop plot. \code{NA}
calculates it with \code{\link[=facet_nrow_ncol]{facet_nrow_ncol()}}.}

\item{scatter_position}{Optional numeric vector \code{c(x, y)} giving the
correlation-statistic annotation position in the one-to-one scatter plot. If
\code{NULL}, the position is calculated from the feature ranges.}

\item{heatmap_sig_height}{Numeric vertical offset, in normalized parent
coordinates, separating correlation values and significance stars in heatmap
cells.}
}
\value{
If \code{plot = FALSE}, a named list containing:
\describe{
\item{\code{cor_pvalue_sig_data}}{A long data frame with common and separate
feature names, correlation estimates, p-values, and significance labels.}
\item{\code{cor_matrix}}{A separate-feature by common-feature numeric
correlation matrix.}
\item{\code{pvalue_matrix}}{A matrix of corresponding p-values.}
\item{\code{significance_matrix}}{A character matrix of significance stars.}
}

If \code{plot = TRUE}, these objects are nested under \code{data}; a second
element, \code{plot}, contains:
\describe{
\item{\code{lollipop_list}}{A named list of separate lollipop ggplots.}
\item{\code{lollipop_facet}}{A faceted lollipop ggplot.}
\item{\code{heatmap}}{A drawn ComplexHeatmap object for one-to-many or
many-to-many input, otherwise \code{NULL}.}
\item{\code{scatter}}{A scatter plot with marginal densities for one-to-one
input, otherwise \code{NULL}.}
}
}
\description{
Test every pair formed by common and separate feature sets, attach
significance stars, reshape the statistics into matrices, and optionally
create lollipop, faceted, scatter, or heatmap visualizations.
}
\details{
Every combination of \code{common_features} and \code{separate_features} is
tested independently with \code{\link[stats:cor.test]{stats::cor.test()}}.
Requested feature names must all occur in \code{cor_df}; otherwise, the
function prints the missing names and stops. Input features should be numeric
and contain enough usable observations for the selected test.

Missing correlation estimates are replaced by zero and missing p-values by
one. Significance labels are \code{***} for p below 0.001, \code{**} below
0.01, \code{*} below 0.05, and an empty string otherwise. No multiple-testing
adjustment is performed.

With plotting enabled, individual and faceted lollipop plots are prepared. A
single common feature paired with a single separate feature additionally
produces a linear-fit scatter plot with marginal density distributions. A
one-to-many or many-to-many input instead produces a clustered ComplexHeatmap
showing correlation estimates and significance stars. The heatmap legend is
labelled \code{Spearman correlation} by the current implementation regardless
of \code{cor_method}.
}
\examples{
\dontrun{
result <- cor_test_star_plot(
  cor_df = iris,
  common_features = c("Sepal.Length", "Sepal.Width"),
  separate_features = c("Petal.Length", "Petal.Width"),
  separate_features_levles = c("Petal.Length", "Petal.Width"),
  cor_method = "spearman",
  plot = FALSE
)

result$cor_matrix
result$pvalue_matrix
result$significance_matrix
}
}
\seealso{
\code{\link[stats:cor.test]{stats::cor.test()}},
\code{\link[=facet_nrow_ncol]{facet_nrow_ncol()}},
\code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}},
\code{\link[ggExtra:ggMarginal]{ggExtra::ggMarginal()}}
}
