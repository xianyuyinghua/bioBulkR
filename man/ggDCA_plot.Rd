\name{ggDCA_plot}
\alias{ggDCA_plot}
\title{Perform decision curve analysis for a logistic nomogram and predictors}
\usage{
ggDCA_plot(
  input_data = df_train,
  group_colname = "Group",
  key_features = key_genes,
  plot_type = "all"
)
}
\arguments{
\item{input_data}{A data frame containing the binary outcome and all predictors
named in \code{key_features}.}

\item{group_colname}{Name of the binary outcome column. The current
implementation converts it with \code{as.numeric(x) - 1}, so a two-level
factor is expected. Its factor-level order determines which class is coded
zero or one.}

\item{key_features}{A character vector naming predictors used for the combined
\code{Nomogram} model and the individual predictor models.}

\item{plot_type}{Models to retain in the decision curve: \code{"all"},
\code{"Nomogram"}, one or more feature names, or a combination of these.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{The decision-curve object returned by
\code{\link[ggDCA:dca]{ggDCA::dca()}}.}
\item{\code{plot}}{A ggplot showing threshold probability against net benefit
for the selected models.}
}
}
\description{
Fit an rms logistic model containing all selected predictors plus a
single-predictor model for each feature, calculate decision curves with ggDCA,
and draw net benefit across threshold probabilities.
}
\details{
Predictor names containing characters other than letters, numbers, and
underscores are cleaned before formulas are built, while model labels retain
the original feature names. An \code{\link[rms:lrm]{rms::lrm()}} model is
fitted for all selected features and for each feature separately. These models
are passed to \code{\link[ggDCA:dca]{ggDCA::dca()}} with names consisting of
\code{Nomogram} followed by the original feature names.

If \code{plot_type} contains feature names, \code{key_features} is first
restricted to their intersection. Consequently, the \code{Nomogram} model is
also rebuilt from that restricted feature set. Use \code{plot_type = "all"}
to retain every model. Unknown requested model names cause the current
implementation to call \code{quit(save = "no", runLast = FALSE)}.

The outcome conversion expects a factor. Supplying an already numeric zero/one
outcome does not preserve those values under the current conversion. The
function also assigns the processed input data to \code{input_data_global} in
the global environment with \code{<<-}; this object is not removed on return.

The plot title is fixed as
\code{Decision Curve Analysis (Disease Prediction)} and colors are obtained
from \code{basicR::get_colors()}.
}
\examples{
\dontrun{
dca_data <- data.frame(
  Group = factor(rep(c("Control", "Disease"), each = 50)),
  GeneA = rnorm(100),
  GeneB = rnorm(100)
)

dca_result <- ggDCA_plot(
  input_data = dca_data,
  group_colname = "Group",
  key_features = c("GeneA", "GeneB"),
  plot_type = c("Nomogram", "GeneA", "GeneB")
)

dca_result$plot
dca_result$result
}
}
\seealso{
\code{\link[ggDCA:dca]{ggDCA::dca()}},
\code{\link[rms:lrm]{rms::lrm()}}
}
