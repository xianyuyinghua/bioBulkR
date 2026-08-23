\name{pROC_plot}
\alias{pROC_plot}
\title{Plot ROC curves and confidence bands for candidate features}
\usage{
pROC_plot(
  data,
  group_col,
  genes,
  data_type = "Train",
  legend_position = c(0.75, 0.2),
  legend_text_size = 12,
  digits = 3
)
}
\arguments{
\item{data}{A data frame containing a binary group column and numeric predictor
columns.}

\item{group_col}{Name of the binary outcome or group column.}

\item{genes}{A character vector naming predictor columns for which ROC curves
are calculated.}

\item{data_type}{A character label appended to the plot title, such as
\code{"Train"} or \code{"Test"}. An empty string omits the parenthetical
label.}

\item{legend_position}{Legend position passed to ggplot2, commonly a
two-element numeric coordinate or a position name such as \code{"right"}.}

\item{legend_text_size}{Numeric legend text size.}

\item{digits}{Number of decimal places used for AUC and confidence interval
values in legend labels.}
}
\value{
A ggplot containing one ROC step curve and translucent 95 percent sensitivity
confidence ribbon per feature, plus the chance diagonal. Each legend label
contains the feature name, AUC, and AUC 95 percent confidence interval.
}
\description{
Calculate a pROC curve, AUC, AUC confidence interval, and sensitivity
confidence band for each requested feature, then combine the curves in a
single ggplot with AUC information in the legend.
}
\details{
The group column is converted to a factor. If a group named \code{control},
\code{normal}, \code{norm}, or \code{low} is present, ignoring case, its first
occurrence is used as the first level; otherwise, levels are sorted
alphabetically. The function is intended for exactly two outcome groups.

ROC direction is selected independently for each feature by comparing its mean
in the first and second groups. The selected level order and direction control
the interpretation of sensitivity and specificity.

AUC confidence intervals are obtained with
\code{\link[pROC:ci.auc]{pROC::ci.auc()}}. Sensitivity confidence intervals are
obtained with \code{\link[pROC:ci.se]{pROC::ci.se()}} over specificities from
zero to one in increments of 0.01, using 200 bootstrap replicates. Curve
coordinates are converted to false-positive rate, ordered, and made
monotonically non-decreasing with \code{cummax()}. Explicit endpoints
\code{(0, 0)} and \code{(1, 1)} are added before plotting.

The function returns only the plot; ROC objects and intermediate confidence
tables are not included in the return value. Confidence calculations can
consume random numbers, so callers can set a seed when reproducibility is
required.
}
\examples{
\dontrun{
roc_data <- data.frame(
  Group = factor(rep(c("Control", "Disease"), each = 30)),
  GeneA = c(rnorm(30, 0), rnorm(30, 1)),
  GeneB = c(rnorm(30, 0), rnorm(30, 0.6))
)

roc_plot <- pROC_plot(
  data = roc_data,
  group_col = "Group",
  genes = c("GeneA", "GeneB"),
  data_type = "Train",
  digits = 3
)

roc_plot
}
}
\seealso{
\code{\link[pROC:roc]{pROC::roc()}},
\code{\link[pROC:auc]{pROC::auc()}},
\code{\link[pROC:ci.auc]{pROC::ci.auc()}},
\code{\link[pROC:ci.se]{pROC::ci.se()}}
}
