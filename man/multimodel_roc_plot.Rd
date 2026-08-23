\name{multimodel_roc_plot}
\alias{multimodel_roc_plot}
\title{Plot pooled ROC curves for multiple models or features}
\usage{
multimodel_roc_plot(
  pred_all,
  features_col = "AlgorithmName",
  observation = "obs",
  probability = "prob",
  pos_class,
  digits = 3,
  sp_grid = seq(0, 1, 0.01),
  title_text = "ROC Curves with 95\% CI",
  legend_position = "right",
  base_size = 16,
  colors = basicR::get_colors(number = 50.1)
)
}
\arguments{
\item{pred_all}{A data frame of pooled predictions, typically out-of-fold
cross-validation predictions, containing model, observed-class, and
positive-class probability columns.}

\item{features_col}{Name of the model or feature identifier column.}

\item{observation}{Name of the observed binary class column.}

\item{probability}{Name of the numeric predicted probability column for the
positive class. Values are normally between zero and one.}

\item{pos_class}{A value in \code{observation} identifying the positive class.}

\item{digits}{Number of decimal places used for AUC and confidence interval
values in legend labels.}

\item{sp_grid}{A numeric specificity grid passed to
\code{\link[pROC:ci.se]{pROC::ci.se()}}. The default evaluates values from zero
to one in increments of 0.01.}

\item{title_text}{Plot title.}

\item{legend_position}{Legend position passed to ggplot2, such as
\code{"right"} or a two-element numeric coordinate.}

\item{base_size}{Base theme text size.}

\item{colors}{A color vector for model curves and confidence ribbons. Supply
at least as many colors as distinct models.}
}
\value{
A ggplot containing one step ROC curve and translucent 95 percent sensitivity
confidence ribbon per model, plus a chance diagonal. Each legend label includes
the model name, AUC, and AUC 95 percent confidence interval. Intermediate ROC
objects and confidence tables are not returned.
}
\description{
Calculate an ROC curve, AUC, AUC confidence interval, and sensitivity
confidence ribbon for each model represented in a prediction table, then
combine all model curves in a single ggplot.
}
\details{
The function first verifies that \code{features_col}, \code{observation}, and
\code{probability} occur in \code{pred_all}. The observed class is converted to
a factor. The negative class is the first factor level other than
\code{pos_class}; the function is therefore intended for binary outcomes.
Each model is evaluated with levels ordered as negative, positive and
\code{direction = "<"}, meaning larger probabilities indicate the positive
class.

AUC confidence intervals are calculated with
\code{\link[pROC:ci.auc]{pROC::ci.auc()}} and sensitivity intervals with
\code{\link[pROC:ci.se]{pROC::ci.se()}}, both receiving \code{boot.n = 200}.
False-positive-rate coordinates are sorted, sensitivities are made
monotonically non-decreasing with \code{cummax()}, and explicit endpoints
\code{(0, 0)} and \code{(1, 1)} are added.

Model order follows first appearance in \code{pred_all}, and that same order is
used for legend labels and colors. Confidence calculations can consume random
numbers, so callers can set a seed when reproducibility is required.
}
\examples{
\dontrun{
predictions <- data.frame(
  AlgorithmName = rep(c("Logistic", "RandomForest"), each = 100),
  obs = rep(rep(c("Control", "Disease"), each = 50), 2),
  prob = c(
    runif(50, 0, 0.5), runif(50, 0.5, 1),
    runif(50, 0, 0.6), runif(50, 0.4, 1)
  )
)

roc_plot <- multimodel_roc_plot(
  pred_all = predictions,
  features_col = "AlgorithmName",
  observation = "obs",
  probability = "prob",
  pos_class = "Disease"
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
