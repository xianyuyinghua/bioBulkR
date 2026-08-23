\name{xgboost_diagnosis}
\alias{xgboost_diagnosis}
\title{Fit an XGBoost diagnostic classifier and rank features}
\usage{
xgboost_diagnosis(
  data = df_train,
  feature_colnames,
  response_colnames,
  n_folds = 5,
  seed = 1
)
}
\arguments{
\item{data}{A data frame containing predictor columns and one categorical
response column.}

\item{feature_colnames}{A character vector naming predictor features used by
XGBoost.}

\item{response_colnames}{Name of the categorical response column. Although the
argument name is plural, the current implementation expects one name.}

\item{n_folds}{Number of cross-validation folds passed to
\code{\link[xgboost:xgb.cv]{xgboost::xgb.cv()}}.}

\item{seed}{Integer random seed set before cross-validation and model fitting.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{An XGBoost feature-importance table sorted by decreasing
\code{Gain}.}
\item{\code{plot}}{A horizontal ggplot bar chart of feature gain.}
\item{\code{best_eta}}{The learning rate selected from the candidate grid.}
\item{\code{best_cv_logloss}}{The smallest recorded cross-validated test log
loss for the selected learning rate.}
}
}
\description{
Select a learning rate using cross-validated binary log loss, fit a final
XGBoost classifier, and create a feature-gain importance plot.
}
\details{
Response levels are ordered with \code{group_new_levels()} and converted to
binary labels: the first level becomes 0 and every other level becomes 1. The
function is intended for binary classification. Predictor data are converted
to a double matrix, and missing or other non-finite matrix values are replaced
by zero.

Four learning rates (\code{0.001}, \code{0.01}, \code{0.1}, and \code{0.3})
are evaluated with 30-round cross-validation. For each rate, the minimum test
log loss across rounds is used as its score. Other tuning parameters are fixed:
maximum tree depth 6, minimum child weight 1, subsample 0.5, column subsample 1,
and gamma 0. The final binary logistic model uses the selected learning rate
and 200 boosting rounds.

Feature importance is calculated with
\code{\link[xgboost:xgb.importance]{xgboost::xgb.importance()}} and sorted by
\code{Gain}. Periods in the feature names prepared for the importance table are
displayed as hyphens. Predictor columns should therefore be numeric and their
model-matrix representation should correspond to the matrix used to train the
model.
}
\examples{
\dontrun{
diagnosis_data <- data.frame(
  GeneA = rnorm(80),
  GeneB = rnorm(80),
  GeneC = rnorm(80),
  Group = factor(rep(c("Control", "Disease"), each = 40))
)

xgb_result <- xgboost_diagnosis(
  data = diagnosis_data,
  feature_colnames = c("GeneA", "GeneB", "GeneC"),
  response_colnames = "Group",
  n_folds = 5,
  seed = 123
)

xgb_result$best_eta
xgb_result$best_cv_logloss
xgb_result$result
xgb_result$plot
}
}
\seealso{
\code{\link[xgboost:xgb.cv]{xgboost::xgb.cv()}},
\code{\link[xgboost:xgb.train]{xgboost::xgb.train()}},
\code{\link[xgboost:xgb.importance]{xgboost::xgb.importance()}}
}
