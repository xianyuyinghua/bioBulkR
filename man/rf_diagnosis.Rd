\name{rf_diagnosis}
\alias{rf_diagnosis}
\title{Tune and fit a random-forest diagnostic classifier}
\usage{
rf_diagnosis(
  data = df_train,
  feature_colnames = gene_name,
  response_colnames = "Group",
  n_folds = nfolds.num,
  seed = 1
)
}
\arguments{
\item{data}{A data frame containing predictor columns and one response column.}

\item{feature_colnames}{A character vector naming predictor features used by
the classifier.}

\item{response_colnames}{Name of the categorical response column. Although the
argument name is plural, the current implementation expects one column name.}

\item{n_folds}{Number of cross-validation folds passed to
\code{\link[caret:createDataPartition]{caret::createFolds()}}. Each class
should have enough observations to be represented appropriately across folds.}

\item{seed}{Integer random seed set before fold creation and model fitting.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{A list with \code{optimal_ntree}, the candidate tree count
with the smallest mean cross-validation error, and \code{importance_df}, a
two-column feature-importance table sorted decreasingly by \code{Importance}.}
\item{\code{plot}}{A list with \code{p_error}, the number-of-trees versus
cross-validation error ggplot, and \code{p_importance}, the horizontal feature
importance ggplot.}
}
The fitted random-forest model and fold-level predictions are not returned.
}
\description{
Use cross-validation to select the number of trees for a random-forest
classifier, fit the final model on all samples, calculate feature importance,
and create error and importance plots.
}
\details{
The function retains only \code{feature_colnames} and
\code{response_colnames}, then creates stratified folds from the response.
Candidate forest sizes are \code{seq(1, 500, by = 5)}, which evaluates 1, 6,
and subsequent increments through 496 trees. For every candidate, a model is
fitted on each training fold and its held-out misclassification rate is
averaged. \code{which.min()} selects the first candidate attaining the minimum
error.

The response is converted to a factor, so the final random forest is a
classifier. The final model is fitted on all observations with
\code{importance = TRUE}. The returned importance value is the first column of
\code{\link[randomForest:importance]{randomForest::importance()}}, which for a
classification forest is the first class-specific importance column rather
than necessarily the overall mean-decrease-accuracy or mean-decrease-Gini
column.

Input predictors should be suitable for
\code{\link[randomForest:randomForest]{randomForest::randomForest()}} and should
not contain unhandled missing values. The search fits one forest for every
candidate/fold combination and can therefore be computationally expensive.
}
\examples{
\dontrun{
diagnosis_data <- data.frame(
  GeneA = rnorm(60),
  GeneB = rnorm(60),
  GeneC = rnorm(60),
  Group = factor(rep(c("Control", "Disease"), each = 30))
)

rf_result <- rf_diagnosis(
  data = diagnosis_data,
  feature_colnames = c("GeneA", "GeneB", "GeneC"),
  response_colnames = "Group",
  n_folds = 5,
  seed = 123
)

rf_result$result$optimal_ntree
rf_result$result$importance_df
rf_result$plot$p_error
rf_result$plot$p_importance
}
}
\seealso{
\code{\link[randomForest:randomForest]{randomForest::randomForest()}},
\code{\link[randomForest:importance]{randomForest::importance()}},
\code{\link[caret:createDataPartition]{caret::createFolds()}}
}
