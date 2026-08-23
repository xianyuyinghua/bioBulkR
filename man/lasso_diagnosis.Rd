\name{lasso_diagnosis}
\alias{lasso_diagnosis}
\title{Fit a cross-validated LASSO diagnostic model}
\usage{
lasso_diagnosis(
  data = df_train,
  feature_colnames,
  response_colnames,
  n_folds = 5,
  seed = seed,
  alpha = 1,
  lambda_type = "lambda.1se",
  feature_colors = basicR::get_colors(number = 100.1)
)
}
\arguments{
\item{data}{A data frame containing predictor columns and one categorical
response column.}

\item{feature_colnames}{A character vector naming the predictor columns.}

\item{response_colnames}{Name of the response column. Although the argument
name is plural, the current implementation expects one column name.}

\item{n_folds}{Number of folds passed to
\code{\link[glmnet:cv.glmnet]{glmnet::cv.glmnet()}}.}

\item{seed}{Integer random seed used before cross-validation. The default
expression refers to an object named \code{seed}, so this argument should
normally be supplied explicitly.}

\item{alpha}{Elastic-net mixing parameter passed to
\code{\link[glmnet:cv.glmnet]{glmnet::cv.glmnet()}}. The default \code{1}
gives LASSO cross-validation.}

\item{lambda_type}{Character string selecting the cross-validated lambda;
normally \code{"lambda.1se"} or \code{"lambda.min"}.}

\item{feature_colors}{A color vector used for coefficient paths and numbers of
non-zero variables in the two plots.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{A list with \code{cvfit}, the cross-validated glmnet fit;
\code{lasso_gene_df}, a two-column data frame containing selected features and
their non-zero coefficients; and \code{lambda_type}, the selected lambda name.}
\item{\code{plot}}{A list with \code{p_path}, the coefficient-path ggplot, and
\code{p_curve}, the cross-validation deviance ggplot.}
}
}
\description{
Fit a binomial model with cross-validated regularization, extract non-zero
coefficients at a selected lambda, and draw coefficient-path and
cross-validation diagnostic plots.
}
\details{
A model matrix is constructed from the selected response and feature columns,
with its intercept column removed. The response is converted to a factor and a
binomial cross-validated glmnet model is fitted using deviance as the assessment
measure. The selected lambda value is printed together with its logarithm.

\code{alpha} is passed to the cross-validation fit. In the current
implementation the separate full coefficient-path fit calls
\code{\link[glmnet:glmnet]{glmnet::glmnet()}} without an explicit \code{alpha},
so that path uses glmnet's default \code{alpha = 1}. Thus, when \code{alpha} is
not \code{1}, the plotted path and cross-validated fit may represent different
mixing parameters.

The response should contain two classes, predictors should be usable in a
numeric model matrix, and the number of observations should support the
requested folds. The supplied color vector should be long enough for all
coefficient curves and non-zero-count categories displayed in the plots.
}
\examples{
\dontrun{
diagnosis_data <- data.frame(
  GeneA = rnorm(60),
  GeneB = rnorm(60),
  GeneC = rnorm(60),
  Group = factor(rep(c("Control", "Disease"), each = 30))
)

lasso_result <- lasso_diagnosis(
  data = diagnosis_data,
  feature_colnames = c("GeneA", "GeneB", "GeneC"),
  response_colnames = "Group",
  n_folds = 5,
  seed = 123,
  alpha = 1,
  lambda_type = "lambda.1se"
)

lasso_result$result$lasso_gene_df
lasso_result$plot$p_path
lasso_result$plot$p_curve
}
}
\seealso{
\code{\link[glmnet:cv.glmnet]{glmnet::cv.glmnet()}},
\code{\link[glmnet:glmnet]{glmnet::glmnet()}},
\code{\link[glmnet:coef.glmnet]{glmnet::coef.glmnet()}}
}
