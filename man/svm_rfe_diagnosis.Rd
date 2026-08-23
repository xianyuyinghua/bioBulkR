\name{svm_rfe_diagnosis}
\alias{svm_rfe_diagnosis}
\title{Perform SVM recursive feature elimination for diagnosis}
\usage{
svm_rfe_diagnosis(
  data,
  feature_colnames,
  response_colnames,
  n_folds = 5,
  seed = 1,
  method = "svmRadial",
  n_cores = 4,
  feature_colors = basicR::get_colors(number = 100.1)
)
}
\arguments{
\item{data}{A data frame containing predictor features and one categorical
response column. For the current subset-size calculation, it should contain
only \code{feature_colnames} plus \code{response_colnames}.}

\item{feature_colnames}{A character vector naming predictor columns.}

\item{response_colnames}{Name of the categorical response column. Although the
argument name is plural, the current implementation expects one column name.
The feature-importance summary expects a binary response.}

\item{n_folds}{Number of cross-validation folds used by
\code{\link[caret:rfeControl]{caret::rfeControl()}}. Each response class should
have enough observations for the requested number of folds.}

\item{seed}{Integer seed used for RFE and parallel resampling seeds.}

\item{method}{Caret model identifier passed to
\code{\link[caret:rfe]{caret::rfe()}}, such as the default
\code{"svmRadial"}. The required model package must be installed.}

\item{n_cores}{Number of PSOCK workers used through doParallel.}

\item{feature_colors}{A color vector associated with
\code{feature_colnames} in the accuracy plot.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{The caret RFE object. It is augmented with
\code{rfe_variables}, a feature table containing first-class, second-class,
and overall mean importance, and with a \code{Feature} column in its
\code{results} table. Standard caret fields such as \code{optVariables} remain
available.}
\item{\code{plot}}{A ggplot of cross-validated accuracy against the number of
selected features, with points colored by the associated feature.}
}
}
\description{
Use caret recursive feature elimination with cross-validated SVM models,
summarize class-specific feature importance, and plot accuracy against the
number of selected features.
}
\details{
The response is converted to a factor and the requested predictors are passed
to caret RFE with \code{caretFuncs}, cross-validation, and parallel execution.
Candidate subset sizes are currently computed as
\code{c(1:ncol(data) - 1)}, producing zero through
\code{ncol(data) - 1}. This assumes that \code{data} consists only of the
requested feature columns and the response; additional columns can produce
subset sizes larger than the predictor count.

Control-like response labels (\code{control}, \code{normal}, \code{norm}, or
\code{low}, ignoring case) are placed first; otherwise labels are sorted
alphabetically. The RFE variable table is summarized into mean importance for
the first class, second class, and overall importance, then sorted decreasingly
by overall importance.

For more than one worker, the function constructs a caret seed list so that
parallel resampling is reproducible. A PSOCK cluster is created and registered
with doParallel. On successful completion it is stopped and the sequential
backend is restored. If RFE raises an error before cleanup, the current
implementation does not use an on-exit handler to guarantee cluster shutdown.
}
\examples{
\dontrun{
svm_data <- data.frame(
  GeneA = rnorm(60),
  GeneB = rnorm(60),
  GeneC = rnorm(60),
  Group = factor(rep(c("Control", "Disease"), each = 30))
)

svm_result <- svm_rfe_diagnosis(
  data = svm_data,
  feature_colnames = c("GeneA", "GeneB", "GeneC"),
  response_colnames = "Group",
  n_folds = 5,
  seed = 123,
  method = "svmRadial",
  n_cores = 2
)

svm_result$result$optVariables
svm_result$result$rfe_variables
svm_result$plot
}
}
\seealso{
\code{\link[caret:rfe]{caret::rfe()}},
\code{\link[caret:rfeControl]{caret::rfeControl()}},
\code{\link[doParallel:registerDoParallel]{doParallel::registerDoParallel()}}
}
