\name{boruta_diagnosis}
\alias{boruta_diagnosis}
\title{Select diagnostic features with the Boruta algorithm}
\usage{
boruta_diagnosis(
  data = df_train,
  feature_colnames = gene_name,
  response_colnames = "Group",
  pvalue = 0.01,
  seed = 1
)
}
\arguments{
\item{data}{A data frame containing predictor columns and one categorical
response column.}

\item{feature_colnames}{A character vector naming predictor features supplied
to Boruta.}

\item{response_colnames}{Name of the categorical response column. Although the
argument name is plural, the current implementation expects one name.}

\item{pvalue}{Numeric confidence-level threshold passed to the \code{pValue}
argument of \code{\link[Boruta:Boruta]{Boruta::Boruta()}}.}

\item{seed}{Integer random seed set before running Boruta.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{A long-format data frame containing \code{Iteration},
\code{Variable}, \code{Importance}, and \code{Decision}. Each row represents
one available variable-importance value from one Boruta iteration.}
\item{\code{plot}}{A ggplot boxplot of importance Z-scores for predictors and
shadow attributes, colored by Boruta decision.}
}
}
\description{
Run Boruta feature selection for a categorical response, reshape the importance
history into a long-format table, and visualize importance distributions by
Boruta decision.
}
\details{
The input is restricted to the response and requested feature columns. The
response is converted to a factor whose levels follow their order of first
appearance. Boruta is run with Bonferroni adjustment, at most 100 runs, stored
importance history, and \code{\link[Boruta:getImpRfZ]{Boruta::getImpRfZ()}}
importance scores.

Importance history is converted to long format, missing and non-finite values
are removed, and original Boruta decisions are joined to the variables. Shadow
attributes are labeled \code{"Shadow"}. Variables are ordered by median
importance, with confirmed variables placed after all other decision groups.

The function calls
\code{\link[Boruta:TentativeRoughFix]{Boruta::TentativeRoughFix()}} and
\code{\link[Boruta:getSelectedAttributes]{Boruta::getSelectedAttributes()}},
but the current return table and plot use decisions from the original Boruta
fit. Consequently, tentatively resolved decisions and the vector of confirmed
features are not returned separately.
}
\examples{
\dontrun{
diagnosis_data <- data.frame(
  GeneA = rnorm(60),
  GeneB = rnorm(60),
  GeneC = rnorm(60),
  Group = factor(rep(c("Control", "Disease"), each = 30))
)

boruta_result <- boruta_diagnosis(
  data = diagnosis_data,
  feature_colnames = c("GeneA", "GeneB", "GeneC"),
  response_colnames = "Group",
  pvalue = 0.01,
  seed = 123
)

head(boruta_result$result)
boruta_result$plot
}
}
\seealso{
\code{\link[Boruta:Boruta]{Boruta::Boruta()}},
\code{\link[Boruta:TentativeRoughFix]{Boruta::TentativeRoughFix()}},
\code{\link[Boruta:getSelectedAttributes]{Boruta::getSelectedAttributes()}}
}
