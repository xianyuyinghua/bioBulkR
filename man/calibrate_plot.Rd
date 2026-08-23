\name{calibrate_plot}
\alias{calibrate_plot}
\title{Fit and plot a logistic-regression calibration curve}
\usage{
calibrate_plot(
  input_data = df_train,
  group_colname = "Group",
  key_features = key_genes,
  calibrate_method = "boot",
  calibrate_B = 1000,
  number_of_groups = 5,
  line_color = basicR::get_colors(number = 2.1, package = "ggsci", name = "jco")
)
}
\arguments{
\item{input_data}{A data frame containing the binary outcome and all predictors
named in \code{key_features}. The current Hosmer-Lemeshow step also requires a
factor column named \code{Group}.}

\item{group_colname}{Name of the binary outcome column used in the logistic
model formula.}

\item{key_features}{A character vector naming predictors included in the
logistic regression.}

\item{calibrate_method}{Calibration method passed to
\code{\link[rms:calibrate]{rms::calibrate()}}, such as \code{"boot"},
\code{"crossvalidation"}, or \code{".632"}.}

\item{calibrate_B}{Number of bootstrap repetitions or resampling iterations
passed as \code{B} to \code{\link[rms:calibrate]{rms::calibrate()}}.}

\item{number_of_groups}{Number of groups used by
\code{\link[ResourceSelection:hoslem.test]{ResourceSelection::hoslem.test()}}.}

\item{line_color}{A color vector of length at least two for the apparent and
bias-corrected curves.}
}
\value{
A named list containing:
\describe{
\item{\code{result}}{A long data frame with \code{Predicted}, \code{Type}, and
\code{Observed} columns. \code{Type} distinguishes the \code{Apparent} and
\code{Bias_corrected} curves.}
\item{\code{plot}}{The calibration ggplot, including the identity line and an
annotation containing the Hosmer-Lemeshow p-value.}
}
}
\description{
Fit an rms logistic regression from selected predictors, estimate apparent and
bias-corrected calibration, perform a Hosmer-Lemeshow goodness-of-fit test, and
draw the calibration curves against the identity line.
}
\details{
The logistic model is fitted with \code{\link[rms:lrm]{rms::lrm()}} using
\code{group_colname ~ key_features}. Calibration output columns \code{predy},
\code{calibrated.orig}, and \code{calibrated.corrected} are converted to a
tidy table with \code{Apparent} and \code{Bias_corrected} curve types.

For the Hosmer-Lemeshow test, the current implementation specifically converts
\code{input_data$Group} with \code{as.numeric(Group) - 1}, regardless of
\code{group_colname}. Therefore, normal use requires a two-level factor named
\code{Group}; its factor-level order determines which class is coded zero or
one. Supplying an already numeric zero/one column does not preserve those
values under this conversion.

The plot title is fixed as \code{Calibration Curve (Pain Prediction Nomogram)}.
The identity line represents perfect calibration. Resampling-based calibration
consumes random numbers, so callers should set a seed when reproducibility is
required.
}
\examples{
\dontrun{
calibration_data <- data.frame(
  Group = factor(rep(c("Control", "Disease"), each = 50)),
  GeneA = rnorm(100),
  GeneB = rnorm(100)
)

set.seed(123)
calibration <- calibrate_plot(
  input_data = calibration_data,
  group_colname = "Group",
  key_features = c("GeneA", "GeneB"),
  calibrate_method = "boot",
  calibrate_B = 200,
  number_of_groups = 5
)

calibration$plot
calibration$result
}
}
\seealso{
\code{\link[rms:lrm]{rms::lrm()}},
\code{\link[rms:calibrate]{rms::calibrate()}},
\code{\link[ResourceSelection:hoslem.test]{ResourceSelection::hoslem.test()}}
}
