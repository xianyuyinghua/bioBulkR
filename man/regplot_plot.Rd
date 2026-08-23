\name{regplot_plot}
\alias{regplot_plot}
\title{Fit a logistic model and draw a regplot nomogram}
\usage{
regplot_plot(
  data = data,
  key_features = key_genes,
  showP = FALSE,
  sample_select = 1,
  formats = c("pdf", "png", "tiff")
)
}
\arguments{
\item{data}{A data frame containing an outcome column named \code{Group} and
all columns listed in \code{key_features}.}

\item{key_features}{A character vector naming predictors included in the
logistic regression model.}

\item{showP}{Logical; pass \code{TRUE} to display predictor p-values in the
regplot output.}

\item{sample_select}{A row index selecting the observation highlighted by
droplines in the nomogram.}

\item{formats}{A character vector of requested output formats. This argument
is retained by the interface but is not used by the current implementation.}
}
\value{
If no selected feature names contain special characters, the recorded
nomogram plot returned by
\code{\link[grDevices:recordplot]{grDevices::recordPlot()}}. If names are
cleaned, a magick image object read from the corrected PDF is returned.
Consequently, the return class depends on the predictor names.
}
\description{
Fit an rms logistic regression of \code{Group} on selected features, draw a
regplot nomogram for a selected observation, and handle feature names that
contain characters unsuitable for model formulas.
}
\details{
Characters other than letters, numbers, and underscores in predictor names
are replaced with underscores before model fitting. The function fits
\code{Group ~ predictors} with \code{\link[rms:lrm]{rms::lrm()}} and passes the
result to \code{\link[regplot:regplot]{regplot::regplot()}} with violin and box
displays. \code{sample_select} must identify an existing row in \code{data}.

The function assigns an rms \code{datadist} object named \code{dd} in
\code{.GlobalEnv}, sets the global \code{datadist} and \code{na.action}
options, writes \code{logistic_regplot.pdf} in the current working directory,
and deletes that PDF before returning. These global objects and options are not
restored. The PDF is written twice during plot capture and replay.

When special feature names are present, the function invokes qpdf at
\code{/data/home/zuoanjian/miniconda3/envs/RNA_Base/bin/qpdf} to replace
cleaned labels with their original names, then uses
\code{\link[magick:editing]{magick::image_read_pdf()}} at 600 DPI. This branch
requires that executable and its supporting shell tools to be available. The
current implementation does not write persistent PNG or TIFF files, regardless
of \code{formats}.
}
\examples{
\dontrun{
model_data <- data.frame(
  Group = factor(rep(c("Control", "Disease"), each = 25)),
  GeneA = rnorm(50),
  GeneB = rnorm(50)
)

nomogram <- regplot_plot(
  data = model_data,
  key_features = c("GeneA", "GeneB"),
  showP = TRUE,
  sample_select = 1
)

replayPlot(nomogram)
}
}
\seealso{
\code{\link[rms:lrm]{rms::lrm()}},
\code{\link[rms:datadist]{rms::datadist()}},
\code{\link[regplot:regplot]{regplot::regplot()}},
\code{\link[magick:editing]{magick::image_read_pdf()}}
}
