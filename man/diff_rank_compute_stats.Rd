\name{diff_rank_compute_stats}
\alias{diff_rank_compute_stats}
\title{Rank numeric features after splitting samples by target-gene expression}
\usage{
diff_rank_compute_stats(
  exp_df,
  target_gene,
  split_method = c("median", "quantile"),
  q = 0.3,
  min_n = 5,
  rank_stat = c("logFC", "t", "signed_log10p"),
  return_full_table = FALSE
)
}
\arguments{
\item{exp_df}{A data frame with samples in rows and expression features in
columns. Only numeric columns are analyzed; numeric metadata columns are
therefore also included.}

\item{target_gene}{A character vector naming one or more numeric columns used
independently to split samples.}

\item{split_method}{Sample-splitting method: \code{"median"} or
\code{"quantile"}. The first value is selected by default through
\code{\link[base:match.arg]{match.arg()}}.}

\item{q}{For quantile splitting, the lower and upper tail proportion. Samples
at or below the \code{q} quantile and at or above the \code{1 - q} quantile
are kept. Values should normally be between zero and 0.5.}

\item{min_n}{Minimum number of retained samples required in both the Low and
High groups.}

\item{rank_stat}{Statistic used to rank features: \code{"logFC"}, \code{"t"},
or \code{"signed_log10p"}. The first value is selected by default.}

\item{return_full_table}{Logical; return full limma results and split metadata
instead of only combined ranking data.}
}
\value{
With \code{return_full_table = FALSE}, a one-element list named after
\code{rank_stat}. Its data frame contains:
\describe{
\item{\code{common_features}}{The ranked numeric feature name.}
\item{\code{separate_features}}{The target gene used for the sample split.}
\item{ranking column}{A numeric column named \code{logFC}, \code{t}, or
\code{signed_log10p}, according to \code{rank_stat}.}
}

With \code{return_full_table = TRUE}, a list named by \code{target_gene}. Each
element contains:
\describe{
\item{\code{ranks}}{A decreasing named vector of ranking statistics.}
\item{\code{deg}}{A limma table with gene, logFC, t statistic, raw and adjusted
p-values, average expression, and B statistic.}
\item{\code{groups}}{Low and High sample counts.}
\item{\code{kept_samples}}{Row names of samples retained in the analysis.}
\item{\code{target_gene}}{The target gene used for this split.}
\item{\code{split_method}}{The selected splitting method.}
\item{\code{q}}{The quantile proportion, or \code{NA} for median splitting.}
}
}
\description{
Split samples into low and high expression groups for one or more target genes,
fit a limma model for every numeric feature, and return ranked statistics or
complete differential-analysis tables.
}
\details{
Median splitting assigns target values greater than or equal to the median to
High and retains every sample. Quantile splitting retains only the two tails.
Both groups must contain at least \code{min_n} samples. Missing target values
can interfere with group construction or model fitting and should be handled
before calling the function.

For each target, the limma design is \code{~ grp}; its second coefficient is
\code{High - Low}, so positive log fold changes indicate higher values in the
High target-expression group. The workflow uses
\code{\link[limma:lmFit]{limma::lmFit()}},
\code{\link[limma:ebayes]{limma::eBayes()}}, and
\code{\link[limma:toptable]{limma::topTable()}} without reordering the initial
limma result.

\code{signed_log10p} is \code{-log10(P.Value)} multiplied by the sign of
\code{logFC}; p-values are bounded below by the smallest positive floating-point
value before transformation. Rankings are sorted from largest to smallest.
The target column itself and every other numeric column in \code{exp_df} are
included in the fitted feature matrix.
}
\examples{
\dontrun{
expression <- as.data.frame(matrix(
  rnorm(40 * 20),
  nrow = 40,
  dimnames = list(
    paste0("Sample", seq_len(40)),
    paste0("Gene", seq_len(20))
  )
))

ranks <- diff_rank_compute_stats(
  exp_df = expression,
  target_gene = c("Gene1", "Gene2"),
  split_method = "quantile",
  q = 0.3,
  min_n = 5,
  rank_stat = "signed_log10p"
)

ranks$signed_log10p
}
}
\seealso{
\code{\link[limma:lmFit]{limma::lmFit()}},
\code{\link[limma:ebayes]{limma::eBayes()}},
\code{\link[limma:toptable]{limma::topTable()}}
}
