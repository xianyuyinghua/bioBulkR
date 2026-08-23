\name{Array_DEGs}
\alias{Array_DEGs}
\title{Identify differentially expressed genes in array data}
\usage{
Array_DEGs(
  expr = data$expr,
  condition = data$condition,
  group_order = data$group_order,
  output_dir = "02_DEGs",
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05,
  top_gene = 10,
  change = c("Up", "Down", "Not")
)
}
\arguments{
\item{expr}{A numeric gene-by-sample expression matrix or data frame. Gene
identifiers should be stored in row names and sample identifiers in column
names.}

\item{condition}{A data frame containing sample group information. It must
have exactly one column named \code{group} or \code{groups}, ignoring case.
Sample identifiers may be supplied in exactly one column named
\code{sample} or \code{samples}, ignoring case, or in row names.}

\item{group_order}{A character vector giving the group order used to build
the design and contrast. The differential comparison is
\code{group_order[1] - group_order[2]}.}

\item{output_dir}{A character string giving the directory in which the
differential-expression result is written. The directory is created
recursively when it does not exist.}

\item{logFC}{A non-negative numeric fold-change threshold on the log2 scale.
A gene is labelled upregulated when \code{logFC} is strictly greater than
this value and downregulated when it is strictly less than its negative.}

\item{p_value}{The significance measure used to classify genes. Use
\code{"p.adj"} for the FDR-adjusted p-value; any other value selects the raw
p-value.}

\item{p}{A numeric significance threshold. Genes must have a p-value strictly
less than this value to be classified as differentially expressed.}

\item{top_gene}{An integer indicating a requested number of top genes. This
argument is retained by the function interface but is not used in the
current calculation.}

\item{change}{A character vector of three labels used, in order, for
upregulated, downregulated, and non-significant genes.}
}
\value{
A data frame containing the complete
\code{\link[limma:toptable]{limma::topTable()}} result after rows with missing
values have been removed. The p-value columns are named \code{p.value} and
\code{p.adj}, and an additional \code{sig} column contains the labels from
\code{change}. The same table is written to \code{01_limma_Res.csv} under
\code{output_dir}.
}
\description{
Perform a two-group differential-expression analysis for microarray data
with the limma workflow. The function aligns expression samples with their
condition records, constructs a no-intercept design matrix, fits the chosen
contrast, classifies genes by fold change and statistical significance, and
saves the complete result.
}
\details{
If all expression values are non-negative and the maximum is at least
10,000, the expression matrix is transformed with \code{log2(expr + 1)}. If
all values are non-negative and every sample's 95th percentile is below 25,
the values are used unchanged. For other ranges, the function displays the
range and asks interactively whether the current data should be used; an
answer other than \code{"yes"} stops the analysis.

Only samples shared by the expression column names and condition row names
are retained. The condition table is arranged by its group labels before the
expression columns are reordered. Differential expression is calculated
with \code{\link[limma:lmFit]{limma::lmFit()}},
\code{\link[limma:contrasts.fit]{limma::contrasts.fit()}}, and
\code{\link[limma:ebayes]{limma::eBayes()}}. Multiple-testing correction in
\code{topTable()} uses the FDR method.
}
\examples{
\dontrun{
expr <- matrix(
  rnorm(100 * 6, mean = 8),
  nrow = 100,
  dimnames = list(
    paste0("Gene", seq_len(100)),
    paste0("Sample", seq_len(6))
  )
)

condition <- data.frame(
  Group = rep(c("Control", "Disease"), each = 3),
  row.names = colnames(expr)
)

deg_result <- Array_DEGs(
  expr = expr,
  condition = condition,
  group_order = c("Control", "Disease"),
  output_dir = tempfile("array_degs_"),
  logFC = 0.5,
  p_value = "p.adj",
  p = 0.05
)
}
}
\seealso{
\code{\link[limma:lmFit]{limma::lmFit()}},
\code{\link[limma:makeContrasts]{limma::makeContrasts()}},
\code{\link[limma:toptable]{limma::topTable()}}
}
