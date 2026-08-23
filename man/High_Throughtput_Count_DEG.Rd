\name{High_Throughtput_Count_DEG}
\alias{High_Throughtput_Count_DEG}
\title{Identify differentially expressed genes from high-throughput count data}
\usage{
High_Throughtput_Count_DEG(
  expr = data$expr,
  condition = data$condition,
  group_order = data$group_order,
  Exp_output_dir = "01_Expression_data",
  DEG_output_dir = "02_DEGs",
  logFC = 0.5,
  p_value = "p.adj",
  expr_cutoff = 10,
  filter_sample_num = nrow(data$condition) * 0.85,
  p = 0.05,
  change = c("Up", "Down", "Not")
)
}
\arguments{
\item{expr}{A gene-by-sample matrix or data frame of raw, non-negative integer
counts. Gene identifiers should be stored in row names and sample identifiers
in column names.}

\item{condition}{A data frame containing a \code{Group} column. Samples must
be identified either by row names or by a \code{sample} column. An optional
\code{Raw_Group} column is removed before the DESeq2 object is constructed.}

\item{group_order}{A character vector defining group factor levels. The
differential contrast is \code{group_order[1]} versus
\code{group_order[2]}.}

\item{Exp_output_dir}{A character string giving the directory in which the
variance-stabilized expression matrix is written. It must already exist.}

\item{DEG_output_dir}{A character string giving the directory in which the
differential-expression result is written. It is created recursively when
necessary.}

\item{logFC}{A non-negative absolute log2 fold-change threshold.}

\item{p_value}{The significance measure used for classification. Use
\code{"p.adj"} for adjusted p-values; any other value selects raw p-values.}

\item{expr_cutoff}{A non-negative count threshold used to filter genes.}

\item{filter_sample_num}{The minimum sample-count criterion for filtering. A
gene is retained only when the number of qualifying samples is strictly
greater than this value. The default is 85 percent of the number of condition
rows.}

\item{p}{A numeric significance threshold.}

\item{change}{Three labels used, in order, for upregulated, downregulated, and
non-significant genes.}
}
\value{
A data frame containing DESeq2 statistics, log2 normalized counts for each
sample, renamed \code{logFC}, \code{p.value}, and \code{p.adj} columns, and a
factor column \code{sig} containing the labels from \code{change}. Rows with
missing result values are removed.
}
\description{
Filter low-expression genes, perform a two-group DESeq2 analysis, classify
genes by fold change and statistical significance, and save both the full
differential-expression result and a variance-stabilized expression matrix.
}
\details{
When \code{expr_cutoff = 0}, genes are retained when counts are greater than
zero in more than \code{filter_sample_num} samples. For a nonzero cutoff,
counts must be greater than or equal to \code{expr_cutoff}. DESeq2 is run with
design \code{~ Group}, independent filtering is disabled in
\code{\link[DESeq2:results]{DESeq2::results()}}, and the contrast compares the
first two values of \code{group_order}.

The condition table is arranged by its factor-valued \code{Group} column and
the expression columns are reordered to match. The function therefore expects
the expression matrix and condition table to describe the same samples.

The returned table and \code{01_Deseq2_Res.csv} contain log2 normalized counts
merged with the DESeq2 result. The separate file
\code{01_Expression_throughput_vst_Norm.csv} contains the matrix produced by
\code{\link[DESeq2:vst]{DESeq2::vst()}} with \code{blind = FALSE}.

The historical function name contains the spelling \code{Throughtput} and is
retained for compatibility.
}
\examples{
\dontrun{
counts <- matrix(
  rnbinom(100 * 6, mu = 100, size = 2),
  nrow = 100,
  dimnames = list(
    paste0("Gene", seq_len(100)),
    paste0("Sample", seq_len(6))
  )
)

condition <- data.frame(
  Group = rep(c("Control", "Disease"), each = 3),
  row.names = colnames(counts)
)

result <- High_Throughtput_Count_DEG(
  expr = counts,
  condition = condition,
  group_order = c("Control", "Disease"),
  Exp_output_dir = "01_Expression_data",
  DEG_output_dir = "02_DEGs",
  logFC = 0.5,
  p_value = "p.adj",
  expr_cutoff = 10,
  filter_sample_num = 2,
  p = 0.05
)
}
}
\seealso{
\code{\link[DESeq2:DESeqDataSet]{DESeq2::DESeqDataSetFromMatrix()}},
\code{\link[DESeq2:DESeq]{DESeq2::DESeq()}},
\code{\link[DESeq2:results]{DESeq2::results()}},
\code{\link[DESeq2:vst]{DESeq2::vst()}}
}
