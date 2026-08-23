\name{exp_metadata_data_bind}
\alias{exp_metadata_data_bind}
\title{Combine sample metadata with expression data}
\usage{
exp_metadata_data_bind(
  exp_data,
  metadata_data,
  metadata_select_col = c("Group", "group"),
  sample_common_special_chars = NA
)
}
\arguments{
\item{exp_data}{An expression matrix or data frame. Sample identifiers must
occur either in row names or column names, and feature identifiers should
occur on the opposite dimension.}

\item{metadata_data}{A metadata data frame. Sample identifiers may be stored
in row names or in a column whose name contains \code{sample} or
\code{Sample}. The function expects no more than one matching sample column.}

\item{metadata_select_col}{A character vector naming metadata columns to
retain in the merged result. Every supplied column must exist in
\code{metadata_data}. The default requests both \code{Group} and
\code{group}; specify the actual column or columns when the metadata contains
only one of these names.}

\item{sample_common_special_chars}{By default, \code{NA}, all matching samples
are considered. Otherwise, a character vector of regular-expression patterns
used to select metadata sample identifiers before merging.}
}
\value{
A data frame with one row per sample shared by the expression and metadata
tables. Row names contain sample identifiers, the selected metadata columns
appear first, and all expression features follow.
}
\description{
Detect whether samples are stored in the rows or columns of an expression
table, orient the table with samples in rows, optionally select samples by
name patterns, and merge selected metadata columns with all expression
features.
}
\details{
The function counts metadata sample identifiers found in the expression row
names and column names. If more matches occur in columns, \code{exp_data} is
transposed; otherwise, its existing orientation is retained. Metadata sample
identifiers are moved to row names before the two tables are merged. Ties in
the orientation check leave \code{exp_data} unchanged.

When \code{sample_common_special_chars} is supplied, its values are collapsed
with \code{"|"} and passed to \code{\link[base:grep]{base::grep()}}, so they
are interpreted as regular expressions. The filtering is applied to metadata
sample identifiers before expression orientation and merging.

If no samples are shared after processing, the function prints a message and
stops. Otherwise, it prints the expression, metadata, and final sample counts
plus a contingency table of the selected metadata columns. The merge uses
sample identifiers and may reorder the returned rows.
}
\examples{
expression <- matrix(
  c(10, 20, 30, 40, 50, 60),
  nrow = 2,
  dimnames = list(
    c("Gene1", "Gene2"),
    c("Sample1", "Sample2", "Sample3")
  )
)

metadata <- data.frame(
  Sample = c("Sample1", "Sample2", "Sample3"),
  Group = c("Control", "Disease", "Disease")
)

combined <- exp_metadata_data_bind(
  exp_data = expression,
  metadata_data = metadata,
  metadata_select_col = "Group"
)
}
\seealso{
\code{\link[base:merge]{base::merge()}},
\code{\link[tibble:rownames]{tibble::column_to_rownames()}}
}
