\name{read_exp_data}
\alias{read_exp_data}
\title{Read an expression matrix and sample condition table}
\usage{
read_exp_data(
  expression_dir = "./01_Expression_data",
  expression_pattern = "Expression.*\\\\.csv$",
  condition_pattern = "*Condition*.csv$"
)
}
\arguments{
\item{expression_dir}{A character string giving the directory searched for
expression and condition files. Subdirectories are searched recursively.}

\item{expression_pattern}{A regular-expression pattern passed to
\code{\link[base:list.files]{list.files()}} to select the expression-data
file. The filename is also used to infer the expression-data type.}

\item{condition_pattern}{A regular-expression pattern passed to
\code{\link[base:list.files]{list.files()}} to select the sample condition
file.}
}
\value{
A named list containing:
\describe{
\item{\code{expr}}{The expression table. Its columns are converted to numeric
when necessary, and rows containing missing values are removed.}
\item{\code{condition}}{The sample condition table.}
\item{\code{exp_data_type}}{The inferred expression-data class:
\code{"array"} or \code{"throughput"}.}
\item{\code{throughput_type}}{For high-throughput data, one of
\code{"Count"}, \code{"FPKM"}, \code{"TPM"}, \code{"RPM"}, or
\code{"unknown"}; \code{NULL} for array data.}
}
}
\description{
Locate and read an expression matrix and its sample condition table from an
expression-data directory. The function standardizes expression values to
numeric data and infers the expression platform and measurement type from
the expression filename.
}
\details{
Both files are located recursively with
\code{\link[base:list.files]{list.files()}} and imported with
\code{\link[=read_data]{read_data()}}. Patterns should select exactly one
expression file and one condition file.

Expression filenames containing \code{"array"} are classified as array
data. Filenames containing \code{"throughput"} are classified as
high-throughput data; the function then searches the filename for
\code{"Count"}, \code{"FPKM"}, \code{"TPM"}, or \code{"RPM"} to determine
the measurement type. These filename checks are case-sensitive.

The patterns use regular-expression syntax. In particular, callers can use
\code{".*Condition.*\\\\.csv$"} as a conventional condition-file pattern.
If no recognized \code{"array"} or \code{"throughput"} marker is present in
the expression filename, the function cannot determine
\code{exp_data_type}.
}
\examples{
\dontrun{
data <- read_exp_data(
  expression_dir = "03_Expression_DEG/GSE28829/01_Expression_data",
  expression_pattern = "01_Expression_array.*\\\\.csv$",
  condition_pattern = ".*Condition.*\\\\.csv$"
)

expression_matrix <- data$expr
sample_groups <- data$condition
}
}
\seealso{
\code{\link[=read_data]{read_data()}},
\code{\link[base:list.files]{list.files()}}
}
