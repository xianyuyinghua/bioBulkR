\name{geodata_download_and_process}
\alias{geodata_download_and_process}
\title{Download and process bulk expression data from GEO}
\usage{
geodata_download_and_process(
  workdir = "/project/project_analysis/26_26YLM199F/02_Result/",
  Project_ID = "Expression_DEG",
  Data_ID = "GSE28829",
  num = 1,
  gencode_annotation = "/reference/genome/human/GRCh38/GENCODE/gencode.v36.annotation.csv",
  gene_type = "protein_coding",
  duplicate_remove_method = "mean-max",
  supplementary_sample_model = "^GSM",
  process_rawdata = NULL,
  gpl_column = NULL,
  group_levels = NULL,
  throughput_data_type = NULL,
  count_file_keywords = NULL,
  count_gene_column = NULL,
  sample_column_keywords = NULL,
  sample_rename_mode = NULL
)
}
\arguments{
\item{workdir}{A character string giving the parent directory in which the
project directory is created or reused. The function changes the R working
directory to this path.}

\item{Project_ID}{A character string used to identify and name the analysis
project directory.}

\item{Data_ID}{A GEO Series accession, such as \code{"GSE28829"}.}

\item{num}{A positive integer selecting the \code{ExpressionSet} to process
when the downloaded GEO Series contains multiple platforms.}

\item{gencode_annotation}{A path to a GENCODE annotation CSV file used by
\code{\link[=trans_gene_ID_type]{trans_gene_ID_type()}} when gene identifiers
need to be converted.}

\item{gene_type}{The gene type retained during identifier conversion. For
example, \code{"protein_coding"}; use \code{"all"} to avoid filtering by gene
type.}

\item{duplicate_remove_method}{The method passed to
\code{\link[=aggregate_duplicate_gene_remove]{aggregate_duplicate_gene_remove()}}
when multiple rows map to the same gene symbol.}

\item{supplementary_sample_model}{A regular expression passed as
\code{sample_modle} to
\code{\link[=supplementary_data_process]{supplementary_data_process()}} for
selecting sample files in a supplementary archive.}

\item{process_rawdata}{For array data with negative expression values,
\code{"yes"} downloads and processes raw data; another value declines. When
\code{NULL}, the function asks interactively.}

\item{gpl_column}{The GPL identifier mapping in the form
\code{"TYPE:COLUMN"}, for example \code{"SYMBOL:Gene Symbol"} or
\code{"ENTREZID:ENTREZ_GENE_ID"}. When \code{NULL} and automatic GPL-column
detection fails, the function asks interactively.}

\item{group_levels}{A character vector giving the desired order of sample
groups. It is used when no group named \code{normal} or \code{control} is
detected; if \code{NULL}, the order is requested interactively.}

\item{throughput_data_type}{The high-throughput expression source to process:
\code{"count"} or \code{"supplementary"}. When \code{NULL}, the function asks
interactively.}

\item{count_file_keywords}{A character string used to identify count files
when more than one downloaded file is available. When \code{NULL}, the
keywords are requested interactively.}

\item{count_gene_column}{The name of the gene-identifier column in a count
table. When \code{NULL} and it cannot be selected automatically, the name is
requested interactively.}

\item{sample_column_keywords}{One or more character strings used to identify
sample-expression columns in a count table. When \code{NULL}, the keywords
are requested interactively.}

\item{sample_rename_mode}{Optional sample-name mappings used when expression
columns cannot be matched to metadata. Supply comma-separated
\code{"original::new"} pairs; the original names can be patterns. When
\code{NULL}, mappings are requested interactively if needed.}
}
\value{
An invisible named list containing:
\describe{
\item{\code{project_folder}}{The project-directory name.}
\item{\code{data_dir}}{The path to the GEO data directory within the project.}
\item{\code{expression_dir}}{The path to the expression-data directory.}
\item{\code{data}}{A list containing the downloaded GEO object, initial
expression matrix, sample metadata, data type, and species taxonomy ID.}
\item{\code{expr}}{The processed gene-by-sample expression matrix, or
\code{NULL} if no matrix was produced.}
\item{\code{condition}}{The processed sample grouping table, or \code{NULL}
if no grouping table was produced.}
}
}
\description{
Download a GEO Series and prepare its expression matrix and sample metadata
for downstream bulk-expression analysis. The workflow handles array data and
high-throughput count or supplementary data, converts gene identifiers to
symbols, resolves duplicate genes, orders samples by group, and saves the
processed results in a project directory.
}
\details{
The function creates (or reuses) a numbered project directory under
\code{workdir}. It saves the initially downloaded data as \code{01_data.rds}
and writes processed expression and condition files under
\code{01_Expression_data}. Depending on the input type, raw and supplementary
files can also be downloaded and extracted. Existing output files may be
overwritten.

Several decisions require user input unless their corresponding arguments are
provided. Supplying these arguments is useful in a non-interactive session,
although supplementary-data processing can itself contain interactive steps.
The function supports human (taxonomy ID 9606) and mouse (taxonomy ID 10090)
gene-annotation settings.
}
\examples{
\dontrun{
result <- geodata_download_and_process(
  workdir = tempdir(),
  Project_ID = "GEO_example",
  Data_ID = "GSE28829",
  process_rawdata = "no",
  group_levels = c("Control", "Disease")
)

result$expr
result$condition
}
}
\seealso{
\code{\link[=download_gset]{download_gset()}},
\code{\link[=metadata_process]{metadata_process()}},
\code{\link[=supplementary_data_process]{supplementary_data_process()}},
\code{\link[=trans_gene_ID_type]{trans_gene_ID_type()}},
\code{\link[=aggregate_duplicate_gene_remove]{aggregate_duplicate_gene_remove()}}
}
