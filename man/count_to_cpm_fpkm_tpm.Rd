\name{count_to_cpm_fpkm_tpm}
\alias{count_to_cpm_fpkm_tpm}
\title{Convert a raw count matrix to CPM, FPKM, or TPM}
\usage{
count_to_cpm_fpkm_tpm(
  df_count,
  gtf_path = "/data/home/zuoanjian/bioinfo_documents/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz",
  core_number = 10,
  trans_to_obj = "fpkm",
  species = 9606,
  gene_type = "Sybmol"
)
}
\arguments{
\item{df_count}{A numeric gene-by-sample count matrix or data frame. Gene
identifiers must be stored in row names.}

\item{gtf_path}{A path to a GTF or gzip-compressed GTF annotation readable by
\code{\link[rtracklayer:import]{rtracklayer::import()}}. If \code{NULL}, a
built-in human or mouse path is selected from \code{species}.}

\item{core_number}{A positive integer giving the number of parallel workers
used to calculate gene lengths.}

\item{trans_to_obj}{The output type: \code{"cpm"}, \code{"fpkm"}, or
\code{"tpm"}. Any other value returns all three matrices in a list. Matching
is case-sensitive.}

\item{species}{An NCBI taxonomy identifier used only to select a fallback GTF
path when \code{gtf_path = NULL}. Supported values are \code{9606} for human
and \code{10090} for mouse. An unsupported value uses the human fallback path
after printing a message.}

\item{gene_type}{The count-row identifier type. Use \code{"Symbol"} to match
GTF \code{gene_name} values or \code{"Ensembl"} to match version-stripped GTF
\code{gene_id} values. Matching is case-sensitive. The historical default
\code{"Sybmol"} does not select either identifier branch, so callers should
specify this argument explicitly.}
}
\value{
A numeric matrix containing the requested CPM, FPKM, or TPM values. If
\code{trans_to_obj} is not one of the three recognized lowercase values, a
named list containing \code{cpm}, \code{fpkm}, and \code{tpm} matrices is
returned. Only genes shared by \code{df_count} and the selected GTF identifier
are retained.
}
\description{
Import exon coordinates from a GTF annotation, calculate non-overlapping exon
length for each gene, retain genes shared by the count matrix and the
annotation, and calculate CPM, FPKM, and TPM expression values.
}
\details{
Gene length is calculated as the number of unique genomic positions covered
by all exons assigned to a gene, so overlapping exon segments are counted
once. Ensembl version suffixes following a period are removed before
matching. CPM is calculated with
\code{\link[edgeR:cpm]{edgeR::cpm()}}. FPKM is calculated from each sample's
library total and gene lengths in kilobases. TPM first divides counts by gene
length and then scales each sample to one million.

The function validates that \code{df_count} is a matrix or data frame with row
names. All three normalized matrices are calculated before
\code{trans_to_obj} selects the return value.

The function creates a PSOCK cluster with
\code{\link[parallel:makeCluster]{parallel::makeCluster()}} for the gene-length
calculation. The current implementation does not explicitly stop that cluster
before returning.
}
\examples{
\dontrun{
counts <- matrix(
  c(100, 120, 80, 95, 200, 240),
  nrow = 2,
  dimnames = list(
    c("TP53", "EGFR"),
    c("Sample1", "Sample2", "Sample3")
  )
)

tpm <- count_to_cpm_fpkm_tpm(
  df_count = counts,
  gtf_path = "annotation/Homo_sapiens.GRCh38.gtf.gz",
  core_number = 2,
  trans_to_obj = "tpm",
  gene_type = "Symbol"
)
}
}
\seealso{
\code{\link[edgeR:cpm]{edgeR::cpm()}},
\code{\link[rtracklayer:import]{rtracklayer::import()}},
\code{\link[parallel:makeCluster]{parallel::makeCluster()}}
}
