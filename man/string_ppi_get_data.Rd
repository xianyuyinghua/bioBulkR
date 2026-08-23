\name{string_ppi_get_data}
\alias{string_ppi_get_data}
\title{Retrieve protein-protein interactions from STRING}
\usage{
string_ppi_get_data(genes, required_score = 400, species = 9606)
}
\arguments{
\item{genes}{A character vector of gene or protein identifiers accepted by the
STRING \code{get_string_ids} API endpoint, typically gene symbols.}

\item{required_score}{A numeric STRING combined-score threshold, usually
between 0 and 1000. Only interactions meeting the requested confidence
threshold are returned by STRING.}

\item{species}{An NCBI taxonomy identifier passed to STRING, such as
\code{9606} for human or \code{10090} for mouse.}
}
\value{
A tibble containing the STRING network interactions. Columns are supplied by
the STRING TSV network endpoint and typically include STRING identifiers,
\code{preferredName_A}, \code{preferredName_B}, the taxonomy ID, \code{score},
and evidence-specific scores. Only rows whose two preferred names both occur
in \code{genes} are retained.
}
\description{
Map a vector of gene identifiers to preferred STRING names, query the STRING
network endpoint, and retain interactions for which both preferred endpoint
names occur in the original input vector.
}
\details{
The function first calls the STRING \code{get_string_ids} endpoint and writes
its TSV response to a temporary file. Preferred names from that response are
then submitted to the STRING \code{network} endpoint with
\code{required_score} and \code{species}. The network response is downloaded
into memory and parsed with \code{\link[readr:read_delim]{readr::read_tsv()}}.
The temporary identifier-mapping file is removed after it is read.

Because both API requests use \code{https://string-db.org}, an active internet
connection is required. Network, HTTP, identifier-mapping, and parsing errors
are propagated to the caller. STRING synonyms that map to a preferred name
not present in the original \code{genes} vector can be removed by the final
endpoint filter.
}
\examples{
\dontrun{
interactions <- string_ppi_get_data(
  genes = c("TP53", "EGFR", "MYC", "AKT1"),
  required_score = 400,
  species = 9606
)

head(interactions)
}
}
\seealso{
\url{https://string-db.org/help/api/}
}
