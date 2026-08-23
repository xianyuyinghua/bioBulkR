\name{facet_nrow_ncol}
\alias{facet_nrow_ncol}
\title{Calculate row and column counts for a faceted plot}
\usage{
facet_nrow_ncol(facet_number)
}
\arguments{
\item{facet_number}{A positive numeric value giving the number of facets to
arrange.}
}
\value{
A numeric vector of length two. The first element is the number of rows and
the second is the number of columns.
}
\description{
Estimate a compact rectangular layout from the requested number of facets.
Three facets are handled specially and arranged in a single row.
}
\details{
For values other than three, the number of rows is
\code{round(sqrt(facet_number))} and the number of columns is
\code{ceiling(sqrt(facet_number))}. When \code{facet_number = 3}, the function
returns \code{c(1, 3)}.
}
\examples{
facet_nrow_ncol(3)
facet_nrow_ncol(6)
facet_nrow_ncol(10)
}
