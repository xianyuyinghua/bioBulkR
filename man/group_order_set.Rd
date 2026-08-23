\name{group_order_set}
\alias{group_order_set}
\title{Set the order of sample groups}
\usage{
group_order_set(vector, group_levels = NULL)
}
\arguments{
\item{vector}{A vector containing sample group labels. Duplicate values are
removed before the group order is determined, and the remaining values are
converted to character strings.}

\item{group_levels}{An optional vector specifying the desired group order.
Every unique value in \code{vector} must occur in \code{group_levels}. When
\code{NULL} or empty, the order is determined automatically.}
}
\value{
A character vector containing the ordered group labels. When
\code{group_levels} does not contain every unique label from \code{vector},
the function prints a message and returns \code{NULL} invisibly.
}
\description{
Determine a consistent order for sample groups. A control-like group is
placed first when one is detected; otherwise, the unique group labels are
sorted alphabetically. A user-supplied order can be used instead.
}
\details{
The input is reduced to unique values in their original order and converted
to character. Automatic ordering recognizes \code{"control"},
\code{"normal"}, and \code{"norm"} without regard to letter case. If more
than one control-like label is present, the first one found is placed first
and all other labels retain their original relative order.

If no control-like label is found, labels are returned in alphabetical
order. When \code{group_levels} is supplied and contains every input label,
it is returned as provided; consequently, it may contain additional labels
that were not present in \code{vector}.
}
\examples{
group_order_set(c("Disease", "Control", "Disease"))

group_order_set(
  c("Treatment", "Control"),
  group_levels = c("Control", "Treatment")
)

group_order_set(c("GroupB", "GroupA"))
}
