#' Calculate row and column counts for a faceted plot
#'
#' Estimate a compact rectangular layout from the requested number of facets.
#' Three facets are handled specially and arranged in a single row.
#'
#' @param facet_number A positive numeric value giving the number of facets to
#'   arrange.
#'
#' @return A numeric vector of length two. The first element is the number of
#'   rows and the second is the number of columns.
#'
#' @details
#' For values other than three, the number of rows is
#' `round(sqrt(facet_number))` and the number of columns is
#' `ceiling(sqrt(facet_number))`. When `facet_number = 3`, the function returns
#' `c(1, 3)`.
#'
#' @examples
#' facet_nrow_ncol(3)
#' facet_nrow_ncol(6)
#' facet_nrow_ncol(10)
#'
#' @export
facet_nrow_ncol <- function(facet_number){
    # 根据分面数量计算行列数
    nrow = round(sqrt(facet_number))
    ncol = ceiling(sqrt(facet_number))
    if(facet_number == 3){
        nrow = 1
        ncol = 3
    }
    return(c(nrow,ncol))
}
