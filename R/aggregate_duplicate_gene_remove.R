#' Resolve duplicated genes in an expression table
#'
#' Produce one expression profile per gene by either retaining the row with
#' the highest mean expression or aggregating duplicated rows with a summary
#' function.
#'
#' @param data A `data.frame` containing a gene identifier column and one or
#'   more expression sample columns.
#' @param gene_col A single character string giving the name of the gene
#'   identifier column in `data`.
#' @param sample_col A character vector giving the names of expression columns
#'   in `data`.
#' @param duplicate_remove_method A character string specifying how duplicated
#'   genes are resolved. The default, `"mean-max"`, retains the row having the
#'   largest mean across `sample_col`. Other values are passed to
#'   [base::match.fun()], for example, `"mean"`, `"median"`, or `"max"`, and
#'   are used to aggregate each sample column by gene.
#'
#' @return A `data.frame` with one row per gene. The first column has the name
#'   supplied in `gene_col`, followed by the columns supplied in `sample_col`.
#'
#' @details Only `gene_col` and `sample_col` are retained in the output. Sample
#'   columns are converted to a numeric matrix before duplicate processing, so
#'   their values must be coercible to numeric.
#'
#'   With `duplicate_remove_method = "mean-max"`, a mean is calculated across
#'   the selected sample columns for every input row using `na.rm = TRUE`.
#'   Within each gene, the row with the largest mean is retained; its original
#'   expression values are returned unchanged after numeric conversion.
#'
#'   For any other method, [stats::aggregate()] applies the selected function
#'   independently to every sample column for each gene, also using
#'   `na.rm = TRUE`.
#'
#' @examples
#' expression_data <- data.frame(
#'   gene = c("TP53", "TP53", "BRCA1"),
#'   sample_1 = c(4, 8, 5),
#'   sample_2 = c(6, 10, 7)
#' )
#'
#' aggregate_duplicate_gene_remove(
#'   expression_data,
#'   gene_col = "gene",
#'   sample_col = c("sample_1", "sample_2")
#' )
#'
#' aggregate_duplicate_gene_remove(
#'   expression_data,
#'   gene_col = "gene",
#'   sample_col = c("sample_1", "sample_2"),
#'   duplicate_remove_method = "mean"
#' )
#'
#' @export
aggregate_duplicate_gene_remove <- function(data, gene_col, sample_col,
                                            duplicate_remove_method = "mean-max"
                                           ) {
    
    stopifnot(gene_col %in% names(data))
    stopifnot(all(sample_col %in% names(data)))
    
    # 1) 只保留必要列，避免 dplyr 处理奇怪类型列
    df <- data[, c(gene_col, sample_col), drop = FALSE]
    
    # 2) 确保样本列是纯 numeric（强制转成 matrix 再转回）
    expr_mat <- as.matrix(df[, sample_col, drop = FALSE])
    storage.mode(expr_mat) <- "double"
    df[, sample_col] <- as.data.frame(expr_mat)
    
    if (duplicate_remove_method != "mean-max") {
        # 按 duplicate_remove_method 聚合（如 mean / median / max）
        agg_fun <- match.fun(duplicate_remove_method)
        out <- aggregate(df[, sample_col, drop = FALSE],
                         by = list(df[[gene_col]]),
                         FUN = agg_fun, na.rm = TRUE)
        colnames(out)[1] <- gene_col
        return(out)
    } else {
        # mean-max：先算每行均值，再每个基因取 rowMean 最大的一条
        df$rowMean <- rowMeans(expr_mat, na.rm = TRUE)
        
        suppressPackageStartupMessages(library(dplyr))
        out <- df %>%
          dplyr::group_by(.data[[gene_col]]) %>%
          dplyr::arrange(dplyr::desc(.data$rowMean), .by_group = TRUE) %>%
          dplyr::slice_head(n = 1) %>%   # 比 slice(1) 更稳
          dplyr::ungroup() %>%
          dplyr::select(dplyr::all_of(c(gene_col, sample_col)))
        
        return(as.data.frame(out))
    }
}
