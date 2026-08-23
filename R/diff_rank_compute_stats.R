#' Rank numeric features after splitting samples by target-gene expression
#'
#' Split samples into low and high expression groups for one or more target
#' genes, fit a limma model for every numeric feature, and return ranked
#' statistics or complete differential-analysis tables.
#'
#' @param exp_df A data frame with samples in rows and expression features in
#'   columns. Only numeric columns are analyzed; numeric metadata columns are
#'   therefore also included.
#' @param target_gene A character vector naming one or more numeric columns used
#'   independently to split samples.
#' @param split_method Sample-splitting method: `"median"` or `"quantile"`.
#' @param q For quantile splitting, the lower and upper tail proportion. Samples
#'   at or below the `q` quantile and at or above the `1 - q` quantile are kept.
#' @param min_n Minimum number of retained samples required in both the Low and
#'   High groups.
#' @param rank_stat Statistic used to rank features: `"logFC"`, `"t"`, or
#'   `"signed_log10p"`.
#' @param return_full_table Logical; return full limma results and split
#'   metadata instead of only combined ranking data.
#'
#' @return With `return_full_table = FALSE`, a one-element list named after
#'   `rank_stat`, containing a combined data frame with feature names, target
#'   genes, and ranking values. With `TRUE`, a list named by `target_gene`; each
#'   element contains `ranks`, `deg`, `groups`, `kept_samples`, `target_gene`,
#'   `split_method`, and `q`.
#'
#' @details
#' Median splitting assigns values greater than or equal to the median to High
#' and retains every sample. Quantile splitting retains only the two tails.
#' For each target, the limma design is `~ grp`; its second coefficient is
#' `High - Low`, so positive log fold changes indicate higher values in the
#' High target-expression group.
#'
#' `signed_log10p` is `-log10(P.Value)` multiplied by the sign of `logFC`.
#' Rankings are sorted from largest to smallest. The target column itself and
#' every other numeric column in `exp_df` are included in the fitted feature
#' matrix.
#'
#' @examples
#' \dontrun{
#' expression <- as.data.frame(matrix(
#'   rnorm(40 * 20),
#'   nrow = 40,
#'   dimnames = list(
#'     paste0("Sample", seq_len(40)),
#'     paste0("Gene", seq_len(20))
#'   )
#' ))
#' ranks <- diff_rank_compute_stats(
#'   exp_df = expression,
#'   target_gene = c("Gene1", "Gene2"),
#'   split_method = "quantile",
#'   q = 0.3,
#'   min_n = 5,
#'   rank_stat = "signed_log10p"
#' )
#' ranks$signed_log10p
#' }
#'
#' @seealso [limma::lmFit()], [limma::eBayes()], [limma::topTable()]
#' @export
diff_rank_compute_stats <- function(exp_df,                                   # 行=样本；列=基因+可能的元信息
                                   target_gene,                               # 可为单个字符 or 向量 c("G1","G2",...)
                                   split_method = c("median","quantile"),
                                   q = 0.3,                                   # quantile 模式：top q vs bottom q
                                   min_n = 5,                                 # 每组最少样本数
                                   rank_stat = c("logFC","t","signed_log10p"),
                                   return_full_table = FALSE
                                  ){
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    split_method <- match.arg(split_method)
    rank_stat    <- match.arg(rank_stat)
    
    # 仅保留数值型表达列
    numeric_cols <- vapply(exp_df, is.numeric, logical(1))
    expr_all <- exp_df[, numeric_cols, drop = FALSE]
    if (ncol(expr_all) == 0){stop("没有数值型表达列!")}
    
    # 统一样本ID
    if (is.null(rownames(expr_all))) {
    rownames(expr_all) <- make.names(seq_len(nrow(expr_all)), unique = TRUE)
    }
    
    # 内部：单个基因的完整流程 -----------------------------------------
    .one_gene <- function(gene){
        if (!(gene %in% colnames(expr_all))) {
            stop(cat("\n",red,"找不到目标基因列：",reset,sprintf("%s", gene),"\n"))
        }
        expr <- expr_all
        gvals <- expr[[gene]]
        
        # 分组
        if (split_method == "median") {
            thr <- stats::median(gvals, na.rm = TRUE)
            grp <- factor(ifelse(gvals >= thr, "High", "Low"),levels = c("Low","High"))
            keep <- seq_len(nrow(expr))
        }else{
            lo_thr <- stats::quantile(gvals, probs = q, na.rm = TRUE)
            hi_thr <- stats::quantile(gvals, probs = 1 - q, na.rm = TRUE)
            keep   <- which(gvals <= lo_thr | gvals >= hi_thr)
            expr   <- expr[keep, , drop = FALSE]
            grp    <- factor(ifelse(gvals[keep] >= hi_thr, "High", "Low"),levels = c("Low","High"))
        }
        
        # 样本量检查
        tab <- table(grp)
        if (any(tab < min_n)) {
            stop(cat("\n",red,"分组样本不足：",reset,sprintf("Low=%d, High=%d",as.integer(tab["Low"]), as.integer(tab["High"])),red,"请降低 q 或使用 median 分组!",reset,"\n"))
        }
        
        # limma 差异
        suppressPackageStartupMessages({ library(limma) })
        design <- model.matrix(~ grp)    # 截距 + grpHigh（High - Low）
        fit <- lmFit(t(expr), design)
        fit <- eBayes(fit)
        tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
        
        # ranks
        ranks <- switch(rank_stat,
                        "logFC" = tt$logFC,
                        "t" = tt$t,
                        "signed_log10p" = (-log10(pmax(tt$P.Value, .Machine$double.xmin))) * sign(tt$logFC)
                       )
        names(ranks) <- rownames(tt)
        ranks <- sort(ranks, decreasing = TRUE)
        
        if(return_full_table){
            tt$gene <- rownames(tt)
            tt <- tt[, c("gene","logFC","t","P.Value","adj.P.Val","AveExpr","B")]
            return(list(
                ranks = ranks,
                deg   = tt,
                groups= tab,
                kept_samples = rownames(expr),
                target_gene  = gene,
                split_method = split_method,
                q = if (split_method == "quantile") q else NA
            ))
        }else{
            ranks_df <- data.frame(common_features = names(ranks),
                                   separate_features = gene,
                                   rank_value = as.numeric(ranks),
                                   stringsAsFactors = FALSE
                                  )
            names(ranks_df)[3] <- rank_stat
            return(ranks_df)
        }
    }
                
    # ------------------------------------------------------------------
    # 返回以基因名命名的列表
    res_list <- lapply(target_gene, .one_gene)
    names(res_list) <- target_gene
                
    if(!return_full_table){
        res_data <- purrr::reduce(res_list, rbind)
        return(setNames(list(res_data), rank_stat))
    }else{
        return(res_list)
    }
}
