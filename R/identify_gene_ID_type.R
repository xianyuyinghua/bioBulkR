#' Identify the predominant gene identifier type
#'
#' Infer the identifier type used by a vector of gene IDs by comparing the IDs
#' with a set of regular-expression patterns.
#'
#' @param gene_IDs A non-empty character vector of gene identifiers without
#'   missing values.
#'
#' @return A single character string indicating the inferred identifier type:
#'   `"ENSEMBL"`, `"ENTREZID"`, `"UCSC"`, `"REFSEQ"`, `"SYMBOL"`, or
#'   `"UNKNOWN"` when none of the supported patterns exceeds the matching
#'   threshold.
#'
#' @details The function checks identifier types in the following order:
#'   Ensembl gene IDs beginning with `ENSG`, numeric Entrez Gene IDs, UCSC IDs
#'   beginning with `uc`, RefSeq transcript IDs beginning with `NM_` or `NR_`,
#'   and gene symbols containing letters, numbers, periods, or hyphens.
#'
#'   A type is selected when more than 80 percent of `gene_IDs` match its
#'   pattern. Because checking stops after the first successful match, the
#'   order above determines the result if more than one pattern exceeds the
#'   threshold. The function recognizes identifier format only; it does not
#'   verify that identifiers exist in an annotation database.
#'
#' @examples
#' identify_gene_ID_type(c("ENSG00000141510", "ENSG00000171862"))
#' identify_gene_ID_type(c("7157", "672", "1956"))
#' identify_gene_ID_type(c("TP53", "BRCA1", "EGFR"))
#' identify_gene_ID_type(c("gene one", "gene two"))
#'
#' @export
identify_gene_ID_type <- function(gene_IDs) {
    # 定义基因 ID 类型及其对应的正则表达式
    gene_id_patterns <- list(
        ENSEMBL = "^ENSG",
        ENTREZID = "^[0-9]+$",
        UCSC = "^uc",
        REFSEQ = "^N[MR]_\\d+",
        SYMBOL = "^[A-Za-z0-9.-]+$"
    )
    
    # 初始化
    gene_ID_type <- "UNKNOWN"
    
    # 遍历计算匹配比例
    for (type in names(gene_id_patterns)) {
    pattern <- gene_id_patterns[[type]]
    match_ratio <- mean(grepl(pattern, gene_IDs))
    if (match_ratio > 0.8) {
        gene_ID_type <- type
        break  # 找到匹配类型后退出循环
    }
    }
    
    return(gene_ID_type)
}
