#' Combine sample metadata with expression data
#'
#' Detect whether samples are stored in the rows or columns of an expression
#' table, orient the table with samples in rows, optionally select samples by
#' name patterns, and merge selected metadata columns with all expression
#' features.
#'
#' @param exp_data An expression matrix or data frame. Sample identifiers must
#'   occur either in row names or column names, and feature identifiers should
#'   occur on the opposite dimension.
#' @param metadata_data A metadata data frame. Sample identifiers may be stored
#'   in row names or in a column whose name contains `sample` or `Sample`.
#' @param metadata_select_col A character vector naming metadata columns to
#'   retain in the merged result. Every supplied column must exist in
#'   `metadata_data`.
#' @param sample_common_special_chars By default, `NA`, all matching samples
#'   are considered. Otherwise, a character vector of regular-expression
#'   patterns used to select metadata sample identifiers before merging.
#'
#' @return A data frame with one row per sample shared by the expression and
#'   metadata tables. Row names contain sample identifiers, the selected
#'   metadata columns appear first, and all expression features follow.
#'
#' @details
#' The function counts metadata sample identifiers found in the expression row
#' names and column names. If more matches occur in columns, `exp_data` is
#' transposed; otherwise, its existing orientation is retained. Metadata
#' sample identifiers are moved to row names before the two tables are merged.
#'
#' When `sample_common_special_chars` is supplied, its values are collapsed
#' with `"|"` and passed to [base::grep()], so they are interpreted as regular
#' expressions. If no samples are shared after processing, the function prints
#' a message and stops. Otherwise, it prints input and final sample counts plus
#' a table of the selected metadata columns.
#'
#' @examples
#' expression <- matrix(
#'   c(10, 20, 30, 40, 50, 60),
#'   nrow = 2,
#'   dimnames = list(c("Gene1", "Gene2"), c("Sample1", "Sample2", "Sample3"))
#' )
#' metadata <- data.frame(
#'   Sample = c("Sample1", "Sample2", "Sample3"),
#'   Group = c("Control", "Disease", "Disease")
#' )
#'
#' combined <- exp_metadata_data_bind(
#'   exp_data = expression,
#'   metadata_data = metadata,
#'   metadata_select_col = "Group"
#' )
#'
#' @seealso [base::merge()], [tibble::column_to_rownames()]
#' @export
exp_metadata_data_bind <- function(exp_data,
                                   metadata_data,
                                   metadata_select_col = c("Group","group"),
                                   sample_common_special_chars = NA
                                  ){
    # 将表达数据和元数据合并，合并元数据中选择的列和所有表达数据
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    
    if( "sample" %in% tolower(colnames(metadata_data)) ){
        # 当元数据中有 sample / Sample
        sample_col_name <- grep("sample|Sample",colnames(metadata_data),value = T)
        meta_samples <- metadata_data[[sample_col_name]]
    }else{
        # 当元数据中无 sample / Sample
        sample_col_name <- "rownames"
        meta_samples <- rownames(metadata_data)
    }

    
    # 检查 表达数据 和 元数据 中的样本所在位置（行/列）
    if( length(sample_common_special_chars) == 1 & any(is.na(sample_common_special_chars)) ){
        # 当不指定样本中包含的特殊字符。筛选出 表达数据 和 元数据 共同的样本
        if( sample_col_name != "rownames" ){
            metadata_df <- metadata_data %>%  `rownames<-`(NULL) %>%  tibble::column_to_rownames(sample_col_name)
        }else{
            metadata_df <- metadata_data
        }
        
        # 获取共同samples
        exp_row_sample_check <-  meta_samples %in% rownames(exp_data) %>% sum()
        exp_col_sample_check <-  meta_samples %in% colnames(exp_data) %>% sum()
        if(exp_col_sample_check > exp_row_sample_check){
            exp_df <- exp_data %>% t()
        }else{
            exp_df <- exp_data
        }
    }else{
        # 指定了样本中包含的特殊字符
        common_special_chars <- paste(sample_common_special_chars,collapse = "|")
        meta_samples_select <- grep(common_special_chars,meta_samples,value = T)
        
        if( sample_col_name != "rownames" ){
            metadata_df <- metadata_data %>% dplyr::filter(!!sym(sample_col_name) %in%  meta_samples_select) %>% 
                        `rownames<-`(NULL) %>%  tibble::column_to_rownames(sample_col_name)
        }else{
            metadata_df <- metadata_data %>% .[rownames(.) %in% meta_samples_select, ]
        }
        
        
        # 获取共同samples
        exp_row_sample_check <-  meta_samples_select %in% rownames(exp_data) %>% sum()
        exp_col_sample_check <-  meta_samples_select %in% colnames(exp_data) %>% sum()
        if(exp_col_sample_check > exp_row_sample_check){
            exp_df <- exp_data %>% t()
        }else{
            exp_df <- exp_data
        }
    }

    # 组合
    exp_meta_df <- merge(metadata_df %>% dplyr::select(all_of(metadata_select_col)),exp_df,by = "row.names") %>% tibble::column_to_rownames("Row.names")
    
    # 检查样本数量
    if(nrow(exp_meta_df) == 0){
        cat("\n\n",red,"Expression matrix and Metadata have no common samples.",reset,"\n\n")
        flush.console()
        stop()
    }else{
        cat("\n\n",yellow,"Expression matrix samples: ",reset,nrow(exp_df),"\n")
        cat("\n",yellow,"Metadata samples: ",reset,nrow(metadata_df),"\n")
        cat("\n",yellow,"Final number of samples: ",reset,"\n\n")
        print(table(exp_meta_df[,metadata_select_col]))
        flush.console()
        return(exp_meta_df)
    }
}
