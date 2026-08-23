#' Interactively process GEO sample metadata
#'
#' Inspect, filter, and simplify GEO sample metadata into a sample grouping
#' table. The function displays frequency tables for likely metadata columns,
#' then interactively asks the user to filter samples, choose a grouping
#' column, and optionally rename its group labels.
#'
#' @param metadata A `data.frame` containing sample metadata. It must include a
#'   `geo_accession` column and the column selected interactively as the sample
#'   grouping variable. Metadata returned by `Biobase::pData()` is a typical
#'   input.
#'
#' @return A `data.frame` with three columns:
#'   \describe{
#'   \item{`Group`}{The final group label after optional renaming.}
#'   \item{`Sample`}{The sample accession taken from `geo_accession`.}
#'   \item{`Raw_Group`}{The original group label before renaming.}
#'   }
#'
#' @details Candidate columns whose names begin with `title`, `geo`, `source`,
#'   `characteristics`, `platform`, `scan`, `data`, `grade`, `batch`, `tissue`,
#'   or `contact` are summarized at the start of the function.
#'
#'   Filtering instructions use the form
#'   `column::value1,value2`. Multiple column filters can be separated with
#'   semicolons, for example,
#'   `column1::value1,value2;column2::value3,value4`. Values are matched
#'   exactly, and filters for multiple columns are combined sequentially.
#'   Press Enter without supplying a filter to retain all samples.
#'
#'   Group renaming instructions use the form
#'   `original_pattern::new_name`, with multiple mappings separated by commas.
#'   Original names are interpreted as regular-expression patterns by
#'   [base::grepl()]. Press Enter to retain the original group labels.
#'
#'   This function calls [base::readline()] and is therefore intended for use
#'   in an interactive R session.
#'
#' @examples
#' \dontrun{
#' metadata <- data.frame(
#'   geo_accession = c("GSM1", "GSM2", "GSM3"),
#'   tissue = c("tumor", "tumor", "normal"),
#'   stringsAsFactors = FALSE
#' )
#'
#' groups <- metadata_process(metadata)
#' }
#'
#' @export
metadata_process <- function(metadata = metadata){
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 原metadata处理为简洁的样本分组信息
    # 展示元数据关键信息
    selected_columns <- metadata[, grepl("^(title|geo|source|characteristics|platform|scan|data|grade|batch|tissue|contact)", names(metadata), ignore.case = TRUE)]
    #selected_columns <- metadata
    metadata_summary <- lapply(selected_columns, function(column) { 
        table(column, useNA = "ifany")  # 包含 NA 值的统计
    })
    for (col_name in names(metadata_summary)) {
        cat("\n",yellow, "Column:", col_name, reset,"\n")
        print(metadata_summary[[col_name]])
    }
    flush.console()  # 强制刷新控制台输出
    # 选择分组列和要筛选的组
    cat("\n",yellow,"Specify the group to filter by '::'\n\tIf filtering, enter a string unique to the group or common within the group,Use ',' to separate multiple groups.\n\tPress Enter to confirm without filtering.",reset)
    cat("\n",yellow,"If multiple groups need to be filtered, separate the groups with a ';'.",reset,"\n")
    cat("\n\t",yellow,"eg1. Single group filter: ",cyan,"characteristics_ch1.1::thompson grade: I,thompson grade: II",reset,"\n")
    cat("\n\t",yellow,"eg2. Multiple group filter: ",cyan,"characteristics_ch1.1::thompson grade: I,thompson grade: II;characteristics_ch1.2::tissue: Nucleus pulposus,tissue: Annulus fibrosus",reset,"\n")
    flush.console()  # 强制刷新控制台输出
    data_filter_by_group <- readline(prompt = "")
    if(data_filter_by_group == ""){
        # 1. 不进行过滤
        cat("\n\n",yellow,"No filtering on metadata",reset,"\n")
    }else {
        # 2. 过滤
        group_select_names <- strsplit(data_filter_by_group, ";") %>% unlist() %>% strsplit(.,"::") %>% sapply(.,"[",1)
        group_select_values <- strsplit(data_filter_by_group, ";") %>% unlist() %>% strsplit(.,"::") %>% sapply(.,"[",2) %>% strsplit(.,",") 
        names(group_select_values) <- group_select_names
    
        metadata <- Reduce(function(d, column_name) {
            filter_values <- group_select_values[[column_name]]
            d[d[[column_name]] %in% filter_values, , drop = FALSE]
        }, names(group_select_values), init = metadata)
    
        cat("\n",yellow,"Number of filtered packets",reset,"\n")
        lapply(names(group_select_values),function(column){
            print(table(metadata[[column]]))
        })
    }
    
    
    cat(paste0("\n",yellow, "Please select a grouping column: ", reset, "\n"))
    flush.console()  # 强制刷新控制台输出
    group_col <- readline(prompt = "")
    
    # 筛选分组信息列
    group_df <- metadata %>% dplyr::select(geo_accession, all_of(group_col))
    group_df <- group_df %>% rename('Group' = all_of(group_col),'Sample' = 'geo_accession')
    group_df$Raw_Group <- group_df$Group
    print(group_df)
    # 更改分组标签
    cat("\n",yellow,"Please set group renaming. Enter the original and new group names, mapped by '::', separated by ','",reset,blue,"\n\t the original name can be a pattern.",reset,magenta,"\n\t eg. recurrent miscarriage::RM,elective termination::Control",reset)
    flush.console()  # 强制刷新控制台输出
    group_rename_mode <- readline(prompt = "")
    if(group_rename_mode != ""){
        rename_pairs  <- strsplit(group_rename_mode, ",") %>% unlist() %>% lapply(function(pair) strsplit(pair,"::")[[1]])
        group_newname_list <- sapply(rename_pairs,"[",2)                                                                     
        names(group_newname_list) <- sapply(rename_pairs,"[",1)
    }
    for(pattern in names(group_newname_list)){group_df$Group[grepl(pattern, group_df$Group)] <- group_newname_list[pattern]}
    cat("\n",yellow,"Revise the grouped data：",reset,"\n")
    print(group_df)
    flush.console()  # 强制刷新控制台输出
    return(group_df)                                                                           
}
