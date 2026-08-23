#' Interactively process GEO supplementary expression data
#'
#' Process either a supplementary-data archive containing one expression file
#' per sample or a single supplementary table. The function interactively asks
#' the user to identify gene, expression, and optional annotation columns and
#' returns standardized expression data.
#'
#' @param file_path A single character string giving the path to a
#'   supplementary-data file. Files ending in `.tar` are handled as multi-file
#'   archives; other supported formats are passed to [read_data()].
#' @param read_data_function A function used by `check_column_consistency()` to
#'   read individual files extracted from a `.tar` archive. It must accept a
#'   file path and return a tabular object. The default is [read_data()].
#' @param sample_modle A regular-expression pattern used to select sample files
#'   recursively from an extracted `.tar` archive. The default, `"^GSM"`,
#'   selects filenames beginning with a GEO Sample accession. The argument name
#'   is retained as `sample_modle` for compatibility.
#'
#' @return A named list containing:
#'   \describe{
#'   \item{`gene_exp_type`}{The expression-data type entered by the user, such
#'     as `"Count"`, `"TPM"`, or `"FPKM"`.}
#'   \item{`gene_expression_data`}{A table containing gene identifiers and
#'     expression values. For a multi-file archive, sample tables are merged by
#'     their common gene identifiers and the gene column is named `Gene_ID`.}
#'   \item{`gene_another_data`}{Optional additional columns selected by the
#'     user, or `NULL` when no additional information is requested.}
#'   }
#'
#' @details For a `.tar` input, the archive is extracted into a directory next
#'   to `file_path` whose name is formed by removing the `.tar` suffix. Nested
#'   `.gz` sample files are decompressed in place with [R.utils::gunzip()]. The
#'   selected sample files are then passed to `check_column_consistency()`. If
#'   their columns differ, the function displays the consistent and
#'   inconsistent groups and stops so that the files can be reviewed.
#'
#'   When archive columns are consistent, the function proposes character
#'   columns as gene-ID candidates and integer columns as expression-value
#'   candidates. The user confirms the desired columns and expression type.
#'   Only genes shared by every sample file are included in the merged
#'   expression table.
#'
#'   For a non-`.tar` input, [read_data()] imports the table. The user supplies
#'   comma-separated patterns for the gene-ID and expression columns. Patterns
#'   are matched against column names with [base::grep()], so a pattern can
#'   select more than one column. Optional additional-column patterns use the
#'   same comma-separated format.
#'
#'   The function is intended for an interactive R session because it uses
#'   [base::readline()]. Archive extraction creates files and directories next
#'   to the input path, and gzip decompression can remove the compressed sample
#'   files after successful extraction.
#'
#' @seealso [read_data()], [R.utils::gunzip()], [utils::untar()]
#'
#' @examples
#' \dontrun{
#' processed <- supplementary_data_process(
#'   file_path = "GSE00000_RAW.tar",
#'   sample_modle = "^GSM"
#' )
#'
#' expression_data <- processed$gene_expression_data
#' }
#'
#' @export
supplementary_data_process <- function(file_path, read_data_function = read_data, sample_modle = "^GSM"){
    # 处理supplementary文件
    file_ext <- tools::file_ext(file_path)
    if ( file_ext == "tar" ){ 
        # 1.tar 文件
        uncompressed_file <- gsub(".tar$", "", file_path)
        untar(file_path, exdir = uncompressed_file)
        files_in_tar <- list.files(path = uncompressed_file, pattern = sample_modle, full.names = T,recursive = T)
        files_in_tar_ext <- tools::file_ext(files_in_tar)
        if( any(files_in_tar_ext == "gz")){
            # 如果tar解压后是gz文件，则解压gz文件
            files_in_tar_gz <- files_in_tar[which(files_in_tar_ext == "gz")]
            lapply(files_in_tar_gz,function(file){
                uncompressed_file_gz <- gsub(".gz$", "", file)
                R.utils::gunzip(file, destname = uncompressed_file_gz, overwrite = TRUE)
            })
        }
        files_in_tar <- list.files(path = uncompressed_file, pattern = sample_modle, full.names = T,recursive = T)
        # 读取tar 文件夹中所有文件并统计列的一致性
        result <- check_column_consistency(files = files_in_tar, read_data_function = read_data_function, sample_modle = sample_modle)
        if( length(result$inconsistent_files) == 0 ){
            # 1.1 当tar文件夹中的文件的列名没有不一致的情况
            cat("\n",yellow,"All files have consistent column names in: \n\t",reset,cyan,uncompressed_file,reset,"\n\n")
            cat("\n",yellow,"The first data head:",reset,"\n\n")
            print(head(result$data_list[[1]]))
            flush.console()  # 强制刷新控制台输出
            # 输入提取的geneID列和表达值列
            integer_columns <- names(result$data_list[[1]])[sapply(result$data_list[[1]], class) == "integer"]
            gene_columns <- names(result$data_list[[1]])[sapply(result$data_list[[1]], class) == "character"]
            if(length(integer_columns) > 1){
                # 判断count列
                count_column <- readline(prompt = "Please select gene expression column name: " )
            }else{
                if( min(result$data_list[[1]][[integer_columns]] ) >= 0 & max(result$data_list[[1]][[integer_columns]]) >= 10000 ){
                    count_column <- integer_columns
                }else{
                    count_column <- readline(prompt = "Please select gene expression column name: " )
                }
            }
            if(length(gene_columns) > 1){
                # 判断gene列
                gene_column <- readline(prompt = "Please select gene ID column name: ")
            }else{
                gene_column <- gene_columns
            }
            
            # 输入表达值得类型
            gene_exp_type <- readline(prompt = "Please input the gene expression type: (Count/TPM/FPKM/RPM/···)")
            
            # 输入提取的其它列
            cat("\n",yellow,"In addition to gene expression data, is there any other information that needs to be preserved?",reset,"\n")
            cat("\n",yellow,"To save other information, enter the column name, otherwise enter nothing.",reset,"\n")
            flush.console()  # 强制刷新控制台输出
            another_columns <- readline(prompt = "Please input the columns name' (Separated by ‘,’)." )
            if(another_columns == ""){
                another_columns_select <- NULL
            }else{
                another_columns_select <- another_columns %>% strsplit(.,"[,]") %>% unlist()
                another_columns_select <- c(gene_column,another_columns_select) %>% unique()
            }
            another_columns_list <- lapply(another_columns_select,function(another_modle){
                another_col_select  <- grep(another_modle,colnames(result$data_list[[1]]),value = T)
                return(another_col_select)
            })
            another_columns_select_list <- Reduce(function(x, y) c(x,y), another_columns_list) %>% unique() 
            
            # 筛选共同 gene 并 提取 gene 和 表达值 列
            common_genes <- Reduce(intersect, lapply(result$data_list, function(df) { df[[gene_column]] }) )                     
            result$exp_select_list <- lapply(seq_along(result$data_list), function(num){
                df <- result$data_list[[num]]
                sample_name <- names(result$data_list)[num]
                data <- df[df[[gene_column]]  %in% common_genes,]    
                expression_data <- data.frame(Gene = data[[gene_column]], Exp_Value = data[[count_column]])
                colnames(expression_data) <- c("Gene_ID", sample_name)
                return(expression_data)
            })                                    
            gene_expression_data <- Reduce(function(x, y) merge(x, y, by = "Gene_ID", all = TRUE), result$exp_select_list)
            
            # 提取其他列
            result$another_select_list <- lapply(seq_along(result$data_list), function(num){
                df <- result$data_list[[num]]
                sample_name <- names(result$data_list)[num]
                data <- df[df[[gene_column]]  %in% common_genes,]
                if(!is.null(another_columns_select_list)){
                    another_data <- data[,another_columns_select_list]
                    colnames(another_data) <- c("Gene_ID",paste0(sample_name,"_",another_columns_select_list[2:length(another_columns_select_list)]))
                }else{
                    another_data <- NULL
                }
                return(another_data)
            }) 
            if( !is.null(result$another_select_list[[1]]) ){
                gene_another_data <- Reduce(function(x, y) merge(x, y, by = "Gene_ID", all = TRUE), result$another_select_list)  
            }else{
                gene_another_data <- NULL
            }
        
                 
        }else{
            # 1.2 当tar文件夹中的文件的列名不一致的情况
            cat("\n",yellow,"All files have inconsistent column names in: \n\t",reset,cyan,uncompressed_file,reset,"\n\n")
            cat("\n",yellow,"The first data: ",names(result$consistent_files)[1],reset,"\n\n")
            print(head(result$data_list[[1]]))
            cat("\n",yellow,"The samples consistent with the first data are: ",length(names(result$consistent_files)),"samples",reset,"\n")
            print(names(result$consistent_files))
            cat("\n",yellow,"The another data head:",reset,"\n\n\n")
            cat("\n",yellow,"The first data of inconsistent_files: ",names(result$inconsistent_files)[1],reset,"\n\n")
            print(head(result$data_list[[names(result$inconsistent_files)[1]]]))
            cat("\n",yellow,"The samples inconsistent with the first data are: ",length(names(result$inconsistent_files)),"samples",reset,"\n")
            print(names(result$inconsistent_files))
            cat("\n",magenta,"Please verify!",reset,"\n")  
            flush.console()  # 强制刷新控制台输出
            stop()
        }
    
    }else{
        # 2. 非 tar
        data <- read_data(file_path = file_path)
        data_show <- data[1:2,]
        for (col_name in names(data_show)) {
            cat("\n",magenta, "Column:",reset,yellow,col_name, reset,"\n")
            print(data_show[[col_name]])
        }
        cat("\n\n",yellow,"The Supplemently Data columns: ",reset,"\n")
        print(colnames(data))
        flush.console()  # 强制刷新控制台输出
        
        gene_exp_columns <- readline(prompt = "Please input the 'gene ID column' and 'gene expression value columns' (Separated by ‘,’)." )
        gene_exp_type <- readline(prompt = "Please input the gene expression type: (Count/TPM/FPKM/RPM/···)")
        gene_exp_columns_modle <- gene_exp_columns %>% strsplit(.,"[,]") %>% unlist()
        gene_exp_columns_select <- lapply(gene_exp_columns_modle,function(gene_modle){
            gene_exp_columns_select  <- grep(gene_modle,colnames(data),value = T)
            return(gene_exp_columns_select)
        })
        gene_expression_columns_select <- Reduce(function(x, y) c(x,y),gene_exp_columns_select)
        gene_expression_data <- data[,gene_expression_columns_select]
        
        # 保存其它信息                                         
        cat("\n",yellow,"In addition to gene expression data, is there any other information that needs to be preserved?",reset,"\n")
        cat("\n",yellow,"To save other information, enter the column name, otherwise enter nothing.",reset,"\n")
        flush.console()  # 强制刷新控制台输出
        another_columns <- readline(prompt = "Please input the columns name' (Separated by ‘,’)." )
        if(another_columns == ""){
            another_columns_select <- NULL
        }else{
            another_columns_select <- another_columns %>% strsplit(.,"[,]") %>% unlist()
        }
        another_columns_list <- lapply(another_columns_select,function(another_modle){
            another_col_select  <- grep(another_modle,colnames(data),value = T)
            return(another_col_select)
        })
        another_columns_select_list <- Reduce(function(x, y) c(x,y),another_columns_list) %>% unique() 
        gene_another_data <- data[,another_columns_select_list]
        
        if(is.null(another_columns_select)){gene_another_data <- NULL}
                
    }
    return(list(gene_exp_type = gene_exp_type, gene_expression_data = gene_expression_data, gene_another_data = gene_another_data )) 
}
