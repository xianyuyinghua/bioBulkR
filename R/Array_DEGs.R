#' Identify differentially expressed genes in array data
#'
#' Perform a two-group differential-expression analysis for microarray data
#' with the limma workflow. The function aligns expression samples with their
#' condition records, constructs a no-intercept design matrix, fits the chosen
#' contrast, classifies genes by fold change and statistical significance, and
#' saves the complete result.
#'
#' @param expr A numeric gene-by-sample expression matrix or data frame. Gene
#'   identifiers should be stored in row names and sample identifiers in
#'   column names.
#' @param condition A data frame containing sample group information. It must
#'   have exactly one column named `group` or `groups`, ignoring case. Sample
#'   identifiers may be supplied in exactly one column named `sample` or
#'   `samples`, ignoring case, or in row names.
#' @param group_order A character vector giving the group order used to build
#'   the design and contrast. The comparison is
#'   `group_order[1] - group_order[2]`.
#' @param output_dir A character string giving the directory in which the
#'   differential-expression result is written. The directory is created
#'   recursively when it does not exist.
#' @param logFC A non-negative numeric fold-change threshold on the log2 scale.
#' @param p_value The significance measure used to classify genes. Use
#'   `"p.adj"` for the FDR-adjusted p-value; any other value selects the raw
#'   p-value.
#' @param p A numeric significance threshold. Genes must have a p-value
#'   strictly less than this value to be classified as differentially
#'   expressed.
#' @param top_gene An integer indicating a requested number of top genes. This
#'   argument is retained by the function interface but is not used in the
#'   current calculation.
#' @param change A character vector of three labels used, in order, for
#'   upregulated, downregulated, and non-significant genes.
#'
#' @return A data frame containing the complete [limma::topTable()] result
#'   after rows with missing values have been removed. The p-value columns are
#'   named `p.value` and `p.adj`, and an additional `sig` column contains the
#'   labels from `change`. The same table is written to `01_limma_Res.csv`
#'   under `output_dir`.
#'
#' @details
#' If all expression values are non-negative and the maximum is at least
#' 10,000, the expression matrix is transformed with `log2(expr + 1)`. If all
#' values are non-negative and every sample's 95th percentile is below 25, the
#' values are used unchanged. For other ranges, the function displays the
#' range and asks interactively whether the current data should be used; an
#' answer other than `"yes"` stops the analysis.
#'
#' Only samples shared by the expression column names and condition row names
#' are retained. Differential expression is calculated with [limma::lmFit()],
#' [limma::contrasts.fit()], and [limma::eBayes()]. Multiple-testing correction
#' in [limma::topTable()] uses the FDR method.
#'
#' @examples
#' \dontrun{
#' expr <- matrix(
#'   rnorm(100 * 6, mean = 8),
#'   nrow = 100,
#'   dimnames = list(
#'     paste0("Gene", seq_len(100)),
#'     paste0("Sample", seq_len(6))
#'   )
#' )
#' condition <- data.frame(
#'   Group = rep(c("Control", "Disease"), each = 3),
#'   row.names = colnames(expr)
#' )
#' deg_result <- Array_DEGs(
#'   expr = expr,
#'   condition = condition,
#'   group_order = c("Control", "Disease"),
#'   output_dir = tempfile("array_degs_"),
#'   logFC = 0.5,
#'   p_value = "p.adj",
#'   p = 0.05
#' )
#' }
#'
#' @seealso [limma::lmFit()], [limma::makeContrasts()], [limma::topTable()]
#' @export
Array_DEGs  <- function(expr = data$expr, 
                        condition = data$condition, 
                        group_order = data$group_order,
                        output_dir = "02_DEGs",
                        logFC = 0.5, 
                        p_value = "p.adj",
                        p = 0.05,
                        top_gene = 10,
                        change = c("Up","Down","Not")){

        # colors
        red <- "\033[31m"
        green <- "\033[32m"
        yellow <- "\033[33m"
        blue <- "\033[34m"
        magenta <- "\033[35m"
        cyan <- "\033[36m"
        reset <- "\033[0m"
        
        # 芯片数据差异分析
        if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
        
        if(p_value == "p.adj"){
            pvalue <- "p.adj"
        }else{
            pvalue <- "p.value"
        }
    
        if( !("data.frame" %in%  class(expr))){
            expr <- as.data.frame(expr,check.names = FALSE)
        }
        
      # 检查表达值 
      if( min(expr) >= 0 & max(expr) >= 10000 ){
        expr <- log2(expr + 1)
      }else if( min(expr) >= 0 & all(sapply(expr, function(x) quantile(x, 0.95)) < 25) ){
        expr <- expr
      }else{
        cat("\n",red,"Please verify The Expression Data!\n\nThe Expression value range: ",reset,"\n")
        cat(paste0("\t",range(expr)))
        flush.console()
        use_verify <- readline(prompt = "Do you want to use the current data? (yes/no): ")  
        if(use_verify != "yes"){
           stop()
        } 
        
      }
    
      # 过滤低表达基因：保留在至少 30% 样本中表达量 >5 的基因
       #expr <- expr[rowSums(expr > 5) >= floor(ncol(expr) * 0.3), ]   
                                            
      # 检查condition格式
      group_col_name <- colnames(condition)[grepl("^group$|^groups$",tolower(colnames(condition)))]
      sample_col_name <- colnames(condition)[grepl("^sample$|^samples$",tolower(colnames(condition)))]                                      
      if(length(group_col_name) == 1){
        condition <- condition %>% rename("Group" = !!sym(group_col_name))  # 修改列名为 Group
      }else{
        cat("\n","\033[33m","Please set the column name 'condition' to 'group/groups (case-insensitive)', and there should be only one column name that matches this criterion.","\033[0m","\n")
        print(head(condition))
        stop()    
      }
      if(length(sample_col_name) == 1){
        condition <- condition %>% rename("sample" = !!sym(sample_col_name))  # 修改列名为 sample
      }else if(length(sample_col_name) == 0){
        if( rownames(condition) %in% colnames(expr) %>% any() ){
          condition <- condition
        }else{
          cat("\n","\033[33m","Sample information not found.","\033[0m","\n")
          print(head(condition))
          stop()
        }
      }else{
        cat("\n","\033[33m","Please set the column name 'condition' to 'sample/samples (case-insensitive)', and there should be only one column name that matches this criterion.","\033[0m","\n")
        print(head(condition))
        stop() 
      }
      if(ncol(condition) != 1){
        # 将sample列放置到行名
        condition <- condition %>% select(c("sample","Group"))
        rownames(condition) <- condition$sample
        condition$sample <- NULL
      }
      
      # 根据condition和expr筛选共同样本
      sample_intersect <- intersect(colnames(expr),rownames(condition))
      expr <- expr %>% dplyr::select(all_of(sample_intersect))
      condition <- condition %>% filter(rownames(condition) %in% sample_intersect)
                                            
      # 按照指定的分组顺序排序分组信息样本和表达数据样本
      condition <- condition %>% arrange(Group)
      expr <- expr %>% dplyr::select(rownames(condition))
                                            
      # 设置比较顺序                                      
      list <- c(condition$Group) %>% factor(., levels = group_order, ordered = F)
      cat("\n","\033[33m","The group level: ","\033[0m","\n")
      print(list)
      
      list <- model.matrix(~factor(list)+0)  #把Group设置成一个model matrix
      colnames(list) <- group_order
      
      # 差异分析
      df.fit <- lmFit(expr, list)  ## 数据与list进行匹配
      df.matrix <- makeContrasts(contrasts = paste0(group_order[1], "-", group_order[2]), levels = list)
      fit <- contrasts.fit(df.fit, df.matrix)
      fit <- eBayes(fit)
      Res <- topTable(fit,n = Inf, adjust = "fdr") %>% na.omit() ##所有的差异结果
      Res <- Res %>% dplyr::rename( "p.value" = "P.Value" , 'p.adj' = 'adj.P.Val')
      
      # 标签
      Res$sig <- ifelse(Res[[pvalue]] < p & Res$logFC >  logFC,  change[1],  # up
                        ifelse(Res[[pvalue]] < p & Res$logFC < -logFC, 
                               change[2], # down
                               change[3]  # Not
                              )
                       )
      cat("\n","\033[33m","Difference conditions: logFC >",logFC,"&",pvalue,"< ",p,"\033[0m")
      cat("\n","\033[33m","Number of differential genes：","\033[0m","\n")
      print(table(Res$sig))
      flush.console()  # 强制刷新控制台输出
      
      write.csv(Res,file.path(output_dir,"01_limma_Res.csv"))
      return(Res)
}
