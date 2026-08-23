#' Identify differentially expressed genes from high-throughput count data
#'
#' Filter low-expression genes, perform a two-group DESeq2 analysis, classify
#' genes by fold change and statistical significance, and save both the full
#' differential-expression result and a variance-stabilized expression matrix.
#'
#' @param expr A gene-by-sample matrix or data frame of raw, non-negative
#'   integer counts. Gene identifiers should be stored in row names and sample
#'   identifiers in column names.
#' @param condition A data frame containing a `Group` column. Samples must be
#'   identified either by row names or by a `sample` column. An optional
#'   `Raw_Group` column is removed before the DESeq2 object is constructed.
#' @param group_order A character vector defining group factor levels. The
#'   differential contrast is `group_order[1]` versus `group_order[2]`.
#' @param Exp_output_dir Directory in which the variance-stabilized expression
#'   matrix is written. It must already exist.
#' @param DEG_output_dir Directory in which the differential-expression result
#'   is written. It is created recursively when necessary.
#' @param logFC A non-negative absolute log2 fold-change threshold.
#' @param p_value The significance measure used for classification. Use
#'   `"p.adj"` for adjusted p-values; any other value selects raw p-values.
#' @param expr_cutoff A non-negative count threshold used to filter genes.
#' @param filter_sample_num The minimum sample-count criterion for filtering.
#'   A gene is retained only when the number of qualifying samples is strictly
#'   greater than this value.
#' @param p A numeric significance threshold.
#' @param change Three labels used, in order, for upregulated, downregulated,
#'   and non-significant genes.
#'
#' @return A data frame containing DESeq2 statistics, log2 normalized counts
#'   for each sample, renamed `logFC`, `p.value`, and `p.adj` columns, and a
#'   factor column `sig` containing the labels from `change`. Rows with missing
#'   result values are removed.
#'
#' @details
#' When `expr_cutoff = 0`, genes are retained when counts are greater than zero
#' in more than `filter_sample_num` samples. For a nonzero cutoff, counts must
#' be greater than or equal to `expr_cutoff`. DESeq2 is run with design
#' `~ Group`, independent filtering is disabled in [DESeq2::results()], and
#' the contrast compares the first two values of `group_order`.
#'
#' The returned table and the file `01_Deseq2_Res.csv` contain log2 normalized
#' counts merged with the DESeq2 result. The separate file
#' `01_Expression_throughput_vst_Norm.csv` contains the matrix produced by
#' [DESeq2::vst()] with `blind = FALSE`.
#'
#' The historical function name contains the spelling `Throughtput` and is
#' retained for compatibility.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(
#'   rnbinom(100 * 6, mu = 100, size = 2),
#'   nrow = 100,
#'   dimnames = list(
#'     paste0("Gene", seq_len(100)),
#'     paste0("Sample", seq_len(6))
#'   )
#' )
#' condition <- data.frame(
#'   Group = rep(c("Control", "Disease"), each = 3),
#'   row.names = colnames(counts)
#' )
#'
#' result <- High_Throughtput_Count_DEG(
#'   expr = counts,
#'   condition = condition,
#'   group_order = c("Control", "Disease"),
#'   Exp_output_dir = "01_Expression_data",
#'   DEG_output_dir = "02_DEGs",
#'   logFC = 0.5,
#'   p_value = "p.adj",
#'   expr_cutoff = 10,
#'   filter_sample_num = 2,
#'   p = 0.05
#' )
#' }
#'
#' @seealso [DESeq2::DESeqDataSetFromMatrix()], [DESeq2::DESeq()],
#'   [DESeq2::results()], [DESeq2::vst()]
#' @export
High_Throughtput_Count_DEG <- function(expr = data$expr,                        
                                       condition = data$condition, 
                                       group_order = data$group_order,
                                       Exp_output_dir = "01_Expression_data",
                                       DEG_output_dir = "02_DEGs",
                                       logFC = 0.5,
                                       p_value = "p.adj",
                                       expr_cutoff = 10,
                                       filter_sample_num = nrow(data$condition) * 0.85,  # 当阈值为0时，设置0.5
                                       p = 0.05, 
                                       change = c("Up","Down","Not")
                                      ){
      # 高通量数据差异分析
      if(!dir.exists(DEG_output_dir)){dir.create(DEG_output_dir,recursive = T)}
      if(p_value == "p.adj"){
        pvalue <- "p.adj"
      }else{
        pvalue <- "p.value"
      }
      
      # DEseq2
      if(any(grepl("Raw_Group",colnames(condition)))){condition$Raw_Group <- NULL}
      condition$Group <- factor(condition$Group,levels = group_order)
      
      # 按照分组顺序排序 分组样本顺序和表达数据样本顺序  
      condition <- condition %>% arrange(Group)
      if("sample" %in% colnames(condition)) {
          expr <- expr %>% dplyr::select(condition$sample)
      }else{
          expr <- expr %>% dplyr::select(rownames(condition))
      } 
      
    
      dds <- DESeqDataSetFromMatrix(countData = expr, # 表达矩阵
                                    colData = condition,        # 表达矩阵列名和分组信息的对应关系
                                    design = ~ Group)         # group为colData中的group，也就是分组信息 
      # 低表达基因过滤
      if(expr_cutoff == 0){
          filter_dds <- rowSums(counts(dds) > expr_cutoff) > filter_sample_num
      }else{
          filter_dds <- rowSums(counts(dds) >= expr_cutoff) > filter_sample_num
      }
      
      
    
      dds <- dds[filter_dds,]  
      dds <- DESeq(dds)
      res = results(dds, independentFiltering = FALSE, c("Group",group_order[1],group_order[2]), alpha = 0.1)
      cat("\n",yellow,"The Summary Res: ",reset,"\n")
      print(summary(res)) 
      flush.console()  # 强制刷新控制台输出
      
      resdata <- merge(as.data.frame(res), as.data.frame(log2(counts(dds,normalized=TRUE) + 1) ), by="row.names", sort=FALSE) 
      rownames(resdata) <- resdata$Row.names
      resdata$Row.names <- NULL
      #norm_exp <- as.data.frame(counts(dds,normalized=TRUE)) 
      norm_exp <- as.data.frame(assay(vst(dds, blind = FALSE)))
      resOrdered <- resdata[order(resdata$padj), ]
      RES <- as.data.frame(resOrdered) %>% na.omit(DEG) 
      RES <- RES %>% dplyr::rename('logFC' = 'log2FoldChange', 'p.adj' = 'padj', 'p.value' = 'pvalue')
      RES$sig = as.factor(ifelse(RES[[pvalue]] < p & abs(RES$logFC) > logFC,
                                 ifelse(RES$logFC > logFC,
                                        change[1],
                                        change[2]),
                                 change[3]))
      write.csv(RES, file = file.path(DEG_output_dir,"01_Deseq2_Res.csv"))
      write.csv(norm_exp, file = file.path(Exp_output_dir,"01_Expression_throughput_vst_Norm.csv"))
    
      cat("\n","\033[33m","Difference conditions: logFC >",logFC,"&",pvalue,"< ",p,"\033[0m")
      cat("\n","\033[33m","Number of differential genes：","\033[0m","\n")
      print(table(RES$sig))
      flush.console()  # 强制刷新控制台输出
      return(RES)
} 
