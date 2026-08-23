#' Select diagnostic features with the Boruta algorithm
#'
#' Run Boruta feature selection for a categorical response, reshape the
#' importance history into a long-format table, and visualize importance
#' distributions by Boruta decision.
#'
#' @param data A data frame containing predictor columns and one categorical
#'   response column.
#' @param feature_colnames Character vector naming predictor features supplied
#'   to Boruta.
#' @param response_colnames Name of the categorical response column. Although
#'   the argument name is plural, the current implementation expects one name.
#' @param pvalue Numeric confidence-level threshold passed to the `pValue`
#'   argument of [Boruta::Boruta()].
#' @param seed Integer random seed set before running Boruta.
#'
#' @return A named list with `result` and `plot`. `result` is a long-format data
#'   frame containing `Iteration`, `Variable`, `Importance`, and `Decision`.
#'   `plot` is a ggplot boxplot of importance Z-scores for predictors and shadow
#'   attributes, colored by decision.
#'
#' @details
#' The input is restricted to the response and requested feature columns. The
#' response is converted to a factor whose levels follow their order of first
#' appearance. Boruta is run with Bonferroni adjustment, at most 100 runs,
#' stored importance history, and [Boruta::getImpRfZ()] importance scores.
#'
#' Importance history is converted to long format, non-finite values are
#' removed, and original Boruta decisions are joined to the variables. Shadow
#' attributes are labeled `"Shadow"`. Variables are ordered by median
#' importance, with confirmed variables placed after all other decision groups.
#'
#' The function calls [Boruta::TentativeRoughFix()] and
#' [Boruta::getSelectedAttributes()], but the current return table and plot use
#' decisions from the original Boruta fit. Consequently, tentatively resolved
#' decisions and the vector of confirmed features are not returned separately.
#'
#' @examples
#' \dontrun{
#' diagnosis_data <- data.frame(
#'   GeneA = rnorm(60),
#'   GeneB = rnorm(60),
#'   GeneC = rnorm(60),
#'   Group = factor(rep(c("Control", "Disease"), each = 30))
#' )
#'
#' boruta_result <- boruta_diagnosis(
#'   data = diagnosis_data,
#'   feature_colnames = c("GeneA", "GeneB", "GeneC"),
#'   response_colnames = "Group",
#'   pvalue = 0.01,
#'   seed = 123
#' )
#'
#' head(boruta_result$result)
#' boruta_result$plot
#' }
#'
#' @seealso [Boruta::Boruta()], [Boruta::TentativeRoughFix()],
#'   [Boruta::getSelectedAttributes()]
#' @export
boruta_diagnosis <- function(data = df_train,feature_colnames = gene_name, response_colnames = "Group",pvalue = 0.01,seed = 1){
    
    suppressPackageStartupMessages({
        library(readr)
        library(Boruta)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
    })

    # 设置数据
    data <- data %>% select(response_colnames,feature_colnames)
    data[[response_colnames]] <- factor(data[[response_colnames]],levels = unique(data[[response_colnames]]))

    set.seed(seed)
    Var.Selec <- Boruta(as.formula(paste(response_colnames, "~ .")),   #构建公式
                        data = data,  #特征筛选所用数据集
                        pValue = pvalue,  #置信水平，默认值0.01
                        mcAdj = TRUE,  #使用Bonferroni方法的多重比较调整
                        maxRuns = 100,   #最大迭代次
                        doTrace = 0,  #运行报告冗长等级
                        holdHistory = TRUE,   #储存重要历史记录
                        getImp = getImpRfZ  #获取属性重要性
                       )   
    final <- TentativeRoughFix(Var.Selec)
    getSelectedAttributes(final, withTentative = FALSE)
    
    # Step 1: 提取长格式
    boruta_long <- as.data.frame(Var.Selec$ImpHistory) %>%
      tibble::rownames_to_column("Iteration") %>%
      pivot_longer(cols = -Iteration, names_to = "Variable", values_to = "Importance") %>%
      drop_na() %>%
      filter(is.finite(Importance))  # 清除非数值
    
    # Step 2: 添加决策标签（含 Shadow）
    decision_map <- data.frame(
      Variable = names(Var.Selec$finalDecision),
      Decision = as.character(Var.Selec$finalDecision)
    )
    
    boruta_long <- left_join(boruta_long, decision_map, by = "Variable")
    boruta_long$Decision <- as.character(boruta_long$Decision)
    boruta_long$Decision[grepl("^shadow", boruta_long$Variable)] <- "Shadow"
    
    # Step 3: 再次设为因子并排序
    boruta_long$Decision <- factor(boruta_long$Decision, levels = c("Confirmed", "Tentative", "Rejected", "Shadow"))
    boruta_long1 <- boruta_long %>% dplyr::filter(Decision != "Confirmed")
    boruta_long2 <- boruta_long %>% dplyr::filter(Decision == "Confirmed")
    
    avg_order1 <- boruta_long1 %>%
      group_by(Variable) %>%
      summarise(meanImp = median(Importance, na.rm = TRUE)) %>%
      arrange(meanImp) %>%
      pull(Variable)
    
    avg_order2 <- boruta_long2 %>%
      group_by(Variable) %>%
      summarise(meanImp = median(Importance, na.rm = TRUE)) %>%
      arrange(meanImp) %>%
      pull(Variable)    

    boruta_long$Variable <- factor(boruta_long$Variable, levels = c(avg_order1,avg_order2))

    boruta_genes <- boruta_long$Variable[boruta_long$Decision == "Confirmed"] %>% unique()   
    
    # Step 4: 颜色映射
    decision_colors <- c(
      "Confirmed" = "#868686FF",   
      "Tentative" = "#EFC000FF",  
      "Rejected" = "#CD534CFF",
      "Shadow" = "#0073C2FF"
    )

    # Step 5: 绘图
   p <-  ggplot(boruta_long, aes(x = Variable, y = Importance, color = Decision)) +
      geom_boxplot(outlier.size = 1, lwd = 0.3) +
      scale_color_manual(values = decision_colors, na.value = "grey80") +
      labs(title = "Boruta Feature Importance",
           x = "",
           y = "Importance (Z-score)"
          ) +
      theme_classic(base_size = 16) +
      theme(
            axis.text.x = element_text(angle = 60, hjust = 1,size = 15,face = "bold"),
            axis.text.y = element_text(size = 15),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
            axis.title  =   element_text(size = 15,face = "bold"),
            legend.title = element_text(size = 16,face = "bold"),
            legend.text = element_text(size = 15),
            plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "mm") 
           )
    
    return(list(result = boruta_long,plot = p))
}
