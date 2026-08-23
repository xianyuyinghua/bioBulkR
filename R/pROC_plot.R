#' Plot ROC curves and confidence bands for candidate features
#'
#' Calculate a pROC curve, AUC, AUC confidence interval, and sensitivity
#' confidence band for each requested feature, then combine the curves in a
#' single ggplot with AUC information in the legend.
#'
#' @param data A data frame containing a binary group column and numeric
#'   predictor columns.
#' @param group_col Name of the binary outcome or group column.
#' @param genes Character vector naming predictor columns for which ROC curves
#'   are calculated.
#' @param data_type Character label appended to the plot title, such as
#'   `"Train"` or `"Test"`. An empty string omits the parenthetical label.
#' @param legend_position Legend position passed to ggplot2, commonly a
#'   two-element numeric coordinate or a position name.
#' @param legend_text_size Numeric legend text size.
#' @param digits Number of decimal places used for AUC and confidence interval
#'   values in legend labels.
#'
#' @return A ggplot containing one ROC step curve and translucent 95 percent
#' sensitivity confidence ribbon per feature, plus the chance diagonal.
#'
#' @details
#' The group column is converted to a factor. If a group named `control`,
#' `normal`, `norm`, or `low` is present, ignoring case, its first occurrence is
#' used as the first level; otherwise, levels are sorted alphabetically. The
#' function is intended for exactly two outcome groups.
#'
#' ROC direction is selected independently for each feature by comparing its
#' mean in the first and second groups. AUC confidence intervals are obtained
#' with [pROC::ci.auc()] and sensitivity confidence intervals with
#' [pROC::ci.se()] over specificities from zero to one in increments of 0.01.
#' The latter uses 200 bootstrap replicates. Curve coordinates are ordered by
#' false-positive rate and made monotonically non-decreasing with `cummax()`.
#'
#' @examples
#' \dontrun{
#' roc_data <- data.frame(
#'   Group = factor(rep(c("Control", "Disease"), each = 30)),
#'   GeneA = c(rnorm(30, 0), rnorm(30, 1)),
#'   GeneB = c(rnorm(30, 0), rnorm(30, 0.6))
#' )
#' roc_plot <- pROC_plot(
#'   data = roc_data,
#'   group_col = "Group",
#'   genes = c("GeneA", "GeneB"),
#'   data_type = "Train",
#'   digits = 3
#' )
#' roc_plot
#' }
#'
#' @seealso [pROC::roc()], [pROC::auc()], [pROC::ci.auc()], [pROC::ci.se()]
#' @export
pROC_plot <- function(data,group_col,genes,data_type = "Train",legend_position = c(0.75, 0.2),legend_text_size = 12,digits = 3){
    # 特征的ROC曲线
    
    # 加载必要包
    library(pROC)
    library(ggplot2)
    library(dplyr)

    group_new_levels <- function(vector,group_levels = NULL){
        
        # colors
        red <- "\033[31m"
        green <- "\033[32m"
        yellow <- "\033[33m"
        blue <- "\033[34m"
        magenta <- "\033[35m"
        cyan <- "\033[36m"
        reset <- "\033[0m"
        
        # 获取列名
        vector <- vector %>% as.character() %>% unique()
        if(length(group_levels) == 0){
            if(tolower(vector) %in% c("control", "normal","norm","low","Low") %>% any()){
    
                # 找出小写化后等于 "control" 或 "normal" 的列名
                target_col <- vector[tolower(vector) %in% c("control", "normal","norm","low","Low")]
                
                # 如果存在符合条件的列
                if (length(target_col) > 0) {
                  # 将目标列放在第一位，其余列保持原顺序但去掉目标列
                  new_order <- c(target_col[1], setdiff(vector, target_col[1]))
                  return(new_order)  
                }
            }else{
                return(sort(vector))  
            }
        }else{
            if(all(vector %in% group_levels)){
                return(group_levels)
            }else{
                cat("\n",yellow,"group_levels Not Match the vector!",reset,"\n")
            }
        }
    }

    # 确保分组为 factor
    data[[group_col]] <- factor(data[[group_col]], levels = group_new_levels(vector = as.character(data[[group_col]]) ))

    # 计算每组基因的平均表达，判断方向
    mean_expr <- data %>% dplyr::group_by(!!sym(group_col)) %>% dplyr::summarise(across(all_of(genes), mean, na.rm = TRUE), .groups = "drop")
    
    
    # 初始化空表
    plot_df <- data.frame()
    ci_ribbon_df <- data.frame()  # 用于保存置信区间
    
    # 循环计算每个基因的 ROC
    for (gene in genes) {
        Direction <- ifelse(mean_expr[[gene]][1] > mean_expr[[gene]][2], ">", "<")
        
        roc_obj <- roc(data[[group_col]], data[[gene]],
                     levels = group_new_levels(vector = as.character(data[[group_col]])),
                     direction = Direction,
                     ci = TRUE)
        
        auc_val <- sprintf(paste0("%.", digits, "f"), as.numeric(auc(roc_obj)))

        # 提取 AUC 和其 95% CI
        auc_ci <- pROC::ci.auc(roc_obj, boot.n = 200)  # boot.n 可按需调高
        auc_ci <- sprintf(paste0("%.", digits, "f"), as.numeric(auc_ci))
        gene_auc_label <- paste0(gene, "\n(AUC=", auc_val, ", 95%CI: ", auc_ci[1], "-", auc_ci[3], ")")
        
        # 坐标提取
        coords_df <- data.frame(
            fpr = 1 - roc_obj$specificities,
            tpr = roc_obj$sensitivities
        ) %>%
            dplyr::arrange(fpr) %>%
            dplyr::mutate(
              tpr = cummax(tpr),
              gene_auc = gene_auc_label
        )
    
      coords_df <- bind_rows(
        data.frame(fpr = 0, tpr = 0, gene_auc = unique(coords_df$gene_auc)),
        coords_df,
        data.frame(fpr = 1, tpr = 1, gene_auc = unique(coords_df$gene_auc))
      ) %>% distinct(fpr, tpr, gene_auc, .keep_all = TRUE)
    
      plot_df <- bind_rows(plot_df, coords_df)


        
      # 提取 CI（按 sensitivity 的 CI，固定横轴 fpr）
      ci_obj <- pROC::ci.se(roc_obj, specificities = seq(0, 1, 0.01), boot.n = 200)
      ci_df <- data.frame(
        fpr = 1 - as.numeric(rownames(ci_obj)),
        tpr_low = ci_obj[, 1],
        tpr_med = ci_obj[, 2],
        tpr_high = ci_obj[, 3],
        gene_auc = unique(coords_df$gene_auc)
      )
      ci_ribbon_df <- bind_rows(ci_ribbon_df, ci_df)
    } 
    
    plot_df$gene_auc <- factor(plot_df$gene_auc,levels = unique(plot_df$gene_auc))
    ci_ribbon_df$gene_auc <- factor(ci_ribbon_df$gene_auc,levels = unique(plot_df$gene_auc))
    
    if(data_type != ""){
        title_text <- paste0("ROC Curves of Candidate Features \n(", data_type, ")")
    }else{
        title_text <- "ROC Curves of Candidate Features"
    }
    
    # 绘图
    ggplot(plot_df, aes(x = fpr, y = tpr, color = gene_auc)) +
      geom_ribbon(data = ci_ribbon_df,
                  aes(x = fpr, ymin = tpr_low, ymax = tpr_high, fill = gene_auc, group = gene_auc),
                  alpha = 0.1,
                  inherit.aes = FALSE
                 ) +
      geom_step(size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey", size = 1) +
      scale_color_manual(values = basicR::get_colors(number = 50.1)) +
      scale_fill_manual(values = basicR::get_colors(number = 50.1), guide = "none") +
      theme_classic(base_size = 16) +
      labs(title = title_text,  x = "1 - Specificity", y = "Sensitivity", color = "", fill = "") +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text_size),
        plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15, face = "bold"),
        plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm"),
        legend.background = element_blank(),   # 图例整体背景（透明）
        legend.box.background = element_blank(),# 图例外框背景（透明）
        legend.key = element_blank(),           # 每个图例key的小方块背景（透明） 
        legend.position = legend_position
      ) -> p
    
    return(p)

}
