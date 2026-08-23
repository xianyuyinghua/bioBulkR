#' Plot pooled ROC curves for multiple models or features
#'
#' Calculate an ROC curve, AUC, AUC confidence interval, and sensitivity
#' confidence ribbon for each model represented in a prediction table, then
#' combine all model curves in a single ggplot.
#'
#' @param pred_all A data frame of pooled predictions, typically out-of-fold
#'   cross-validation predictions, containing model, observed-class, and
#'   positive-class probability columns.
#' @param features_col Name of the model or feature identifier column.
#' @param observation Name of the observed binary class column.
#' @param probability Name of the numeric predicted probability column for the
#'   positive class.
#' @param pos_class Value in `observation` identifying the positive class.
#' @param digits Number of decimal places used for AUC and confidence interval
#'   values in legend labels.
#' @param sp_grid Numeric specificity grid passed to [pROC::ci.se()].
#' @param title_text Plot title.
#' @param legend_position Legend position passed to ggplot2, such as
#'   `"right"` or a two-element numeric coordinate.
#' @param base_size Base theme text size.
#' @param colors Color vector for model curves and confidence ribbons.
#'
#' @return A ggplot containing one step ROC curve and translucent 95 percent
#' sensitivity confidence ribbon per model, plus a chance diagonal. Legend
#' labels include model name, AUC, and AUC 95 percent confidence interval.
#'
#' @details
#' The observed class is converted to a factor. The negative class is the first
#' factor level other than `pos_class`; the function is therefore intended for
#' binary outcomes. Each model is evaluated with levels ordered as negative,
#' positive and `direction = "<"`, meaning larger probabilities indicate the
#' positive class.
#'
#' AUC confidence intervals are calculated with [pROC::ci.auc()] and
#' sensitivity intervals with [pROC::ci.se()], both receiving `boot.n = 200`.
#' False-positive-rate coordinates are sorted, sensitivities are made
#' monotonically non-decreasing, and explicit `(0, 0)` and `(1, 1)` endpoints
#' are added. Only the combined plot is returned.
#'
#' @examples
#' \dontrun{
#' predictions <- data.frame(
#'   AlgorithmName = rep(c("Logistic", "RandomForest"), each = 100),
#'   obs = rep(rep(c("Control", "Disease"), each = 50), 2),
#'   prob = c(runif(50, 0, 0.5), runif(50, 0.5, 1),
#'            runif(50, 0, 0.6), runif(50, 0.4, 1))
#' )
#' roc_plot <- multimodel_roc_plot(
#'   pred_all = predictions,
#'   features_col = "AlgorithmName",
#'   observation = "obs",
#'   probability = "prob",
#'   pos_class = "Disease"
#' )
#' roc_plot
#' }
#'
#' @seealso [pROC::roc()], [pROC::auc()], [pROC::ci.auc()], [pROC::ci.se()]
#' @export
multimodel_roc_plot <- function(pred_all,
                                features_col = "AlgorithmName",  # 特征列列名
                                observation = "obs",  # 观察列列名（分组）
                                probability = "prob",  # 正类预测概率
                                pos_class,  # 正类名，通常是疾病的组名
                                digits = 3,
                                sp_grid = seq(0, 1, 0.01),
                                title_text = "ROC Curves with 95% CI",
                                legend_position = "right",  # left/right/bottom/top/或者是具体的坐标c(0.5,0.5)
                                base_size = 16,
                                colors = basicR::get_colors(number = 50.1)
                               ){
    
    ### 用于基于多模型（或多特征）的预测结果绘制ROC曲线
    ### 输入必须是一个多个模型（特征）在交叉验证（CV）中得到的“折外预测” 数据据框，并且至少包含以下三类列：模型/特征标识列，真实标签列（分组），正类概率列
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(pROC)
    })

  stopifnot(all(c(features_col, observation, probability) %in% names(pred_all)))

  # 统一obs为factor，并保证pos/neg顺序稳定
  pred_all <- pred_all %>% dplyr::mutate(!!sym(observation) := factor(!!sym(observation)))

  if (!pos_class %in% levels(pred_all[[observation]])){
      stop(paste0("pos_class not found in ",observation,"."))
  }
  neg_class <- setdiff(levels(pred_all[[observation]]), pos_class)[1]

  features <- unique(pred_all[[features_col]])
  
  plot_df <- data.frame()
  ci_ribbon_df <- data.frame()
  label_order <- character()

  for (m in features) {
    dfm <- pred_all %>% dplyr::filter(!!sym(features_col) == m)

    # 方向：prob越大越倾向pos，所以 direction="<" + levels=c(neg,pos) 最常用
    roc_obj <- pROC::roc(
      response = dfm[[observation]],
      predictor = dfm[[probability]],
      levels = c(neg_class, pos_class),
      direction = "<",
      ci = TRUE,
      quiet = TRUE
    )

    auc_val <- sprintf(paste0("%.", digits, "f"), as.numeric(pROC::auc(roc_obj)))

    # AUC 95%CI（bootstrap）
    auc_ci <- pROC::ci.auc(roc_obj, boot.n = 200)
    auc_ci <- sprintf(paste0("%.", digits, "f"), as.numeric(auc_ci))
      
    # 你的gene_auc_label风格：用于图例
    model_auc_label <- paste0(m, "\n(AUC=", auc_val, ", 95%CI: ", auc_ci[1], "-", auc_ci[3], ")")
    label_order <- c(label_order, model_auc_label)

    # ROC坐标
    coords_df <- data.frame(fpr = 1 - roc_obj$specificities,
                            tpr = roc_obj$sensitivities ) %>% 
      arrange(fpr) %>% 
      dplyr::mutate(tpr = cummax(tpr),          # 确保单调
                    model_auc = model_auc_label
                   )

    # 补(0,0)与(1,1)，并去重
    coords_df <- bind_rows(data.frame(fpr = 0, tpr = 0, model_auc = model_auc_label),
                           coords_df,
                           data.frame(fpr = 1, tpr = 1, model_auc = model_auc_label)
                          ) %>%
      distinct(fpr, tpr, model_auc, .keep_all = TRUE)

    plot_df <- bind_rows(plot_df, coords_df)

    # 计算敏感度CI：ci.se 输入 specificities 网格（0..1）
    # 输出矩阵行名是specificity，列为low/median/high
    ci_obj <- pROC::ci.se(roc_obj,specificities = sp_grid,boot.n = 200)

    ci_df <- data.frame(fpr = 1 - as.numeric(rownames(ci_obj)),  # rownames是specificities
                        tpr_low = ci_obj[, 1],
                        tpr_med = ci_obj[, 2],
                        tpr_high = ci_obj[, 3],
                        model_auc = model_auc_label
                       ) %>% arrange(fpr)

    ci_ribbon_df <- bind_rows(ci_ribbon_df, ci_df)
  }

  # 保持图例顺序与颜色一致
  plot_df$model_auc <- factor(plot_df$model_auc, levels = label_order)
  ci_ribbon_df$model_auc <- factor(ci_ribbon_df$model_auc, levels = label_order)

  # 绘图
  p <- ggplot(plot_df, aes(x = fpr, y = tpr, color = model_auc)) +
    geom_ribbon(data = ci_ribbon_df,
                aes(x = fpr, ymin = tpr_low, ymax = tpr_high, fill = model_auc),
                alpha = 0.05,
                colour = NA,
                inherit.aes = FALSE
               ) +
    geom_step(linewidth = 1) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors, guide = "none") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    theme_classic(base_size = base_size) +
    labs(title = title_text, x = "1 - Specificity", y = "Sensitivity", color = "", fill = "") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.background = element_blank(),   # 图例整体背景（透明）
      legend.box.background = element_blank(),# 图例外框背景（透明）
      legend.key = element_blank(),           # 每个图例key的小方块背景（透明）  
      legend.position = legend_position
    ) +
    guides(fill = "none")  # 只保留color图例；如你想fill也显示就删掉这行

  return(p)
}
