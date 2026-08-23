#' Fit and plot a logistic-regression calibration curve
#'
#' Fit an rms logistic regression from selected predictors, estimate apparent
#' and bias-corrected calibration, perform a Hosmer-Lemeshow goodness-of-fit
#' test, and draw the calibration curves against the identity line.
#'
#' @param input_data A data frame containing the binary outcome and all
#'   predictors named in `key_features`. The current Hosmer-Lemeshow step also
#'   requires a factor column named `Group`.
#' @param group_colname Name of the binary outcome column used in the logistic
#'   model formula.
#' @param key_features Character vector naming predictors included in the
#'   logistic regression.
#' @param calibrate_method Calibration method passed to [rms::calibrate()],
#'   such as `"boot"`, `"crossvalidation"`, or `".632"`.
#' @param calibrate_B Number of bootstrap repetitions or resampling iterations
#'   passed as `B` to [rms::calibrate()].
#' @param number_of_groups Number of groups used by
#'   [ResourceSelection::hoslem.test()].
#' @param line_color A color vector of length at least two for the apparent and
#'   bias-corrected curves.
#'
#' @return A named list with `result`, a long data frame containing `Predicted`,
#'   `Type`, and `Observed` calibration values, and `plot`, the calibration
#'   ggplot annotated with the Hosmer-Lemeshow p-value.
#'
#' @details
#' The logistic model is fitted with [rms::lrm()] using
#' `group_colname ~ key_features`. Calibration output columns `predy`,
#' `calibrated.orig`, and `calibrated.corrected` are converted to a tidy table
#' with `Apparent` and `Bias_corrected` curve types.
#'
#' For the Hosmer-Lemeshow test, the current implementation specifically
#' converts `input_data$Group` with `as.numeric(Group) - 1`, regardless of
#' `group_colname`. Therefore, normal use requires a two-level factor named
#' `Group`; its factor-level order determines which class is coded zero or one.
#' Resampling-based calibration consumes random numbers, so set a seed before
#' calling when reproducibility is required.
#'
#' @examples
#' \dontrun{
#' calibration_data <- data.frame(
#'   Group = factor(rep(c("Control", "Disease"), each = 50)),
#'   GeneA = rnorm(100),
#'   GeneB = rnorm(100)
#' )
#' set.seed(123)
#' calibration <- calibrate_plot(
#'   input_data = calibration_data,
#'   group_colname = "Group",
#'   key_features = c("GeneA", "GeneB"),
#'   calibrate_method = "boot",
#'   calibrate_B = 200,
#'   number_of_groups = 5
#' )
#' calibration$plot
#' calibration$result
#' }
#'
#' @seealso [rms::lrm()], [rms::calibrate()],
#'   [ResourceSelection::hoslem.test()]
#' @export
calibrate_plot <- function(input_data = df_train,
                           group_colname = "Group",
                           key_features = key_genes,
                           calibrate_method = "boot",  # boot/crossvalidation/.632
                           calibrate_B = 1000, # 1000/10/1000
                           number_of_groups = 5,
                           line_color = basicR::get_colors(number = 2.1,package = "ggsci",name = "jco")
                          ){
    
    library(rms)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ResourceSelection)

    foumular <- as.formula( paste0(group_colname," ~ ", paste(sprintf("`%s`", key_features), collapse = " + ")))
    # 拟合 Logistic 回归模型
    fit_logit <- lrm(foumular,  data = input_data, x = TRUE, y = TRUE)

    cal <- calibrate(fit_logit, method = calibrate_method, B = calibrate_B)
    
    # 构造 data.frame
    cal_df <- data.frame(
        Predicted = cal[, "predy"],
        Apparent = cal[, "calibrated.orig"],
        Bias_corrected = cal[, "calibrated.corrected"]
    )
    
    # 转为长数据格式
    cal_long <- cal_df %>%  pivot_longer(cols = c("Apparent", "Bias_corrected"),
                                         names_to = "Type", 
                                         values_to = "Observed"
                                        )

    # p
    # 确保 Group 是数值型二分类（0/1）
    input_data$Group <- as.numeric(input_data$Group) - 1
    
    # 使用 predict 获取概率预测
    probs <- predict(fit_logit, type = "fitted")
    
    # 进行 HL 检验，g = 10 表示分为10组
    hl <- hoslem.test(input_data$Group, probs, g = number_of_groups)
    
    p_text <- paste0("p = ", format.pval(hl$p.value, digits = 3, eps = 0.001)) 
    
    # ggplot 校准曲线
    ggplot(cal_long, aes(x = Predicted, y = Observed, color = Type, linetype = Type)) +
      geom_line(size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1, color = "black") +
      scale_color_manual(values = c("Apparent" = line_color[1], "Bias_corrected" = line_color[2])) +
      scale_linetype_manual(values = c("Apparent" = "solid", "Bias_corrected" = "solid")) +
      coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) + 
      labs(title = "Calibration Curve (Pain Prediction Nomogram)",
           x = "Nomogram-Predicted Probability",
           y = "Actual Probability",
           color = "Calibration Type",
           linetype = "Calibration Type") +
      annotate("text", x = 0.05, y = 0.95, label = p_text, hjust = 0, size = 6)+
      theme_classic(base_size = 16) +
      theme(text = element_text(size = 15),
            plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 14),
            axis.line = element_line(linewidth = 0.7),  
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 15),
            legend.position = c(0.7, 0.2)
      ) -> p_cal

    return(list(result = cal_long,plot = p_cal)) 
}
