#' Tune and fit a random-forest diagnostic classifier
#'
#' Use cross-validation to select the number of trees for a random-forest
#' classifier, fit the final model on all samples, calculate feature importance,
#' and create error and importance plots.
#'
#' @param data A data frame containing predictor columns and one response
#'   column.
#' @param feature_colnames Character vector naming predictor features used by
#'   the classifier.
#' @param response_colnames Name of the categorical response column.
#' @param n_folds Number of cross-validation folds passed to
#'   [caret::createFolds()].
#' @param seed Integer random seed set before fold creation and model fitting.
#'
#' @return A named list with `result` and `plot`. `result` contains
#'   `optimal_ntree`, the tree count with the smallest mean cross-validation
#'   error, and `importance_df`, a decreasing feature-importance table. `plot`
#'   contains `p_error`, the tree-count error curve, and `p_importance`, the
#'   horizontal importance plot.
#'
#' @details
#' The function retains only `feature_colnames` and `response_colnames`, then
#' creates stratified folds from the response. Candidate forest sizes are
#' `seq(1, 500, by = 5)`, which evaluates 1, 6, ..., 496 trees. For every
#' candidate, a model is fitted on each training fold and its held-out
#' misclassification rate is averaged. The first minimum determines the final
#' tree count.
#'
#' The response is converted to a factor, so the final random forest is a
#' classifier. Feature importance is enabled in the final fit; the returned
#' value is the first column of [randomForest::importance()]. Input predictors
#' should be suitable for `randomForest` and should not contain unhandled
#' missing values.
#'
#' @examples
#' \dontrun{
#' diagnosis_data <- data.frame(
#'   GeneA = rnorm(60),
#'   GeneB = rnorm(60),
#'   GeneC = rnorm(60),
#'   Group = factor(rep(c("Control", "Disease"), each = 30))
#' )
#' rf_result <- rf_diagnosis(
#'   data = diagnosis_data,
#'   feature_colnames = c("GeneA", "GeneB", "GeneC"),
#'   response_colnames = "Group",
#'   n_folds = 5,
#'   seed = 123
#' )
#' rf_result$result$optimal_ntree
#' rf_result$plot$p_error
#' rf_result$plot$p_importance
#' }
#'
#' @seealso [randomForest::randomForest()], [randomForest::importance()],
#'   [caret::createFolds()]
#' @export
rf_diagnosis <- function(data = df_train,feature_colnames = gene_name, response_colnames = "Group",n_folds = nfolds.num,seed = 1){
    
    suppressPackageStartupMessages({
        library(randomForest)
        library(caret)
    })
    
    # check_stable_min_error <- function(cv_results, window_size = 50) {
    #     # 获取最小 error 值,判断设定之后连续多少树都是最低点
    #     min_error <- min(cv_results$error)
        
    #     # 找到第一个达到最小值的位置
    #     first_min_index <- which(cv_results$error == min_error)[1]
        
    #     # 检查数据是否足够长
    #     if (first_min_index + window_size - 1 > nrow(cv_results)) {
    #     warning("数据不足以检查连续 ", window_size, " 个值。")
    #     return(NULL)
    #     }
        
    #     # 提取从该点起的连续 error 值
    #     error_window <- cv_results$error[first_min_index:(first_min_index + window_size - 1)]
        
    #     # 判断是否全部等于最小值
    #     all_equal <- all(error_window == min_error)
        
    #     # 返回结果列表
    #     return(list(
    #     min_error = min_error,
    #     first_min_index = first_min_index,
    #     first_min_ntree = cv_results$ntree[first_min_index],
    #     all_equal = all_equal
    #     ))
    # }

    # Data filter
    data <- data %>% dplyr::select(all_of(c(feature_colnames,response_colnames)))

    #--main-----------------------------------------------
    set.seed(seed)
    cv_folds <- createFolds(data[[response_colnames]], k = n_folds, list = TRUE)
    
    # 设定 ntree 搜索范围（从 1 到 500，每次增加 5）
    ntree_range <- seq(1, 500, by = 5)
    cv_error <- numeric(length(ntree_range))
    
    # 遍历不同的 ntree 进行交叉验证
    for (i in seq_along(ntree_range)) {
        ntree <- ntree_range[i]
        fold_errors <- numeric(length(cv_folds))
        
        for (j in seq_along(cv_folds)) {
            fold <- cv_folds[[j]]
            
            # 划分训练集和测试集
            train_data <- data %>% dplyr::slice(-fold) %>% dplyr::select(-response_colnames) 
            train_label <- as.factor(data[[response_colnames]][-fold])
            test_data <- data %>% dplyr::slice(fold) %>% dplyr::select(-response_colnames)
            test_label <- as.factor(data[[response_colnames]][fold])
            
            # 训练随机森林模型
            rf_model <- randomForest(x = train_data, y = train_label, ntree = ntree)
            
            # 预测测试集并计算错误率
            predictions <- predict(rf_model, newdata = test_data)
            fold_errors[j] <- mean(predictions != test_label)
        }
        
        # 计算当前 ntree 的平均交叉验证误差
        cv_error[i] <- mean(fold_errors)
    }
    cv_results <- data.frame(ntree = ntree_range, error = cv_error)

    # 选择最优 ntree（使错误率最小）
    optimal_ntree <- ntree_range[which.min(cv_error)]

    # 训练最终模型
    rf_model <- randomForest(x = data %>% dplyr::select(-response_colnames), y = as.factor(data[[response_colnames]]), ntree = optimal_ntree, importance = TRUE)
    # 获取基因重要性
    importance_df <- data.frame(Gene = rownames(importance(rf_model)),Importance = importance(rf_model)[,1])
    # 按重要性排序
    importance_df <- importance_df %>% dplyr::arrange(desc(Importance))

    # 绘制 ntree vs. 交叉验证误差曲线
    p1 <- ggplot(cv_results, aes(x = ntree, y = error)) +
        geom_line(color = "#3576AD", linewidth = 1) +
        geom_point(color = "red", size = 2) +
        labs(title = "Error vs. Number of Trees (Cross-Validation)",
             x = "Number of Trees",
             y = "Cross-Validation Error Rate") +
        theme_classic()+
        theme(axis.text = element_text(size = 12,face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            plot.title = element_text(hjust = 0.5),
            text = element_text(size = 20, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.margin = ggplot2::margin(t = 5,r = 5,b = 5,l = 5,unit = "mm")
           )

    # 重要度 出图
    p2 <- ggplot(importance_df, aes(x = reorder(Gene, Importance), y = Importance)) +
        geom_point(aes(color = Gene), size = 5) +
        scale_color_manual(values = basicR::get_colors(number = 100.1))+
        coord_flip() +  # 旋转坐标轴，横向显示
        labs(title = "",
           x = "",
           y = "Importance",
           color = ""
          ) +
        theme_classic() +
        theme(axis.text = element_text(size = 15),
            axis.text.y = element_text(size = 15,face = "bold"),
            axis.title = element_text(size = 16,face = "bold"),
            axis.line = element_line(linewidth = 0.7),
            legend.position =  "none",
            plot.margin = ggplot2::margin(t = 5,r = 5,b = 5,l = 5,unit = "mm")
           )

    return(list(result = list(optimal_ntree = optimal_ntree, importance_df = importance_df), plot = list(p_error = p1,p_importance = p2)))
}
