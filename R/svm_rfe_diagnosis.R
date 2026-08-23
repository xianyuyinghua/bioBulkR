#' Perform SVM recursive feature elimination for diagnosis
#'
#' Use caret recursive feature elimination with cross-validated SVM models,
#' summarize class-specific feature importance, and plot accuracy against the
#' number of selected features.
#'
#' @param data A data frame containing predictor features and one categorical
#'   response column. For the current subset-size calculation, it should contain
#'   only `feature_colnames` plus `response_colnames`.
#' @param feature_colnames Character vector naming predictor columns.
#' @param response_colnames Name of the categorical response column. The current
#'   importance summary expects a binary response.
#' @param n_folds Number of cross-validation folds used by
#'   [caret::rfeControl()].
#' @param seed Integer seed used for RFE and parallel resampling seeds.
#' @param method Caret model identifier passed to [caret::rfe()], such as the
#'   default `"svmRadial"`.
#' @param n_cores Number of PSOCK workers used through doParallel.
#' @param feature_colors Color vector associated with `feature_colnames` in the
#'   accuracy plot.
#'
#' @return A named list with `result`, the caret RFE object augmented with an
#'   `rfe_variables` feature-importance summary and a `Feature` column in its
#'   results table, and `plot`, the cross-validated accuracy ggplot.
#'
#' @details
#' The response is converted to a factor and the requested predictors are
#' passed to caret RFE with `caretFuncs`, cross-validation, and parallel
#' execution. Candidate subset sizes are currently computed as
#' `c(1:ncol(data) - 1)`, producing zero through `ncol(data) - 1`; this assumes
#' `data` consists only of the feature and response columns.
#'
#' Control-like response labels (`control`, `normal`, `norm`, or `low`, ignoring
#' case) are placed first; otherwise labels are alphabetically ordered. The RFE
#' variable table is summarized into mean importance for the first class,
#' second class, and overall importance, then sorted by overall importance.
#'
#' A PSOCK cluster is created and registered with doParallel. On successful
#' completion it is stopped and the sequential backend is restored. If RFE
#' raises an error before cleanup, the current implementation does not use an
#' on-exit handler to guarantee cluster shutdown.
#'
#' @examples
#' \dontrun{
#' svm_data <- data.frame(
#'   GeneA = rnorm(60),
#'   GeneB = rnorm(60),
#'   GeneC = rnorm(60),
#'   Group = factor(rep(c("Control", "Disease"), each = 30))
#' )
#' svm_result <- svm_rfe_diagnosis(
#'   data = svm_data,
#'   feature_colnames = c("GeneA", "GeneB", "GeneC"),
#'   response_colnames = "Group",
#'   n_folds = 5,
#'   seed = 123,
#'   method = "svmRadial",
#'   n_cores = 2
#' )
#' svm_result$result$optVariables
#' svm_result$plot
#' }
#'
#' @seealso [caret::rfe()], [caret::rfeControl()], [doParallel::registerDoParallel()]
#' @export
svm_rfe_diagnosis <- function(data,feature_colnames, response_colnames, n_folds = 5, seed = 1, method = "svmRadial", n_cores = 4,feature_colors = basicR::get_colors(number = 100.1)) {
    
    suppressPackageStartupMessages({
        require(parallel)
        require(doParallel)
        require(caret)
        require(e1071)
    })

    make_seeds <- function(n_folds, sizes, seed = 123) {
        set.seed(seed)
        seeds <- vector("list", n_folds + 1)
        for (i in 1:n_folds) {
        # 每个折对应长度 = sizes 数量
        seeds[[i]] <- sample.int(1000000, length(sizes))
        }
        # 最后一个元素是单个整数
        seeds[[n_folds + 1]] <- sample.int(1000000, 1)
        return(seeds)
    }
    
    group_new_levels <- function(vector,group_levels = NULL){
        # 设置分组水平
        
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

    # colors
    feature_colors_use <- feature_colors[1:length(feature_colnames)]
    feature_colors_use <- feature_colors_use %>% setNames(feature_colnames)
    
    
    # 准备数据
    response <- as.factor(data[[response_colnames]])  # 响应变量
    features <- data %>% dplyr::select(all_of(feature_colnames)) 
    sizes <- c(1:ncol(data)-1)
    
    # 创建并行集群
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    # 设定随机种子（可复现）
    set.seed(seed)

    if (n_cores > 1) {
        seeds_obj <- make_seeds(n_folds = n_folds, sizes = sizes, seed = seed)
    } else {
        seeds_obj <- NA
    }
    
    # 设置 rfe 控制参数
    rfe_control <- rfeControl(
        functions = caretFuncs,
        method = "cv",
        number = n_folds,
        seeds  = seeds_obj,
        allowParallel = TRUE
    )

    # 执行 RFE
    rfe_results <- suppressMessages(suppressWarnings(
        rfe(
          x = features, # features
          y = response,  # response
          sizes = sizes,  # seq(1, ncol(features), by = 1)
          rfeControl = rfe_control,
          method = method
        )
    ))

    stopCluster(cl)
    registerDoSEQ()
    
    rfe_results$rfe_variables <- rfe_results$variables %>%
        group_by(var) %>%
        summarise(
            meanControl = mean(!!sym(group_new_levels(as.character(data[[response_colnames]]))[1]), na.rm = TRUE),
            meanCase = mean(!!sym(group_new_levels(as.character(data[[response_colnames]]))[2]), na.rm = TRUE),
            meanOverall = mean(Overall, na.rm = TRUE)
        ) %>%
        arrange(desc(meanOverall)) 

    rfe_results$results$Feature <- c(NA,rfe_results$rfe_variables$var)
    rfe_results$results$Feature <- factor(rfe_results$results$Feature,levels = unique(rfe_results$results$Feature))

    df_plot <- rfe_results$results %>% dplyr::filter(Variables > 0 )
    
    # plot
    ggplot(df_plot, aes(x = Variables, y = Accuracy)) +
        geom_line(color = "#8FC9E2",linewidth  =1.5) +
        geom_point(aes(color = Feature), size = 5) +
        scale_color_manual(values = feature_colors_use[unique(df_plot$Feature)]) +
        labs( title = "SVM–RFE Analysis",
             x = "Number of Selected Features",
             y = "Accuracy (Cross-Validation)"
            ) +
        theme_classic(base_size = 16) +
        theme(
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 16,face = "bold"),
            legend.text = element_text(size = 15),  
            plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm")
        )-> p
    
    return(list(result = rfe_results,plot = p))
}
