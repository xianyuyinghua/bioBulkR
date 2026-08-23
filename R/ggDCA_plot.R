#' Perform decision curve analysis for a logistic nomogram and predictors
#'
#' Fit an rms logistic model containing all selected predictors plus a
#' single-predictor model for each feature, calculate decision curves with
#' ggDCA, and draw net benefit across threshold probabilities.
#'
#' @param input_data A data frame containing the binary outcome and all
#'   predictors named in `key_features`.
#' @param group_colname Name of the binary outcome column. The current
#'   implementation converts it with `as.numeric(x) - 1`, so a two-level factor
#'   is expected.
#' @param key_features Character vector naming predictors used for the combined
#'   `Nomogram` model and the individual predictor models.
#' @param plot_type Models to retain in the decision curve: `"all"`,
#'   `"Nomogram"`, one or more feature names, or a combination of these.
#'
#' @return A named list with `result`, the object returned by [ggDCA::dca()],
#'   and `plot`, a ggplot showing threshold probability against net benefit.
#'
#' @details
#' Predictor names containing characters other than letters, numbers, and
#' underscores are cleaned before formulas are built, while model labels retain
#' the original feature names. An [rms::lrm()] model is fitted for all selected
#' features and for each feature separately. These models are passed to
#' [ggDCA::dca()] with their names.
#'
#' If `plot_type` contains feature names, `key_features` is first restricted to
#' their intersection. Consequently, the `Nomogram` model is also rebuilt from
#' that restricted feature set. Unknown requested model names cause the current
#' implementation to call `quit(save = "no", runLast = FALSE)`.
#'
#' The function assigns the processed input data to `input_data_global` in the
#' global environment with `<<-`. This object is not removed on return.
#'
#' @examples
#' \dontrun{
#' dca_data <- data.frame(
#'   Group = factor(rep(c("Control", "Disease"), each = 50)),
#'   GeneA = rnorm(100),
#'   GeneB = rnorm(100)
#' )
#' dca_result <- ggDCA_plot(
#'   input_data = dca_data,
#'   group_colname = "Group",
#'   key_features = c("GeneA", "GeneB"),
#'   plot_type = c("Nomogram", "GeneA", "GeneB")
#' )
#' dca_result$plot
#' dca_result$result
#' }
#'
#' @seealso [ggDCA::dca()], [rms::lrm()]
#' @export
ggDCA_plot <- function(input_data = df_train,
                       group_colname = "Group",
                       key_features = key_genes,
                       plot_type = "all"  # 'all'/'Nomogram'/another features
                      ){
    
    library(dplyr)
    library(rms)

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    if(length(intersect(plot_type,key_features)) >0){
        key_features <- intersect(plot_type,key_features)
    }

    # 确保 Group 是数值型二分类（0/1）
    input_data[[group_colname]] <- as.numeric(input_data[[group_colname]]) - 1
    
    # 检查关键特征中是否有特殊字符
    get_name_map_with_cleaning <- function(genes, pattern = "[^[:alnum:]_]", replacement = "_"){
        # 将有特殊字符的生成替换前后的对应关系
        cleaned <- gsub(pattern, replacement, genes)
        changed <- cleaned != genes
        setNames(genes[changed], cleaned[changed])
    }
    check_result <- get_name_map_with_cleaning(key_features)

    if(length(check_result) > 0){
        # 存在包含特殊字符的特征，替换
        cat("\n",blue,"Features with special characters:",reset,check_result,"\n")
        flush.console()

        # 更改数据的列名
        input_data <- input_data %>% rename_with(~ ifelse(.x %in% check_result,
                                  names(check_result)[match(.x, check_result)],
                                  .x))
        # 更改特征名
        key_features_updated <- ifelse( key_features %in% check_result,
                                       names(check_result)[match(key_features, check_result)],
                                       key_features
                                      )
    }else{
        # 不存在特殊字符的特征
        key_features_updated = key_features
    }

    # 创建 lrm list（支持含特殊字符的列名）
    lrm_list <- lapply(0:length(key_features_updated), function(i) {
        if (i == 0) {
            # 全部变量一起
            foumular <- as.formula(paste0(group_colname," ~ ", paste(sprintf("`%s`", key_features_updated), collapse = " + ")))
        } else {
            # 单个变量
            foumular <- as.formula(paste0(group_colname," ~ `", key_features_updated[i], "`"))
        }
        
        input_data_global <<- input_data
        lrm(foumular, data = input_data_global, x = TRUE, y = TRUE)
    })
    lrm_list <- setNames(lrm_list, c("Nomogram", key_features))

    # 选择画图的特征
    if(length(plot_type)== 1){
        if(plot_type == "all"){
            lrm_list <- lrm_list
        }else{
            if(all(plot_type %in% names(lrm_list) )){
                lrm_list <- lrm_list[names(lrm_list) %in% plot_type]
            }else{
                cat("\n",yellow,setdiff(plot_type,names(lrm_list)),"Not found!",reset,"\n")
                quit(save = "no", runLast = FALSE)
            }            
        }
        
    }else{
        if(all(plot_type %in% names(lrm_list) )){
            lrm_list <- lrm_list[names(lrm_list) %in% plot_type]
        }else{
            cat("\n",yellow,setdiff(plot_type,names(lrm_list)),"Not found!",reset,"\n")
            quit(save = "no", runLast = FALSE)
        }
    }

    dca_lrm <- do.call(ggDCA::dca, c(unname(lrm_list), list(model.names = names(lrm_list))))

    # Step 4: 绘图
    p <- ggplot(dca_lrm,linetype = F) +
      scale_color_manual(values = basicR::get_colors(number = 11.1)) +
      labs(title = "Decision Curve Analysis (Disease Prediction)",
           x = "Threshold Probability",
           y = "Net Benefit") +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        legend.position = "right",
        legend.background = element_blank(),
        axis.title = element_text(face = "bold")  
      )
    
    return(list(result = dca_lrm,plot = p))
}
