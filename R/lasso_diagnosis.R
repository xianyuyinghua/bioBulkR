#' Fit a cross-validated LASSO diagnostic model
#'
#' Fit a binomial model with cross-validated regularization, extract non-zero
#' coefficients at a selected lambda, and draw coefficient-path and
#' cross-validation diagnostic plots.
#'
#' @param data A data frame containing predictor columns and one categorical
#'   response column.
#' @param feature_colnames Character vector naming the predictor columns.
#' @param response_colnames Name of the response column. Although the argument
#'   name is plural, the current implementation expects one column name.
#' @param n_folds Number of folds passed to [glmnet::cv.glmnet()].
#' @param seed Integer random seed used before cross-validation. The default
#'   expression refers to an object named `seed`, so this argument should
#'   normally be supplied explicitly.
#' @param alpha Elastic-net mixing parameter passed to [glmnet::cv.glmnet()].
#'   The default `1` gives LASSO cross-validation.
#' @param lambda_type Character string selecting the cross-validated lambda;
#'   normally `"lambda.1se"` or `"lambda.min"`.
#' @param feature_colors Color vector used for coefficient paths and numbers of
#'   non-zero variables in the two plots.
#'
#' @return A named list with `result` and `plot`. `result` contains `cvfit`, the
#'   cross-validated glmnet fit; `lasso_gene_df`, a data frame of selected
#'   features and their non-zero coefficients; and `lambda_type`. `plot`
#'   contains `p_path`, the coefficient-path plot, and `p_curve`, the
#'   cross-validation deviance plot.
#'
#' @details
#' A model matrix is constructed from the selected response and feature
#' columns, with its intercept column removed. The response is converted to a
#' factor and a binomial cross-validated glmnet model is fitted using deviance
#' as the assessment measure. The selected lambda value is printed together
#' with its logarithm.
#'
#' `alpha` is passed to the cross-validation fit. In the current implementation
#' the separate full coefficient-path fit calls [glmnet::glmnet()] without an
#' explicit `alpha`, so that path uses glmnet's default `alpha = 1`. Thus, when
#' `alpha` is not `1`, the plotted path and cross-validated fit may represent
#' different mixing parameters.
#'
#' The response should contain two classes, predictors should be usable in a
#' numeric model matrix, and the number of observations should support the
#' requested folds. The supplied color vector should be long enough for all
#' coefficient curves and non-zero-count categories displayed in the plots.
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
#' lasso_result <- lasso_diagnosis(
#'   data = diagnosis_data,
#'   feature_colnames = c("GeneA", "GeneB", "GeneC"),
#'   response_colnames = "Group",
#'   n_folds = 5,
#'   seed = 123,
#'   alpha = 1,
#'   lambda_type = "lambda.1se"
#' )
#'
#' lasso_result$result$lasso_gene_df
#' lasso_result$plot$p_path
#' lasso_result$plot$p_curve
#' }
#'
#' @seealso [glmnet::cv.glmnet()], [glmnet::glmnet()], [glmnet::coef.glmnet()]
#' @export
lasso_diagnosis <- function(data = df_train,feature_colnames, response_colnames,n_folds = 5,seed = seed,alpha = 1,
                            lambda_type = "lambda.1se", # lambda.1se、lambda.min
                            feature_colors = basicR::get_colors(number = 100.1)
                           ){
    # LASSO
    suppressPackageStartupMessages({
        library(dplyr)
        library(rstatix)   # 提供 t_test / wilcox_test 等
        library(tibble)
        library(knitr)
        library(tidyr)
        library(DescTools)
        library(glmnet) 
    })
    
    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 构建数据矩阵
    features <- model.matrix(as.formula(paste(response_colnames, "~ .")),data = data[, c(response_colnames, feature_colnames)])[, -1]  # 去掉 Intercept
    response <- as.factor(data[[response_colnames]])  # 响应变量


    # LASSO 回归 + 交叉验证（10折）
    set.seed(seed)
    cvfit <- cv.glmnet(features, # features
                       response, # response
                       family = "binomial", alpha = alpha, nfolds = n_folds, type.measure = "deviance")
    fit <- glmnet(features, # features
                  response, # response
                  family = "binomial")
    
    # 输出系数
    coef.min = coef(cvfit, s = lambda_type) ## lambda.min & lambda.1se 取一个
    coef_df <- as.data.frame(as.matrix(coef.min))
    coef_df <- data.frame(coef_df) %>% tibble::rownames_to_column("Feature")
    coef_df <- coef_df %>% filter(Feature != "(Intercept)") %>% filter(!!sym(lambda_type) != 0)
    colnames(coef_df) <- c('Feature','Coef')
    coef_df$Feature <- gsub("`", "", coef_df$Feature)

    cat("\n",yellow,paste0("LASSO ",lambda_type,": "),reset,cvfit[[lambda_type]],"\n\t")
    flush.console()
    cat("\n",yellow,paste0("Log(LASSO ",lambda_type,"): "),reset,cvfit[[lambda_type]] %>% log(),"\n\t")
    flush.console()
    
    # lasso Gene and coef
    lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
    lasso_genecof <- as.data.frame(coef.min@x)
    lasso_gene<- as.data.frame(cbind(lasso_geneids,lasso_genecof))  
    
    # 轨迹
    x <- as.matrix(coef(fit))  # 确保 x 是矩阵
    tmp <- as.data.frame(as.matrix(x)) 
    tmp$coef <- row.names(tmp) 
    tmp$coef <- gsub("`([^`]*)`", "\\1  ", tmp$coef)
    tmp <- reshape::melt(tmp, id = "coef") 
    tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
    tmp$coef <- gsub('_','-',tmp$coef) 
    tmp$lambda <- fit$lambda[tmp$variable+1] 
    # extract the lambda values 
    x_sub <- x[-1, , drop = FALSE]
    tmp$norm <- apply(abs(x_sub), 2, sum)[tmp$variable + 1]
    tmp <- tmp %>% filter(coef != "(Intercept)")

    # 偏似然曲线图
    xx <- data.frame(lambda=cvfit[["lambda"]],
                     cvm=cvfit[["cvm"]],
                     cvsd=cvfit[["cvsd"]], 
                     cvup=cvfit[["cvup"]],
                     cvlo=cvfit[["cvlo"]],
                     nozezo=cvfit[["nzero"]]) 
    xx$ll<- log(xx$lambda) 
    xx$NZERO<- paste0(xx$nozezo,' vars')
    xx$NZERO <- factor(xx$NZERO, levels = paste0(0:50, " vars"))


    # 轨迹图
    p_path <- ggplot(tmp,aes(log(lambda),value,color = coef)) + 
      # geom_vline(xintercept = log(cvfit$lambda.min),
      #            linewidth=0.8,color='grey60',
      #            alpha=0.8,linetype=2)+
      geom_vline(xintercept = log(cvfit[[lambda_type]]),
             linewidth=0.8,color='grey60',
             alpha=0.8,linetype=2)+
      geom_line(size=1) + 
      labs(x = "Log(Lambda)",
           y = "Coefficients",
           subtitle = paste0(lambda_type, " : ", cvfit[[lambda_type]])
          )+
      guides(col=guide_legend(ncol = ceiling(length(unique(tmp$coef)) / 20)))+
      scale_color_manual(values= feature_colors)+ 
      scale_y_continuous(expand = c(0.01,0.01))+ 
      theme_classic(base_rect_size = 2,base_size = 16)+ 
      theme(panel.grid = element_blank(), 
            axis.title = element_text(size=16,color='black',face = "bold"), 
            axis.text = element_text(size=14,color='black'), 
            legend.title = element_blank(), 
            legend.text = element_text(size=14,color='black'), 
            legend.position = 'right',
            plot.subtitle = element_text(face = "italic",size = 12,hjust = 0),
            plot.margin = ggplot2::margin(t = 5,r = 5,b = 5,l = 5,unit = "mm")
           )
      

    # 偏似然曲线图
    p_curve <- ggplot(xx,aes(ll,cvm,color = NZERO))+ 
          geom_errorbar(aes(x = ll,ymin = cvlo,ymax = cvup),
                        width = 0.05,size = 1)+ 
          # geom_vline(xintercept = log(cvfit$lambda.min),
          #            size = 0.8, color = 'grey60', alpha = 0.8, 
          #            linetype = 2)+
          geom_vline(xintercept = log(cvfit[[lambda_type]]),
                     size = 0.8, color = 'grey60', alpha = 0.8, 
                     linetype = 2)+
          geom_point(size=2)+
          labs(x = "Log(Lambda)",
               y = "Partial Likelihood Deviance",
               subtitle = paste0(lambda_type, " : ", cvfit[[lambda_type]])
              )+
          guides(col=guide_legend(ncol = ceiling(length(unique(xx$NZERO)) / 20)))+
          scale_color_manual(values= feature_colors)+ 
          scale_x_continuous(expand = c(0.05,0.02))+ 
          scale_y_continuous(expand = c(0.02,0.02))+ 
          theme_classic(base_rect_size = 1.5,base_size = 16)+ 
          theme(panel.grid = element_blank(),
                plot.subtitle = element_text(face = "italic",size = 12,hjust = 0),
                axis.title = element_text(size=16,color='black',face = "bold"), 
                axis.text = element_text(size=14,color='black'), 
                legend.title = element_blank(), 
                legend.text = element_text(size=14,color='black'), 
                legend.position = 'right'
          )
          
    return(list(result = list( cvfit = cvfit, lasso_gene_df = coef_df,lambda_type = lambda_type) ,plot = list(p_path = p_path,p_curve = p_curve) ))
    
}
