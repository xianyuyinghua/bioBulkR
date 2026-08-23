#' Fit a logistic model and draw a regplot nomogram
#'
#' Fit an rms logistic regression of `Group` on selected features, draw a
#' regplot nomogram for a selected observation, and handle feature names that
#' contain characters unsuitable for model formulas.
#'
#' @param data A data frame containing an outcome column named `Group` and all
#'   columns listed in `key_features`.
#' @param key_features A character vector naming predictors included in the
#'   logistic regression model.
#' @param showP Logical; pass `TRUE` to display predictor p-values in the
#'   regplot output.
#' @param sample_select A row index selecting the observation highlighted by
#'   droplines in the nomogram.
#' @param formats A character vector of requested output formats. This argument
#'   is retained by the interface but is not used by the current implementation.
#'
#' @return If no selected feature names contain special characters, the
#'   recorded nomogram plot returned by [grDevices::recordPlot()]. If names are
#'   cleaned, a magick image object read from the corrected PDF is returned.
#'
#' @details
#' Characters other than letters, numbers, and underscores in predictor names
#' are replaced with underscores before model fitting. The function fits
#' `Group ~ predictors` with [rms::lrm()] and passes the result to
#' [regplot::regplot()] with violin and box displays.
#'
#' The function assigns an rms `datadist` object named `dd` in `.GlobalEnv`,
#' sets global `datadist` and `na.action` options, writes
#' `logistic_regplot.pdf` in the current working directory, and deletes that
#' PDF before returning. These global objects and options are not restored.
#'
#' When special feature names are present, the function invokes a qpdf binary
#' at a fixed path to replace cleaned labels with their original names, then
#' uses [magick::image_read_pdf()] at 600 DPI. This branch requires that qpdf
#' executable and its supporting shell tools to be available.
#'
#' @examples
#' \dontrun{
#' model_data <- data.frame(
#'   Group = factor(rep(c("Control", "Disease"), each = 25)),
#'   GeneA = rnorm(50),
#'   GeneB = rnorm(50)
#' )
#' nomogram <- regplot_plot(
#'   data = model_data,
#'   key_features = c("GeneA", "GeneB"),
#'   showP = TRUE,
#'   sample_select = 1
#' )
#' replayPlot(nomogram)
#' }
#'
#' @seealso [rms::lrm()], [rms::datadist()], [regplot::regplot()],
#'   [magick::image_read_pdf()]
#' @export
regplot_plot <- function(data = data,key_features = key_genes,showP = FALSE,sample_select = 1,formats = c("pdf","png","tiff")){
    
    library(dplyr)
    library(rms)
    library(regplot)
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
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
        data <- data %>% rename_with(~ ifelse(.x %in% check_result,
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

    # 构建 datadist
    assign("dd", datadist(data), envir = .GlobalEnv)
    options(datadist = "dd")
    options(na.action = "na.delete")
    # 构建公式
    foumular <- as.formula( paste0("Group ~ ", paste(sprintf("`%s`", key_features_updated), collapse = " + ")))
    # 拟合 Logistic 回归模型
    fit_logit <- lrm(foumular,  data = data, x = TRUE, y = TRUE)
    
    # 绘图保存 PDF
    pdf(paste0("logistic_regplot.pdf"), height = 6, width = 9)
    regplot(fit_logit,
            plots=c("violin", "boxes"), 
            showP = showP,
            droplines = T,
            clickable = F,
            points=T,
            subticks = T,
            failtime = NULL,
            title = 'Nomogram',
            observation = data[sample_select,]
           )
    p <- recordPlot()
    dev.off()
    
    # 重新绘制保存
    pdf(paste0("logistic_regplot.pdf"), height=6, width=9)
    replayPlot(p)
    dev.off()

    
    # 有特殊字符的话更改特征名为原本
    if( length(check_result) > 0 ){
        # 存在包含特殊字符的特征，替换为原名
        fix_pdf_gene_name_replace <- function(pdf_file, replacements, qpdf_path = "/data/home/zuoanjian/miniconda3/envs/RNA_Base/bin/qpdf") {
                library(tools)
                
                if (!file.exists(qpdf_path)) {
                stop("❌ qpdf binary not found at: ", qpdf_path)
                }
                
                qdf_file <- sub("\\.pdf$", "_uncompressed.pdf", pdf_file)
                #final_file <- sub("\\.pdf$", "_fixed.pdf", pdf_file)
                
                # Step 1: 解压
                cmd1 <- paste(shQuote(qpdf_path), "--qdf --object-streams=disable",shQuote(pdf_file), shQuote(qdf_file))
                message("🔧 解压命令：", cmd1)
                if (system(cmd1) != 0) stop("❌ 解压失败，请检查 qpdf 是否正常运行")
                
                # Step 2: 替换文本
                for (i in seq_along(replacements)) {
                    old <- names(replacements)[i]
                    new <- replacements[i]
                    cmd2 <- sprintf('sed -i "s/%s/%s/g" %s', old, new, shQuote(qdf_file))
                    message("🔁 替换命令：", cmd2)
                    system(cmd2)
                }
                
                # Step 3: 重新压缩
                cmd3 <- paste(shQuote(qpdf_path), shQuote(qdf_file), shQuote(pdf_file))
                message("📦 压缩命令：", cmd3)
                if (system(cmd3) != 0) stop("❌ 重新压缩失败，请检查")
            
                file.remove(qdf_file)
                message("✅ 替换完成：", pdf_file)
        }
        fix_pdf_gene_name_replace(pdf_file = paste0("logistic_regplot.pdf"),replacements = check_result,qpdf_path = "/data/home/zuoanjian/miniconda3/envs/RNA_Base/bin/qpdf")
        
        # pdf 转png
        img <- magick::image_read_pdf(paste0("logistic_regplot.pdf"), density = 600)
        file.remove(paste0("logistic_regplot.pdf"))
        
        return(img)
    }else{
        # 不存在包含特殊字符的特征，生成png、tiff
        file.remove(paste0("logistic_regplot.pdf"))
        return(p)
    }
  
}
