#' Calculate feature correlations and optionally create significance plots
#'
#' Test every pair formed by common and separate feature sets, attach
#' significance stars, reshape the statistics into matrices, and optionally
#' create lollipop, faceted, scatter, or heatmap visualizations.
#'
#' @param cor_df A data frame containing observations in rows and all requested
#'   numeric features in columns. Row names are used to align reshaped data.
#' @param common_features Character vector of features each compared with every
#'   member of `separate_features`.
#' @param separate_features Character vector defining the other side of each
#'   pairwise correlation.
#' @param separate_features_levles Factor levels controlling the order of
#'   separate features. The historical argument spelling is retained.
#' @param cor_method Correlation method passed to [stats::cor.test()]:
#'   `"spearman"`, `"kendall"`, or `"pearson"`.
#' @param plot Logical; create plots in addition to result tables.
#' @param plot_title Optional title for lollipop plots. `NULL` uses
#'   `"Correlation"`.
#' @param axis_x_angle Numeric angle for common-feature labels on lollipop
#'   x-axes.
#' @param facet_ncol Number of columns in the faceted lollipop plot. `NA`
#'   calculates it with [facet_nrow_ncol()].
#' @param scatter_position Optional numeric vector `c(x, y)` giving the
#'   correlation-statistic annotation position in the one-to-one scatter plot.
#' @param heatmap_sig_height Numeric vertical offset, in normalized parent
#'   coordinates, separating correlation values and significance stars in
#'   heatmap cells.
#'
#' @return If `plot = FALSE`, a list containing `cor_pvalue_sig_data`,
#'   `cor_matrix`, `pvalue_matrix`, and `significance_matrix`. If `plot = TRUE`,
#'   a list with `data`, containing those four objects, and `plot`, containing
#'   `lollipop_list`, `lollipop_facet`, `heatmap`, and `scatter`.
#'
#' @details
#' Each feature pair is tested independently with [stats::cor.test()]. Missing
#' correlation estimates are replaced by zero and missing p-values by one.
#' Significance labels are `***` for p below 0.001, `**` below 0.01, `*` below
#' 0.05, and an empty string otherwise.
#'
#' With plotting enabled, individual and faceted lollipop plots are always
#' prepared. A single common feature paired with a single separate feature
#' additionally produces a linear-fit scatter plot with marginal density
#' distributions. One-to-many or many-to-many input instead produces a
#' clustered ComplexHeatmap showing correlations and significance stars.
#'
#' @examples
#' \dontrun{
#' result <- cor_test_star_plot(
#'   cor_df = iris,
#'   common_features = c("Sepal.Length", "Sepal.Width"),
#'   separate_features = c("Petal.Length", "Petal.Width"),
#'   separate_features_levles = c("Petal.Length", "Petal.Width"),
#'   cor_method = "spearman",
#'   plot = FALSE
#' )
#' result$cor_matrix
#' }
#'
#' @seealso [stats::cor.test()], [facet_nrow_ncol()],
#'   [ComplexHeatmap::Heatmap()], [ggExtra::ggMarginal()]
#' @export
cor_test_star_plot <- function(cor_df,
                               common_features,
                               separate_features,
                               separate_features_levles,
                               cor_method = "spearman",  # spearman/kendall/pearson
                               plot = FALSE,
                               plot_title = NULL,
                               axis_x_angle = 45,
                               facet_ncol = NA,
                               scatter_position = NULL,
                               heatmap_sig_height = 0.015
                         ){
    # common_features: 一起和其它的特征计算相关性
    # separate_features： 分别和共同特征计算相关性
    # cor_df：包含了common_features和separate_features特征列的数据

    # library
    library(ggplot2)

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 检查 特征列是否都在数据中
    if( !all(c(common_features,separate_features) %in% colnames(cor_df)) ){
        not_exist_feature <- setdiff(c(common_features,separate_features),colnames(cor_df))
        cat("\n\n",red,not_exist_feature,reset,"Not found.","\n\n")
        flush.console()
        stop()
    }

    # 宽数据转为长数据
    separate_features_df <- cor_df %>% dplyr::select(all_of(separate_features)) %>%  tibble::rownames_to_column("RowName")
    # 1. 选取相关基因并转换数据格式
    cor_long <- cor_df %>%
         dplyr::select(all_of(common_features)) %>%
         tibble::rownames_to_column("RowName") %>%
         tidyr::pivot_longer(cols = common_features, names_to =  "common_features", values_to = "Value_com")
    cor_long <- merge(cor_long,separate_features_df,by = "RowName")
    cor_long <- cor_long %>% tidyr::pivot_longer(cols = separate_features, names_to = "separate_features", values_to = "Value_separate")
    
    # 2. 计算相关性
    cor_results <- cor_long %>%
          dplyr::group_by(common_features, separate_features) %>%
          dplyr::summarize(
            cor = cor.test(Value_com, Value_separate, method = cor_method)$estimate,
            p_value = cor.test(Value_com, Value_separate, method = cor_method)$p.value,
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            cor = ifelse(is.na(cor), 0, cor),
            p_value = ifelse(is.na(p_value), 1, p_value)  # p 值 NA 的设为 1
          )
    
    # 3. 添加显著性标记
    cor_results <- cor_results %>%
         dplyr::mutate(significance = case_when(p_value < 0.001  ~ "***",
                                         p_value < 0.01   ~ "**",
                                         p_value < 0.05   ~ "*",
                                         TRUE  ~ ""
                                        ))

    
    cor_results$common_features <- as.factor(cor_results$common_features)
    cor_results$separate_features <- factor(cor_results$separate_features,separate_features_levles)

    # 转为3个宽数据矩阵
    cor_matrix <- cor_results %>% dplyr::select(-c(p_value,significance)) %>%
      tidyr::pivot_wider(names_from = common_features, values_from = cor) %>%
      tibble::column_to_rownames("separate_features") %>%
      as.matrix()
    
    p_value_matrix <- cor_results %>% dplyr::select(-c(cor,significance)) %>%
      tidyr::pivot_wider(names_from = common_features, values_from = p_value) %>%
      tibble::column_to_rownames("separate_features") %>%
      as.matrix()
    
    significance_matrix <- cor_results %>% dplyr::select(-c(p_value,cor)) %>%
      tidyr::pivot_wider(names_from = common_features, values_from = significance) %>%
      tibble::column_to_rownames("separate_features") %>%
      as.matrix()

   # plot
   if(plot){
        # 根据字符长度计算左侧x距离
        if(axis_x_angle == 90){
          x_left_distence = 0.03 
        }else{
          x_left_distence =  0.03 + max(nchar(as.character(cor_results$common_features)))*0.001
        }
        # 设置x text 角度
        if(axis_x_angle == 0){
            vjust_value = 0.5   
            hjust_value = 0.5 
        }else if(axis_x_angle == 90){
            vjust_value = 0.5
            hjust_value = 1
        }else{
            vjust_value = 1
            hjust_value = 1 
        }


       
       # 计算分面列
       if(is.na(facet_ncol)){
           facet_ncol = facet_nrow_ncol(facet_number = length(unique(separate_features)))[2]
       }
       if(is.null(plot_title)){
           title_character <- "Correlation"
       }else{
           title_character <- plot_title
       }
       
       # 1.1 单个特征棒棒糖图
       lollipop_list <- lapply(unique(cor_results$separate_features),function(feature){
                     # 单独每个特征出图
                    cor_result <- cor_results %>% dplyr::filter(separate_features == feature)
                     # 画棒棒糖图（以 gene_risk 分组）
                    p1 <- ggplot(cor_result, aes(x = reorder(common_features, cor), y = cor, color = -log10(p_value))) +
                            geom_segment(aes(x = common_features, xend = common_features, y = 0, yend = cor), linewidth = 0.5) + # 连接线
                            geom_point(size = 3) +  # 圆点
                            geom_text(aes(label = significance, vjust = ifelse(cor > 0, -0.5, 2)), size = 4) + # 添加显著性标记
                            scale_color_viridis_c(option = "D", direction = 1) +  # 颜色根据 p 值显著性变化
                            theme_classic(base_size = 20) +
                            ylim(min(min(cor_results$cor) - 0.1,-0.01), max(max(cor_results$cor) + 0.1,0.01)) +
                            scale_x_discrete(expand = expansion(mult = c(x_left_distence, 0.04)))+
                            labs(x = "", 
                                 y = "Correlation Coefficient", 
                                 color = "-log10(p)", 
                                 title = paste(title_character,"\n", feature)) +  # 添加每个图的标题
                            geom_hline(yintercept = c(0.3, -0.3), linetype = "dashed", color = "grey50") +  # 添加虚线
                            theme(axis.text.x = element_text(angle = axis_x_angle, hjust = hjust_value,vjust = vjust_value,face = "bold"),
                                  strip.background = element_blank(),
                                  plot.title = element_text(face = "bold"),
                                  axis.title = element_text(face = "bold"),
                                  plot.margin = ggplot2::margin(t = 5, r = 5, b =5, l = 20, unit = "mm")
                                 )
                    return(p1)
            })
       lollipop_list <- setNames(lollipop_list,levels(cor_results$separate_features))

       
       # 1.2 分面棒棒糖图
       lollipop_facet <- ggplot(cor_results, aes(x = reorder(common_features, cor), y = cor, color = -log10(p_value))) +
                geom_segment(aes(x = common_features, xend = common_features, y = 0, yend = cor), size = 0.5) +  # 连接线
                geom_point(size = 3) +  # 圆点
                geom_text(aes(label = significance, vjust = ifelse(cor > 0, -0.5, 2)), size = 4) +  # 添加显著性标记
                scale_color_viridis_c(option = "D", direction = 1) +  # 颜色根据 p 值显著性变化
                theme_classic(base_size = 20) +
                ylim(min(min(cor_results$cor) - 0.1,-0.01), max(max(cor_results$cor) + 0.1,0.01)) +
                scale_x_discrete(expand = expansion(mult = c(x_left_distence, 0.04)))+
                labs(x = "", y = "Correlation Coefficient", color = "-log10(p)", 
                     title = title_character) +  # 添加每个图的标题
                geom_hline(yintercept = c(0.3, -0.3), linetype = "dashed", color = "grey50") +  # 添加虚线
                theme(axis.text.x = element_text(angle = axis_x_angle, hjust = hjust_value,vjust = vjust_value,,face = "bold"),
                      strip.background = element_blank(),
                      strip.text = element_text(size = 21, face = "bold"),
                      plot.title = element_text(face = "bold"),
                      axis.title = element_text(face = "bold"),
                      plot.margin = ggplot2::margin(t = 5, r = 5, b =5, l = 20,unit = "mm")
                     )+
             facet_wrap(~ separate_features, ncol = facet_ncol)

       # 2.1 散点图或者热图
       if(length(separate_features) == 1 & length(common_features) == 1){
           if(is.null(scatter_position)){
               x_position <- min(cor_df[[separate_features]]) + ((max(cor_df[[separate_features]]) - min(cor_df[[separate_features]]))/2)
               y_position <- max(cor_df[[common_features]]) - ((max(cor_df[[common_features]]) - min(cor_df[[common_features]]))/15) 
           }else{
               x_position <- scatter_position[1]
               y_position <- scatter_position[2]
           }
            # 1对1 相关性散点图
            r_value <- round(cor_results$cor, 3)
            p_value <- signif(cor_results$p_value, 3)

            # 绘制散点图
            p_scatter <- ggplot(cor_df, aes(x = !!sym(separate_features), y = !!sym(common_features))) +
                geom_point(color = "#4682B4", alpha = 0.7) +  # 蓝色散点
                geom_smooth(method = "lm", color = "blue", fill = "lightblue", se = TRUE) +  # 拟合线
                annotate("text", 
                         x = x_position, 
                         y = y_position,
                         label = paste("r :", r_value, " p :", p_value), size = 9, color = "black") +  # **手动标注 R 和 P**
                theme_classic(base_size = 20) +
                labs(x = gsub("_"," ",separate_features), y = gsub("_"," ",common_features)) +  # Y 轴设置为当前药物名称
                theme(
                  axis.text = element_text(size = 20,face = "bold"),  # 坐标轴文字
                  axis.title = element_text(size = 21,face = "bold"),  # 轴标题
                  plot.margin = margin(10, 20, 10, 10)  # **增加右侧边距，防止裁剪**
                )
            # 添加边缘密度图
            p_scatter_marginal <- ggExtra::ggMarginal(p_scatter, type = "density", fill = "#4682B4", color = NA, alpha = 0.5)
           
            p_heatmap <- NULL
       }else if(sum(length(separate_features),length(common_features)) > 2){
           # 1对多或多对多相关性 相关性热图
            suppressPackageStartupMessages({
                library(viridis)
                library(ComplexHeatmap)
                library(circlize)
            })   
    
            cor_mat <- cor_matrix
            sig_mat <- significance_matrix
            viridis_colors <- viridis::viridis(n = 100)
    
           # 绘图
            p_heatmap <- ComplexHeatmap::Heatmap(
                matrix = cor_mat,
                name = "Spearman\ncorrelation",
                col = colorRamp2(seq(min(cor_mat), max(cor_mat), length.out = 100), viridis_colors),
                cell_fun = function(j, i, x, y, width, height, fill) {
                    # R值
                    r_value <- sprintf("%.2f", cor_mat[i, j])
                    # 显著性
                    star <- sig_mat[i, j]
                    # 上面显示R
                    grid.text(r_value,x = x,y = y*0.99 + unit(heatmap_sig_height, "npc"),gp = gpar(fontsize = 14))
                    # 下面显示星号
                    grid.text(star,x = x,y = y*0.99 - unit(heatmap_sig_height, "npc"),gp = gpar(fontsize = 14))},
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_names_side = "left",
                column_names_rot = 45,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 10),
                                            title_position = "topleft",
                                            legend_width = unit(3, "cm")
                                           ),
                row_names_gp = gpar(fontsize = 12,fontface = "bold"),
                column_names_gp = gpar(fontsize = 12,fontface = "bold"),
                row_names_max_width = unit(10, "cm")
            )
           p_heatmap <- draw(p_heatmap,padding = unit(c(5, 10, 5, 10), "mm"))
        
           p_scatter_marginal <- NULL
       }
       
       # 输出数据+图
       return(list(data = list(cor_pvalue_sig_data = cor_results,
                               cor_matrix = cor_matrix,
                               pvalue_matrix = p_value_matrix,
                               significance_matrix = significance_matrix
                              ),
                   plot = list(lollipop_list = lollipop_list,
                               lollipop_facet = lollipop_facet,
                               heatmap = p_heatmap,
                               scatter = p_scatter_marginal
                              )
                  )
             )
       
   }else{

    # 输出数据
    return(list(cor_pvalue_sig_data = cor_results,
                cor_matrix = cor_matrix,
                pvalue_matrix = p_value_matrix,
                significance_matrix = significance_matrix
               )
          )
    }
}
