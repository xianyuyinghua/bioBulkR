#' Draw grouped box or violin plots with pairwise significance annotations
#'
#' Standardize a long-format measurement table, perform all pairwise group
#' comparisons within each feature, and draw a box or violin plot with optional
#' sample points and raw or adjusted p-value annotations.
#'
#' @param data_long A long-format data frame containing feature, group, and
#'   numeric value columns.
#' @param type_col Name of the feature column. Set to `NA` for data containing
#'   only groups and values, without a feature dimension.
#' @param group_col Name of the group column.
#' @param value_col Name of the numeric measurement column.
#' @param test_method Pairwise test: a value beginning with `"wilcox"` uses
#'   [rstatix::pairwise_wilcox_test()]; `"t"` uses
#'   [rstatix::pairwise_t_test()].
#' @param type_level Optional feature order. `NA` uses observed factor order.
#' @param group_level Optional group order. `NA` places `Normal` or `Control`
#'   first when present, then keeps the remaining observed order.
#' @param x_text_angle Numeric angle for x-axis labels.
#' @param add_point Logical; overlay jittered observations.
#' @param point_size,point_alpha,point_width Point size, transparency, and
#'   jitter width.
#' @param group_color Color vector for groups.
#' @param plot_type Plot family. Values beginning with `"box"` select a box
#'   plot; values beginning with `"vi"` select a violin plot.
#' @param diff_show_type Annotation form containing `"sig"` for stars or
#'   `"value"` for numeric p-values.
#' @param diff_show_name Significance source containing `"value"` for raw
#'   p-values or `"adj"` for Bonferroni-adjusted p-values.
#' @param sig_bracket_vjust Vertical justification of annotation text.
#' @param step_increase Incremental vertical spacing between comparison
#'   brackets.
#' @param hide_ns Logical; hide non-significant comparisons.
#' @param bracket_size Bracket line width.
#' @param bracket_y Vertical bracket nudge.
#' @param axis_y_title,axis_x_title,legend_title Axis and color-legend titles.
#' @param title_name Optional plot title.
#'
#' @return A list with `plot`, the annotated ggplot object, and `stat`, the
#' pairwise-test table including group pairs, p-values, adjusted p-values,
#' significance labels, and plotting coordinates.
#'
#' @details
#' Input columns are renamed internally to `Type`, `Group`, and `Value`. All
#' pairwise comparisons among factor levels of `Group` are tested separately
#' for each `Type`, with Bonferroni adjustment. Features whose values are all
#' zero are excluded from testing and receive non-significant placeholder
#' rows. For certain constant-group cases, the implementation adds 0.0001 to
#' one value before testing.
#'
#' When no feature column is requested, groups are drawn on the x-axis and the
#' legend is hidden. Otherwise, feature types are placed on the x-axis and
#' groups control color and dodging. `add_point = TRUE` suppresses boxplot
#' outliers and overlays the individual observations.
#'
#' @examples
#' \dontrun{
#' plot_data <- data.frame(
#'   Gene = rep(c("GeneA", "GeneB"), each = 12),
#'   Group = rep(rep(c("Control", "Disease"), each = 6), 2),
#'   Value = c(rnorm(6, 5), rnorm(6, 7), rnorm(6, 4), rnorm(6, 4.5))
#' )
#' result <- box_viloin_sig_plot(
#'   data_long = plot_data,
#'   group_level = c("Control", "Disease"),
#'   test_method = "wilcox",
#'   plot_type = "box",
#'   diff_show_type = "signif",
#'   diff_show_name = "p.adj",
#'   add_point = TRUE
#' )
#' result$plot
#' result$stat
#' }
#'
#' @seealso [ggpubr::ggboxplot()], [ggpubr::ggviolin()],
#'   [ggpubr::stat_pvalue_manual()], [rstatix::pairwise_wilcox_test()]
#' @export
box_viloin_sig_plot <- function(data_long, type_col = "Gene",group_col = "Group",value_col = "Value",test_method = "wilcox",
                                type_level = NA, group_level = NA,x_text_angle = 0,
                                add_point = FALSE,point_size = 1.8,point_alpha = 0.7,point_width = 0.15,
                                group_color = basicR::get_colors(package = "ggsci",name = "jco",number = 5.1),
                                plot_type = "box",diff_show_type = "signif",diff_show_name = "p.value",
                                sig_bracket_vjust = 0.7,step_increase = 0,hide_ns = TRUE,bracket_size = 0.8,bracket_y = 0,
                                axis_y_title = "Expression", axis_x_title = "Gene",legend_title = "Group",
                                title_name = NULL
                               ){
    library(ggplot2)
    ## 组间差异盒子图或小提琴图
    # rename
    if(is.na(type_col)){
        # no Type
        another_cols <- setdiff(colnames(data_long),c(group_col,value_col))
        if(c("Group","Value") %in% another_cols %>% any()){
            data_long <- data_long %>%
                              rename_with(
                                .cols = all_of(another_cols[another_cols %in% c("Group", "Value")]),
                                .fn = ~ paste0(.x, "_raw")
                              )
        }
        data_long <- data_long %>%  dplyr::rename("Group" = group_col,"Value"  = value_col)
        data_long$Type <- ""
    }else{
        # Type
        another_cols <- setdiff(colnames(data_long),c(group_col,value_col,value_col))
        if(c("Type","Group","Value") %in% another_cols %>% any()){
            data_long <- data_long %>%
                  rename_with(
                    .cols = all_of(another_cols[another_cols %in% c("Type","Group", "Value")]),
                    .fn = ~ paste0(.x, "_raw")
                  )
        }
        data_long <- data_long %>%  dplyr::rename("Type" = type_col,"Group" = group_col,"Value"  = value_col)
        
        # factor
        if(any(is.na(type_level))){
            data_long$Type <- as.factor(data_long$Type)
        }else{
            if( all(type_level %in% unique(data_long$Type)) ){
                data_long$Type <- factor(data_long$Type,levels = type_level)
            }else{
                stop("type_level not match data!")
            } 
        }
    }
    
    # group levels check
    if( length(group_level) == 1 && is.na(group_level) ){
        # no group_level
        if(any(c("Normal", "Control") %in% unique(data_long$Group))){
            normal_name <- intersect(c("Normal", "Control"),unique(data_long$Group))
            group_levels <- c(normal_name,setdiff(unique(data_long$Group),normal_name))
            data_long$Group <- factor(data_long$Group,levels = group_levels)
        }else{
            data_long$Group <- as.factor(data_long$Group)
        }
        
    }else{
        # group_level yes
        if( all(group_level %in% unique(data_long$Group)) ){
            data_long$Group <- factor(data_long$Group,levels = group_level)
        }else{
            stop("group_level not match data!")
        } 
    }
        
    # compare group
    my_comparisons <- utils::combn(levels(data_long$Group), 2, simplify = FALSE)

    deleted_types <- data_long %>%
      dplyr::group_by(Type) %>%
      dplyr::summarise(all_zero = all(Value == 0), .groups = "drop") %>%
      dplyr::filter(all_zero) %>%
      dplyr::pull(Type)
    
    filtered_types <- setdiff(unique(data_long$Type),deleted_types)
    
    # 检查数据中是否任意2个分组的值都为0
    zero_status <- data_long %>%
      dplyr::group_by(Type, Group) %>%
      dplyr::summarise(all_same = n_distinct(Value) == 1, .groups = "drop") %>%
      dplyr::group_by(Type) %>%
      dplyr::summarise(any_two_groups_same = any(combn(all_same, 2, FUN = function(x) all(x))), .groups = "drop")
                          
    deleted_types_special <- setdiff(zero_status %>% dplyr::filter(any_two_groups_same) %>% dplyr::pull(Type),deleted_types)
    if( length(deleted_types_special) != 0 ){                                            
        add_location <- lapply(deleted_types_special,function(feature){  which(data_long$Type %in% feature)[1] }) %>% unlist()
        raw_value <- data_long$Value[add_location]
        data_long$Value[add_location] <- raw_value + 0.0001 
    }                                            
    # 获取过滤后的数据
    filtered_data <- data_long %>% filter(Type %in% filtered_types)                                            

    df_mean_deleted <- data_long %>%
      dplyr::filter(Type %in% deleted_types) %>%
      dplyr::group_by(Type) %>%
      dplyr::summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
                                                
    # sig stat 
    if( grepl("^wilcox",test_method) ){
        stat.test <- filtered_data %>%
                    group_by(Type) %>%
                    rstatix::pairwise_wilcox_test(Value ~ Group,
                                                  comparisons = my_comparisons,
                                                  p.adjust.method = "bonferroni",
                                                  detailed = FALSE
                                                 )
        if(length(deleted_types) != 0){
            add_stat_test <- data.frame("Type" = rep(deleted_types,length(my_comparisons)),
                                        ".y." = rep("Value",length(my_comparisons)),
                                        "group1" = sapply(my_comparisons, `[`, 1),
                                        "group2" = sapply(my_comparisons, `[`, 2),
                                        "n1" = rep(unique(stat.test$n1),length(my_comparisons)),
                                        "n2" = rep(unique(stat.test$n2),length(my_comparisons)),
                                        "statistic" = df_mean_deleted$mean_value,
                                        "p" = rep(NA,length(my_comparisons)),
                                        "p.adj" = rep(NA,length(my_comparisons)),
                                        "p.adj.signif" = rep("ns",length(my_comparisons))
                                       )
            stat.test <- rbind(stat.test,add_stat_test)
        }
    }else if(grepl("^t$",test_method)){
        stat.test <- filtered_data %>%
                    group_by(Type) %>%
                    rstatix::pairwise_t_test(Value ~ Group,
                                             comparisons = my_comparisons,
                                             p.adjust.method = "bonferroni")
        if(length(deleted_types) != 0){
            add_stat_test <- data.frame("Type" = rep(deleted_types,length(my_comparisons)),
                            ".y." = rep("Value",length(my_comparisons)),
                            "group1" = sapply(my_comparisons, `[`, 1),
                            "group2" = sapply(my_comparisons, `[`, 2),
                            "n1" = rep(unique(stat.test$n1),length(my_comparisons)),
                            "n2" = rep(unique(stat.test$n2),length(my_comparisons)),
                            "p" = rep(NA,length(my_comparisons)),
                            "p.signif" = rep("ns",length(my_comparisons)),       
                            "p.adj" = rep(NA,length(my_comparisons)),
                            "p.adj.signif" = rep("ns",length(my_comparisons))
                           )
            stat.test <- rbind(stat.test,add_stat_test)
        }
    }else{
        stop("method no matched!")
    }

    # 添加 p.signif
    stat.test <- stat.test %>%
              mutate(p.signif = case_when( 
                  #p < 0.0001 ~ "****",
                  p < 0.001 ~ "***",
                  p < 0.01 ~ "**",
                  p < 0.05 ~ "*",
                  TRUE ~ "ns" 
              ),
                     p.adj.signif =case_when( p.adj.signif == "****" ~ "***",
                                             TRUE ~ p.adj.signif
                                            )
                    )
    
    
    # 添加 x y  position
    if( length(stat.test$Type) == 1 ){
        # 只有1类
        x_name <- "Group"
        stat.test$y.position = max(filtered_data$Value) * 1.1
        stat.test$xmin = 1
        stat.test$xmax = 2

        stat.test  <- as.data.frame(stat.test)
        stat.test$groups <- paste0(stat.test$group2,"_vs_",stat.test$group1)
    }else{
        # 有多类
        x_name <- "Type"
        stat.test <- stat.test %>% rstatix::add_xy_position(x = "Type") 

        stat.test  <- as.data.frame(stat.test)
        stat.test$groups <- sapply(stat.test$groups, function(x) paste(x, collapse = "_vs_"))
    }
    

                               
    # plot
    if(add_point){outlier_shape = NA}else{outlier_shape = 19}
    # 根据 plot_type 设置绘图函数
    plot_func <- if (grepl("^box", plot_type)) {
      ggpubr::ggboxplot
    } else if (grepl("^vi", plot_type)) {
      ggpubr::ggviolin
    } else {
      stop("不支持的 plot_type")
    }
    
    # 调用绘图函数，并统一设置
    p1 <- plot_func(data_long, x = x_name, y = "Value", width = 0.6,
                    color = "Group",
                    palette = "jco",
                    bxp.errorbar = FALSE,
                    bxp.errorbar.width = 0.5,
                    size = 0.5,
                    outlier.shape = outlier_shape,
                    legend = "right"  ) +
            labs(y = axis_y_title, x = axis_x_title,color = legend_title)
    # 是否加散点                               
    if (add_point) {
        if (x_name == "Group") {
            p1 <- p1 +
                geom_jitter(
                    aes(x = Group, y = Value, color = Group),
                    width = point_width,
                    size = point_size,
                    alpha = point_alpha
                )
        } else {
            p1 <- p1 +
                geom_jitter(
                    aes(x = Type, y = Value, color = Group),
                    position = position_jitterdodge(
                        jitter.width = point_width,
                        dodge.width = 0.8
                    ),
                    size = point_size,
                    alpha = point_alpha
                )
        }
    }
    # 是否显示图例
    if(x_name == "Group"){
        p1 <- p1 + theme(legend.position = "none")
    }
    
    # 设置x text 角度
    if(x_text_angle == 0){
        vjust_value = 0.5   
        hjust_value = 0.5 
    }else if(x_text_angle == 90){
        vjust_value = 0.5
        hjust_value = 1
    }else{
        vjust_value = 1
        hjust_value = 1 
    }
        
  p1 <- p1 + labs(title = title_name)  
  p1 <- p1 + theme(axis.text.x = element_text(angle = x_text_angle, vjust = vjust_value, hjust = hjust_value,face = "bold",size = 17),
                   axis.title = element_text(face = "bold", size = 20),
                   axis.text.y = element_text(size = 17),
                   legend.title = element_text(face = "bold", size = 17),
                   legend.text = element_text(size = 17),
                   plot.title = element_text(hjust = 0.5,vjust = 1,face = "bold", size = 18),
                   plot.margin = ggplot2::margin(t = 5,r = 5,b = 5,l = 5,unit = "mm")
                  )
    if(is.na(type_col)){
        p1 <- p1 + theme(axis.ticks.x = element_blank())
    }
                                                
    if(length(group_color) > 1){
        # 有指定颜色
        p1 <- p1 +  scale_color_manual(values = group_color)
    }

    # 显示显著性
    if(grepl("sig",diff_show_type)){
        lable_suffix <- ".signif"
    }else if(grepl("value",diff_show_type)){
        lable_suffix <- ""
    }

    if(grepl("value",diff_show_name)){
        lable_pre <- "p"
    }else if(grepl("adj",diff_show_name)){
        lable_pre <- "p.adj"
    }

    label_show <- paste0(lable_pre,lable_suffix)
        
    p1 + ggpubr::stat_pvalue_manual(
      stat.test,
      label = label_show,
      step.increase = step_increase,
      hide.ns = hide_ns,
      tip.length = 0,  # 调整括号末端长度
      vjust = sig_bracket_vjust, # 调整 * 线y位置
      bracket.size = bracket_size,
      bracket.nudge.y = bracket_y, # 调整 线y位置
      size = 5
    ) -> p1
    return(list(plot = p1, stat = stat.test))
}
