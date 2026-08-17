enrich_go_kegg <- function(genes,
                            species, 
                            p_name = "p.value",  # p.value/p.adj
                            filter_pvalue = 0.05,
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 1,
                            kegg_analysis_method = "online", # online/local
                            plot_type = NULL, # NULL,barplot,dotplot,sankeyplot,treeplot,circularplot
                            scankeyplot_mode = "sankey_buble", # buble_sankey/sankey_buble
                            fill_colors = basicR::get_colors(number = 4.8),
                            go_kegg_topn = c(10,10),
                            description_wrap_width = 100
                            ){
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(grid)
        library(viridis)
    	library(stringr)
    })
    ###################################################################################################
    ########################################### Function ##############################################
    ###################################################################################################
    # Barplot
    barplot <- function(data,p_name = "pvalue",colors = basicR::get_colors(number = 4.8)){
        # 输入 富集后的数据
        if( all(unique(data$ONTOLOGY) %in% c('BP', 'CC', 'MF', 'KEGG')) ){
            data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = rev(c('BP', 'CC', 'MF', 'KEGG'))))
        }else if( all(unique(data$ONTOLOGY) %in% c('BP', 'CC', 'MF')) ){
            data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = rev(c('BP', 'CC', 'MF'))))
        }else{
           data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = ONTOLOGY )) 
        }
        data <- data %>% dplyr::arrange(ONTOLOGY,desc(!!sym(p_name)), Count) %>% dplyr::mutate(Description = factor(Description, levels = unique(Description)))
        data <- data %>% tibble::rowid_to_column('index')

        # 将 geneID 分割并插入换行符
        data <- data %>% dplyr::mutate(geneID_trunc = purrr::map_chr(str_split(geneID, "/"), ~ if (length(.) > 12){ paste(.[1:12], collapse = "/")}else{paste(., collapse = "/")}))
              
        # 左侧分类标签和基因数量点图的宽度
        # x 轴长度
        xaxis_max <- max(-log10(data[[p_name]])) + 1
        width <- max(-log10(data[[p_name]])) * 0.05
        # 左侧分类标签数据
        rect.data <- group_by(data, ONTOLOGY) %>% reframe(n = n()) %>% ungroup() %>% dplyr::mutate(xmin = -3 * width,xmax = -2 * width,ymax = cumsum(n),ymin = lag(ymax, default = 0) + 0.6,ymax = ymax + 0.4)
        
        ggplot2::ggplot(data = data,aes(-log10(!!sym(p_name)), y = index, fill = ONTOLOGY)) +
            gground::geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
            ggplot2::geom_text(aes(x = 0.05, label = Description),hjust = 0, size = 6,lineheight = 0.6
                               #fontface = "bold"
                              ) +
            ggplot2::geom_text(aes(x = 0.1, label = geneID_trunc, colour = ONTOLOGY), hjust = 0, vjust = 2.6, size = 4, fontface = 'italic',show.legend = FALSE) +
            # 基因数量
            ggplot2::geom_point(aes(x = -width, size = Count),shape = 21) +
            ggplot2::geom_text(aes(x = -width, label = Count),size = 5) +
            scale_size_continuous(name = 'Count', range = c(5, 8)) +
            # 分类标签
            gground::geom_round_rect(data = rect.data,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill = ONTOLOGY),radius = unit(2, 'mm'),inherit.aes = FALSE ) +
            ggplot2::geom_text(data = rect.data, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),inherit.aes = FALSE,fontface = "bold",angle = 90,size = 5) +
            ggplot2::geom_segment(aes(x = 0, y = 0, xend = xaxis_max, yend = 0),linewidth = 1.5,inherit.aes = FALSE ) +
            labs(y = NULL) +
            scale_fill_manual(name = 'Category', values = colors,guide = guide_legend(reverse = TRUE)) +
            scale_colour_manual(values =  colors) +
            scale_x_continuous(breaks = seq(0, xaxis_max, 4),expand = expansion(c(0, 0))) +
            ggprism::theme_prism(base_size = 18, base_family = "") +
            theme(axis.text.y = element_blank(),
                  axis.line = element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.title = element_text()
                 ) -> p
        
            return(p)
    }

    # Dotplot
    dotplot <- function(data, 
                        x_col = "GeneRatio",
                        size_col = "Count",
                        fill_col = "pvalue",
                        fill_colors = fill_colors,
                        base_size = 18
                       ){
        # plot
        ggplot(data,aes(x = !!sym(x_col), y = Description)) +
            geom_point(aes(size = !!sym(size_col),color = !!sym(fill_col))) +
            facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
            scale_color_gradientn(colours = fill_colors,guide = guide_colorbar(reverse = TRUE))+
            scale_size(range = c(2, 8)) +
            labs(x = x_col,y = NULL) +
            theme_bw(base_size = base_size) +
            theme(panel.grid.major = element_blank(),   # 去掉主网格线
                  panel.grid.minor = element_blank(),   # 去掉次网格线
                  strip.text = element_text(face = "bold"),
                  strip.background = element_rect(fill = "grey95"),
                  legend.position = "right",
                  legend.title = element_text(face = "bold"),
                  axis.text.x = element_text(size = base_size),
                  axis.title.x = element_text(face = "bold"),
                  axis.text.y = element_text(size = base_size + 2,lineheight = 0.6)
               )
    }

    # sankeyplot
    sankeyplot <- function(data,
                           plot_mode = "sankey_buble",  # buble_sankey/sankey_buble
                           flow_colors = rep(basicR::get_colors(number = 100.1),5),
                           fill_colors = basicR::get_colors(number = 10.1,package = "viridis",name = "D") %>% rev(),
                           description_width = 100,
                           buble_panel_border = FALSE
                          ){
        library(ggplot2)
        library(ggsankey)
        library(scales)
        # 检查数据
        
        # 出图模式
        if(plot_mode == "sankey_buble"){
            gene_location = 0.02
            pathway_location = 2.40
            bubble_x_min <- 2.5
            bubble_x_max <- 3.3
            hjust_left = 1
            hjust_right = 1
            nudge_x_left = 0
            nudge_x_right = -0.15
            fontface_left = "plain"
            fontface_right = "bold"
            size_left = 4
            size_right = 6
            xlim = c(-0.25, 3.25)
            legend_position = "right"
            plot_margin = ggplot2::margin(b = 10, l = 0, t = 10, r = 20)
            
            
        }else if(plot_mode == "buble_sankey"){
            gene_location = 2.30
            pathway_location = 0.70
            bubble_x_min <- 0.00
            bubble_x_max <- 0.6
            hjust_left = 0
            hjust_right = 0
            nudge_x_left = 0.08
            nudge_x_right = 0.04  
            fontface_left = "bold"
            fontface_right = "plain"
            size_left = 6
            size_right = 4
            xlim = c(0.00, 2.60)
            legend_position = "left"
            plot_margin =  ggplot2::margin(b = 10, l = 30, t = 10, r = 10)
        }
        
        # 1. 富集通路基因流数据
        data <- data %>% as.data.frame(check.names = F)   
        df_long <- data %>% dplyr::select(Description, geneID) %>% tidyr::separate_rows(geneID, sep = "/") %>% dplyr::rename(pathway = Description, gene = geneID)
        sankey_df <- df_long %>% ggsankey::make_long(pathway, gene)
        
        # 2. 把 sankey 的 x 改成数值 
        sankey_num <- sankey_df %>% dplyr::mutate(x_num = dplyr::case_when(x == "pathway" ~ pathway_location,x == "gene" ~ gene_location,TRUE ~ NA_real_),
                                                  next_x_num = dplyr::case_when(next_x == "pathway" ~ pathway_location,next_x == "gene" ~ gene_location,TRUE ~ NA_real_),
                                                  pathway_fill = node
                                                 )
    
        # 3. 临时 sankey：提取节点坐标 =============================================================================
        p_tmp <- ggplot(sankey_num,aes(x = x_num,next_x = next_x_num,node = node,next_node = next_node,fill = pathway_fill)) +
                    geom_sankey(flow.alpha = 0.65,smooth = 8, width = 0.08, show.legend = FALSE) +
                    scale_fill_manual(values = flow_colors) +
                    theme_void()
        pb <- ggplot_build(p_tmp)
        node_df <- pb$data[[2]] %>% dplyr::distinct(node, x, xmin, xmax, ymin, ymax, fill)
    
        if(plot_mode == "sankey_buble"){
            # gene 节点（最左）
            left_nodes <- node_df %>% dplyr::filter(x == min(x)) %>% dplyr::mutate(x_text = xmin - 0.03,y = (ymin + ymax) / 2,
                                                                                   label_wrap = ifelse(nchar(node) > 10, str_wrap(node, width = 10), node)
                                                                                  )
            # pathway 节点（中间）
            right_nodes <- node_df %>% dplyr::filter(x == max(x)) %>% dplyr::mutate(x_text = xmax + 0.05,y = (ymin + ymax) / 2,
                                                                                  label_wrap = ifelse(nchar(node) > description_width, str_wrap(node, width = description_width), node)
                                                                                 )
        }else if(plot_mode == "buble_sankey"){
            # 左侧 pathway 节点
            left_nodes <- node_df %>% dplyr::filter(x == min(x)) %>% dplyr::mutate(x_text = xmin,y = (ymin + ymax) / 2,
                                                                                   label_wrap = ifelse(nchar(node) > description_width,str_wrap(node, width = description_width),node)
                                                                                  )
            # 右侧 gene 节点
            right_nodes <- node_df %>% dplyr::filter(x == max(x)) %>% dplyr::mutate(x_text = (xmin + xmax) / 2,y = (ymin + ymax) / 2,
                                                                                    label_wrap = ifelse(nchar(node) > 10,str_wrap(node, width = 10),node)
                                                                                   )
        }
    
        # 4. 富集结果 join 到 pathway 的真实 y 坐标 ===================================================================
        if(plot_mode == "sankey_buble"){    
            bubble_df <- data %>% dplyr::inner_join(right_nodes %>% dplyr::transmute(Description = node,y,ymin,ymax,pathway_fill = fill,label_wrap),by = "Description")
        }else if(plot_mode == "buble_sankey"){
            bubble_df <- data %>% dplyr::inner_join(left_nodes %>% dplyr::transmute(Description = node,y,ymin,ymax,pathway_fill = fill,label_wrap),by = "Description")
        }
    
        # 5. 定义布局位置 ============================================================================================
        # Sankey 的 pathway 列是 x=1，gene 列是 x=2
        # bubble 放在左边，pathway 标签放在 bubble 和 sankey 之间
        # 把 RichFactor 映射到 bubble 区域
        rescale_to <- function(x, from, to) {(x - from[1]) / (from[2] - from[1]) * (to[2] - to[1]) + to[1]}
        bubble_df <- bubble_df %>% dplyr::mutate(bubble_x = rescale_to(RichFactor,from = c(min(RichFactor), max(RichFactor)),to = c(bubble_x_min, bubble_x_max)))
    
        # 6. 气泡图区边框和 x 轴参数 ===================================================================================
        bubble_box_xmin <- bubble_x_min - 0.05
        bubble_box_xmax <- bubble_x_max + 0.05
        
        if(plot_mode == "sankey_buble"){
            bubble_box_ymin <- min(right_nodes$ymin)
            bubble_box_ymax <- max(right_nodes$ymax) * 1.02               
        }else if(plot_mode == "buble_sankey"){
            bubble_box_ymin <- min(left_nodes$ymin)
            bubble_box_ymax <- max(left_nodes$ymax) * 1.02   
        }
    
        # x轴位置
        axis_y <- bubble_box_ymin * 1.03
        
        # 原始 RichFactor 刻度
        x_breaks_raw <- pretty(data$RichFactor, n = 3)
        x_breaks_raw <- x_breaks_raw[x_breaks_raw >= min(data$RichFactor) & x_breaks_raw <= max(data$RichFactor)]
        
        # 映射到当前 bubble 区域的 x 坐标
        x_breaks_plot <- rescale_to(x_breaks_raw,from = c(min(data$RichFactor), max(data$RichFactor)),to = c(bubble_x_min, bubble_x_max))
        
        axis_df <- data.frame(x_raw  = x_breaks_raw,x_plot = x_breaks_plot)
    
        # 7. 正式作图 =====================================================================================================
        p_final <- ggplot() +
            geom_sankey(data = sankey_num,aes(x = x_num,next_x = next_x_num,node = node,next_node = next_node,fill = pathway_fill),flow.alpha = 0.6,smooth = 9,width = 0.05,show.legend = FALSE) + # Sankey 主体 
            scale_fill_manual(values = flow_colors) +
            ggnewscale::new_scale_color() +  # 为 bubble 开启新的颜色标尺
            geom_point(data = bubble_df,aes(x = bubble_x,y = y,size = Count,color = p.adjust),alpha = 0.95) + # 左侧 bubble
            scale_color_gradientn(colours = fill_colors,guide = guide_colorbar(reverse = TRUE,order = 1,theme = theme(legend.key.height = grid::unit(5, "lines"),legend.key.width = grid::unit(1.2, "lines")))) +
            scale_size_continuous(range = c(4, 8),name = "Size",guide = guide_legend(order = 2,override.aes = list(color = "grey40",alpha = 0.95),theme = theme(legend.key.height = grid::unit(1, "lines"),legend.key.width = grid::unit(1, "lines")))) +
            annotate("segment",x = bubble_box_xmin,xend = bubble_box_xmax,y = axis_y,yend = axis_y,linewidth = 0.8,color = "black") + # 气泡图 x 轴 
            geom_segment(data = axis_df,aes(x = x_plot,xend = x_plot,y = axis_y,yend = axis_y - 0.5),inherit.aes = FALSE,linewidth = 0.8,color = "black") + # x轴刻度线
            geom_text(data = axis_df,aes(x = x_plot,y = axis_y - 0.7,label = label_number(accuracy = 0.01)(x_raw)),inherit.aes = FALSE,size = 6,vjust = 1) + # x轴刻度文字
            annotate("text",x = (bubble_box_xmin + bubble_box_xmax) / 2,y = axis_y - (2.5+(nrow(data)*0.1)), label = "Ratio",size = 8) + # x轴标题
            geom_text(data = left_nodes,aes(x = x_text,y = y,label = label_wrap),hjust = hjust_left,nudge_x = nudge_x_left,vjust = 0.5,size = size_left,lineheight = 0.6,fontface = fontface_left) + # 中间 pathway 标签
            geom_text(data = right_nodes,aes(x = x_text,y = y,label = label_wrap),hjust = hjust_right,nudge_x = nudge_x_right,vjust = 0.5,size = size_right,lineheight = 0.6,fontface = fontface_right) + # 右侧 gene 标签
            coord_cartesian(xlim = xlim,ylim = c(min(node_df$ymin) - 14, max(node_df$ymax)),clip = "off") + # 坐标范围
            theme_void(base_size = 20) +
            theme(legend.position = legend_position,legend.box = "vertical",legend.key.height = grid::unit(0.45, "lines"),
                  legend.spacing.y = grid::unit(0.05, "lines"),legend.margin = ggplot2::margin(0, 0, 0, 0),plot.margin = plot_margin
                 )
    
        if(buble_panel_border){
            p_final <- p_final + 
                annotate("rect",xmin = bubble_box_xmin,xmax = bubble_box_xmax,ymin = axis_y,ymax = bubble_box_ymax,fill = NA,color = "black",linewidth = 0.8) # 气泡图区四边框
        }
        
        return(p_final)
    }

    # circosEnrichmentPlot
    circosEnrichmentPlot <- function(df,
                                     topN = 6,
                                     classCol = c("#f7cb16", "#65c3fc", "#bfe046", "#bdb5e3", "#fbccca", "#54beaa"),
                                     ifLog = FALSE,
                                     type = "base"
                                    ){
    
        # 代码来源：https://www.r2omics.cn/
        # 加载必要库
        library(tidyverse)
        library(circlize)
        library(RColorBrewer)
        library(ComplexHeatmap)
        
        df.top  <- df %>% group_by(Class) %>% slice_min(PValue, n = topN, with_ties = FALSE) %>% dplyr::arrange(desc(Ratio), .by_group = TRUE) %>% ungroup()
    
        main.col  <- classCol[as.numeric(as.factor(df.top$Class))]
    
        # 创建四个数据层
        df1 <- data.frame(TermID  = df.top$TermID,  start = 0, end = max(df.top$bg_pro_num,  na.rm  = TRUE))
        
        # 颜色映射
        p_max <- round(max(df.top$Enrichment))  + 1
        color_assign <- colorRamp2(0:p_max, colorRampPalette(c('#fee5d9', '#fb6a4a'))(p_max + 1))
        
        df2 <- data.frame(TermID = df.top$TermID,
                          start = 0,
                          end = df.top$bg_term_num,
                          bg_term_num = df.top$bg_term_num,
                          bg_term_num_col = color_assign(df.top$Enrichment)
                         )
    
        # 根据类型处理前景数据
        if (type == "base") {
            df3 <- data.frame(TermID = df.top$TermID,
                              start = 0,
                              end = df.top$fg_term_num,
                              fg_term_num = df.top$fg_term_num,
                              color = "#ba55d3"
                             )
        } else if (type == "sig") {
            tempLong <- if (ifLog) log10(max(df.top$bg_pro_num,  na.rm  = TRUE) + 1) else max(df.top$bg_pro_num,  na.rm  = TRUE)
        
            df3 <- bind_rows(
                data.frame(TermID = df.top$TermID,
                         start = 0,
                         end = df.top$Up  / (df.top$Up  + df.top$Down)  * tempLong,
                         count = df.top$Up,
                         color = "#69115dA1"
                        ),
                data.frame(TermID = df.top$TermID,
                           start = df.top$Up  / (df.top$Up  + df.top$Down)  * tempLong,
                           end = tempLong,
                           count = df.top$Down,
                           color = "#838bc5"
                          )
                ) %>% dplyr::mutate(count = ifelse(count == 0, "", count))
        }else{
            stop("type must be 'base' or 'sig'")
        }
    
      # 比例数据
      df4 <- data.frame(
        TermID = df.top$TermID,
        start = 0,
        end = max(df.top$bg_pro_num,  na.rm  = TRUE),
        ratio = df.top$Ratio  / max(df.top$Ratio,  na.rm  = TRUE) * 10,
        col = main.col
      )
    
      # 坐标轴对数转换
      if (ifLog) {
        df1 <- dplyr::mutate(df1, end = log10(end + 1))
        df2 <- dplyr::mutate(df2, end = log10(end + 1))
        if (type == "base") df3 <- dplyr::mutate(df3, end = log10(end + 1))
        df4 <- dplyr::mutate(df4, end = log10(end + 1))
      }
    
      # 绘图设置
      par(omi = c(0.1, 0.1, 0.1, 0.1), xpd = NA) # xpd = NA 允许最外圈坐标文字画到默认绘图区之外；较大的
      circos.par(track.margin  = c(0.01, 0.01),
                 circle.margin = c(0.05, 0.05, 0.05, 0.05), # 为扇区首、末端的 0/100/10000 标签预留空间。
                 gap.degree = 0.1
                )
    
      # 初始化轨道
      circos.genomicInitialize(df1,  plotType = "none")
    
      # 第一轨道：分类标签
      circos.trackPlotRegion(
        ylim = c(0, 1),
        panel.fun  = function(x, y) {
          sector.index  <- get.cell.meta.data("sector.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          #circos.text(mean(xlim),  mean(ylim), sector.index,cex = 0.6, facing = "bending.inside",  niceFacing = TRUE)
          if (grepl("^GO:?", sector.index)) {
            # bending.inside 对字符串中的 \n 支持不好，因此将 GO 和编号
            # 分别绘制成两行；增大 cex 时也不容易与相邻扇区重叠。
            circos.text(mean(xlim), 0.8, "GO",cex = 0.8, facing = "inside", niceFacing = TRUE)
            circos.text(mean(xlim), 0.2, sub("^GO:?", "", sector.index),cex = 0.8, facing = "inside", niceFacing = TRUE)
          }else{
            circos.text(mean(xlim), mean(ylim), sector.index,cex = 0.8, facing = "inside", niceFacing = TRUE)
          }          
        },
        track.height  = 0.1,
        bg.border  = NA,
        bg.col  = main.col
      )
    
      # 坐标轴
      if (!ifLog) {
        for (si in get.all.sector.index())  {
          circos.axis(
                h = "top",
                labels.cex  = 1,
                sector.index  = si,
                track.index  = 1,
                major.at  = pretty(c(0, max(df1$end, na.rm  = TRUE)), n = 3),
                labels.facing  = "clockwise",
                labels.pos.adjust = TRUE
          )
        }
      } else {
        for (si in get.all.sector.index())  {
          circos.axis(
                h = "top",
                labels.cex  = 1,
                sector.index  = si,
                track.index  = 1,
                major.at  = c(0, 2, 4),
                labels = c(0,100, 10000),
                labels.facing  = "clockwise",
                labels.pos.adjust = TRUE
          )
        }
      }
    
      # 第二轨道：背景基因数
      circos.genomicTrack(
        df2,
        ylim = c(0, 1),
        track.height  = 0.1,
        bg.border  = "white",
        panel.fun  = function(region, value, ...) {
          circos.genomicRect(region,  value, ytop = 0, ybottom = 1,col = value[, 2], border = NA, ...)
          circos.genomicText(region,  value, y = 0.5, labels = value[, 1],adj = c(0.5, 0.5), cex = 1, ...)
        }
      )
    
      # 第三轨道：前景基因数
      circos.genomicTrack(
        df3,
        ylim = c(0, 1),
        track.height  = 0.1,
        bg.border  = "white",
        panel.fun  = function(region, value, ...) {
          circos.genomicRect(region,  value, ytop = 0, ybottom = 1,col = value[, 2], border = NA, ...)
          circos.genomicText(region,  value, y = 0.5, labels = value[, 1],cex = 1, adj = c(0.5, 0.5), ...)
        }
      )
    
      # 第四轨道：富集比例
      circos.genomicTrack(
        df4,
        ylim = c(0, 10),
        track.height  = 0.35,
        bg.border  = "white",
        bg.col  = "grey90",
        panel.fun  = function(region, value, ...) {
          cell.xlim  <- get.cell.meta.data("cell.xlim")
          cell.ylim  <- get.cell.meta.data("cell.ylim")
          for (j in 1:9) {
            y <- cell.ylim[1]  + (cell.ylim[2]  - cell.ylim[1])  / 10 * j
            circos.lines(cell.xlim,  c(y, y), col = "#FFFFFF", lwd = 0.3)
          }
          circos.genomicRect(region,  value, ytop = 0, ybottom = value[, 1],
                             col = value[, 2], border = NA, ...)
        }
      )
    
      circos.clear()
    
      # 添加图例
      if (type == "base") {
        middle.legend  <- Legend(
            labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
            type = "points",
            pch = c(15, 15, 17),
            legend_gp = gpar(col = c('pink', '#BA55D3', main.col[1])),
            labels_gp = gpar(fontsize = 12),
            title = "",
            nrow = 3,
            size = unit(5, "mm"),
            row_gap = unit(3, "mm")
        )
      } else {
        middle.legend  <- Legend(
            labels = c('Number of Genes', 'Number of Up', 'Number of Down', 'Rich Factor(0-1)'),
            type = "points",
            pch = c(15, 15, 15, 17),
            legend_gp = gpar(col = c('pink', '#69115dA1', '#838bc5', main.col[1])),
            labels_gp = gpar(fontsize = 12),
            title = "",
            nrow = 4,
            size = unit(5, "mm"),
            row_gap = unit(3, "mm")
        )
      }
    
      circle_size <- unit(1, "snpc")
      draw(middle.legend,  x = circle_size * 0.5)
    
      # 主图例
      main.legend  <- Legend(
            labels = unique(df.top$Class),
            type = "points",
            pch = 15,
            legend_gp = gpar(col = classCol),
            labels_gp = gpar(fontsize = 10),
            title = "Class",
            title_position = 'topcenter',
            nrow = 4,
            size = unit(8, "mm"),
            grid_height = unit(6, "mm"),
            grid_width = unit(6, "mm")
      )
    if(identical(kdata_plot_process$pvalue,kdata_plot_process$PValue)){
        legend_pname = "p.value"
    }else if(identical(kdata_plot_process$p.adjust,kdata_plot_process$PValue)){
        legend_pname = "p.adj"
    }
      # P值图例
      logp.legend  <- Legend(
        col_fun = colorRamp2(
          round(seq(0, p_max, length.out  = 6), 0),
          colorRampPalette(c('#fee5d9', '#fb6a4a'))(6)
        ),
        legend_height = unit(3, 'cm'),
        labels_gp = gpar(fontsize = 10),
        title_position = 'topcenter',
        title = paste0("-Log10(",legend_pname,")")
      )
    
      #lgd <- packLegend(main.legend,logp.legend)
      draw(main.legend, x = circle_size * 0.06, y = circle_size * 0.88, just = "left")
      draw(logp.legend, x = circle_size * 0.9, y = circle_size * 0.85, just = "left")  
    }

    
    # circularplot 需要的富集结果数据处理函数
    process_enrichment_data_base <- function(data,pname = "pvalue") {
        
        # 检查必要的列是否存在
        required_cols <- c("ONTOLOGY", "ID", "Description", "GeneRatio", pname, "BgRatio", "GeneRatio_raw")
        missing_cols <- required_cols[!required_cols %in% colnames(data)]
        if (length(missing_cols) > 0) {
            stop("数据缺少以下必需列: ", paste(missing_cols, collapse = ", "))
        }
        
        # 重命名列
        result <- data %>% dplyr::rename("Class" = "ONTOLOGY","TermID" = "ID","Term" = "Description","Ratio" = "GeneRatio")
        
        # 使用基础 R 处理字符串分割
        result$bg_term_num <- sapply(strsplit(result$BgRatio, "/"), `[`, 1) %>% as.numeric()
        result$bg_pro_num <- sapply(strsplit(result$BgRatio, "/"), `[`, 2) %>% as.numeric()
        result$fg_term_num <- sapply(strsplit(result$GeneRatio_raw, "/"), `[`, 1) %>% as.numeric()
        result$fg_pro_num <- sapply(strsplit(result$GeneRatio_raw, "/"), `[`, 2) %>% as.numeric()
        result$PValue <- result[[pname]]
        result$Enrichment <- -log10(result[[pname]])
        
        return(result)
    }
    ###################################################################################################
    ############################################ Enrich ###############################################
    ###################################################################################################
    # GO、KEGG富集分析
    if(p_name == "p.value"){p_name = "pvalue"}else if(p_name == "p.adj"){p_name = "p.adjust"}
    if(species == "9606"){ species_OrgDb <- "org.Hs.eg.db";organism = "hsa"  }else if(species == "10090"){ species_OrgDb <- "org.Mm.eg.db";organism = "mmu" }
    gene_ids <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType ="ENTREZID", OrgDb = species_OrgDb)
  
    #--KEGG--------------------------------------------------------------------------------
    # install.packages('R.utils')
    R.utils::setOption("clusterProfiler.download.method", "auto") # KEGG富集有的时候没法下载数据，需要修改参数
    
    if(kegg_analysis_method == "online"){
        kegg <- clusterProfiler::enrichKEGG(gene = gene_ids$ENTREZID,
                                            organism = organism,
                                            keyType = "kegg",
                                            pvalueCutoff = pvalueCutoff,
                                            qvalueCutoff = qvalueCutoff,
                                            pAdjustMethod = 'BH')
    }else if(kegg_analysis_method == "local"){
        library(KEGG.db)
        kegg <- clusterProfiler::enrichKEGG(gene = gene_ids$ENTREZID,
                                            organism = organism,
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = pvalueCutoff,
                                            qvalueCutoff = qvalueCutoff,
                                            use_internal_data =T)
    }
        
    kegg <- clusterProfiler::setReadable(kegg, species_OrgDb, "ENTREZID")
    sig_kegg <- dplyr::filter(kegg,!!sym(p_name) < filter_pvalue)
    
    kdata <- sig_kegg@result %>% dplyr::mutate(ONTOLOGY = "KEGG")
    kdata$GeneRatio_raw <- kdata$GeneRatio
    kdata$GeneRatio <- sapply(strsplit(kdata$GeneRatio, "/"),function(x) as.numeric(x[1]) / as.numeric(x[2]))
    if(any(grepl(" - Mus musculus \\(house mouse\\)", kdata$Description))){
        kdata$Description <- gsub(" - Mus musculus \\(house mouse\\)","",kdata$Description)
    }                          
    
    #--GO----------------------------------------------------------------------------------------
    ego <- clusterProfiler::enrichGO(
                    gene = gene_ids$ENTREZID,
                    OrgDb = species_OrgDb,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = 'BH',
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    readable = TRUE) 
                     
    sig_ego <- dplyr::filter(ego,!!sym(p_name) < filter_pvalue)
    
    pdata <- sig_ego@result
    pdata$GeneRatio_raw <- pdata$GeneRatio
    pdata$GeneRatio <- sapply(strsplit(pdata$GeneRatio, "/"),function(x) as.numeric(x[1]) / as.numeric(x[2]))
    pdata <- pdata %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")))  # 确保 ONTOLOGY 是因子，并指定排序顺序
    pdata <- pdata %>% dplyr::mutate(Description = factor(Description, levels = unique(Description)))               

    # GO & KEGG ---------------------------------------------------------------------------------------
    pdata_select <- pdata %>% dplyr::select(ONTOLOGY,ID,Description,GeneRatio_raw,GeneRatio,BgRatio,RichFactor,FoldEnrichment,zScore,pvalue,p.adjust,qvalue,geneID,Count)
    kdata_select <- kdata %>% dplyr::select(ONTOLOGY,ID,Description,GeneRatio_raw,GeneRatio,BgRatio,RichFactor,FoldEnrichment,zScore,pvalue,p.adjust,qvalue,geneID,Count)
    pkdata <- rbind(pdata_select,kdata_select)
    pkdata$ONTOLOGY <- factor(pkdata$ONTOLOGY,levels = c("BP","CC","MF","KEGG"))                     
    pkdata <- pkdata %>% dplyr::mutate(Description = factor(Description, levels = unique(Description)))
    pkdata <- pkdata %>% dplyr::arrange(ONTOLOGY,pvalue,desc(GeneRatio),desc(Count))

    #########################################################################################################################                         
    ########################################################## Plot #########################################################
    #########################################################################################################################                          
    if(length(plot_type) != 0){
        
        # 获取Top 数据
        if(length(go_kegg_topn) == 1){
            go_topn = go_kegg_topn
            kegg_topn = go_kegg_topn
        }else if(length(go_kegg_topn) == 2){
            go_topn = go_kegg_topn[1]
            kegg_topn = go_kegg_topn[2]
        }
            
        pdata_plot <- pkdata %>% dplyr::filter(ONTOLOGY %in% c("BP","CC","MF")) %>% group_by(ONTOLOGY) %>% dplyr::arrange(ONTOLOGY,pvalue,desc(GeneRatio),desc(Count)) %>% slice_head(n = go_topn) %>%   
                dplyr::mutate(Description_raw = as.character(Description),Description = stringr::str_wrap(Description_raw, width = description_wrap_width))
        kdata_plot <- pkdata %>% dplyr::filter(ONTOLOGY %in% c("KEGG")) %>% group_by(ONTOLOGY) %>% dplyr::arrange(ONTOLOGY,pvalue,desc(GeneRatio),desc(Count)) %>% slice_head(n = kegg_topn) %>% 
                dplyr::mutate(Description_raw = as.character(Description),Description = stringr::str_wrap(Description_raw, width = description_wrap_width))
        pkdata_plot <- rbind(pdata_plot,kdata_plot) %>% 
                dplyr::mutate(Description_raw = as.character(Description),Description = stringr::str_wrap(Description_raw, width = description_wrap_width))

        # 输出图片列表    
        plot_list <- list(barplot = NULL,dotplot = NULL,dotplot = NULL,sankeyplot = NULL,treeplot = NULL,circularplot = NULL)
        
        # 出图
        if("barplot" %in% plot_type){
            # 条形图
            p_go <- barplot(data = pdata_plot, p_name = p_name, colors = fill_colors)
            p_kegg <- barplot(data = kdata_plot, p_name = p_name, colors = fill_colors)
            p_gokegg <- barplot(data = pkdata_plot, p_name = p_name, colors = fill_colors)

            plot_list$barplot$go = p_go
            plot_list$barplot$kegg = p_kegg
            plot_list$barplot$gokegg = p_gokegg
        }
        if("dotplot" %in% plot_type){
            # 气泡图
            p_go <- dotplot(data = pdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)
            p_kegg <- dotplot(data = kdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)
            p_gokegg <- dotplot(data = pkdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)

            plot_list$dotplot$go = p_go
            plot_list$dotplot$kegg = p_kegg
            plot_list$dotplot$gokegg = p_gokegg
        }
        if("sankeyplot" %in% plot_type){
            # 桑基图
            p_go <- sankeyplot(data = pdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            p_kegg <- sankeyplot(data = kdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            p_gokegg <- sankeyplot(data = pkdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            
            plot_list$sankeyplot$go = p_go
            plot_list$sankeyplot$kegg = p_kegg
            plot_list$sankeyplot$gokegg = p_gokegg           
        }
        # if("treeplot" %in% plot_type){
        #     # 树图
            
            
        #     plot_list$treeplot$go = p_go
        #     plot_list$treeplot$kegg = p_kegg
        #     plot_list$treeplot$gokegg = p_gokegg           
        # }
        if("circularplot" %in% plot_type){
            # 环形图
            pdata_plot_process <- process_enrichment_data_base(pdata_plot,pname = p_name)
            kdata_plot_process <- process_enrichment_data_base(kdata_plot,pname = p_name)
            pkdata_plot_process <- process_enrichment_data_base(pkdata_plot,pname = p_name)
            p_go <- function(){
               circosEnrichmentPlot(pdata_plot_process,topN = go_topn,classCol = fill_colors,
                                     ifLog = T,    # 是否将坐标轴取个log10
                                     type = "base"  # 2中类型，base和sig。sig相比base是前景基因显示上下调
                                    )             
            }

            #p_go <- recordPlot()
            p_kegg <- function(){
                circosEnrichmentPlot(kdata_plot_process,topN = kegg_topn,classCol = fill_colors,
                                     ifLog = T,    # 是否将坐标轴取个log10
                                     type = "base"  # 2中类型，base和sig。sig相比base是前景基因显示上下调
                                    )              
            }

            #p_kegg <- recordPlot()
            p_gokegg <- function(){
                circosEnrichmentPlot(pkdata_plot_process,topN = go_topn,classCol = fill_colors,
                                     ifLog = T,    # 是否将坐标轴取个log10
                                     type = "base"  # 2中类型，base和sig。sig相比base是前景基因显示上下调
                                    )              
            }

            #p_gokegg <- recordPlot()

            plot_list$circularplot$go = p_go
            plot_list$circularplot$kegg = p_kegg
            plot_list$circularplot$gokegg = p_gokegg              
        }
    }else{
        # 不出图
        plot_list <- NULL 
        
    }                          

   return(list(result = pkdata, plot = plot_list))                           
                        
}
