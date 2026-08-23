#' Create volcano, rectangular heatmap, and circular heatmap plots
#'
#' Reclassify a differential-expression result using specified fold-change
#' and p-value thresholds, label top genes, and create a volcano plot plus
#' rectangular and circular heatmaps for the selected genes.
#'
#' @param Res A differential-expression data frame. Row names must contain gene
#'   symbols, and columns must include `logFC_name`, `p.value`, and `p.adj`.
#' @param expr A numeric gene-by-sample expression matrix or data frame whose
#'   row names match the gene symbols in `Res`.
#' @param condition A sample annotation data frame containing a `group` or
#'   `groups` column, ignoring case. Sample names may be in row names or a
#'   `sample` or `samples` column.
#' @param sample_color A character vector of sample colors. Retained by the
#'   interface but not used in the current plotting code.
#' @param logFC_name The name of the fold-change column in `Res`.
#' @param logFC A non-negative absolute fold-change threshold on the log2
#'   scale.
#' @param p_value The significance column used for classification. Use
#'   `"p.adj"` for adjusted p-values; any other value selects `"p.value"`.
#' @param p A numeric significance threshold.
#' @param top_gene Number of top upregulated and top downregulated genes used
#'   for labels and heatmaps.
#' @param heatmap_scale_range A numeric vector of length two giving lower and
#'   upper truncation limits for the scaled rectangular heatmap matrix.
#' @param group_order A character vector specifying the sample-group order and
#'   the group names displayed in plot annotations.
#' @param group_order_color A color vector corresponding to `group_order`.
#' @param change Three labels used, in order, for upregulated, downregulated,
#'   and non-significant genes.
#' @param change_colors Three colors corresponding to `change`.
#' @param heatmap_value_color A color vector used by the rectangular and
#'   circular heatmaps.
#' @param prominent_genes An optional character vector of genes to emphasize
#'   in the volcano plot and include in the heatmaps.
#' @param prominent_genes_color Colors for the regulation classes represented
#'   by `prominent_genes`.
#' @param prominent_genes_size Numeric point-size control for emphasized genes;
#'   their text size is twice this value.
#' @param track_height Numeric height of the circular heatmap track.
#' @param dend_track_height Numeric height of its dendrogram track.
#'
#' @return A named list with `result`, the updated differential-expression
#'   table, and `plot`, a list containing `valcano` (the volcano ggplot),
#'   `rHeatmap` (the ComplexHeatmap object), and `cHeatmap` (a recorded circular
#'   heatmap plot). The historical name `valcano` is retained by the function.
#'
#' @details
#' Genes are classified using strict p-value and absolute fold-change
#' thresholds. Within each direction, genes are ordered first by fold change
#' and then by the selected p-value. Up to `top_gene` genes from each direction
#' are labelled and included in the heatmaps.
#'
#' The rectangular heatmap scales each selected gene across samples and clips
#' values to `heatmap_scale_range`. The circular heatmap clusters the unscaled
#' selected expression matrix, draws directly on the active graphics device,
#' and is captured with [grDevices::recordPlot()].
#'
#' @examples
#' \dontrun{
#' plots <- Volcano_Heatmap_plot(
#'   Res = deg_result,
#'   expr = expression_matrix,
#'   condition = sample_condition,
#'   group_order = c("Control", "Disease"),
#'   group_order_color = c("#FF7F00", "#1177B0"),
#'   logFC = 0.5,
#'   p_value = "p.adj",
#'   p = 0.05,
#'   top_gene = 10
#' )
#' plots$plot$valcano
#' ComplexHeatmap::draw(plots$plot$rHeatmap)
#' replayPlot(plots$plot$cHeatmap)
#' }
#'
#' @seealso [Array_DEGs()], [ggplot2::ggplot()], [ComplexHeatmap::Heatmap()],
#'   [circlize::circos.heatmap()]
#' @export
Volcano_Heatmap_plot <- function(Res = Res,
                                 expr = data$expr,
                                 condition = data$condition,
                                 sample_color = colors,
                                 logFC_name = "logFC",
                                 logFC = 0.5, 
                                 p_value = "p.adj",
                                 p = 0.05,
                                 top_gene = 10,
                                 heatmap_scale_range = c(-2,2),
                                 group_order = data$group_order,   
                                 group_order_color = c("#FF7F00", "#1177b0"),
                                 change = c("Up","Down","Not"),
                                 change_colors = c("orange","#1177b0", "#7b7c7d"),
                                 heatmap_value_color = viridis(10),
                                 prominent_genes = NULL,
                                 prominent_genes_color = c("red"),
                                 prominent_genes_size = 2,
                                 track_height = 0.3,  # 环形热图的热图本体那一圈有多“厚”
                                 dend_track_height = 0.2  # 聚类树那一圈的厚度
                                ){

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
  
    # 差异判断 
    if(p_value == "p.adj"){
        pvalue <- "p.adj"
    }else{
        pvalue <- "p.value"
    }
    names(change_colors) <- change
  
    Res$sig = as.factor(ifelse(Res[[pvalue]] < p & abs(Res[[logFC_name]]) > logFC,
                             ifelse(Res[[logFC_name]] > logFC,
                                    change[1],
                                    change[2]),
                             change[3])) 
    Res$sig <- factor(Res$sig, levels = change)
    Res$Symbols <- rownames(Res)

    #
    up_genes_num <- Res[which(Res$sig == change[1]),] %>% nrow()
    down_genes_num <- Res[which(Res$sig == change[2]),] %>% nrow()
    
    # up/down type     
    up_DEG <- Res[which(Res$sig == change[1]),] #%>% dplyr::filter(!grepl("^ENSG|^LINC|^HP$",Symbols))
    up_DEG <- up_DEG[order(-up_DEG[[logFC_name]] , up_DEG[[pvalue]] ),]
    down_DEG <- Res[which(Res$sig == change[2]),] #%>% dplyr::filter(!grepl("^ENSG|^LINC|^HP$",Symbols))
    down_DEG <- down_DEG[order(down_DEG[[logFC_name]], down_DEG[[pvalue]] ),]
    data_repel_pre <- rbind(up_DEG[1:top_gene,], down_DEG[1:top_gene,]) %>% na.omit()
    
    Res$lable <- ifelse(Res$Symbols %in%  data_repel_pre$Symbols, Res$Symbols,NA)
    if(any(!is.na(prominent_genes))){
    # 检查 prominent_genes_color
    if(all(prominent_genes %in% Res$Symbols )){
      # 特定基因都能匹配到
      prominent_genes_sig_length <-  Res$sig[Res$Symbols %in% prominent_genes] %>% unique() %>% length()
    }else{
      # 部分或全部基因未能匹配到
      cat("\n","\033[31m","The following specific genes could not be matched.","\033[0m","\n")
      print(prominent_genes[!prominent_genes %in% Res$Symbols])
      flush.console()
      stop()
    }
    
    if(length(prominent_genes_color) != prominent_genes_sig_length){
      cat("\n","\033[31m","The number of prominent_genes_color does not match the number of up/down regulation types corresponding to prominent_genes.","\033[0m","\n")
      cat("\n","\033[33m","Please set ",prominent_genes_sig_length," color values.","\033[0m","\n")
      flush.console()
      stop()
    }
    # 凸显特定基因数据处理
    Res$prominent <- ifelse(Res$Symbols %in% prominent_genes,Res$Symbols,NA)
    Res_select <- Res[is.na(Res$prominent),]
    }else{
        Res_select <- Res
    }
  
    data_repel <- Res[data_repel_pre$Symbols,]
  
    # 火山图--------------------------------------------------------------------------------------------------------------------------------
    cat("\n\t",green,"Volcano Plot!",reset,"\n")
    if( nrow(up_DEG) > 0 & nrow(down_DEG) > 0 ){ 
    xintercept <- c(-logFC,logFC) 
    }else if(nrow(up_DEG) > 0){
    xintercept <- c(logFC)
    }else{
    xintercept <- c(-logFC)
    }
  
  ggplot(data = Res_select, aes(x = !!sym(logFC_name), y = -log10(!!sym(pvalue)))) +
    geom_point(aes(color = sig),alpha = 1, size = 1.5,shape = 19) +
    geom_text_repel(aes(label = lable),colour = "black", fontface = 'bold', 
                    size = 2,  max.overlaps = 5000,  position = 'identity')+
    guides(color = guide_legend(title = "",  
                                title.theme = element_text(size = 10, face="bold"), 
                                label.vjust = 0.5, 
                                label.theme = element_text(size = 10), 
                                override.aes = list(size = 3), 
                                ncol = 1
                                ),
    size = "none") +
    labs(x = expression(Log[2]("Fold Change")), y = expression(-Log[10](p_value))) +
    scale_color_manual(values = change_colors)+
    geom_vline(xintercept = xintercept,lty=2,col="black",linewidth=0.4) +
    geom_hline(yintercept = -log10(p),lty=2,col="black",linewidth=0.4) +
    theme_bw(base_size = 12) +
    theme(
      #legend.justification = c(1,1),
      legend.position = "right",
      legend.direction = "vertical",
      legend.background = element_rect(fill = "white", color = "white", linewidth = 0.2),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(face="bold",color="black",size=13),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",size = 20),
      axis.text.x = element_text(face = "bold",color = "black",size = 17),
      axis.text.y = element_text(face = "bold",color = "black",size = 17),
      axis.title.x = element_text(face = "bold",color = "black",size = 20),
      axis.title.y = element_text(face = "bold",color = "black",size = 20),
      plot.subtitle = element_text(hjust = 0.5,size = 14, face = "italic", colour = "black")) +
    labs(x = gsub("_"," ",logFC_name),y = paste0("-log10 (",p_value,")"),
         title = paste0(group_order[1]," vs ",group_order[2]),
         subtitle = paste(sprintf("%s: %s;", p_value, p),
                          sprintf('%s: %s;', gsub("_"," ",logFC_name), logFC),
                          sprintf(paste0(change[1],": %1.0f; ",change[2],": %1.0f;"), up_genes_num, down_genes_num),
                          sprintf('Total: %1.0f', up_genes_num + down_genes_num))) -> p_valcano
  
  if( !is.null(prominent_genes) ){
    # 凸显特定基因
    p_valcano <- p_valcano + ggnewscale::new_scale_color()+
      geom_point(data = Res[!is.na(Res$prominent),],
                 aes(color = sig),alpha = 1, size = prominent_genes_size, shape = 19, show.legend = FALSE) + 
      scale_color_manual(values = prominent_genes_color)+
      geom_text_repel(data = Res[!is.na(Res$prominent),],
                      aes(label = prominent,colour = sig), fontface = 'bold', size = prominent_genes_size*2,  
                      max.overlaps = 5000,  position = 'identity',show.legend = FALSE)+
      scale_colour_manual(values = prominent_genes_color)
    
  }

  # 矩形热图-----------------------------------------------------------------------------------------------------------------------------------------
      cat("\n\t",green,"Rectangular  Heatmap Plot!",reset,"\n")
      if( any( !is.na(prominent_genes) )){
        data_repel_specific <- Res[Res$Symbol %in% prominent_genes,]
        specific_no_Top_genes <- setdiff(data_repel_specific$Symbols, data_repel$Symbols)
        data_psecific_no_top <- Res[Res$Symbol %in% specific_no_Top_genes,]
        gene <- do.call(rbind,list(data_repel,data_psecific_no_top)) 
      }else{
        gene <- data_repel
      }  
      
      select.expr <- expr[rownames(gene),]
      if( !any(grepl("^sample$|^samples$",tolower(colnames(condition)))) ){
        condition$sample <- rownames(condition)
      }else{
        colnames(condition)[grep("^sample$|^samples$",tolower(colnames(condition)))] <- "sample"
      }
      
      colnames(condition)[grep("^group$|^groups$",tolower(colnames(condition)))] <- "Group"
      
      mat <- t(scale(t(select.expr))) #归一化
      mat[mat < heatmap_scale_range[1]] <- heatmap_scale_range[1]
      mat[mat > heatmap_scale_range[2]] <- heatmap_scale_range[2]
      
      condition$Group <- factor(condition$Group, levels = group_order)
      condition <- condition[order(condition$Group), ]
      annotation_col <- data.frame(row.names = condition$sample, Group = condition$Group)
      annotation_row <- data.frame(row.names = rownames(data_repel), sig = data_repel$sig)
      
      group_order_colors <- setNames(group_order_color, group_order)
      
      mat <- mat[,rownames(annotation_col)]  # 对样本按照分组顺序排列
      
      # densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
      #   HeatmapAnnotation(Group = annotation_col$Group, col = list(Group = group_order_colors)) %v%
      #   #Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,show_row_names = T,name = "mat", height = unit(8, "cm"),col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
      #   Heatmap(mat, row_names_gp = gpar(fontsize = 9),show_column_names = F,show_row_names = T,name = "mat", height = unit(8, "cm"), col = heatmap_value_color) -> hp
      p_hp <- densityHeatmap(mat,title = "Distribution as heatmap",ylab = " ",height = unit(6, "cm")
                             # column_names_gp = gpar(fontfamily = "Times"),
                             # heatmap_legend_param = list(title_gp = gpar(fontfamily = "Times"),
                             #                             labels_gp = gpar(fontfamily = "Times")
                             #                            )
                          ) %v% 
            HeatmapAnnotation(Group = annotation_col$Group,
                              col = list(Group = group_order_colors)
                              # annotation_name_gp = gpar(fontfamily = "Times"),
                              # gp = gpar(fontfamily = "Times")
                             ) %v% 
            Heatmap(mat,
                    # row_names_gp = gpar(fontsize = 9, fontfamily = "Times"),
                    # column_names_gp = gpar(fontfamily = "Times"),
                    show_column_names = FALSE,
                    show_row_names = TRUE,
                    name = "mat",
                    height = unit(8, "cm"),
                    col = heatmap_value_color
                    # heatmap_legend_param = list(title_gp = gpar(fontfamily = "Times"),
                    #                             labels_gp = gpar(fontfamily = "Times")
                    #                            )
                   )
    
    # 环形热图-------------------------------------------------------------------------------------------------------------------------------------
    cat("\n\t",green,"Circular Heatmap Plot!",reset,"\n")
    viridis_colors <- heatmap_value_color
    
    mat <- expr[rownames(expr) %in% rownames(gene),]
    group_df <- condition
    group_df$sample <- NULL
    annotation_colors <- list(group = group_order_colors)
    
    pheatmap_result <- pheatmap::pheatmap(mat,silent = T)
    gene_ordered_names <- rownames(mat)[pheatmap_result$tree_row$order]
    sample_ordered_names <- colnames(mat)[pheatmap_result$tree_col$order]
    mat <- mat[pheatmap_result$tree_row$order, pheatmap_result$tree_col$order]
    
    circos_collor_fun <- circlize::colorRamp2(seq(min(mat), max(mat), length.out = 10), heatmap_value_color)
    mat_use <- mat[rownames(mat) %in% rownames(gene),] %>% t()
    
    # plot
    group_df_sample_order <- condition[sample_ordered_names,]
    rownames_colors <- sapply(group_df_sample_order$Group, function(group) { group_order_colors[group] })
    names(rownames_colors) <- sample_ordered_names
    
    heatmap_circos_plot <- function(mat_use, circos_collor_fun, rownames_colors, group_colors){
        #grob_obj <- grid::grid.grabExpr({
            scale_value <- function(x, min_val = 6, max_val = 50, min_scale = 1, max_scale = 0.8) {
              # 定义根据样本数计算样本字体大小函数
              if (x <= min_val) {
                return(min_scale)
              }
              if (x >= max_val) {
                return(max_scale)
              }
              scaled_value <- min_scale + (max_scale - min_scale) * (x - min_val) / (max_val - min_val)
              return(scaled_value)
            }
            sample_text_size <- scale_value(length(rownames_colors))
            
            
            circlize::circos.clear() 
            # 设置 circos 参数
            circos.par(start.degree = 45, gap.after = 90,
                       track.margin = c(0.001,0.001),
                       cell.padding = c(0,0,0,0),
                       canvas.xlim = c(-0.9, 0.9), 
                       canvas.ylim = c(-0.9, 0.9)
                      )
            circos.heatmap(mat_use,
                           col = circos_collor_fun,
                           track.height = track_height,
                           cluster = TRUE,
                           dend.side = "inside",
                           dend.track.height = dend_track_height,
                           dend.callback = function(dend,m,si) { color_branches(dend, k = ceiling(sqrt(length(rownames_colors)^1.7)), 
                                                                                col = 1:ceiling(sqrt(length(rownames_colors)^1.7))) },
                           rownames.side = "outside",
                           rownames.cex = sample_text_size,
                           rownames.col = rownames_colors
            ) 
            
            circos.track(track.index = 2, # 将列名添加在第二个轨道（就是热图所在的环形轨道）
                         panel.fun = function(x, y) {
                           if(CELL_META$sector.numeric.index == 1) { # the last sector
                             cn = colnames(mat_use) %>% rev() # 取得列名 
                             n = length(cn)
                             circos.text(x = rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), # x轴坐标
                                         y = 1:n - convert_y(1.5, "mm"), # y轴坐标
                                         cn, # 输入要展示的列名
                                         cex = 0.6, # 列名(基因名)的大小
                                         adj = c(0, 0.5),
                                         facing = "inside",
                                         gp = gpar(fontsize = 12)
                                        )
                           }
                         }, bg.border = NA)
            
            marker_exp_legend = Legend(title = "",
                                       col_fun = circos_collor_fun,
                                       direction = c( "horizontal"),
                                       grid_height = unit(0.15, "npc"),
                                       legend_width = unit(0.25, "npc"),
                                       labels_gp = gpar(fontsize = 14),
                                       tick_length = unit(0.01, "npc")
            )
            pushViewport(viewport(x = 0.5, y = 0.7, width = 2.8, height = 0.1, just = c("center", "bottom")))
            grid.draw(marker_exp_legend)
            popViewport()
            
            graphics_list <- function(col_list){
              return(
                lapply(1:length(col_list), function(i) {
                  #num_text <- tsne_merge$cluster_new[which(tsne_merge$celltype2 %in% names(col_list))] %>% unique() %>% gsub("C","",.)
                  
                  function(x, y, w, h) {
                    grid.points(x , y, gp = gpar(col = col_list[i],fill = col_list[i]), pch = 21, size = unit(0.2,'npc'))
                    #grid.text(num_text[i], x = x + unit(-6.8, "points"), y = y, just = "left", gp = gpar(col = "black"))
                  }
                })
              )
            }
            
            grouptype_legend <- Legend(title  = "", 
                                       ncol = 2,
                                       labels = group_order, 
                                       legend_gp = gpar(col = group_colors),
                                       labels_gp = gpar(fontsize = 12),
                                       graphics = graphics_list(group_colors),  # 设置类型
                                       title_gp = gpar(fontsize = 16,box_fill = "#C199BA",box_col  = "#C199BA"),
                                       title_position = "leftcenter",
                                       title_gap = unit(0.03, "npc"),
                                       grid_height = unit(0.01, "npc"), # 每一个 legend的高度
                                       grid_width = unit(0.01, "npc"),
                                       legend_height  = unit(0.01,"npc"),
                                       row_gap   = unit(0.03, "npc"),  
                                       column_gap = unit(0.03,"npc")
                                       #pch = 20,
                                       #type = 'points',
                                       #size  = unit(10,"points")
                                       #background = "white"  # grid background
                                      )
            pushViewport(viewport(x = 0.5, y = 0.7, width = 1, height = 0.1, just = c("center", "bottom")))
            grid.draw(grouptype_legend)
            popViewport()
            
            pushViewport(viewport(x = 0.5, y = 0.95, width = 1, height = 0.1, just = c("center", "top")))
            grid.text(paste0("Heatmap Circos Top",top_gene," DEGs"), gp = gpar(fontsize = 18, fontface = "bold"))
            popViewport()
        #})
        
        #return(grob_obj) 
      }

      heatmap_circos_plot(mat_use = mat_use,
                          circos_collor_fun = circos_collor_fun,
                          rownames_colors = rownames_colors,
                          group_colors = group_order_color)
      p_chp <- recordPlot()

      return(list(result = Res,plot = list(valcano = p_valcano,rHeatmap = p_hp,cHeatmap = p_chp)))
}
