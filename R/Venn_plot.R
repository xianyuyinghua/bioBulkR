#' Draw a Venn diagram and return genes shared by all sets
#'
#' Create a Venn diagram from a named collection of gene vectors, optionally
#' format region labels as percentages and counts, and calculate the
#' intersection shared by every supplied gene set.
#'
#' @param target_genes A list of vectors containing the genes or other values
#'   in each set. At least two sets are expected.
#' @param target_genes_names A character vector of display names corresponding
#'   to the elements of `target_genes`.
#' @param fill_colors A character vector of fill colors. It must contain at
#'   least one color for every element of `target_genes`.
#' @param base_size Numeric size passed to `VennDiagram::venn.diagram()` for
#'   set-category labels.
#' @param label_size Numeric size used for region count or percentage labels.
#' @param sigdigs Numeric control retained by the function interface. The
#'   current implementation uses two significant digits internally.
#' @param digits `NULL`, or a non-negative integer specifying the decimal
#'   places used when recalculating percentage labels. This is applied only
#'   when `print_mode` contains `"percent"`.
#' @param print_mode A character vector passed to
#'   [VennDiagram::venn.diagram()] to select region-label content, commonly
#'   `"raw"`, `"percent"`, or `c("percent", "raw")`.
#' @param font A font-family name used for region and set-category labels.
#'
#' @return A named list with `common`, the intersection of all vectors in
#'   `target_genes`, and `plot`, a grid grob containing the Venn diagram.
#'
#' @details
#' Diagrams are drawn without proportional scaling. For two sets, category
#' labels are positioned at 180 degrees; diagrams with more sets use the
#' defaults from `VennDiagram`. When percentage output and `digits` are both
#' requested, the text grobs are rewritten to display the recalculated
#' percentage followed by the raw count on a new line. Label line height is
#' set to 0.7.
#'
#' The function sets the `VennDiagramLogger` threshold in `futile.logger` to
#' `ERROR` to suppress lower-priority VennDiagram log messages.
#'
#' @examples
#' genes <- list(
#'   DEGs = c("TP53", "EGFR", "MYC"),
#'   Targets = c("EGFR", "MYC", "AKT1")
#' )
#'
#' venn <- Venn_plot(
#'   target_genes = genes,
#'   target_genes_names = names(genes),
#'   fill_colors = c("#F89FA8", "#F9E9A4"),
#'   print_mode = c("percent", "raw"),
#'   digits = 1
#' )
#' venn$common
#' grid::grid.draw(venn$plot)
#'
#' @seealso [VennDiagram::venn.diagram()], [grid::grid.draw()]
#' @export
Venn_plot <- function(target_genes = list(data$DEGs,data$target_gene),
                      target_genes_names = c("DEGs",data$target_gene_name),
                      fill_colors = c("#F89FA8","#F9E9A4"),
                      base_size = 2,
                      label_size = 2,
                      sigdigs = 3,
                      digits = NULL,
                      print_mode = c("percent", "raw"),
                      font = "Arial"
                     ){
    
      library(stringr)
    
      if( length(fill_colors) < length(target_genes) ){
          stop("colors less")
      }else{
        fill_color = fill_colors[1:length(target_genes)]
      }
      
      colours = rep("white",length(target_genes))
    
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")  # 抑制日志信息
    
      if(length(target_genes) == 2){
        venn_ploy <- VennDiagram::venn.diagram(x = setNames(target_genes,target_genes_names),
                                             filename = NULL,
                                             scaled = FALSE ,
                                             fill = fill_color,
                                             print.mode = print_mode,
                                             col = fill_color,
                                             sigdigs = 2,
                                             cat.pos = c(180,180),
                                             margin = 0.1,
                                             lwd = 0,
                                             cex = label_size,
                                             cat.cex = base_size,
                                             disable.logging = TRUE,
                                               fontfamily = font,
                                               cat.fontfamily = font
                                              )
      }else{
        venn_ploy <- VennDiagram::venn.diagram(x = setNames(target_genes,target_genes_names),
                                               filename = NULL,
                                                 scaled = FALSE ,
                                                 fill = fill_color,
                                                 print.mode = print_mode,
                                                 col = fill_color,
                                                 margin = 0.1,
                                                 sigdigs = 2,
                                                 lwd = 0,
                                                 cex = label_size,
                                                 cat.cex = base_size,
                                                 disable.logging = TRUE,
                                               fontfamily = font,
                                               cat.fontfamily = font
                                              ) 
      }


      # 根据小数位参数进行设置小数位数
      if("percent" %in% print_mode){
          if(!is.null(digits)){
                # 找出所有标签图层的位置
                label_indices <- which(sapply(venn_ploy, function(x) inherits(x, "text")))
                
                # 查看标签原内容
                original_labels <- sapply(venn_ploy[label_indices], function(x) x$label)
                
                counts <- as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", original_labels))
                counts[is.na(counts)] <- 0
                total <- sum(counts[counts > 0])
                
                # 替换百分比标签
                new_labels <- ifelse(
                    counts == 0,
                    original_labels,  # 保留原始标签（如"A", "B"）
                    paste0(format(round(counts / total * 100, digits), nsmall = digits), "%\n(", counts, ")")
                )
                for (i in seq_along(label_indices)) {
                    venn_ploy[[label_indices[i]]]$label <- new_labels[i]
                }
                
          }
      }
      # 控制换行后的行距
      text_indices <- which(sapply(venn_ploy, function(x) inherits(x, "text")))
      for(i in text_indices){
            venn_ploy[[i]]$gp$lineheight <- 0.7
      }
      venn_grob <- grid::grobTree(children = do.call(grid::gList, venn_ploy))                                    

      return(list(common = Reduce(intersect, target_genes),plot = venn_grob))

}
