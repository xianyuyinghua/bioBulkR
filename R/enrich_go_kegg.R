#' Perform GO and KEGG enrichment analysis and create enrichment plots
#'
#' Converts gene symbols to Entrez identifiers, performs Gene Ontology (GO)
#' and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation
#' analysis, combines the significant results, and optionally creates several
#' publication-oriented visualizations.
#'
#' @param genes A character vector of gene symbols, for example
#'   `c("TP53", "EGFR", "MYC")`. Duplicate and unmapped symbols are handled
#'   by [clusterProfiler::bitr()]; genes that cannot be mapped are omitted from
#'   the enrichment analysis.
#' @param species A single character string containing an NCBI taxonomy ID.
#'   Currently supported values are `"9606"` for human and `"10090"` for
#'   mouse. Human uses `org.Hs.eg.db` with KEGG organism code `"hsa"`; mouse
#'   uses `org.Mm.eg.db` with KEGG organism code `"mmu"`.
#' @param p_name A character string selecting the significance column used to
#'   filter and display results. Use `"p.value"` for the raw p-value or
#'   `"p.adj"` for the Benjamini-Hochberg adjusted p-value. The default is
#'   `"p.value"`.
#' @param filter_pvalue A single numeric value between 0 and 1. After GO and
#'   KEGG enrichment, terms with the column selected by `p_name` strictly less
#'   than this value are retained. The default is `0.05`.
#' @param pvalueCutoff A single numeric value between 0 and 1 passed to
#'   [clusterProfiler::enrichGO()] and [clusterProfiler::enrichKEGG()] as their
#'   initial raw p-value cutoff. It should normally be greater than or equal to
#'   `filter_pvalue`. The default is `0.5`.
#' @param qvalueCutoff A single numeric value between 0 and 1 passed to the
#'   clusterProfiler enrichment functions as the q-value cutoff. The default
#'   is `1`, which does not impose an additional restrictive q-value filter.
#' @param kegg_analysis_method A character string specifying the source of KEGG
#'   annotations. `"online"` (default) queries the current KEGG data through
#'   clusterProfiler. `"local"` uses the legacy internal data supplied by
#'   `KEGG.db` and requires that package to be installed.
#' @param plot_type `NULL` to skip plotting, or a character vector containing
#'   one or more of `"barplot"`, `"dotplot"`, `"sankeyplot"`, and
#'   `"circularplot"`. Multiple plot families can be requested in one call,
#'   for example `c("barplot", "dotplot")`. Each family contains separate
#'   GO, KEGG, and combined GO/KEGG outputs. Circular plot entries are drawing
#'   functions and must be invoked, whereas the other entries are ggplot
#'   objects.
#' @param scankeyplot_mode A character string controlling the Sankey layout.
#'   Use `"sankey_buble"` (default) to place genes on the left, pathways in the
#'   middle, and the enrichment bubbles on the right. Use `"buble_sankey"` to
#'   place bubbles and pathways on the left and genes on the right. The argument
#'   name and values retain the historical spelling for backward compatibility.
#' @param fill_colors A character vector of colors used for ontology fills and
#'   continuous plot scales. By default colors are obtained from
#'   `basicR::get_colors(number = 4.8)`. Supply enough distinct colors for all
#'   requested ontology classes, particularly for circular plots.
#' @param go_kegg_topn A positive numeric vector of length one or two controlling
#'   the number of terms shown in each plot. A length-one value is used for both
#'   analyses. With two values, the first is the number retained separately for
#'   each GO ontology (BP, CC, and MF), and the second is the number retained for
#'   KEGG. This affects plots only; the complete filtered table is returned in
#'   `result`. The default is `c(10, 10)`.
#' @param description_wrap_width A positive number giving the approximate number
#'   of characters at which pathway and ontology descriptions are wrapped in
#'   plots. The default is `100`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{result}{A data frame containing the combined, filtered GO and KEGG
#'   results. It includes ontology, term identifiers and descriptions, gene and
#'   background ratios, enrichment statistics, p-values, mapped gene symbols,
#'   and gene counts.}
#'   \item{plot}{`NULL` when `plot_type = NULL`; otherwise, a nested list keyed
#'   by plot family. Each requested family contains `go`, `kegg`, and `gokegg`.
#'   Bar, dot, and Sankey entries are ggplot objects. Circular entries are
#'   zero-argument functions that draw the plot when called.}
#' }
#'
#' @details
#' Input identifiers must be gene symbols. The function uses Benjamini-Hochberg
#' correction in both enrichment analyses. GO enrichment is performed for BP,
#' CC, and MF together. Only human and mouse annotation packages are currently
#' supported. Online KEGG analysis requires network access at run time.
#'
#' Sankey plots require `ggsankey`, `ggnewscale`, and `gground`. Circular plots
#' require `circlize` and `ComplexHeatmap`. The `basicR` package supplies the
#' default color palette.
#'
#' @examples
#' \dontrun{
#' genes <- c("TP53", "EGFR", "MYC", "CDKN1A", "BAX")
#'
#' # Analysis only
#' res <- enrich_go_kegg(
#'   genes = genes,
#'   species = "9606",
#'   p_name = "p.adj",
#'   filter_pvalue = 0.05
#' )
#' head(res$result)
#'
#' # Analysis with bar and dot plots
#' res <- enrich_go_kegg(
#'   genes = genes,
#'   species = "9606",
#'   plot_type = c("barplot", "dotplot"),
#'   go_kegg_topn = c(10, 10)
#' )
#' res$plot$barplot$gokegg
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @import grid
#' @import stringr
#' @import ggsankey
#' @import scales
#' @import circlize
#' @import ComplexHeatmap
#' @importFrom graphics par
#' @importFrom grDevices colorRampPalette
#' @export
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
    if (!is.character(genes) || length(genes) == 0L || anyNA(genes)) {
        stop("`genes` must be a non-empty character vector without missing values.",
             call. = FALSE)
    }
    if (!is.character(species) || length(species) != 1L ||
        !species %in% c("9606", "10090")) {
        stop("`species` must be \"9606\" (human) or \"10090\" (mouse).",
             call. = FALSE)
    }
    if (!is.character(p_name) || length(p_name) != 1L ||
        !p_name %in% c("p.value", "p.adj")) {
        stop("`p_name` must be \"p.value\" or \"p.adj\".", call. = FALSE)
    }
    if (!is.numeric(filter_pvalue) || length(filter_pvalue) != 1L ||
        is.na(filter_pvalue) || filter_pvalue < 0 || filter_pvalue > 1) {
        stop("`filter_pvalue` must be a single number between 0 and 1.",
             call. = FALSE)
    }
    if (!is.numeric(pvalueCutoff) || length(pvalueCutoff) != 1L ||
        is.na(pvalueCutoff) || pvalueCutoff < 0 || pvalueCutoff > 1) {
        stop("`pvalueCutoff` must be a single number between 0 and 1.",
             call. = FALSE)
    }
    if (!is.numeric(qvalueCutoff) || length(qvalueCutoff) != 1L ||
        is.na(qvalueCutoff) || qvalueCutoff < 0 || qvalueCutoff > 1) {
        stop("`qvalueCutoff` must be a single number between 0 and 1.",
             call. = FALSE)
    }
    if (!is.character(kegg_analysis_method) ||
        length(kegg_analysis_method) != 1L ||
        !kegg_analysis_method %in% c("online", "local")) {
        stop("`kegg_analysis_method` must be \"online\" or \"local\".",
             call. = FALSE)
    }
    allowed_plots <- c("barplot", "dotplot", "sankeyplot", "circularplot")
    if (!is.null(plot_type) &&
        (!is.character(plot_type) || any(!plot_type %in% allowed_plots))) {
        stop("`plot_type` must be NULL or contain only: ",
             paste(allowed_plots, collapse = ", "), ".", call. = FALSE)
    }
    if (!is.character(scankeyplot_mode) || length(scankeyplot_mode) != 1L ||
        !scankeyplot_mode %in% c("sankey_buble", "buble_sankey")) {
        stop("`scankeyplot_mode` must be \"sankey_buble\" or \"buble_sankey\".",
             call. = FALSE)
    }
    if (!is.numeric(go_kegg_topn) || !length(go_kegg_topn) %in% c(1L, 2L) ||
        anyNA(go_kegg_topn) || any(go_kegg_topn <= 0)) {
        stop("`go_kegg_topn` must contain one or two positive numbers.",
             call. = FALSE)
    }
    if (!is.numeric(description_wrap_width) ||
        length(description_wrap_width) != 1L ||
        is.na(description_wrap_width) || description_wrap_width <= 0) {
        stop("`description_wrap_width` must be a single positive number.",
             call. = FALSE)
    }
    
    ###################################################################################################
    ########################################### Function ##############################################
    ###################################################################################################
    # Barplot
    barplot <- function(data,p_name = "pvalue",colors = basicR::get_colors(number = 4.8)){
        if( all(unique(data$ONTOLOGY) %in% c('BP', 'CC', 'MF', 'KEGG')) ){
            data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = rev(c('BP', 'CC', 'MF', 'KEGG'))))
        }else if( all(unique(data$ONTOLOGY) %in% c('BP', 'CC', 'MF')) ){
            data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = rev(c('BP', 'CC', 'MF'))))
        }else{
           data <- data %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY,  levels = ONTOLOGY )) 
        }
        data <- data %>% dplyr::arrange(ONTOLOGY,desc(!!sym(p_name)), Count) %>% dplyr::mutate(Description = factor(Description, levels = unique(Description)))
        data <- data %>% tibble::rowid_to_column('index')

        data <- data %>% dplyr::mutate(geneID_trunc = purrr::map_chr(str_split(geneID, "/"), ~ if (length(.) > 12){ paste(.[1:12], collapse = "/")}else{paste(., collapse = "/")}))
              
        xaxis_max <- max(-log10(data[[p_name]])) + 1
        width <- max(-log10(data[[p_name]])) * 0.05
        rect.data <- group_by(data, ONTOLOGY) %>% reframe(n = n()) %>% ungroup() %>% dplyr::mutate(xmin = -3 * width,xmax = -2 * width,ymax = cumsum(n),ymin = lag(ymax, default = 0) + 0.6,ymax = ymax + 0.4)
        
        ggplot2::ggplot(data = data,aes(-log10(!!sym(p_name)), y = index, fill = ONTOLOGY)) +
            gground::geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
            ggplot2::geom_text(aes(x = 0.05, label = Description),hjust = 0, size = 6,lineheight = 0.6
                               #fontface = "bold"
                              ) +
            ggplot2::geom_text(aes(x = 0.1, label = geneID_trunc, colour = ONTOLOGY), hjust = 0, vjust = 2.6, size = 4, fontface = 'italic',show.legend = FALSE) +
            ggplot2::geom_point(aes(x = -width, size = Count),shape = 21) +
            ggplot2::geom_text(aes(x = -width, label = Count),size = 5) +
            scale_size_continuous(name = 'Count', range = c(5, 8)) +
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
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
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
        
        data <- data %>% as.data.frame(check.names = F)   
        df_long <- data %>% dplyr::select(Description, geneID) %>% tidyr::separate_rows(geneID, sep = "/") %>% dplyr::rename(pathway = Description, gene = geneID)
        sankey_df <- df_long %>% ggsankey::make_long(pathway, gene)
        
        sankey_num <- sankey_df %>% dplyr::mutate(x_num = dplyr::case_when(x == "pathway" ~ pathway_location,x == "gene" ~ gene_location,TRUE ~ NA_real_),
                                                  next_x_num = dplyr::case_when(next_x == "pathway" ~ pathway_location,next_x == "gene" ~ gene_location,TRUE ~ NA_real_),
                                                  pathway_fill = node
                                                 )
    
        p_tmp <- ggplot(sankey_num,aes(x = x_num,next_x = next_x_num,node = node,next_node = next_node,fill = pathway_fill)) +
                    geom_sankey(flow.alpha = 0.65,smooth = 8, width = 0.08, show.legend = FALSE) +
                    scale_fill_manual(values = flow_colors) +
                    theme_void()
        pb <- ggplot_build(p_tmp)
        node_df <- pb$data[[2]] %>% dplyr::distinct(node, x, xmin, xmax, ymin, ymax, fill)
    
        if(plot_mode == "sankey_buble"){
            left_nodes <- node_df %>% dplyr::filter(x == min(x)) %>% dplyr::mutate(x_text = xmin - 0.03,y = (ymin + ymax) / 2,
                                                                                   label_wrap = ifelse(nchar(node) > 10, str_wrap(node, width = 10), node)
                                                                                  )
            right_nodes <- node_df %>% dplyr::filter(x == max(x)) %>% dplyr::mutate(x_text = xmax + 0.05,y = (ymin + ymax) / 2,
                                                                                  label_wrap = ifelse(nchar(node) > description_width, str_wrap(node, width = description_width), node)
                                                                                 )
        }else if(plot_mode == "buble_sankey"){
            left_nodes <- node_df %>% dplyr::filter(x == min(x)) %>% dplyr::mutate(x_text = xmin,y = (ymin + ymax) / 2,
                                                                                   label_wrap = ifelse(nchar(node) > description_width,str_wrap(node, width = description_width),node)
                                                                                  )
            right_nodes <- node_df %>% dplyr::filter(x == max(x)) %>% dplyr::mutate(x_text = (xmin + xmax) / 2,y = (ymin + ymax) / 2,
                                                                                    label_wrap = ifelse(nchar(node) > 10,str_wrap(node, width = 10),node)
                                                                                   )
        }
    
        if(plot_mode == "sankey_buble"){    
            bubble_df <- data %>% dplyr::inner_join(right_nodes %>% dplyr::transmute(Description = node,y,ymin,ymax,pathway_fill = fill,label_wrap),by = "Description")
        }else if(plot_mode == "buble_sankey"){
            bubble_df <- data %>% dplyr::inner_join(left_nodes %>% dplyr::transmute(Description = node,y,ymin,ymax,pathway_fill = fill,label_wrap),by = "Description")
        }
    
        rescale_to <- function(x, from, to) {(x - from[1]) / (from[2] - from[1]) * (to[2] - to[1]) + to[1]}
        bubble_df <- bubble_df %>% dplyr::mutate(bubble_x = rescale_to(RichFactor,from = c(min(RichFactor), max(RichFactor)),to = c(bubble_x_min, bubble_x_max)))
    
        bubble_box_xmin <- bubble_x_min - 0.05
        bubble_box_xmax <- bubble_x_max + 0.05
        
        if(plot_mode == "sankey_buble"){
            bubble_box_ymin <- min(right_nodes$ymin)
            bubble_box_ymax <- max(right_nodes$ymax) * 1.02               
        }else if(plot_mode == "buble_sankey"){
            bubble_box_ymin <- min(left_nodes$ymin)
            bubble_box_ymax <- max(left_nodes$ymax) * 1.02   
        }
    
        axis_y <- bubble_box_ymin * 1.03
        
        x_breaks_raw <- pretty(data$RichFactor, n = 3)
        x_breaks_raw <- x_breaks_raw[x_breaks_raw >= min(data$RichFactor) & x_breaks_raw <= max(data$RichFactor)]
        
        x_breaks_plot <- rescale_to(x_breaks_raw,from = c(min(data$RichFactor), max(data$RichFactor)),to = c(bubble_x_min, bubble_x_max))
        
        axis_df <- data.frame(x_raw  = x_breaks_raw,x_plot = x_breaks_plot)
    
        p_final <- ggplot() +
            geom_sankey(data = sankey_num,aes(x = x_num,next_x = next_x_num,node = node,next_node = next_node,fill = pathway_fill),flow.alpha = 0.6,smooth = 9,width = 0.05,show.legend = FALSE) +
            scale_fill_manual(values = flow_colors) +
            ggnewscale::new_scale_color() +
            geom_point(data = bubble_df,aes(x = bubble_x,y = y,size = Count,color = p.adjust),alpha = 0.95) +
            scale_color_gradientn(colours = fill_colors,guide = guide_colorbar(reverse = TRUE,barheight = grid::unit(5, "lines"),barwidth = grid::unit(1.2, "lines"))) +
            scale_size_continuous(range = c(4, 8),name = "Size") +
            annotate("segment",x = bubble_box_xmin,xend = bubble_box_xmax,y = axis_y,yend = axis_y,linewidth = 0.8,color = "black") +
            geom_segment(data = axis_df,aes(x = x_plot,xend = x_plot,y = axis_y,yend = axis_y - 0.5),inherit.aes = FALSE,linewidth = 0.8,color = "black") +
            geom_text(data = axis_df,aes(x = x_plot,y = axis_y - 0.7,label = label_number(accuracy = 0.01)(x_raw)),inherit.aes = FALSE,size = 6,vjust = 1) +
            annotate("text",x = (bubble_box_xmin + bubble_box_xmax) / 2,y = axis_y - (2.5+(nrow(data)*0.1)), label = "Ratio",size = 8) +
            geom_text(data = left_nodes,aes(x = x_text,y = y,label = label_wrap),hjust = hjust_left,nudge_x = nudge_x_left,vjust = 0.5,size = size_left,lineheight = 0.6,fontface = fontface_left) +
            geom_text(data = right_nodes,aes(x = x_text,y = y,label = label_wrap),hjust = hjust_right,nudge_x = nudge_x_right,vjust = 0.5,size = size_right,lineheight = 0.6,fontface = fontface_right) +
            coord_cartesian(xlim = xlim,ylim = c(min(node_df$ymin) - 14, max(node_df$ymax)),clip = "off") +
            theme_void(base_size = 20) +
            theme(legend.position = legend_position,legend.box = "vertical",legend.key.height = grid::unit(0.45, "lines"),
                  legend.spacing.y = grid::unit(0.05, "lines"),legend.margin = ggplot2::margin(0, 0, 0, 0),plot.margin = plot_margin
                 )
    
        if(buble_panel_border){
            p_final <- p_final + 
                annotate("rect",xmin = bubble_box_xmin,xmax = bubble_box_xmax,ymin = axis_y,ymax = bubble_box_ymax,fill = NA,color = "black",linewidth = 0.8)
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
    
        df.top  <- df %>% group_by(Class) %>% slice_min(PValue, n = topN, with_ties = FALSE) %>% dplyr::arrange(desc(Ratio), .by_group = TRUE) %>% ungroup()
    
        main.col  <- classCol[as.numeric(as.factor(df.top$Class))]
    
        df1 <- data.frame(TermID  = df.top$TermID,  start = 0, end = max(df.top$bg_pro_num,  na.rm  = TRUE))
        
        p_max <- round(max(df.top$Enrichment))  + 1
        color_assign <- colorRamp2(0:p_max, colorRampPalette(c('#fee5d9', '#fb6a4a'))(p_max + 1))
        
        df2 <- data.frame(TermID = df.top$TermID,
                          start = 0,
                          end = df.top$bg_term_num,
                          bg_term_num = df.top$bg_term_num,
                          bg_term_num_col = color_assign(df.top$Enrichment)
                         )
    
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
    
      df4 <- data.frame(
        TermID = df.top$TermID,
        start = 0,
        end = max(df.top$bg_pro_num,  na.rm  = TRUE),
        ratio = df.top$Ratio  / max(df.top$Ratio,  na.rm  = TRUE) * 10,
        col = main.col
      )
    
      if (ifLog) {
        df1 <- dplyr::mutate(df1, end = log10(end + 1))
        df2 <- dplyr::mutate(df2, end = log10(end + 1))
        if (type == "base") df3 <- dplyr::mutate(df3, end = log10(end + 1))
        df4 <- dplyr::mutate(df4, end = log10(end + 1))
      }
    
      par(omi = c(0.1, 0.1, 0.1, 0.1), xpd = NA)
      circos.par(track.margin  = c(0.01, 0.01),
                 circle.margin = c(0.05, 0.05, 0.05, 0.05),
                 gap.degree = 0.1
                )
    
      circos.genomicInitialize(df1,  plotType = "none")
    
      circos.trackPlotRegion(
        ylim = c(0, 1),
        panel.fun  = function(x, y) {
          sector.index  <- get.cell.meta.data("sector.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          #circos.text(mean(xlim),  mean(ylim), sector.index,cex = 0.6, facing = "bending.inside",  niceFacing = TRUE)
          if (grepl("^GO:?", sector.index)) {
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

    
    process_enrichment_data_base <- function(data,pname = "pvalue") {
        
        required_cols <- c("ONTOLOGY", "ID", "Description", "GeneRatio", pname, "BgRatio", "GeneRatio_raw")
        missing_cols <- required_cols[!required_cols %in% colnames(data)]
        if (length(missing_cols) > 0) {
            stop("Data are missing required columns: ",
                 paste(missing_cols, collapse = ", "))
        }
        
        result <- data %>% dplyr::rename("Class" = "ONTOLOGY","TermID" = "ID","Term" = "Description","Ratio" = "GeneRatio")
        
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
    if(p_name == "p.value"){p_name = "pvalue"}else if(p_name == "p.adj"){p_name = "p.adjust"}
    if(species == "9606"){
        species_OrgDb <- org.Hs.eg.db::org.Hs.eg.db
        organism <- "hsa"
    }else if(species == "10090"){
        species_OrgDb <- org.Mm.eg.db::org.Mm.eg.db
        organism <- "mmu"
    }
    gene_ids <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType ="ENTREZID", OrgDb = species_OrgDb)
  
    #--KEGG--------------------------------------------------------------------------------
    # install.packages('R.utils')
    R.utils::setOption("clusterProfiler.download.method", "auto")
    
    if(kegg_analysis_method == "online"){
        kegg <- clusterProfiler::enrichKEGG(gene = gene_ids$ENTREZID,
                                            organism = organism,
                                            keyType = "kegg",
                                            pvalueCutoff = pvalueCutoff,
                                            qvalueCutoff = qvalueCutoff,
                                            pAdjustMethod = 'BH')
    }else if(kegg_analysis_method == "local"){
        kegg <- clusterProfiler::enrichKEGG(gene = gene_ids$ENTREZID,
                                            organism = organism,
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = pvalueCutoff,
                                            qvalueCutoff = qvalueCutoff,
                                            use_internal_data = TRUE)
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
    pdata <- pdata %>% dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")))
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

        plot_list <- list(barplot = NULL, dotplot = NULL, sankeyplot = NULL,
                          circularplot = NULL)
        
        if("barplot" %in% plot_type){
            p_go <- barplot(data = pdata_plot, p_name = p_name, colors = fill_colors)
            p_kegg <- barplot(data = kdata_plot, p_name = p_name, colors = fill_colors)
            p_gokegg <- barplot(data = pkdata_plot, p_name = p_name, colors = fill_colors)

            plot_list$barplot$go = p_go
            plot_list$barplot$kegg = p_kegg
            plot_list$barplot$gokegg = p_gokegg
        }
        if("dotplot" %in% plot_type){
            p_go <- dotplot(data = pdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)
            p_kegg <- dotplot(data = kdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)
            p_gokegg <- dotplot(data = pkdata_plot,x_col = "GeneRatio", size_col = "Count", fill_col = p_name, fill_colors = fill_colors, base_size = 20)

            plot_list$dotplot$go = p_go
            plot_list$dotplot$kegg = p_kegg
            plot_list$dotplot$gokegg = p_gokegg
        }
        if("sankeyplot" %in% plot_type){
            p_go <- sankeyplot(data = pdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            p_kegg <- sankeyplot(data = kdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            p_gokegg <- sankeyplot(data = pkdata_plot, plot_mode = scankeyplot_mode, fill_colors = fill_colors,description_width = description_wrap_width, buble_panel_border = FALSE)
            
            plot_list$sankeyplot$go = p_go
            plot_list$sankeyplot$kegg = p_kegg
            plot_list$sankeyplot$gokegg = p_gokegg           
        }
        # if("treeplot" %in% plot_type){
            
            
        #     plot_list$treeplot$go = p_go
        #     plot_list$treeplot$kegg = p_kegg
        #     plot_list$treeplot$gokegg = p_gokegg           
        # }
        if("circularplot" %in% plot_type){
            pdata_plot_process <- process_enrichment_data_base(pdata_plot,pname = p_name)
            kdata_plot_process <- process_enrichment_data_base(kdata_plot,pname = p_name)
            pkdata_plot_process <- process_enrichment_data_base(pkdata_plot,pname = p_name)
            p_go <- function(){
               circosEnrichmentPlot(pdata_plot_process,topN = go_topn,classCol = fill_colors,
                                     ifLog = TRUE,
                                     type = "base"
                                    )             
            }

            #p_go <- recordPlot()
            p_kegg <- function(){
                circosEnrichmentPlot(kdata_plot_process,topN = kegg_topn,classCol = fill_colors,
                                     ifLog = TRUE,
                                     type = "base"
                                    )              
            }

            #p_kegg <- recordPlot()
            p_gokegg <- function(){
                circosEnrichmentPlot(pkdata_plot_process,topN = go_topn,classCol = fill_colors,
                                     ifLog = TRUE,
                                     type = "base"
                                    )              
            }

            #p_gokegg <- recordPlot()

            plot_list$circularplot$go = p_go
            plot_list$circularplot$kegg = p_kegg
            plot_list$circularplot$gokegg = p_gokegg              
        }
    }else{
        plot_list <- NULL 
        
    }                          

   return(list(result = pkdata, plot = plot_list))                           
                        
}
