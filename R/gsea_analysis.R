#' Run GSEA for multiple ranked feature groups and MSigDB collections
#'
#' Build a ranked gene vector for each level of a grouping column, load human
#' or mouse MSigDB GMT collections, run clusterProfiler GSEA for every
#' group-collection combination in parallel, and create enrichment plots for
#' the strongest significant pathways.
#'
#' @param feature_value_df A data frame containing gene identifiers, ranking
#'   values, and a grouping column.
#' @param feature_col Name of the gene or feature identifier column.
#' @param value_col Name of the numeric ranking-statistic column.
#' @param separate_col Name of the column whose factor levels define separate
#'   ranked lists and analyses.
#' @param background_gene_sets Character vector selecting collections. Supported
#'   patterns include `"GO"`, `"BP"`, `"CC"`, `"MF"`, `"KEGG"`,
#'   `"Reactome"`, `"Hallmark"`, and `"Wikipathway"`.
#' @param colors Color vector for enrichment-score curves.
#' @param species NCBI taxonomy identifier: `9606` for human or `10090` for
#'   mouse.
#' @param p_threshold_type Result column used for the final 0.05 filter,
#'   normally `"pvalue"` or `"p.adjust"`.
#' @param topn Maximum number of significant pathways displayed per plot.
#' @param subplots Integer subset of `1:3` selecting the enrichment-score,
#'   gene-hit, and ranked-metric panels.
#' @param base_size Base plotting text size.
#' @param legend_text_size Enrichment-curve legend text size.
#' @param width_break Approximate character width used to wrap pathway labels.
#' @param gmt_MSigDB_dir Root directory containing the expected human and mouse
#'   MSigDB version 2025.1 GMT directory structure.
#'
#' @return A nested list. The first level is named by each level of
#' `separate_col`, and the second by each selected gene-set collection. A
#' successful element contains `result`, the filtered GSEA table; `gsea`, the
#' enrichment plot; and `topn`. A failed combination is returned as a
#' `gsea_error` object containing an `error` message.
#'
#' @details
#' Within each separate group, ranking values are named by `feature_col` and
#' sorted decreasingly. [clusterProfiler::GSEA()] is run with
#' `pvalueCutoff = 1`; results are subsequently filtered at 0.05 using
#' `p_threshold_type` and ordered by decreasing absolute normalized enrichment
#' score.
#'
#' The function expects MSigDB 2025.1 symbol GMT files below
#' `human/v2025.1` or `mouse/v2025.1`. Selecting `"GO"` includes BP, CC, and
#' MF. The internal all-C2 collection is currently named `"CM2All"`; therefore
#' the default label `"C2All"` does not match that collection under the current
#' pattern selection.
#'
#' Analyses use a multisession future plan with up to 15 workers and seeded
#' furrr execution. The function sets `future.globals.maxSize` to 10 GiB and
#' changes the future plan to sequential at the end; prior option and plan
#' values are not restored.
#'
#' @examples
#' \dontrun{
#' gsea_input <- data.frame(
#'   gene = rep(paste0("Gene", seq_len(1000)), 2),
#'   rank = rnorm(2000),
#'   target = rep(c("TargetA", "TargetB"), each = 1000)
#' )
#' gsea_result <- gsea_analysis(
#'   feature_value_df = gsea_input,
#'   feature_col = "gene",
#'   value_col = "rank",
#'   separate_col = "target",
#'   background_gene_sets = c("GO", "KEGG", "Hallmark"),
#'   species = 9606,
#'   p_threshold_type = "p.adjust",
#'   topn = 5,
#'   gmt_MSigDB_dir = "/reference/functional/gene_sets/MSigDB"
#' )
#' }
#'
#' @seealso [clusterProfiler::GSEA()], [clusterProfiler::read.gmt()],
#'   [enrichplot::gseaplot2()], [furrr::future_imap()]
#' @export
gsea_analysis <- function(feature_value_df,
                          feature_col,
                          value_col,
                          separate_col,
                          background_gene_sets = c("GO","KEGG","Reactome","Hallmark","C2All","Wikipathway"),
                          colors = basicR::get_colors(number = 10.1),
                          species = 9606,
                          p_threshold_type = "pvalue", # pvalue/p.adjust
                          topn = 5,
                          subplots = 1:3,
                          base_size = 16,
                          legend_text_size = 16,
                          width_break = 60,
                          gmt_MSigDB_dir = "/reference/functional/gene_sets/MSigDB"
                         ){
    # feature_value_df： 包含特征的相关值得数据
    # feature_col：特征列名
    # value_col：排序的值列名
    # separate_col：目标特征列名
    # background_gene_sets：背景 gmt。 "GO"、"BP"、"CC"、"MF"、"KEGG"、"Reactome"、"Hallmark"、“C2All”

    library(magrittr)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(DOSE)
    library(RColorBrewer)
    library(ggplot2)

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    ## --Function-- ------------------------------------------------------------------------------------
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    gseaplot3 <- function(x, geneSetID, title = "", color = "green", base_size = 16,legend_text_size = 16,
                          size=1.8,rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE,
                          ES_geom = "line",width_break = 60){
        library(grid)
        library(DOSE)

        # 放在 gseaplot3() 里（p.res 创建之后、加 scale 之前）
        wrap50 <- function(x, width = width_break){
            stringr::str_wrap(gsub("_", " ", x), width = width)  # 下划线当空格再换行
        }
        
        gseaScores <- getFromNamespace("gseaScores", "DOSE")
        ES_geom <- match.arg(ES_geom, c("line", "dot"))
        geneList <- position <- NULL
        if (length(geneSetID) == 1) {
            gsdata <- gsInfo(x, geneSetID)
        }else{
            gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
        }
        desc_levels <- unique(gsdata$Description) # 为了设置标签
        gsdata$Description <- factor(gsdata$Description, levels = desc_levels)
        gsdata <- gsdata[order(gsdata$Description, gsdata$x), ]   # 关键：排序对齐
        
        cols <- if (length(color) == 1) {
            setNames(rep(color, length(desc_levels)), desc_levels)
        } else {
            setNames(color[seq_along(desc_levels)], desc_levels)
        }
        
        p <- ggplot(gsdata, aes(x = x)) + xlab(NULL) + theme_classic(base_size) + 
            theme(panel.grid.major = element_line(colour = "grey92"), 
                  panel.grid.minor = element_line(colour = "grey92"), 
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank()) + 
            scale_x_continuous(expand = c(0, 0))
        if (ES_geom == "line") {
            es_layer <- geom_line(aes(y = runningScore, color = Description),size = size)
        }else {
            es_layer <- geom_point(aes(y = runningScore, color = Description), 
                                   size = size, data = subset(gsdata, position == 1))
        }
        p.res <- p + es_layer + 
            labs(y = "Enrichment Score")+
            theme(legend.position = c(0.8, 0.8),
                  legend.title = element_blank(),
                  legend.background = element_rect(fill = "transparent")
                 )
        p.res <- p.res + 
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_text(face = "bold"),
                      axis.title.y = element_text(size = base_size +1,face = "bold"),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank(),
                      plot.background = element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2,unit = "cm")
                     )
        
        i <- 0
        for (term in unique(gsdata$Description)) {
            idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
            gsdata[idx, "ymin"] <- i
            gsdata[idx, "ymax"] <- i + 1
            i <- i + 1
        }
        p2 <- ggplot(gsdata, aes(x = x)) + 
                    geom_linerange(aes(ymin = ymin,ymax = ymax, color = Description)) + 
                    scale_x_continuous(expand = c(0, 0)) + 
                    scale_y_continuous(expand = c(0, 0)) +
                    labs(x = NULL)+
                    theme_classic(base_size = base_size) + 
                    theme(legend.position = "none",
                          plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
                          axis.ticks = element_blank(),
                          axis.text = element_blank(), 
                          axis.line.x = element_blank()
                         )

        if (length(geneSetID) == 1) {
            v <- seq(1, sum(gsdata$position), length.out = 9)
            inv <- findInterval(rev(cumsum(gsdata$position)), v)
            if (min(inv) == 0) 
                inv <- inv + 1
            col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
            ymin <- min(p2$data$ymin)
            yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
            xmin <- which(!duplicated(inv))
            xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
            d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                xmax = xmax, col = col[unique(inv)])
            p2 <- p2 + geom_rect(aes(xmin = xmin, xmax = xmax, 
                ymin = ymin, ymax = ymax, fill = I(col)), data = d, 
                alpha = 0.9, inherit.aes = FALSE)
        }
        df2 <- p$data
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + 
                geom_segment(data = df2, aes(x = x, xend = x,y = y, yend = 0), color = "grey")
        p.pos <- p.pos + 
                labs(y = "Rank Metric",x = "Gene Rank")+
                theme(plot.background = element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      axis.text = element_text(face = "bold"),
                      axis.title = element_text(size = base_size + 1,face = "bold"),
                      plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,  l = 0.2, unit = "cm")
                     )
        if (!is.null(title) && !is.na(title) && title != ""){
            p.res <- p.res + ggtitle(title)+
                theme(plot.title = element_text(size = base_size+2,face = "bold"))
        }
        
        if (length(color) >= length(geneSetID)) {
            p.res <- p.res + scale_color_manual(values = cols,breaks = desc_levels,labels = wrap50)
            if (length(color) == 1) {
                p.res <- p.res + theme(legend.position = "top",
                                       legend.key.size = unit(0.3, "cm"),
                                       legend.text = element_text(size = legend_text_size),  # 设置图例字体大小
                                       legend.title = element_text(size = 15)  # 设置图例标题字体大小
                                      )
                p2 <- p2 + scale_color_manual(values = "black")
            }else {
                p2 <- p2 + scale_color_manual(values = cols, breaks = desc_levels) + 
                theme(legend.position = "none",
                      plot.margin=margin(t = 0.1, r = 0.2, b = 0.1, l = 0.2, unit = "cm")
                     )
            }
        }
        if (pvalue_table) {
            pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
            rownames(pd) <- pd$Description
            pd <- pd[, -1]
            for (i in seq_len(ncol(pd))) {
                pd[, i] <- format(pd[, i], digits = 3)
            }
            tp <- tableGrob2(pd, p.res)
            p.res <- p.res + theme(legend.position = "top",
                                   legend.key.size = unit(0.3, "cm"),
                                   plot.margin=margin(t = 0.8, r = 0.2, b = 0.2, l = 0.2,unit = "cm")
                                  )+
            guides(color = guide_legend(ncol = 1))
        }

        if(setequal(1:2, subplots)){
            p2 <- p2 + labs(x = "Gene Rank")
        }else if(setequal(1:3, subplots)){
            p2 <- p2 + theme(axis.title = element_blank())
        }
        if(setequal(1, subplots)){
            p.res <- p.res + labs(x = "Gene Rank")
        }
        
        plotlist <- list(p.res, p2, p.pos)[subplots]
        
        n <- length(plotlist)
        plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
            axis.ticks.x = element_line(), axis.text.x = element_text())
        if (length(subplots) == 1){ 
            return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2,r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
        }
        if (length(rel_heights) > length(subplots)){
            rel_heights <- rel_heights[subplots]
        }

        return(aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights))
    }
    
    gsInfo <- function(object, geneSetID) {
        geneList <- object@geneList
    
        if (is.numeric(geneSetID))
            geneSetID <- object@result[geneSetID, "ID"]
    
        geneSet <- object@geneSets[[geneSetID]]
        exponent <- object@params[["exponent"]]
        df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
        df$ymin <- 0
        df$ymax <- 0
        pos <- df$position == 1
        h <- diff(range(df$runningScore))/20
        df$ymin[pos] <- -h
        df$ymax[pos] <- h
        df$geneList <- geneList
    
        df$Description <- object@result[geneSetID, "Description"]
        return(df)
    }
    
    tableGrob2 <- function(d, p = NULL) {
        # has_package("gridExtra")
        d <- d[order(rownames(d)),]
        tp <- gridExtra::tableGrob(d)
        if (is.null(p)) {
            return(tp)
        }
    
        # Fix bug: The 'group' order of lines and dots/path is different
        p_data <- ggplot_build(p)$data[[1]]
        # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
        p_data <- p_data[order(p_data[["group"]]), ]
        pcol <- unique(p_data[["colour"]])
        ## This is fine too
        ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
        j <- which(tp$layout$name == "rowhead-fg")
    
        for (i in seq_along(pcol)) {
            tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
        }
        return(tp)
    }
    

    # -- Main ------------------------------------------------------------------------------------------------------------------------------------------------
    # 提取数据并排序
    feature_value_df[[separate_col]] <- as.factor(feature_value_df[[separate_col]])

    feature_df <- feature_value_df %>% arrange(!!sym(separate_col),desc(!!sym(value_col)))
    featurelist_list <- lapply(levels(feature_value_df[[separate_col]]),function(feature){
        feature_df_sub <- feature_df %>% dplyr::filter(!!sym(separate_col) == feature)
        featurelist <- setNames(feature_df_sub[[value_col]], feature_df_sub[[feature_col]])
        return(featurelist)
    })
    featurelist_list <- setNames(featurelist_list,levels(feature_value_df[[separate_col]]))

    # gmt_list
    if(species == 9606){
        gmt_list <- lapply(c("c5.go.bp","c5.go.cc","c5.go.mf","c2.cp.kegg_medicus","c2.cp.reactome","h.all","c2.all","c2.cp.wikipathways"),function(set){
            gene_set <- clusterProfiler::read.gmt(file.path(gmt_MSigDB_dir,"human/v2025.1",paste0(set,".v2025.1.Hs.symbols.gmt")))
            return(gene_set)
        })
    }else if(species == 10090){
        gmt_list <- lapply(c("m5.go.bp","m5.go.cc","m5.go.mf","m2.cp","m2.cp.reactome","mh.all","m2.all","m2.cp.wikipathways"),function(set){
            gene_set <- clusterProfiler::read.gmt(file.path(gmt_MSigDB_dir,"mouse/v2025.1",paste0(set,".v2025.1.Mm.symbols.gmt")))
            return(gene_set)
        })      
    }
    gmt_list <- setNames(gmt_list,c("GO_BP","GO_CC","GO_MF","KEGG","Reactome","Hallmark","CM2All","Wikipathway"))

    set_select_name <- grep(paste(tolower(background_gene_sets),collapse = "|"),
                            c("GO_BP","GO_CC","GO_MF","KEGG","Reactome","Hallmark","CM2All","Wikipathway"),
                            value = TRUE,
                            ignore.case = TRUE
                           )
    # 使用 lapply 从 gmt_list 中提取这些元素
    selected_gmt_sets <- lapply(set_select_name, function(set){ 
        gmt_select <- gmt_list[[set]]
        return(gmt_select)
    })
    selected_gmt_sets <- setNames(selected_gmt_sets,set_select_name)

    # 设置循环集
    combos <- expand.grid(feature = names(featurelist_list),
                          set = set_select_name,
                          stringsAsFactors = FALSE
                         )

    run_one <- function(featurelist_list,selected_gmt_sets,feature,set,p_threshold_type = "pvalue",topn = 5,title = "",subplots = 1:3){
        
        featurelist <- featurelist_list[[feature]] %>% sort(., decreasing = T)

        gene_backgrand_set <- selected_gmt_sets[[set]]
        gene_backgrand_list <- selected_gmt_sets[[set]] %>% split(.$term) %>% lapply( "[[", 2)
        result_gsea <- clusterProfiler::GSEA(geneList = featurelist, TERM2GENE = gene_backgrand_set, pvalueCutoff = 1)
        result <- result_gsea@result
        result <- result[result[[p_threshold_type]] < 0.05,]
        result <- result[order(abs(result$NES), decreasing = T),]
        result.write <- result[,-c(8,9,10)]
        
        gseaplot3(x = result_gsea,
                 geneSetID = result.write$ID[1:min(length(result_gsea$ID),topn)], 
                 title = title,
                 color = colors,
                 size = 1,
                 pvalue_table = T,
                 base_size = base_size,
                 legend_text_size = legend_text_size,
                 rel_heights = c(1.5, 0.8, 0.5), 
                 subplots = subplots, 
                 ES_geom = "line",
                 width_break = width_break
                ) -> gsea
        
        return(list(result = result.write,gsea = gsea,topn = topn))
    }
 
    future::plan(future::multisession, workers = min(15, length(unique(combos$feature))))
    options(future.globals.maxSize = 10 * 1024^3) # 2GB
   results_list <- combos %>%
  split(.$feature) %>%
  furrr::future_imap(
    function(df, feat){

      sets <- as.character(df$set)

      out <- purrr::map(
        sets,
        ~ tryCatch(
          run_one(
            featurelist_list = featurelist_list,
            selected_gmt_sets = selected_gmt_sets,
            feature = feat,
            set = .x,
            p_threshold_type = p_threshold_type,
            topn = topn,
            title = paste(feat, .x, sep = " "),
            subplots = subplots
          ),
          error = function(e)
            structure(
              list(error = conditionMessage(e)),
              class = "gsea_error"
            )
        )
      )

      names(out) <- sets
      out
    },
    .options = furrr::furrr_options(seed = 123)
  ) 
    future::plan(future::sequential)
    
    return(results_list)
}
