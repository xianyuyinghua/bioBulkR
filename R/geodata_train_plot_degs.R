#' Run differential expression, Venn, enrichment, and PPI analyses for GEO data
#'
#' Continue a processed GEO project through array or count-based differential
#' expression, volcano and heatmap plotting, target-gene intersections, GO and
#' KEGG enrichment, and STRING protein-interaction network visualization.
#'
#' @param workdir Parent analysis directory. The function changes the working
#'   directory to this path.
#' @param Project_ID Project directory identifier used to find or create a
#'   numbered project folder.
#' @param Data_ID GEO Series identifier naming the data subdirectory.
#' @param exp_dir Expression-data subdirectory within the GEO data directory.
#' @param group_levels Optional ordering of sample groups passed to
#'   [group_order_set()].
#' @param logFC Absolute log2 fold-change threshold for differential genes.
#' @param p_value Differential-expression significance measure: `"p.value"`
#'   or `"p.adj"`.
#' @param p Differential-expression significance threshold.
#' @param enrich_figure_type Enrichment plot family requested from
#'   [enrich_go_kegg()], such as `"barplot"`, `"dotplot"`, `"sankeyplot"`,
#'   `"treeplot"`, or `"circularplot"`.
#' @param enrich_pname Enrichment significance column, `"p.value"` or
#'   `"p.adj"`.
#' @param go_kegg_topn Number of top GO and KEGG terms displayed in enrichment
#'   plots.
#' @param description_wrap_width Approximate character width for wrapping
#'   enrichment descriptions.
#' @param sample_color Colors used for sample annotations in heatmaps.
#' @param confidence STRING required-score threshold used for the PPI query.
#' @param target_gene_file Response to the prompt asking whether a target-gene
#'   CSV has been placed in `03_Venn`; `NULL` keeps the interactive prompt.
#' @param target_file Path selected when multiple target-gene CSV files are
#'   found; `NULL` keeps the interactive prompt.
#' @param target_gene_column Name of the target file's gene-symbol column when
#'   it cannot be identified automatically; `NULL` prompts interactively.
#' @param target_name Name shown for the target set in Venn plots; `NULL`
#'   prompts when needed.
#' @param target_set Response indicating whether `Target_venn_result.csv` is
#'   available; `NULL` keeps the interactive prompt.
#'
#' @return The differential-expression result table returned in
#'   `Res$result` after volcano and heatmap processing. The workflow also saves
#'   its full intermediate `data` list to `00_data_result/02_data.rds`.
#'
#' @details
#' The function expects output from the GEO download and processing workflow,
#' including `01_data.rds`, an expression CSV, and `02_Condition.csv` below the
#' selected project and GEO directories. It creates `02_DEGs`, `03_Venn`,
#' `04_Enrich`, and `05_PPI` directories as necessary.
#'
#' Array data are analyzed with [Array_DEGs()], whereas high-throughput count
#' data use [High_Throughtput_Count_DEG()]. Both branches create volcano,
#' rectangular heatmap, and circular heatmap outputs. The target-gene file in
#' `03_Venn` is intersected with all, upregulated, and downregulated DEGs.
#' Enrichment is then run separately for those three intersections, followed
#' by an online STRING query and network plot.
#'
#' This is a file-writing and network-enabled workflow. It writes CSV, RDS,
#' PDF, and PNG outputs, can overwrite existing results, changes the
#' process working directory, and can request console input. If genes in an
#' existing `Target_venn_result.csv` disagree with the newly calculated
#' intersection, the current implementation calls `quit(save = "no")`.
#'
#' @examples
#' \dontrun{
#' deg_result <- geodata_train_plot_degs(
#'   workdir = "/project/project_analysis/26_26YLM199F/02_Result",
#'   Project_ID = "03_Expression_DEG",
#'   Data_ID = "GSE28829",
#'   group_levels = c("Control", "Disease"),
#'   logFC = 0.5,
#'   p_value = "p.adj",
#'   p = 0.05,
#'   enrich_figure_type = "treeplot",
#'   enrich_pname = "p.value",
#'   go_kegg_topn = c(10, 10),
#'   confidence = 400,
#'   target_gene_column = "Symbol",
#'   target_name = "Target",
#'   target_set = "no"
#' )
#' }
#'
#' @seealso [geodata_download_and_process()], [Array_DEGs()],
#'   [High_Throughtput_Count_DEG()], [Venn_plot()], [enrich_go_kegg()],
#'   [string_ppi_get_data()], [network_graph_plot()]
#' @export
geodata_train_plot_degs <- function(
    workdir = "/project/project_analysis/26_26YLM199F/02_Result",
    Project_ID = "03_Expression_DEG",
    Data_ID = "GSE28829",
    exp_dir = "01_Expression_data",
    group_levels = NULL,
    logFC = 0.5,
    p_value = "p.adj",
    p = 0.05,
    enrich_figure_type = "barplot",
    enrich_pname = "p.value",
    go_kegg_topn = c(10,10),
    description_wrap_width = 100,
    sample_color = basicR::get_colors(number = 10.1),
    confidence = 150,
    target_gene_file = NULL,
    target_file = NULL,
    target_gene_column = NULL,
    target_name = NULL,
    target_set = NULL
) {
    setwd(workdir)

    ##--设置参数-----------------
    opt <- list()
    opt$Project_ID <- Project_ID
    opt$Data_ID <- Data_ID
    opt$exp_dir <- exp_dir

    suppressPackageStartupMessages({
        library(dplyr)
        library(basicR)
        library(ggpubr)
        library(limma)
        library(optparse)
        library(stringr)
        library(ggplot2)
        library(ggrepel)
        library(viridis)
        library(cowplot)
        library(DESeq2)
        library(circlize)
        library(data.table)
        library(dendextend)
        library(gridtext)
        library(VennDiagram)
        library(ComplexHeatmap)
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(org.Mm.eg.db)
        library(GOplot)
        library(treemap)
        library(tidyr)
        library(ggraph)
        library(ccgraph)
        library(tidygraph)
        library(ggsci)
    })

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 设置帮助信息
    description <- paste0(
      "\n",
      blue, "This script checks and processes Bulk RNA data from GEO.", reset, "\n\n",
      cyan, "Arguments:", reset, "\n\n\t",
      
      yellow, str_pad("-p  --Project_ID",width = 30,side = "right"), reset, str_pad("./Rdata/2_seurat_Norm_Cluster.rds", width = 50, side = "right"), "The analysis project ID\n\n\t",
      yellow, str_pad("-d  --Data_ID",width = 30,side = "right"), reset, str_pad("SingleR", width = 50, side = "right"), "The analysis data ID\n\n"
    )
    
    # 定义命令行
    option_list <-list(
      make_option(
        c("-p","--Project_ID"),
        action = "store",
        default = NA,
        type = "character",
        help = "The analysis project ID"
      ),
      make_option(
        c("-d","--Data_ID"),
        action = "store",
        default = NA,
        type = "character",
        help = "The analysis data ID"
      )
    )
    
    # 检查参数
    if(is.na(opt$Project_ID)){stop(paste0(red,"\n\tPlease set the 'Project ID'.",reset,"\n\n"),call. = FALSE)}
    if(is.na(opt$Data_ID)){stop(paste0(red,"\tPlease set the 'Data ID'.",reset,"\n\n"),call. = FALSE)}
    
    if (length(commandArgs(trailingOnly = TRUE)) == 0) {cat(description)}
    
    cat("\n\n")
    print(str(opt))
    cat("\n\n")
    
    # 检查项目文件夹并创建
    folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE) %>% basename() %>% .[!grepl("^\\.", basename(.))]
    if(grepl(opt$Project_ID,folders) %>% any()){
      cat(paste0("\n",magenta,"Warring  Project ID: ",reset,blue,opt$Project_ID,reset,magenta," already exists in the current path.",reset,"\n\n"))
      project_folder <- grep(opt$Project_ID,folders,value = T)
    }else{
      dir_num <- folders %>% strsplit(.,"_") %>% sapply(.,"[",1) %>% as.numeric() %>% max() %>%  `+`(1) %>% as.character()
      if(str_length(dir_num) != 2){dir_num <- paste0("0",dir_num)}
      project_folder <- paste0(dir_num,"_",opt$Project_ID)
    }
    
    # 创建文件夹
    data_dir_name <- file.path(project_folder,opt$Data_ID)
    save_dirs <- c(opt$exp_dir,"02_DEGs","03_Venn","04_Enrich","05_PPI")
    lapply(save_dirs,function(dirname){
        if (!dir.exists(file.path(data_dir_name,dirname))){dir.create(file.path(data_dir_name,dirname),recursive = T)}
    })
    
    # 处理样本信息
    sample_group <- read.csv(file.path(data_dir_name,save_dirs[1],"02_Condition.csv"),row.names = 1)
    if( "group" %in% colnames(sample_group) ){
        sample_group$group <- stringr::str_to_title(sample_group$group)
        sample_group <- sample_group %>% dplyr::rename( "Group" = "group") %>% tibble::column_to_rownames("sample")
        write.csv(sample_group,file = file.path(data_dir_name,save_dirs[1],"02_Condition.csv"))
    }else{
        head(sample_group)
    }
    
    #--Main-------------------------------------------------------------------------------------------------------------------
    
    # 获取数据信息（物种、是芯片数据还是高通量数据），这些信息由前期数据处理脚本完成
    raw_data_path <- list.files(file.path(data_dir_name,"00_data_result"),pattern = "^01.*.rds$",recursive = T,full.names = T)
    if(length(raw_data_path) != 0){
        raw_data <- readRDS(file = raw_data_path) 
    }else{
        raw_data <- NULL
    }
    
    
    expr_files <- list.files(file.path(data_dir_name,save_dirs[1]),pattern = "^01_Expression_",recursive = T,full.names = T)
    expr_file_pattern  <- grep("count|array",expr_files,value= T,ignore.case = TRUE) %>% basename()
    if(length(expr_file_pattern) > 1){
        expr_file_pattern  <- grep("array",expr_files,value= T,ignore.case = TRUE) %>% basename()
    }
    data <- read_exp_data(expression_dir = file.path(data_dir_name,save_dirs[1]),
                          expression_pattern = expr_file_pattern, 
                          condition_pattern ="*Condition*.csv$"
                         )
    data$group_order <- group_order_set(vector = data$condition$Group,group_levels = group_levels)
    
    # 判断读入的data 数据是否和新读入的表达数据的样本数一致
    if("condition" %in% names(raw_data)){
        if(nrow(raw_data$condition) == nrow(data$condition)){
            data <- c(raw_data,data)
        }else{
            data$species = raw_data$species
            data$throughput_type = raw_data$throughput_type
        }
    }else{
        if(length(raw_data) != 0){
            data$species = raw_data$species
        }else{
            data$species = 9606
        }
    }
    
    if(length(raw_data) != 0){
        if(raw_data$data_type == 'Expression profiling by array' ){
            data$throughput_type  <- NULL
        }else{
            data$throughput_type  <- "Count"
        }
        
    }else{
        data$throughput_type  <- "Count"
    }
    
    if(data$exp_data_type == "array"){
        data$throughput_type  <- NULL
    }
    data$logFC <- logFC
    data$p_value <- p_value
    
    # array数据和throughput数据
    if(data$exp_data_type == "array"){
        # limma 差异分析
        cat("\n\n\n",yellow,"Start Array DEGs Plot!",reset,"\n")
        Res <- Array_DEGs(expr = data$expr,   # 表达矩阵
                          condition = data$condition,  # 样本分组信息
                          group_order = data$group_order,  # 分组顺序
                          output_dir = file.path(data_dir_name,save_dirs[2]),  # 输出文件位置
                          logFC = logFC,   # logFC 阈值
                          p_value = p_value,  # pvalue / padj, 值的话内部固定设置 0.05
                          top_gene = 10,    # 差异出图展示的标签的TOP基因数
                          p = p,
                          change = c("Up","Down","Not")  # 对应上下不调的名字
                         )
        
        # 差异出图（Volcano、Rectangular Heatmap、Circular Heatmap）
        cat("\n\n\n",yellow,"Start Volcano、Rectangular Heatmap、Circular Heatmap Plot!",reset,"\n")
        Res <- Volcano_Heatmap_plot(Res = Res,
                                     expr = data$expr,
                                     condition = data$condition,
                                     sample_color = sample_color,
                                     logFC_name = "logFC",
                                     logFC = logFC,
                                     p_value = p_value,
                                     p = p,
                                     top_gene = 10,
                                     heatmap_scale_range = c(-3,3),
                                     group_order = data$group_order,
                                     group_order_color = c("#CD534CFF", "#1177b0"),
                                     change = c("Up","Down","Not"), 
                                     change_colors = c("orange","#1177b0", "#7b7c7d"),
                                     heatmap_value_color = viridis(10),
                                     prominent_genes = NULL,
                                     prominent_genes_color = c("blue"),
                                     prominent_genes_size = 3,
                                     track_height = 0.3,  # 环形热图的热图本体那一圈有多“厚”
                                     dend_track_height = 0.1  # 聚类树那一圈的厚度
                                    )
        # Save Plot
        save_figure(obj = Res$plot$valcano,
                    filename =  file.path(file.path(data_dir_name,save_dirs[2]),"04_Top10_Volcano"),
                    width = 6.3,
                    height = 6,
                    res = 600,
                    formats = c("pdf","png") # "pdf","tiff","png","svg"
                   )
        save_figure(obj = Res$plot$rHeatmap,
                    filename =  file.path(file.path(data_dir_name,save_dirs[2]),"05_Top10_Heatmap"),
                    width = 6,
                    height = 7,
                    res = 600,
                    formats = c("pdf","png") # "pdf","tiff","png","svg"
                   )
        save_figure(obj = Res$plot$cHeatmap,
                    filename =  file.path(file.path(data_dir_name,save_dirs[2]),"05_Top10_Heatmap_circos"),
                    width = 7.5,
                    height = 7.5,
                    res = 600,
                    formats = c("pdf","png") # "pdf","tiff","png","svg"
                   )
        
        # save Data
        data$Res <- Res$result
        data$DEGs <- Res$result[Res$result$sig != "Not",]  %>% rownames()  # 获取差异基因
        data$DEGs_Up <- Res$result[Res$result$sig == "Up",]  %>% rownames()
        data$DEGs_Down <- Res$result[Res$result$sig == "Down",]  %>% rownames()
                
    }else if(data$exp_data_type == "throughput"){
        # 判断数据为高通量数据
        # DEseq2 差异分析
        if(data$throughput_type == "Count"){
            # DEseq2
             Res <- High_Throughtput_Count_DEG(expr = data$expr,
                                               condition = data$condition,
                                               group_order = data$group_order,
                                               Exp_output_dir = file.path(data_dir_name,opt$exp_dir),  # 输出标准化后的表达矩阵
                                               DEG_output_dir = file.path(data_dir_name,"02_DEGs"),   # 输出差异结果
                                               logFC = logFC,   # logFC 阈值
                                               p = p,
                                               p_value = p_value,  # pvalue / padj, 值的话内部固定设置 0.05
                                               expr_cutoff = 0,  # 低表达基因的过滤阈值
                                               filter_sample_num = nrow(data$condition) * 0.5,  # 表达量超过阈值的样本数要求
                                               change = c("Up","Down","Not")
                                              )
        
             data$Norm_expr <- read.csv(file = file.path(file.path(data_dir_name,opt$exp_dir),"01_Expression_throughput_vst_Norm.csv"),row.names = 1,check.names = F)
             data$log_expr  <- log2(data$Norm_expr + 1)
            
             # 差异出图（Volcano、Rectangular Heatmap、Circular Heatmap）
             cat("\n\n\n",yellow,"Start Volcano、Rectangular Heatmap、Circular Heatmap Plot!",reset,"\n")
             Res <- Volcano_Heatmap_plot(Res = Res,
                                         expr = data$log_expr,
                                         condition = data$condition,
                                         sample_color = sample_color,
                                         logFC_name = "logFC",
                                         logFC = logFC,
                                         p_value = p_value,
                                         p = p,
                                         top_gene = 10,
                                         heatmap_scale_range = c(-3,3),
                                         group_order = data$group_order,
                                         group_order_color = c("#FF7F00", "#1177b0"),
                                         change = c("Up","Down","Not"), 
                                         change_colors = c("orange","#1177b0", "#7b7c7d"),
                                         heatmap_value_color = viridis(10),
                                         prominent_genes = NULL,
                                         prominent_genes_color = c("blue"),
                                         prominent_genes_size = 3,
                                         track_height = 0.3,  # 环形热图的热图本体那一圈有多“厚”
                                         dend_track_height = 0.1  # 聚类树那一圈的厚度
                                        )
    
            # Save Plot
            save_figure(obj = Res$plot$valcano,
                        filename =  file.path(file.path(data_dir_name,save_dirs[2]),"04_Top10_Volcano"),
                        width = 6.3,
                        height = 6,
                        res = 600,
                        formats = c("pdf","png") # "pdf","tiff","png","svg"
                       )
            save_figure(obj = Res$plot$rHeatmap,
                        filename =  file.path(file.path(data_dir_name,save_dirs[2]),"05_Top10_Heatmap"),
                        width = 6,
                        height = 7,
                        res = 600,
                        formats = c("pdf","png") # "pdf","tiff","png","svg"
                       )
            save_figure(obj = Res$plot$cHeatmap,
                        filename =  file.path(file.path(data_dir_name,save_dirs[2]),"05_Top10_Heatmap_circos"),
                        width = 7.5,
                        height = 7.5,
                        res = 600,
                        formats = c("pdf","png") # "pdf","tiff","png","svg"
                       )        
            data$Res <- Res$result
            data$DEGs <- Res$result[Res$result$sig != "Not",]  %>% rownames()
            data$DEGs_Up <- Res$result[Res$result$sig == "Up",]  %>% rownames()
            data$DEGs_Down <- Res$result[Res$result$sig == "Down",]  %>% rownames()
            
            
        }else if(data$throughput_type == "FPKM"){
            # 转为 TPM后，limma
            
        }else if(data$throughput_type == "TPM"){
            # limma
            
        }
    
            
    }
    
    
    ############################################# Venn Plot #############################################################
    
    cat("\n\n\n",yellow,"Start Venn Plot!",reset,"\n")
    target_files <- list.files(path = file.path(data_dir_name,save_dirs[3]),pattern = ".csv$",recursive = T,full.names = T)
    target_files <- target_files[!grepl("intersect_genes.csv|Target_venn_result.csv",target_files)]
    
    # 获取 Venn 文件夹下有没有目标基因csv文件
    if(length(target_files) == 0){
        # 没有获取到任何csv文件的情况，提示手动放入文件值文件夹
        cat("\n",yellow,"No target gene CSV file was detected. Please verify and place the CSV file in the '03_Venn' folder.",reset,"\n")
        flush.console()  # 强制刷新控制台输出
        target_gene_file <- if (is.null(target_gene_file)) readline(prompt = "Is there already a target gene CSV file in the '03_Venn' folder? (yes/no): ") else target_gene_file
        if(target_gene_file == "yes"){
            target_files <- list.files(path = file.path(data_dir_name,save_dirs[3]),pattern = ".csv$",recursive = T,full.names = T)
        }
    }else if(length(target_files) > 1){
        # 当获取到多个csv文件，提示手动选择
        cat("\n",yellow,"Multiple CSV files matched. Please select the target file.",reset,"\n")
        print(target_files)
        flush.console()  # 强制刷新控制台输出
        target_files <- if (is.null(target_file)) readline(prompt = "Please enter the selected CSV file: ") else target_file
    }
    
    target <- read_data(file_path = target_files)
    target_symbol_name <- names(sapply(target,identify_gene_ID_type) == "SYMBOL")
    if(length(target_symbol_name) != 1){
        cat("\n",yellow,"The Head target gene file:",reset,"\n")
        print(head(target))
        flush.console()  # 强制刷新控制台输出
        select_genecol <- if (is.null(target_gene_column)) readline(prompt = "Please select the target gene data column containing the SYMBOL: ") else target_gene_column
    }else{
        select_genecol <- target_symbol_name
    }
    
    data$target_gene <- target[[select_genecol]] %>% unique()
    if( "gene" %in% tolower(select_genecol) ){
        cat("\n",yellow,"The Head target gene file:",reset,"\n")
        print(head(target))
        flush.console()  # 强制刷新控制台输出
        target_name <- if (is.null(target_name)) readline(prompt = "Please input target name(Show on the Venn plot): ") else target_name
    }else{
        target_name <- select_genecol
    }
    data$target_gene_name <- target_name
    
    # ALL DEGS
    Intersection_gene_list <- list(data$DEGs, data$target_gene)
    Intersection_gene_list <- setNames(Intersection_gene_list,c("DEGs",data$target_gene_name))
    
    com_genes <- Venn_plot(target_genes = Intersection_gene_list,
                           target_genes_names = names(Intersection_gene_list),
                           fill_colors = basicR::get_colors(number = 2.1,package = "ggsci",name = "jco"),
                           base_size = 1.7,
                           label_size = 1.5,
                           sigdigs = 3,  # 设置有效位数
                           digits = 2,    # 设置小数位数
                           print_mode = c("percent","raw") # "percent", "raw"
                          )
    # Save Plot
    save_figure(obj = com_genes$plot,
            filename =  file.path(file.path(data_dir_name,save_dirs[3]),paste0("01_Venn_",paste(names(Intersection_gene_list),collapse = "_"))),
            width = 6,
            height = 6,
            res = 600,
            formats = c("pdf","png") # "pdf","tiff","png","svg"
           )
    
    # Up DEGs
    Intersection_gene_list <- list(data$DEGs_Up, data$target_gene)
    Intersection_gene_list <- setNames(Intersection_gene_list,c("DEGs_Up",data$target_gene_name))
    
    com_genes_up <- Venn_plot(target_genes = Intersection_gene_list,
                           target_genes_names = names(Intersection_gene_list),
                           fill_colors = basicR::get_colors(number = 2.1,package = "ggsci",name = "jco"),
                           base_size = 1.7,
                           label_size = 1.5,
                           sigdigs = 3,  # 设置有效位数
                           digits = 2,    # 设置小数位数
                           print_mode = c("percent","raw") # "percent", "raw"
                          )
    # Save Plot
    save_figure(obj = com_genes_up$plot,
            filename =  file.path(file.path(data_dir_name,save_dirs[3]),paste0("02_Venn_",paste(names(Intersection_gene_list),collapse = "_"))),
            width = 6,
            height = 6,
            res = 600,
            formats = c("pdf","png") # "pdf","tiff","png","svg"
           )
    
    
    
    # Down DEGs
    Intersection_gene_list <- list(data$DEGs_Down, data$target_gene)
    Intersection_gene_list <- setNames(Intersection_gene_list,c("DEGs_Down",data$target_gene_name))
    
    com_genes_down <- Venn_plot(target_genes = Intersection_gene_list,
                           target_genes_names = names(Intersection_gene_list),
                           fill_colors = basicR::get_colors(number = 2.1,package = "ggsci",name = "jco"),
                           base_size = 1.7,
                           label_size = 1.5,
                           sigdigs = 3,  # 设置有效位数
                           digits = 2,    # 设置小数位数
                           print_mode = c("percent","raw") # "percent", "raw"
                          )
    # Save Plot
    save_figure(obj = com_genes_down$plot,
            filename =  file.path(file.path(data_dir_name,save_dirs[3]),paste0("03_Venn_",paste(names(Intersection_gene_list),collapse = "_"))),
            width = 6,
            height = 6,
            res = 600,
            formats = c("pdf","png") # "pdf","tiff","png","svg"
           )
    
    # 排序
    if(file.exists(file.path(data_dir_name,"03_Venn","Target_venn_result.csv"))){
        df_Venn_target <- read.csv(file.path(data_dir_name,"03_Venn","Target_venn_result.csv"))
    }else{
        cat("\n","Target_venn_result.csv Not in '03_Venn' dir, please set the file.")
        flush.console()
        target_set <- if (is.null(target_set)) readline(prompt = "Target_venn_result.csv in '03_Venn'?(yes/no): ") else target_set
        if(tolower(target_set) == "yes"){
            df_Venn_target <- read.csv(file.path(data_dir_name,"03_Venn","Target_venn_result.csv"))
        }else{
            df_Venn_target <- NULL
        }
    }
    if(is.null(df_Venn_target)){
        # 无比较文件，直接使用共有基因
        com_genes = com_genes$common
        com_up_genes = com_genes_up$common
        com_down_genes = com_genes_down$common
    }else{
        # 和已有文件对比
        if(p_value == "p.value"){pnameset = "PValue"}else if(p_value == "p.adj"){pnameset = "adj_PVal"}
        select_colname <- paste0("Target_venn_logfc_",logFC,"_",pnameset,"_",p)
        target_Venn_genes <- df_Venn_target[[select_colname]] %>% na.omit()
    
        if(setequal(target_Venn_genes, com_genes$common)){
            com_genes = target_Venn_genes
        }else{
            cat("\n",blue,"Target_Venn_genes Not at all in ComVennGenes!",reset,"\n")
            cat("\n",yellow,"Target_Venn_genes: ",sort(target_Venn_genes),reset,"\n")
            cat("\n",yellow,"ComVennGenes: ",sort(com_genes$common),reset,"\n")
            flush.console()
            quit(save = "no")
        }
    }
    
    # Com DEGs
    write.csv(data.frame(Symbol = com_genes),file = file.path(data_dir_name,"03_Venn","01_intersect_genes.csv"))
    data$common_genes_df <- data.frame(Gene = com_genes,Group = data$target_gene_name)
    # Com Up DEGs
    write.csv(data.frame(Symbol = com_up_genes),file = file.path(data_dir_name,"03_Venn","02_intersect_up_genes.csv"))
    data$common_Up_genes_df <- data.frame(Gene = com_up_genes,Group = data$target_gene_name)
    # Com Down DEGs
    write.csv(data.frame(Symbol = com_down_genes),file = file.path(data_dir_name,"03_Venn","03_intersect_down_genes.csv"))
    data$common_Down_genes_df <- data.frame(Gene = com_down_genes,Group = data$target_gene_name)
    
    cat("\n",yellow,"Intersection of DEGs and Target gene: ",nrow(data$common_genes_df),"genes",reset,"\n")
    print(data$common_genes_df)
    flush.console()  # 强制刷新控制台输出
    
                 
    ##################################################################################################################    
    ############################################# Enrich #############################################################
    ##################################################################################################################
    cat("\n\n\n",yellow,"Start Enrich Analysis and Plot!",reset,"\n")
        # Bar Plot ALL        
        enrich_go_kegg(genes = data$common_genes_df$Gene,
                        species = data$species, 
                        p_name = enrich_pname,  # p.value/p.adj
                        filter_pvalue = 0.05,
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 1,
                        kegg_analysis_method = "online", # online/local
                        plot_type = enrich_figure_type, # NULL,barplot,dotplot,sankeyplot,circularplot
                        scankeyplot_mode = "sankey_buble", # buble_sankey/sankey_buble
                        fill_colors = basicR::get_colors(number = 4.1),
                        go_kegg_topn = go_kegg_topn,
                        description_wrap_width = description_wrap_width
                       ) -> result_enrich
        
        write.csv(result_enrich$result,file =file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_GO_KEGG_Enrich_result.csv")))
        
        lapply(names(result_enrich$plot[[enrich_figure_type]]),function(type){
            if(enrich_figure_type == "barplot"){
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12; height = 1.5 + rownum * 0.5
            }
            if(enrich_figure_type == "dotplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]@data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12 + max_lenth*0.07; height = 1.5 + rownum * 0.4
            }
            if(enrich_figure_type == "sankeyplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data)
                width = 12 + max_lenth*0.07; height = 4 + rownum * 0.4
            }    
            if(enrich_figure_type == "circularplot"){
                width = 8; height = 8 
            }
            if(enrich_figure_type == "treeplot"){
                width = 8; height = 8 
            }
            basicR::save_figure(result_enrich$plot[[enrich_figure_type]][[type]],
                                filename = file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_",stringr::str_to_title(type),"_",enrich_figure_type)),
                                width = width, 
                                height = height, 
                                res = 600, 
                                formats = c("pdf","png")
                               )
        })

        # Bar Plot Up
        enrich_go_kegg(genes = data$common_Up_genes_df$Gene,
                        species = data$species, 
                        p_name = enrich_pname,  # p.value/p.adj
                        filter_pvalue = 0.05,
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 1,
                        kegg_analysis_method = "online", # online/local
                        plot_type = enrich_figure_type, # NULL,barplot,dotplot,sankeyplot,circularplot
                        scankeyplot_mode = "sankey_buble", # buble_sankey/sankey_buble
                        fill_colors = basicR::get_colors(number = 4.1),
                        go_kegg_topn = go_kegg_topn,
                        description_wrap_width = description_wrap_width
                       ) -> result_enrich
        
        write.csv(result_enrich$result,file =file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_GO_KEGG_Enrich_result_Up.csv")))
        
        lapply(names(result_enrich$plot[[enrich_figure_type]]),function(type){
            if(enrich_figure_type == "barplot"){
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12; height = 1.5 + rownum * 0.5
            }
            if(enrich_figure_type == "dotplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]@data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12 + max_lenth*0.07; height = 1.5 + rownum * 0.4
            }
            if(enrich_figure_type == "sankeyplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data)
                width = 12 + max_lenth*0.07; height = 3 + rownum * 0.4
            }    
            if(enrich_figure_type == "circularplot"){
                width = 8; height = 8 
            }
            if(enrich_figure_type == "treeplot"){
                width = 8; height = 8 
            }
            
            basicR::save_figure(result_enrich$plot[[enrich_figure_type]][[type]],
                                filename = file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_",stringr::str_to_title(type),"_",enrich_figure_type,"_Up")),
                                width = width, 
                                height = height, 
                                res = 600, 
                                formats = c("pdf","png")
                               )
        })
        
        # Bar Plot Down
        enrich_go_kegg(genes = data$common_Down_genes_df$Gene,
                        species = data$species, 
                        p_name = enrich_pname,  # p.value/p.adj
                        filter_pvalue = 0.05,
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 1,
                        kegg_analysis_method = "online", # online/local
                        plot_type = enrich_figure_type, # NULL,barplot,dotplot,sankeyplot,circularplot
                        scankeyplot_mode = "sankey_buble", # buble_sankey/sankey_buble
                        fill_colors = basicR::get_colors(number = 4.1),
                        go_kegg_topn = go_kegg_topn,
                        description_wrap_width = description_wrap_width
                       ) -> result_enrich
        
        write.csv(result_enrich$result,file =file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_GO_KEGG_Enrich_result_Down.csv")))
        
        lapply(names(result_enrich$plot[[enrich_figure_type]]),function(type){
            if(enrich_figure_type == "barplot"){
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12; height = 1.5 + rownum * 0.5
            }
            if(enrich_figure_type == "dotplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]@data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]@data)
                width = 12 + max_lenth*0.07; height = 1.5 + rownum * 0.4
            }
            if(enrich_figure_type == "sankeyplot"){
                max_lenth <- result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data$Description %>% str_length() %>% max()
                rownum <- nrow(result_enrich$plot[[enrich_figure_type]][[type]]$layers[[3]]$data)
                width = 12 + max_lenth*0.07; height = 4 + rownum * 0.4
            }    
            if(enrich_figure_type == "circularplot"){
                width = 8; height = 8 
            }
            if(enrich_figure_type == "treeplot"){
                width = 8; height = 8 
            }
            
            basicR::save_figure(result_enrich$plot[[enrich_figure_type]][[type]],
                                filename = file.path(file.path(data_dir_name,"04_Enrich"),paste0("01_",stringr::str_to_title(type),"_",enrich_figure_type,"_Down")),
                                width = width, 
                                height = height, 
                                res = 600, 
                                formats = c("pdf","png")
                               )
        }) 
      
    ###################################################    PPI  #############################################################
    
    cat("\n\n\n",yellow,"Start PPI analysis and Plot!",reset,"\n")
    # 获区PPI string 数据
    ppi_df <- string_ppi_get_data(genes = data$common_genes_df$Gene, required_score = confidence, species = data$species)
    write.csv(ppi_df,file = file.path(data_dir_name,"05_PPI","00_string_interactions.csv"))
    
    ppi_df_select <- ppi_df %>% dplyr::select(preferredName_A,preferredName_B,score)
    
    node_df <- table(c(ppi_df_select$preferredName_A,ppi_df_select$preferredName_B)) %>% as.data.frame() %>% setNames(c("Gene","Dgree"))
    node_df <- merge(node_df,data$common_genes_df,by = "Gene")
    
    net <- network_graph_plot(edge_df = ppi_df_select,node_df = node_df,
                              edge_from_colname = "preferredName_A",edge_to_colname = "preferredName_B",
                              edge_weight_colname = "score",edge_group_colname = "preferredName_A",edge_group_levels = NULL,
                              node_name_colname = "Gene",node_group_colname = "Group",node_size_colname = "Dgree",
                              node_group_levels = NULL,node_levels = NULL,
                              sources_use = NULL, targets_use = NULL, top_n = 1, # 筛选相互关系的top数据
                              plot_circular = TRUE,layout_type = "linear", 
                              label_size = 4,label_repel = FALSE,base_size = 14,
                              ggraph_node_angle= TRUE,
                              label_radius_ratio = 1.1,label_align_mode = "inner",  # "center","inner","outer"
                              ggraph_center_node_angle = FALSE,
                              label_radius_center_ratio = 1,label_align_center_mode = "center",  # "center","inner","outer"
                              edge_directed = FALSE,
                              edge_arrow = FALSE,
                              edge_arrow_length = 4,
                              edge_arrow_angle = 20,
                              edge_arrow_type =  "closed",
                              edge_radian = -5,
                              edge_width_range = c(0.5,1.5),
                              node_size_range = c(3, 8),
                              node_shape = FALSE,node_center = NULL, # unique(df_miRNA_TF$mRNA),NULL
                              node_group_color = c(basicR::get_colors(number = 10.1,package = "ggsci",name = "jco"),basicR::get_colors(number = 50.1)),
                              edge_group_color = c(basicR::get_colors(number = 100.1),basicR::get_colors(number = 100.1)),
                              legend_position = "right",legend_box_spacing = 15,
                              legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize = c(FALSE,FALSE,FALSE,FALSE),
                              legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize = c(1,1,1,1),
                              legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize = c(1,0.02,1,1),
                              legend_breaknumber_NodeSize_EdgeSize = c(3,3),
                              legend_numberdigits_NodeSize_EdgeSize = c(0,3),
                              legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep("top",4)),
                              legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(0,4)),
                              legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(0.5,4)),
                              legend_layout = NULL,legend_every_position = NULL,
                              title_name = "",
                              plot_margin = c(15,15,10,15)  # 边距 b,l,t,r 
                      )
    
    save_figure(obj = net$ggplotGrob_data,
                filename =  file.path(file.path(data_dir_name,"05_PPI"),paste0("01_Network.pdf")),
                width = 7,
                height = 7,
                res = 600,
                formats = c("pdf","png") # "pdf","tiff","png","svg"
               )
    
    write.csv(node_df,file = file.path(data_dir_name,"05_PPI","01_Node_Data.csv"))
    # save data
    saveRDS(data,file = file.path(data_dir_name,"00_data_result","02_data.rds"))

    return(Res$result)
}
