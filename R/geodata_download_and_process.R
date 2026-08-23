# -*- coding: utf-8 -*-
# GEO 数据下载、表达矩阵与元数据处理函数
#
# 默认情况下，运行到需要人工判断的位置时仍会使用 readline() 交互。
# 对应参数传入非 NULL 值时，可跳过该项交互，用于非交互分析。

geodata_download_and_process <- function(
    workdir = "/project/project_analysis/26_26YLM199F/02_Result/",
    Project_ID = "Expression_DEG",
    Data_ID = "GSE28829",
    num = 1,
    gencode_annotation = "/reference/genome/human/GRCh38/GENCODE/gencode.v36.annotation.csv",
    gene_type = "protein_coding",
    duplicate_remove_method = "mean-max",
    supplementary_sample_model = "^GSM",
    process_rawdata = NULL,
    gpl_column = NULL,
    group_levels = NULL,
    throughput_data_type = NULL,
    count_file_keywords = NULL,
    count_gene_column = NULL,
    sample_column_keywords = NULL,
    sample_rename_mode = NULL
) {
    setwd(workdir)

    opt <- list()
    opt$Project_ID <- Project_ID
    opt$Data_ID <- Data_ID

    suppressPackageStartupMessages({
        library(optparse)
        library(stringr)
        library(bioBulkR)
        library(GEOquery)
        library(oligo)
        library(affy)
        library(limma)
        library(ggplot2)
        library(ggrepel)
        library(pheatmap)
        library(dplyr)
        library(affyPLM)
        library(data.table)
        library(clusterProfiler)
        library(viridis)
        library(cowplot) 
        library(tibble)
        library(ComplexHeatmap)
    })

    # ANSI 颜色代码
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    # colors
    colors <- c('#4DBBD5FF','#E64B35FF','#00A087FF','#3C5488FF','#925E9FFF','#91D1C2FF',
                '#8491B4FF','#7E6148FF','#0072B5FF','#E18727FF','#B09C85FF','#20854EFF',
                '#6F99ADFF','#FFDC91FF','#00468BFF','#FDAF91FF','#B24745FF','#6699FFFF',
                '#99991EFF','#FFCCCCFF','#358000FF','#99CCFFFF','#FFCC00FF','#8C564BFF',
                '#BCBD22FF','#996600FF','#5CB85CFF','#F39B7FFF','#CE3D32FF','#749B58FF',
                '#466983FF','#F0E685FF','#D595A7FF','#924822FF','#7A65A5FF','#C75127FF',
                '#FFA319FF','#8A9045FF','#8F3931FF','#00AF66FF','#748AA6FF','#D0DFE6FF',
                '#C71000FF','#008EA0FF','#8A4198FF','#D5E4A2FF','#5A9599FF','#FF6348FF',
                '#B7E4F9FF','#FF95A8FF','#526E2DFF','#FB6467FF','#E89242FF','#69C8ECFF',
                '#917C5DFF','#FED439FF','#709AE1FF','#D2AF81FF','#FD7446FF','#4DBBD5FF',
                '#E64B35FF','#00A087FF','#3C5488FF','#925E9FFF','#91D1C2FF', '#8491B4FF')

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
    if(is.na(Project_ID)){stop(paste0(red,"\n\tPlease set the 'Project ID'.",reset,"\n\n"),call. = FALSE)}
    if(is.na(Data_ID)){stop(paste0(red,"\tPlease set the 'Data ID'.",reset,"\n\n"),call. = FALSE)}

    if (length(commandArgs(trailingOnly = TRUE)) == 0) {cat(description)}

    cat("\n\n")
    print(str(opt))
    cat("\n\n")

    # 检查项目文件夹并创建
    folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE) %>% basename() %>% .[!grepl("^\\.", basename(.))]
    if(grepl(Project_ID,folders) %>% any()){
      cat(paste0("\n",magenta,"Warring  Project ID: ",reset,blue,Project_ID,reset,magenta," already exists in the current path.",reset,"\n\n"))
      project_folder <- grep(Project_ID,folders,value = T)
    }else{
      if(length(folders) == 0){
          dir_num <- 01
      }else{
          dir_num <- folders %>% strsplit(.,"_") %>% sapply(.,"[",1) %>% as.numeric() %>% max() %>%  `+`(1) %>% as.character()
      }
      if(str_length(dir_num) != 2){dir_num <- paste0("0",dir_num)}
      project_folder <- paste0(dir_num,"_",Project_ID)
    }

    # dir
    data_dir_name <- file.path(project_folder,Data_ID)
    exp_dir_name <- file.path(data_dir_name,"01_Expression_data")
    if (!dir.exists(exp_dir_name)){dir.create(exp_dir_name,recursive = T)}

    # 获取表达矩阵
    gset <- download_gset(data_id = Data_ID, destdir = exp_dir_name)
    eset <- exprs(gset[[num]]) %>% na.omit()  # 表达矩阵
    metadata <- pData(gset[[num]]) # 元信息
    data_type <- gset[[num]]@experimentData@other$type
    species <- gset[[num]]@experimentData@other$platform_taxid

    # 保存
    if(!dir.exists( file.path(data_dir_name,"00_data_result"))){dir.create( file.path(data_dir_name,"00_data_result"),recursive = T)}
    data <- list()
    data$gset <- gset
    data$eset <- eset
    data$metadata <- metadata
    data$data_type <- data_type
    data$species <- species
    saveRDS(data, file = file.path(data_dir_name, "00_data_result", "01_data.rds"))

    #---正式数据处理---------------------------------------------------------------------------------------------

    if(species == "9606"){
        species_OrgDb <- "org.Hs.eg.db"
        useMart_dataset <- "hsapiens_gene_ensembl"
    }else if(species == "10090"){
        species_OrgDb <- "org.Mm.eg.db"
        useMart_dataset <- "mmusculus_gene_ensembl"
    }
    
    if( grepl("array",data_type) ){
        #--处理array数据-----------------------------------------------------------------------------------
        cat("\n",magenta,"This is Array data.",reset,"\n")
    
        # 处理metadata数据
        group_df <- metadata_process(metadata = metadata)
    																		
        #--处理表达矩阵---------------------------------------------------------------------------------
        if( min(eset) >= 0 & max(eset) > 1000 ){
            # 要进行log转换的count
            eset <- log2(eset +1 ) 
            cat(paste0("\n",yellow,"The range of the Expression matrix: ",range(eset)),reset)
        }else if( min(eset) < 0 ){
          # 当有表达之小于0的情况，则处理原始数据
          cat(paste0("\n",yellow,"The range of the Expression matrix: ",range(eset)),reset,"\n")
          cat(paste0(blue,"\nDownload the raw data for processing?(yes/no)\n",reset))
          flush.console()  # 强制刷新控制台输出
          process_rawdata <- if (is.null(process_rawdata)) readline(prompt = "") else process_rawdata
          if(tolower(process_rawdata) == "yes"){
            raw_data_dirname <- file.path(data_dir_name,"00_Raw_data")
            if(!dir.exists(raw_data_dirname)){dir.create(raw_data_dirname,recursive = T)}
            cat(paste0(magenta,"\nDownloading raw data, please wait...",reset,"\n"))
            flush.console()  # 强制刷新控制台输出
            gset <- download_array_raw_gset(data_id = Data_ID, destdir = raw_data_dirname)
            CLLgcrma <- gcrma(gset)
            eset <- exprs(CLLgcrma)
        	cat(paste0("\n",yellow,"The range of the Expression matrix: ",range(eset)),reset,"\n")
          }
        }
    
        #--ID转换-------------------------------------------------------------------
        platform_id <- unique(metadata$platform_id)
        cat("\n",yellow,"Downloading GPL file,Please waiting...",reset,"\n")
        flush.console()  # 强制刷新控制台输出
        gpl <- download_gset(data_id = platform_id, destdir = exp_dir_name)
        if(class(gpl) == "GPL"){gpl <- Table(gpl)}
        # 判断是Symbol/assignment/RefSeq
        if(any(grepl("Symbol|assignment",colnames(gpl)))){
            probe_symbol_data <- gpl[, c("ID", grep("Symbol|assignment",colnames(gpl),value = T))]
        }else{
            cat("\n",yellow,"The head gpl:",reset,"\n")
            print(head(gpl))
            flush.console()  # 强制刷新控制台输出
            cat("\n",yellow,"Eg：",reset,"SYMBOL:Name / REFSEQ:GB_ACC / ENTREZID:ENTREZ_GENE_ID / ENSEMBL:ENSEMBL_ID",reset,"\n")
            flush.console()  # 强制刷新控制台输出
            select_gplcol <- if (is.null(gpl_column)) readline(prompt = "Please select the GPL column containing the SYMBOL/RefSeq: ") else gpl_column
            select_col_name <- sapply(strsplit(select_gplcol,":"),"[",2)
            gene_select_type <- sapply(strsplit(select_gplcol,":"),"[",1)
            if(grepl("SYMBOL",gene_select_type)){
                # Symbol
                probe_symbol_data <- gpl[, c("ID",select_col_name)]
            }else if(gene_select_type %in% c("ENSEMBL","ENTREZID","REFSEQ")){
                gple_select <- gpl[, c("ID",select_col_name)]
                geneid <- gple_select[[select_col_name]]
            
                if(data$species == 9606){
                    library(org.Hs.eg.db)
                    gene_transform <- clusterProfiler::bitr(geneid,
                                            fromType = gene_select_type,
                                            toType = "SYMBOL",
                                            OrgDb = org.Hs.eg.db
                                           )
                }else if(data$species == 10090){
                    library(org.Mm.eg.db)
                    gene_transform <- clusterProfiler::bitr(geneid,
                                            fromType = gene_select_type,
                                            toType = "SYMBOL",
                                            OrgDb = org.Mm.eg.db
                                           )
                }
 
                gpl_select_merge <- merge(gple_select,gene_transform,by.x = select_col_name,by.y = gene_select_type)
                probe_symbol_data <- gpl_select_merge[,c("ID","SYMBOL")]
            
            }else{
                cat("还需完善对应的脚本")
                stop()
            }
        }

    
        if(str_length(na.omit(probe_symbol_data[,2])) %>% median() > 20){
           probe_symbol_data$Gene_id <- sapply(strsplit(as.character(probe_symbol_data[,2]), " // "), `[`, 2)
        }else{
           probe_symbol_data$Gene_id <- probe_symbol_data[,2] 
        }
    
        anno1 <- probe_symbol_data[,c("ID", "Gene_id")] %>% na.omit()
        anno1 <- anno1[anno1$Gene_id != "",]
        anno1 <- anno1[anno1$Gene_id != "---",]
        anno1 <- anno1[anno1$Gene_id != "--",]
    
        # ID 转换
        gene_ID_type <- identify_gene_ID_type(gene_IDs = anno1$Gene_id[!grepl("/",anno1$Gene_id)])
        if(gene_ID_type == "UNKNOWN"){
             if(gene_select_type %in% c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")){
                 gene_ID_type <- "SYMBOL"
             }
        }
        gene_trans_df <- trans_gene_ID_type(gene_ID_type = gene_ID_type, useMart_dataset = useMart_dataset, gene_id = anno1$Gene_id, species_OrgDb = species_OrgDb,
                                            gencode_annotation = gencode_annotation,
                                            gene_type = gene_type  # all / protein_coding
                                           )
        anno1 <- merge(anno1, gene_trans_df, by.x = "Gene_id", by.y = "Gene_id")
        anno1 <- anno1[, c("ID","Symbol")]
        cat(paste0(yellow,"\nThe platform annotation infomation:",reset,"\n"))
        print(head(anno1))
    
        # 检查eset 和anno1 探针ID是否能匹配
        if(!any(anno1$ID %in% rownames(eset))){
          stop(paste0(red,"\n\tThe probe IDs in the expression matrix do not match those in the platform annotation file. Please verify.",reset,"\n\n"),call. = FALSE)
        }
    
        # 转换表达矩阵中 探针ID 为 Symbol
        eset <- eset %>% data.frame(check.names = F) %>%  mutate(ID = rownames(.))
        merg <- merge(eset, anno1, by="ID")
        gene <- lapply(merg$Symbol, function(y) { strsplit(as.character(y)," /// ")[[1]][1]}) %>% unlist()
        merg$gene <- gene
    
        # 均值去重
        aggr <- aggregate_duplicate_gene_remove(data = merg ,gene_col = "gene",sample_col = setdiff(colnames(eset),"ID"),duplicate_remove_method = duplicate_remove_method)
        aggr <- aggr %>% tibble::column_to_rownames(var = "gene")

        colnames(aggr) <- regmatches(colnames(aggr), gregexpr("GSM[0-9]+", colnames(aggr))) %>% unlist()
 
        # 匹配表达矩阵和元数据
        condition <- data.frame(row.names = group_df$Sample,Group = group_df$Group)
        expr <- subset(aggr,select = rownames(condition))
        if(identical(rownames(condition),colnames(expr))){
          cat(paste0("\n",yellow,"Sample names are identical between expr and condition.",reset,"\n")) 
          print(condition$Group)
        }
        expr <- na.omit(expr)

        # 按照Control Disease的顺序排列
        group_names <- condition$Group %>% unique()
        if( grepl("normal|control",tolower(group_names)) %>% any() ){
            # 当分组中包含 Normal 或 Control
            number <- tolower(group_names) %in% c("normal","control") %>% which()
            group_levels <- c( group_names[number],setdiff(group_names,group_names[number]) )
        
        }else{
            # 当分组中不包含 Normal 或 Control，人工设置顺序
            cat(paste0("\n",yellow,"The group names: ",group_names,reset,"\n"))
            flush.console()  # 强制刷新控制台输出
            group_levels <- if (is.null(group_levels)) readline(prompt = "Please enter the group levels：") else group_levels
        
        }
    
        condition <- condition %>%  arrange(factor(Group, levels = group_levels), rownames(.))
        expr <- expr[,rownames(condition)]
    
        # 剔除在所有样本中都为0的基因
        expr <- expr[rowSums(expr) != 0, ]
    
        cat("\n",yellow,"The array expression range: ",reset,magenta,print(range(expr)),reset,"\n")
                          
        # 保存至表达文件夹中
        write.csv(expr, file = paste0(exp_dir_name,'/01_Expression_array.csv'))
        write.csv(condition, file = paste0(exp_dir_name,'/02_Condition.csv'))
        save(expr,condition,file = file.path(exp_dir_name,paste0(Data_ID,"_Expression_array.rda")))	
    
    }else if( grepl("throughput",data_type) ){
        # 处理throughput数据
        raw_data_dirname <- file.path(data_dir_name,"00_Raw_data")
        if(!dir.exists(raw_data_dirname)){dir.create(raw_data_dirname,recursive = T)}
    
        cat("\n",magenta,"This is Throughput data.",reset,"\n")
        cat("\n",yellow,"Please download the corresponding high-throughput 'count' data or the data in the 'Supplementary file' to '00_Raw_data' directory.",reset,"\n")
        cat("\n",yellow,"If available, please download the corresponding high-throughput annotation data to the ‘00_Raw_data’ directory.",reset,"\n")
        cat("\n",yellow,"If available, please download the corresponding high-throughput FPKM data to the ‘01_Expression_data’ directory.",reset,"\n")
        flush.console()  # 强制刷新控制台输出
        throughput_data_type <- if (is.null(throughput_data_type)) readline(prompt = "Please enter the Gene Expression data type to download (count/supplementary)：") else throughput_data_type
    
        if(throughput_data_type == "count"){
            # 1.下载的是count数据
            counts_files <- list.files(raw_data_dirname,pattern = "counts", recursive = T,full.names = T)
            if(length(counts_files) < 1){
                cat("\n",yellow,"No files containing 'counts' were matched.",reset,"\n")
                flush.console() # 强制刷新控制台输出
                file_name_key_words <- if (is.null(count_file_keywords)) readline(prompt = "Please input the key words in the count file names: " ) else count_file_keywords
                counts_files <- list.files(raw_data_dirname, pattern = file_name_key_words, recursive = T,full.names = T)
            }
            counts_file <- counts_files[grep(paste0("^",Data_ID),basename(counts_files))]
            counts_exp <- read_data(counts_file) %>% data.frame()
            data_show <- counts_exp[1:2,]
            for (col_name in names(data_show)) {
                cat("\n",magenta, "Column:",reset,yellow,col_name, reset,"\n")
                print(data_show[[col_name]])
            }
            flush.console()  # 强制刷新控制台输出
        
            # 识别 gene ID 列
            if(any(grepl("^GSM",colnames(counts_exp)))){
                # GSM开头的样本名
                if( any(grepl("^GSM",colnames(counts_exp))) & length(colnames(counts_exp)[grepl("^GSM[0-9]+",colnames(counts_exp))])  > 1 ){
                    gene_col <- colnames(counts_exp)[!grepl("^GSM[0-9]+",colnames(counts_exp))]
                    if(length(gene_col) > 1 ){
                        gene_col <- if (is.null(count_gene_column)) readline(prompt = "Please input the name of gene ID columns: " ) else count_gene_column
                    }
                }else{
                    gene_col <- if (is.null(count_gene_column)) readline(prompt = "Please input the name of gene ID columns: " ) else count_gene_column
                }
                sample_col_names <- grep("^GSM",colnames(counts_exp),value = T)
            }else{
                # 非 GSM开头的样本名
                sample_name_key_word <- if (is.null(sample_column_keywords)) readline(prompt = "Please input the key words in the count data columns: " ) else paste(sample_column_keywords, collapse = ",")
                sample_name_key_words <- sample_name_key_word %>% strsplit("[,]") %>% unlist()
                matched_indices <- lapply(sample_name_key_words, function(key_word) { 
                    sample_name <- grep(key_word, colnames(counts_exp),value = T) 
                    return(sample_name)
                })
                sample_col_names <- matched_indices %>% unlist()
            }
        
            # gene ID 转换
            gene_ID_type <- identify_gene_ID_type(gene_IDs = counts_exp[[gene_col]][1:5] )  # 为软件ID转换做基因类型识别
            annot_file <- list.files(raw_data_dirname,pattern = "annot", recursive = T,full.names = T)  # 查看是否有匹配的注释文件
            if(length(annot_file) == 1){
                # 用注注释文件进行ID转换
                cat("\n\n",yellow,"Use Annotation files: ",reset,annot_file,reset,"\n\n")
                flush.console()
            
                annot_df <- data.table::fread(annot_file) %>% data.frame()
                gene_trans_df <- annot_df %>% dplyr::filter(GeneType == "protein-coding") %>% dplyr::select(GeneID,Symbol) %>% dplyr::rename("Gene_id" = "GeneID") 
            
            }else{
                # 用 软件进行ID转换
                gene_trans_df <- trans_gene_ID_type(gene_ID_type = gene_ID_type, useMart_dataset = useMart_dataset, gene_id = counts_exp[[gene_col]], species_OrgDb = species_OrgDb,
                                                    gencode_annotation = gencode_annotation,
                                                    gene_type = gene_type  # all / protein_coding
                                                   )
            }
        
        
            counts_exp[[gene_col]] <- as.character(counts_exp[[gene_col]])
            counts_df <- merge(counts_exp, gene_trans_df, by.x = gene_col, by.y = "Gene_id")
        
            # 均值最大值去重
            count <- aggregate_duplicate_gene_remove(data = counts_df, gene_col = "Symbol", sample_col = sample_col_names, duplicate_remove_method = duplicate_remove_method)
            count <- count %>% tibble::column_to_rownames(var = "Symbol")
        
            # metadata 处理
            group_df <- metadata_process(metadata = metadata)
            condition <- data.frame(row.names = group_df$Sample,Group = group_df$Group, Raw_Group = group_df$Raw_Group)
        
            # 筛选匹配的样本表达数据
            sample_names <- intersect(rownames(condition) , colnames(count))
            expr <- subset(count,select = sample_names)
            condition <- condition[sample_names,]
        
            if( min(range(expr)) < 0 | max(range(expr)) < 10000 ){
                count_type <- "Not_Count" 
            }else if( min(range(expr)) >= 0 & max(range(expr)) > 10000){
                count_type <- "Count" 
            }
            expr <- na.omit(expr)

            # 按照Control Disease的顺序排列
            group_names <- condition$Group %>% unique()
            if( grepl("normal|control",tolower(group_names)) %>% any() ){
                # 当分组中包含 Normal 或 Control
                number <- tolower(group_names) %in% c("normal","control") %>% which()
                group_levels <- c( group_names[number],setdiff(group_names,group_names[number]) )
            
            }else{
                # 当分组中不包含 Normal 或 Control，人工设置顺序
                cat(paste0("\n",yellow,"The group names: ",group_names,reset,"\n"))
                flush.console()  # 强制刷新控制台输出
                group_levels <- if (is.null(group_levels)) readline(prompt = "Please enter the group levels：") else group_levels
            
            }
        
            condition <- condition %>%  arrange(factor(Group, levels = group_levels), rownames(.))
            expr <- expr[,rownames(condition)]
        
            # 剔除在所有样本中都为0的基因
            #expr <- expr[rowSums(expr) != 0, ]

            # FPKM Data
            file_fpkms <- list.files(exp_dir_name,pattern = "FPKM",full.names = T,recursive = T)
            file_fpkm <- file_fpkms[!grepl("log",file_fpkms)]
            if(length(file_fpkm) > 0){
                # 存在 FPKM文件
                if(length(file_fpkm) == 1){
                    df_fpkm <- read_data(file_fpkm)
                    df_fpkm <- merge(df_fpkm, gene_trans_df, by.x = "GeneID", by.y = "Gene_id") %>%  as.data.frame() %>% 
                            aggregate_duplicate_gene_remove(data = ., gene_col = "Symbol", sample_col = rownames(condition), duplicate_remove_method = duplicate_remove_method) %>% 
                            tibble::column_to_rownames("Symbol")
                    if(max(df_fpkm) > 1000){
                        df_fpkm <- log2(df_fpkm + 1)
                    }
                    cat("\n",yellow,"The Throughput FPKM range: ",reset,magenta,print(range(df_fpkm)),reset,"\n")
                    write.csv(df_fpkm, file = paste0(exp_dir_name,"/01_Expression_throughput_","log2_FPKM",".csv"))              
                }else{
                    # 有多个FPKM
                    cat("\n",yellow,"The number of FPKM data exceeds 1，please Check it.",reset,"\n")
                    cat(file_fpkm)
                    stop()
                }

            }else{
                # 不存在FPKM文件
                cat("\n",yellow,"No FPKM data was found!",reset,"\n")
            
            }

        
            # Print range
            cat("\n",yellow,"The Throughput count range: ",reset,magenta,print(range(expr)),reset,"\n")
            flush.console()
        
                                      
            # 保存至表达文件夹中
            write.csv(expr, file = paste0(exp_dir_name,"/01_Expression_throughput_",count_type,".csv"))
            write.csv(condition, file = paste0(exp_dir_name,'/02_Condition.csv'))
            save(expr,condition,file = file.path(exp_dir_name,paste0(Data_ID,"_Expression_throughput_",count_type,".rda")))
        
        }else if(throughput_data_type == "supplementary") {
            # 2.下载的是supplementary数据
            raw_data_dirname <- file.path(data_dir_name,"00_Raw_data")
            if(!dir.exists(raw_data_dirname)){dir.create(raw_data_dirname,recursive = T)}
            supplementary_file <- list.files(raw_data_dirname, pattern = paste0("^",Data_ID), recursive = T,full.names = T)
            supplementary_datalist <- supplementary_data_process(file_path = supplementary_file, read_data_function = read_data, sample_modle = supplementary_sample_model )
            gene_col <- colnames(supplementary_datalist$gene_expression_data)[lapply(supplementary_datalist$gene_expression_data,class) == "character"]
            sample_col_names <- setdiff(colnames(supplementary_datalist$gene_expression_data),gene_col)
            gene_ids <-  supplementary_datalist$gene_expression_data[[gene_col]] %>% .[!startsWith(., "#")]

            #拆分多类型ID （需要完善函数）
            #supplementary_datalist$gene_expression_data[[gene_col]] <- multiple_gene_type_process(gene_ids = gene_ids,special_symbols = c("(", ")", "[", "]", "{", "}", "#", "&", "@", "$", "*"))
        
            # ID转换
            gene_ID_type <- identify_gene_ID_type(gene_IDs = supplementary_datalist$gene_expression_data[[gene_col]])
            gene_trans_df <- trans_gene_ID_type(gene_ID_type = gene_ID_type, 
                                                useMart_dataset = useMart_dataset, 
                                                gene_id = supplementary_datalist$gene_expression_data[[gene_col]], 
                                                species_OrgDb = species_OrgDb,
                                                gencode_annotation = gencode_annotation,
                                                gene_type = gene_type  # all / protein_coding
                                               )
            counts_df <- merge(supplementary_datalist$gene_expression_data, gene_trans_df, by.x = gene_col, by.y = "Gene_id")

            # 均值去重
            count <- aggregate_duplicate_gene_remove(data = counts_df, gene_col = "Symbol", sample_col = sample_col_names, duplicate_remove_method = duplicate_remove_method)
            count <- count %>% tibble::column_to_rownames(var = "Symbol")    
                                      
            # metadata 处理
            group_df <- metadata_process(metadata = metadata)
            condition <- data.frame(row.names = group_df$Sample, Group = group_df$Group, Raw_Group = group_df$Raw_Group)                              

            # 替换count列名中的'-'为'_'; condition 中的Raw_Group列中的'-'为'_'
            colnames(count) <- gsub("-","_",colnames(count))
            condition$Raw_Group <- gsub("-","_",condition$Raw_Group)
        
            # 筛选匹配的样本表达数据
            if( all(rownames(condition) %in% colnames(count)) ){
                expr <- subset(count,select = rownames(condition))
            }else{
                if(all(condition$Raw_Group %in% colnames(count) )){
                    rownames(condition) <- condition$Raw_Group
                    expr <- subset(count,select = rownames(condition))
                }else{
                    com_str <- find_common_string(strings = colnames(count))
                    count_replace_colnames <- gsub(com_str,"",colnames(count))
                    if( all(condition$Raw_Group %in% count_replace_colnames )){
                        colnames(count) <-  count_replace_colnames
                        rownames(condition) <- condition$Raw_Group
                    }else{
                        cat("\n",yellow,"The conditon Raw_Group: ",reset,"\n")
                        lapply(condition$Raw_Group,function(sample){print(sample)})
                        cat("\n\n",yellow,"The conditon sample: ",reset,"\n")
                        lapply(rownames(condition),function(sample){print(sample)})
                        cat("\n\n",yellow,"The expression sample: ",reset,"\n")
                        lapply(colnames(count),function(sample){print(sample)})
                        flush.console()  # 强制刷新面板
                        cat("\n",yellow,"Please set sample renaming. Enter the original and new sample names, mapped by '::', separated by ','",reset,blue,"\n\t the original name can be a pattern.",reset,magenta,"\n\t eg. recurrent miscarriage::sample1,elective termination::sample2",reset)
                        flush.console()  # 强制刷新控制台输出
                        sample_rename_mode <- if (is.null(sample_rename_mode)) readline(prompt = "") else sample_rename_mode
                        if(sample_rename_mode != ""){
                              rename_pairs  <- strsplit(sample_rename_mode, ",") %>% unlist() %>% lapply(function(pair) strsplit(pair,"::")[[1]])
                              sample_newname_list <- sapply(rename_pairs,"[",2)                                                                     
                              names(sample_newname_list) <- sapply(rename_pairs,"[",1)
                              sample_df <- data.frame(new_sample = sample_newname_list)
                              if(all(rownames(sample_df) %in% colnames(count))){
                                  # 修改表达矩阵sample name
                                  for( sample in names(sample_newname_list) ){ colnames(count)[ colnames(count) %in% sample ] <- sample_newname_list[sample] }
                              }else if(all(rownames(sample_df) %in% rownames(condition))){
                                  # 修改元数据 smaple name
                                  for( sample in names(sample_newname_list) ){ rownames(condition)[ rownames(condition) %in% sample ] <- sample_newname_list[sample] }
                              }
                                                                                                 
                        }else{
                             stop() 
                        }
                    }
                    expr <- subset(count,select = rownames(condition))
                }
            }
            expr <- expr[!grepl("^#",rownames(expr)),]
                                      
            # another data gene ID trans
            if(!is.null(supplementary_datalist$gene_another_data)){
                # 还需完善
                another_df <- merge(supplementary_datalist$gene_another_data, gene_trans_df, by.x = gene_col, by.y = "Gene_id") 
            }      
                                     
    
            expr <- na.omit(expr)  
                                                                                                     
            # 按照Control Disease的顺序排列
            group_names <- condition$Group %>% unique()
            if( grepl("normal|control",tolower(group_names)) %>% any() ){
                # 当分组中包含 Normal 或 Control
                number <- tolower(group_names) %in% c("normal","control") %>% which()
                group_levels <- c( group_names[number],setdiff(group_names,group_names[number]) )
            
            }else{
                # 当分组中不包含 Normal 或 Control，人工设置顺序
                cat(paste0("\n",yellow,"The group names: ",group_names,reset,"\n"))
                flush.console()  # 强制刷新控制台输出
                group_levels <- if (is.null(group_levels)) readline(prompt = "Please enter the group levels：") else group_levels
            
            }
        
            condition <- condition %>%  arrange(factor(Group, levels = group_levels), rownames(.))
            expr <- expr[,rownames(condition)]
        
            # 剔除在所有样本中都为0的基因
            expr <- expr[rowSums(expr) != 0, ]
                                                                                                     
            cat( "\n",magenta,"The Throughput",yellow,supplementary_datalist$gene_exp_type,reset,magenta,"count range: ",print(range(expr)),reset,"\n")
        
            # 保存至表达文件夹中
            write.csv(expr, file = paste0(exp_dir_name,"/01_Expression_throughput_",supplementary_datalist$gene_exp_type,".csv"))
            write.csv(condition, file = paste0(exp_dir_name,'/02_Condition.csv'))
            save(expr,condition,file = file.path(exp_dir_name,paste0(Data_ID,"_Expression_throughput_",supplementary_datalist$gene_exp_type,".rda")))
            if(!is.null(supplementary_datalist$gene_another_data)){
                write.csv(another_df, file = paste0(exp_dir_name,'/03_Another_data.csv'),quote=F)          
            }
        }    
    }

    invisible(list(
        project_folder = project_folder,
        data_dir = data_dir_name,
        expression_dir = exp_dir_name,
        data = data,
        expr = if (exists("expr", inherits = FALSE)) expr else NULL,
        condition = if (exists("condition", inherits = FALSE)) condition else NULL
    ))
}
