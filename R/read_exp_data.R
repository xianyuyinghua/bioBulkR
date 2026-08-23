read_exp_data <- function(expression_dir = "./01_Expression_data", 
                          expression_pattern = "Expression.*\\.csv$", 
                          condition_pattern ="*Condition*.csv$"
                         ){
  # 读取表达数据文件夹的内容
  expr_file <- list.files(expression_dir, pattern = expression_pattern, recursive = T, full.names = T)
  condition_file <- list.files(expression_dir, pattern = condition_pattern, recursive = T, full.names = T)
  
  expr <- read_data(expr_file)
  condition <- read_data(condition_file)
  if(all(sapply(expr,class) == "numeric") != TRUE){
    expr[] <- sapply(expr, as.numeric)
  }
  expr <- na.omit(expr)
  
  # expression data type
  if(grepl("array",expr_file)){
    exp_data_type <- "array"
    throughput_type <- NULL
  }else if(grepl("throughput",expr_file)){
    exp_data_type <- "throughput"
    if(grepl("Count",expr_file)){
      throughput_type <- "Count"
    }else if (grepl("FPKM",expr_file)){
      throughput_type  <- "FPKM"
    }else if(grepl("TPM",expr_file)){
      throughput_type  <- "TPM"
    }else if(grepl("RPM",expr_file)){
      throughput_type  <- "RPM"
    }else{
      throughput_type  <- "unknown"
    }
  }
  
  return(list(expr = expr,condition = condition, exp_data_type = exp_data_type,throughput_type = throughput_type))
}