# Group levels set
group_order_set <- function(vector,group_levels = NULL){
    # 设置分组水平
    
    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"
    
    # 获取列名
    vector <- unique(vector) %>% as.character()
    if(length(group_levels) == 0){
        if(tolower(vector) %in% c("control", "normal","norm") %>% any()){

            # 找出小写化后等于 "control" 或 "normal" 的列名
            target_col <- vector[tolower(vector) %in% c("control", "normal","norm")]
            
            # 如果存在符合条件的列
            if (length(target_col) > 0) {
              # 将目标列放在第一位，其余列保持原顺序但去掉目标列
              new_order <- c(target_col[1], setdiff(vector, target_col[1]))
              return(new_order)  
            }
        }else{
            return(sort(vector))  
        }
    }else{
        if(all(vector %in% group_levels)){
            return(group_levels)
        }else{
            cat("\n",yellow,"group_levels Not Match the vector!",reset,"\n")
        }
    }
}