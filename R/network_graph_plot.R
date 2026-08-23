#' Draw a configurable node-edge network graph
#'
#' Filter and standardize edge and node tables, construct an igraph object,
#' and draw a circular or non-circular ggraph network with configurable labels,
#' arrows, node and edge aesthetics, legends, and margins.
#'
#' @param edge_df A data frame containing one row per edge.
#' @param node_df A node data frame. Pass `NULL` to derive node names, degrees,
#'   and a common group from `edge_df`. The historical default is `node_df`.
#' @param edge_from_colname,edge_to_colname Names of edge source and target
#'   columns.
#' @param edge_weight_colname Optional edge-weight column name. If `NULL`, all
#'   weights are set to one; also set `top_n = NULL` in that case.
#' @param edge_group_colname Optional edge-group column name, or the original
#'   source or target column name. `NULL` assigns all edges to `Interaction`.
#' @param edge_group_levels Optional ordering of edge groups.
#' @param node_name_colname Name of the node identifier column.
#' @param node_group_colname Optional node-group column; `NULL` assigns `Node`.
#' @param node_size_colname Optional numeric node-size column; `NULL` uses one.
#' @param node_group_levels Optional ordering of node groups.
#' @param node_levels Optional node order used for a manually constructed
#'   circular layout with `node_center`.
#' @param sources_use,targets_use Optional vectors restricting source or target
#'   nodes.
#' @param top_n Proportion of highest-weight edges to retain. `NULL` disables
#'   weight filtering.
#' @param plot_circular Logical; use a circular layout when `TRUE`.
#' @param layout_type A layout name accepted by [ggraph::create_layout()].
#' @param ggraph_node_angle Logical; rotate outer labels with
#'   [ggraph::node_angle()].
#' @param ggraph_center_node_angle Logical; similarly rotate center labels.
#' @param label_size,label_repel Label size and overlap-repelling control.
#' @param base_size Base theme text size.
#' @param label_radius_ratio,label_radius_center_ratio Radial multipliers for
#'   outer and center labels.
#' @param label_align_mode,label_align_center_mode Label alignment, one of
#'   `"center"`, `"inner"`, or `"outer"`.
#' @param edge_directed Logical; construct a directed igraph when `TRUE`.
#' @param edge_radian Curvature strength for loop and arc edges.
#' @param edge_arrow Logical; draw edge arrows.
#' @param edge_arrow_length,edge_arrow_angle,edge_arrow_type Arrow length in
#'   millimetres, angle, and type (`"open"` or `"closed"`).
#' @param edge_width_range,node_size_range Two-element display ranges for edge
#'   widths and node sizes.
#' @param node_shape Logical; map node groups to shapes as well as colors.
#' @param node_center Optional vector of nodes placed at the center of a manual
#'   circular layout.
#' @param node_group_color,edge_group_color Color vectors for node and edge
#'   groups.
#' @param legend_position,legend_box_spacing Legend position and spacing from
#'   the plot in millimetres.
#' @param legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize Four logical
#'   switches for node-color, node-size, edge-color, and edge-width legends.
#' @param legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize Four legend column
#'   counts in the same order.
#' @param legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize Four legend
#'   key-size multipliers.
#' @param legend_breaknumber_NodeSize_EdgeSize Two break counts for node-size
#'   and edge-width legends.
#' @param legend_numberdigits_NodeSize_EdgeSize Two rounding precisions for
#'   node-size and edge-width legend labels.
#' @param legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize Four
#'   legend-title positions.
#' @param legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize Four title
#'   horizontal justifications.
#' @param legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize Four title
#'   vertical justifications.
#' @param legend_layout Optional legend-box direction, such as `"horizontal"`
#'   or `"vertical"`.
#' @param legend_every_position Optional legend-box justification.
#' @param title_name Plot title, or `NULL`.
#' @param title_margin_bottom Space below the title.
#' @param plot_margin Numeric bottom, left, top, and right margins in
#'   millimetres.
#' @param help Logical; print the built-in parameter guide and return
#'   invisibly without drawing when `TRUE`.
#'
#' @return When `help = FALSE`, a list containing `ggplot_p`, the ggraph/ggplot
#'   object, and `ggplotGrob_data`, a grob whose panel clipping is disabled so
#'   labels outside the panel remain visible. With `help = TRUE`, invisible
#'   `NULL`.
#'
#' @details
#' Source, target, and top-weight filters are applied before columns are
#' standardized. `top_n` uses the `(1 - top_n)` quantile and retains weights at
#' or above that threshold. Supplied node data are reduced after source or
#' target filtering. If `node_df = NULL`, nodes and degree-like sizes are
#' generated from the filtered edge endpoints.
#'
#' Circular plots can use a ggraph circular layout directly or, when
#' `node_center` is supplied, a manual inner/outer-ring layout. Non-circular
#' plots use `layout_type` with `circular = FALSE`. Self-loops and other edges
#' are drawn separately, with optional arrows.
#'
#' @examples
#' \dontrun{
#' edges <- data.frame(
#'   from = c("TP53", "TP53", "EGFR"),
#'   to = c("EGFR", "MYC", "MYC"),
#'   score = c(0.9, 0.8, 0.7)
#' )
#' net <- network_graph_plot(
#'   edge_df = edges,
#'   node_df = NULL,
#'   top_n = 1,
#'   plot_circular = TRUE,
#'   layout_type = "linear"
#' )
#' net$ggplot_p
#' grid::grid.draw(net$ggplotGrob_data)
#' }
#'
#' @seealso [igraph::graph_from_data_frame()], [ggraph::ggraph()],
#'   [ggraph::create_layout()]
#' @export
network_graph_plot  <- function(edge_df, node_df = node_df,
                                edge_from_colname = "from", edge_to_colname = "to",
                                edge_weight_colname = "score", edge_group_colname = NULL,
                                edge_group_levels = NULL,
                                node_name_colname = "name", node_group_colname = NULL, node_size_colname = NULL,
                                node_group_levels = NULL,node_levels = NULL,
                                sources_use = NULL, targets_use = NULL, top_n = 0.1, # 筛选相互关系的top数据
                                plot_circular = TRUE,layout_type = "linear",ggraph_node_angle = FALSE,ggraph_center_node_angle = FALSE,
                                label_size = 3,label_repel = FALSE,base_size = 16,label_radius_ratio = 1.2,label_align_mode = "inner", # "center","inner","outer"
                                label_radius_center_ratio = 1,label_align_center_mode = "center",  # "center","inner","outer"
                                edge_directed = TRUE,edge_radian = 0.5,
                                edge_arrow = TRUE,edge_arrow_length = 3,edge_arrow_angle = 20,edge_arrow_type = "closed",
                                edge_width_range = c(0.2, 2),
                                node_size_range = c(5, 10),
                                node_shape = TRUE,node_center = NULL,
                                node_group_color = basicR::get_colors(number = 100.1),
                                edge_group_color = basicR::get_colors(number = 100.1),
                                legend_position = "right",legend_box_spacing = 5,
                                legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize = c(TRUE,TRUE,TRUE,TRUE),
                                legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize = c(2,2,2,2),
                                legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize = c(1,1,1,1),
                                legend_breaknumber_NodeSize_EdgeSize = c(4,4),
                                legend_numberdigits_NodeSize_EdgeSize = c(0,2),
                                legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep("top",4)),
                                legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(0,4)),
                                legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(1,4)),
                                legend_layout = NULL,legend_every_position = NULL,
                                title_name = NULL,
                                title_margin_bottom = 10,
                                plot_margin = c(5,5,5,5),  # b,l,t,r
                                help = FALSE
                               ){
    
    ########################################## help #########################################################
    help_text <- "
    ===========================  network_graph_plot 参数说明  ===========================
    
    通用说明：用于绘制网络图（环形或普通布局），支持节点和边的多分组、多图例、图例布局等控制。
             如果要设置边颜色，则避免使用'from'或'to'作为列名
    
    【输入数据】
    ------------------------------------------------------------------------------------
    - edge_df: 边数据框，必须包含 from、to、可选 edge weight 和 edge group 信息
    - node_df: 节点数据框，包含节点名、可选分组和大小（默认自动生成，但是没有分类和权重信息）

    【筛选数据】
    - sources_use:  使用指定源数据
    - targets_use： 使用指定目标数据
    - top_n:        使用边权重的前百分比数据出图
    
    【边使用参数】
    ------------------------------------------------------------------------------------
    - edge_from_colname:       边的起点列名（默认 'from'）
    - edge_to_colname:         边的终点列名（默认 'to'）
    - edge_weight_colname:     边的权重列名（默认 'score'；为空时设为 1）
    - edge_group_colname:      边的分组变量名（可为上述起止列之一或数据列）
    - edge_group_levels:       手动设置边分组的 factor 顺序（默认自动）

    【节点使用参数】
    ------------------------------------------------------------------------------------
    - node_name_colname:       节点名称列名（默认 'name'）
    - node_group_colname:      节点分组列名（默认 NULL，统一为 'Node'）
    - node_size_colname:       节点大小列名（默认 NULL，统一为 1）
    - node_group_levels:       节点分组顺序（factor levels）
    - node_levels:             节点顺序（factor levels）

    【图形布局与标签】
    ------------------------------------------------------------------------------------
    - plot_circular:             是否为环形布局（TRUE）或普通布局（FALSE）
    - layout_type:               布局类型（circular:'tree','stress','linear';
                                          no circular:'circle','kk','dh','lgl','graphopt','mds','random','nicely','tree','dendrogram','manual','stress')
    - ggraph_node_angle:         是否自动设置外环标签角度对齐圆环文字
    - ggraph_center_node_angle:  是否自动设置中心节点标签角度对齐圆环文字
    - label_size:                标签字体大小（默认 3）
    - label_repel:               是否使用 ggrepel 方式避免标签重叠
    - label_radius_ratio:        非中心节点标签距中心比例（调整标签环的位置）
    - label_radius_center_ratio: 中心节点标签距中心比例（调整标签环的位置）
    - label_align_mode:          非中心节点标签对齐方式（'center','inner','outer'） 
    - label_align_center_mode:   中心节点标签对齐方式（'center','inner','outer'） 
    
    【视觉样式】
    ------------------------------------------------------------------------------------
    - edge_directed:           是否设置方向，由from到to，否则随机方向
    - edge_arrow:              是否带箭头
    - edge_arrow_length:       箭头长度
    - edge_arrow_angle:        箭头角度
    - edge_arrow_type:         箭头类型（'open','closed'）
    - edge_radian:             边的弧度 (0:直线；其它带弧度，并且正负代表弧度方向性)
    - edge_width_range:        边宽范围，如 c(0.2, 2)
    - node_shape:              是否对节点使用形状映射
    - node_size_range:         节点大小范围，如 c(5, 10)
    - node_center:             是否设置特定节点与中央环(默认NULL：就一个环。此参数是在plot_circular=TRUE前提下)
    - node_group_color:        节点颜色向量（用 basicR::get_colors() 等函数生成）
    - edge_group_color:        边颜色向量 （用 basicR::get_colors() 等函数生成）
    - plot_margin:             4边距
    
    【图例控制】
    ------------------------------------------------------------------------------------
    - legend_position:         图例位置（'bottom', 'right','top','left'）
    - legend_box_spacing:      图例与主图之间间距（单位 mm）
    - legend_layout:           图例布局方向（'horizontal','vertical'）
    - legend_every_position:   图例对齐（'center','left','right','bottom','top'）
    
    【图例显示开关与列数设置（四个图例对应）】
    - legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize:         是否显示对应图例，长度为 4 的逻辑向量
    - legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize:            图例列数，长度为 4 的数值向量
    - legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize:    图例键缩放，长度为 4 的数值向量
    - legend_breaknumber_NodeSize_EdgeSize:                         图例项显示个数，长度为 2 的数值向量
    - legend_numberdigits_NodeSize_EdgeSize:                        图例项显示小数位数，长度为 2 的数值向量
    - legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize:  图例标题位置，长度为 4 的逻辑向量（'top','bottom','left','right'）
    - legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize:     图例水平位置，长度为 4 的数值向量（0-1）
    - legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize:     图例垂直位置，长度为 4 的数值向量（0-1）

    【标题与文本】
    - title_name:          图Title名
    - title_margin_bottom: Title和主图之间的距离

    【输出】
    ------------------------------------------------------------------------------------
    - 返回一个列表：ggplot_p（原始 ggplot 对象），ggplotGrob_data（修复标签裁切的 grob 对象）
    
    ====================================================================================
    "
    # 2. 如果用户传入 help = TRUE，则直接打印帮助并退出
    if (help) {
        cat(help_text)
        return(invisible(NULL))
    }

    ####################################### main #########################################################
    suppressPackageStartupMessages({
        library(igraph)
        library(ggraph)
        library(ggplot2)
    })
    # 输入关联数据和节点信息数据绘制网络图
    # "tree","stress","linear"
    # "circle", "kk", "dh","lgl","graphopt","mds","random","nicely","tree","dendrogram","manual","stress"

    # colors
    red <- "\033[31m"
    green <- "\033[32m"
    yellow <- "\033[33m"
    blue <- "\033[34m"
    magenta <- "\033[35m"
    cyan <- "\033[36m"
    reset <- "\033[0m"

    # 筛选数据
    if(!is.null(sources_use)){
        # 筛选边数据
        edge_df <- edge_df %>% dplyr::filter(!!sym(edge_from_colname) %in% sources_use)

        # 筛选节点数据
        if(!is.null(node_df)){
            node_df <- node_df %>% dplyr::filter(!!sym(node_name_colname) %in% unique(c(edge_df[[edge_from_colname]],edge_df[[edge_to_colname]])))
        }
        
    }
    if(!is.null(targets_use)){
        # 筛选边数据
        edge_df <- edge_df %>% dplyr::filter(!!sym(edge_to_colname) %in% targets_use)

        # 筛选节点数据
        if(!is.null(node_df)){
            node_df <- node_df %>% dplyr::filter(!!sym(node_name_colname) %in% unique(c(edge_df[[edge_from_colname]],edge_df[[edge_to_colname]])))
        }
    }       
        
    if(!is.null(top_n)){
        # 筛选边数据
        threshold <- quantile(edge_df[[edge_weight_colname]], probs = 1 - top_n, na.rm = TRUE)
        edge_df <- edge_df %>% dplyr::filter(!!sym(edge_weight_colname) >= threshold)

        # 筛选节点数据
        if(!is.null(node_df)){
            node_df <- node_df %>% dplyr::filter(!!sym(node_name_colname) %in% unique(c(edge_df[[edge_from_colname]],edge_df[[edge_to_colname]])))
        }  
    }

    # 标准化列名
    # edge_df
    colnames(edge_df)[colnames(edge_df) == edge_from_colname] <- "from"
    colnames(edge_df)[colnames(edge_df) == edge_to_colname] <- "to"

    if (!is.null(edge_weight_colname)) {
        colnames(edge_df)[colnames(edge_df) == edge_weight_colname] <- "edge_weight"
    } else {
        edge_df$edge_weight <- 1
    }
    if (!is.null(edge_group_colname)) {
        if(edge_group_colname %in% colnames(edge_df)){
            colnames(edge_df)[colnames(edge_df) == edge_group_colname] <- "edge_group"
        }else if(edge_group_colname %in% c(edge_from_colname,edge_to_colname)){  
            group_select <- match(edge_group_colname, c(edge_from_colname, edge_to_colname))
            edge_df$edge_group <- edge_df[[c("from","to")[group_select]]]
        }else{
            cat("\n",magenta,"edge_group_colname does not match any column name in edge_df.",reset,"\n")
            stop("Function terminated due to unmatched edge_group_colname.")
        }
    } else {
        edge_df$edge_group <- "Interaction"
    }
    
    # node_df
    if(is.null(node_df)){
        # 未提供节点数据，节点数据从边数据中获取
        node_df <- table(c(edge_df$from, edge_df$to)) %>% as.data.frame() %>% setNames(c("name", "node_size"))
        node_df$node_group <- "Node"
        
    }else{
        # 提供节点数据，标准化节点数据
        colnames(node_df)[colnames(node_df) == node_name_colname] <- "name"
        # 设置顺序
        if (!is.null(node_group_colname)) {
            colnames(node_df)[colnames(node_df) == node_group_colname] <- "node_group"
        } else {
            node_df$node_group <- "Node"
        }
        if (!is.null(node_size_colname)) {
            colnames(node_df)[colnames(node_df) == node_size_colname] <- "node_size"
        } else {
            node_df$node_size <- 1
        }
    }

    # 列排序，
    edge_df <- edge_df %>% dplyr::select(from,to,edge_weight,edge_group,setdiff(colnames(.),c('from','to','edge_weight','edge_group')))

    
    
    # 构建图对象
    if(edge_directed){    
        graph <- graph_from_data_frame(d = edge_df, vertices = node_df, directed = TRUE)
    }else{
        graph <- graph_from_data_frame(d = edge_df, vertices = node_df, directed = FALSE)
    }
    # 因子设置（用于颜色一致性）
    # 节点分组水平设置 
    if(!is.null(node_group_levels)){
        # 给定了节点分组水平，则检查并用给定的节点分组碎屏
        if(is.null(node_group_colname)){
            # 当node group 为NULL，直接用后续添加的
            V(graph)$node_group <- factor(V(graph)$node_group, levels = unique(node_df$node_group))
        }else{
            # 当node group 不为NULL，则检查参数和数据是否匹配
            node_group_levels_check <- unique(node_df$node_group) %in% node_group_levels %>% all() 
            if(node_group_levels_check){
                V(graph)$node_group <- factor(V(graph)$node_group, levels = intersect(node_group_levels,unique(node_df$node_group))) 
            }else{
                cat("\n",magenta,"node_group_levels do not match the node grouping in the node data.",reset,"\n")
                cat("\n",yellow,"Use the unique values of the node group column as the levels.",reset,"\n")
                V(graph)$node_group <- factor(V(graph)$node_group, levels = unique(node_df$node_group))
            }
        }   
    }else{
        # 不给定了节点分组水平，则用节点分组数据的唯一值
        V(graph)$node_group <- factor(V(graph)$node_group, levels = unique(node_df$node_group))
    }
    
    # 边分组水平设置
    if(!is.null(edge_group_levels)){
        # 给定了边分组水平，则检查并用给定的边分组碎屏
        if(is.null(edge_group_colname)){
            # 当node group 为NULL，直接用后续添加的
            E(graph)$edge_group <- factor(E(graph)$edge_group, levels = unique(edge_df$edge_group))
        }else{
            # 当node group 不为NULL，则检查参数和数据是否匹配
            edge_group_levels_check <- unique(edge_df$edge_group) %in% edge_group_levels %>% all()
            if(edge_group_levels_check){
                E(graph)$edge_group <- factor(E(graph)$edge_group, levels = intersect(edge_group_levels,unique(edge_df$edge_group))) 
            }else{
                cat("\n",magenta,"edge_group_levels do not match the node grouping in the node data.",reset,"\n")
                cat("\n",yellow,"Use the unique values of the edge group column as the levels.",reset,"\n")
                E(graph)$edge_group <- factor(E(graph)$edge_group, levels = unique(edge_df$edge_group))
            }
        }  
    }else{
        # 不给定了边分组水平，则用边分组数据的唯一值
        E(graph)$edge_group <- factor(E(graph)$edge_group, levels = unique(edge_df$edge_group))
    }


    
    # 绘图
    # 设置位置
    if(label_align_mode == "center"){
        hjust1 = 0.5;hjust2 = 0.5
    }else if(label_align_mode == "inner"){
        hjust1 = 0;hjust2 = 1
    }else if(label_align_mode == "outer"){
        hjust1 = 1;hjust2 = 0
    }

    if(label_align_center_mode == "center"){
        hjust_center1 = 0.5;hjust_center2 = 0.5
    }else if(label_align_center_mode == "inner"){
        hjust_center1 = 0;hjust_center2 = 1
    }else if(label_align_center_mode == "outer"){
        hjust_center1 = 1;hjust_center2 = 0
    }


        
    # 设置节点标签
    if (ggraph_node_angle) {    
        aes_text <- aes(x = label_radius_ratio * x,
                        y = label_radius_ratio * y,
                        label = name, 
                        angle = ggraph::node_angle(x, y),
                        hjust = ifelse(y < 0 & abs(x) < 1e-6,  # 特判：在正下方（x≈0）
                                                  hjust1,  # 外对齐时靠左，内对齐时靠右
                                                  ifelse(x > 0, hjust1, hjust2)        # 其他情况按原逻辑
                                                 )
                       )

        if(ggraph_center_node_angle){
            # 额外设置中心节点
            aes_text_center <- aes(x = label_radius_center_ratio * x, 
                                   y = label_radius_center_ratio * y, 
                                   label = name, 
                                   angle = ggraph::node_angle(x, y),
                                   hjust = ifelse(y < 0 & abs(x) < 1e-6,  # 特判：在正下方（x≈0）
                                                  hjust_center1,  # 外对齐时靠左，内对齐时靠右
                                                  ifelse(x > 0, hjust_center1, hjust_center2)        # 其他情况按原逻辑
                                                 )
                                  )
        }else{
            aes_text_center <- aes(x = label_radius_center_ratio * x,
                                   y = label_radius_center_ratio * y,
                                   label = name,
                                   hjust = ifelse(y < 0 & abs(x) < 1e-6,  # 特判：在正下方（x≈0）
                                                  hjust_center1,  # 外对齐时靠左，内对齐时靠右
                                                  ifelse(x > 0, hjust_center1, hjust_center2)        # 其他情况按原逻辑
                                                 )
                                  )
        }
    }else{
        aes_text <- aes(x = label_radius_ratio * x,
                        y = label_radius_ratio * y,
                        label = name, 
                        hjust = ifelse(x > 0, hjust1, hjust2)
                       )
        aes_text_center <- aes(x = label_radius_center_ratio * x, 
                               y = label_radius_center_ratio * y,
                               label = name,
                               hjust = ifelse(y < 0 & abs(x) < 1e-6,  # 特判：在正下方（x≈0）
                                              hjust_center1,  # 外对齐时靠左，内对齐时靠右
                                              ifelse(x > 0, hjust_center1, hjust_center2)        # 其他情况按原逻辑
                                             )
                              )
    }

    # edge  Legend
    # 转换为 guide 类型向量
    legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize <- ifelse(legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize,"legend", "none")
    legend_title_position <- legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize
    legend_title_hjust <- legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize
    legend_title_vjust <- legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize 

    # 设置 范围分隔函数并设置节点图例
    get_integer_breaks <- function(range, n = 3, digits = NULL) {
        # 设置图例数据量
        if(is.null(digits)){
            unique(seq(range[1], range[2], length.out = n))
        }else{
            round(unique(seq(range[1], range[2], length.out = n)), digits)
        }
    }

    scale_to_range <- function(vec, new_range = c(0, 1)) {
        # 缩放
        old_min <- min(vec)
        old_max <- max(vec)
        new_min <- new_range[1]
        new_max <- new_range[2]
        scaled <- (vec - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
        return(scaled)
    }
    node_size_breaks <- get_integer_breaks(range(node_df$node_size), n = legend_breaknumber_NodeSize_EdgeSize[1],digits = legend_numberdigits_NodeSize_EdgeSize[1]) # 节点大小图例键数量
    edge_size_breaks <- get_integer_breaks(range(edge_df$edge_weight), n = legend_breaknumber_NodeSize_EdgeSize[2],digits = legend_numberdigits_NodeSize_EdgeSize[2])  # 边大小图例键数量
    
    #  各图例键缩放后大小
    NodeColorLegendSize <- 5 * legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize[1]
    NodeSizeLegendSize <- ((node_size_breaks * legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize[2]) + 1) * node_size_range[1] # 节点大小图例键缩放
    EdgeColorLegendSize <- 2 * legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize[3]
    EdgeSizeLegendSize <- scale_to_range(edge_size_breaks, edge_width_range) * legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize[4]  # 边大小图例键缩放
    
    # 构建 guides 列表
    guides_list <- list()
    if(legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize[1]){
        guides_list$fill <- guide_legend(override.aes = list(size = NodeColorLegendSize),title = "Node Group",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[1],
                                         title.position = legend_title_position[1],title.hjust = legend_title_hjust[1],title.vjust = legend_title_vjust[1])
        guides_list$color <- guide_legend(override.aes = list(size = NodeColorLegendSize),title = "Node Group",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[1],
                                         title.position = legend_title_position[1],title.hjust = legend_title_hjust[1],title.vjust = legend_title_vjust[1])
        guides_list$shape <- guide_legend(override.aes = list(size = NodeColorLegendSize),title = "Node Group",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[1],
                                         title.position = legend_title_position[1],title.hjust = legend_title_hjust[1],title.vjust = legend_title_vjust[1])
    }
    if(legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize[2]){
        guides_list$size <- guide_legend(override.aes = list(size = NodeSizeLegendSize),title = "Node Size",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[2],
                                         title.position = legend_title_position[2],title.hjust = legend_title_hjust[2],title.vjust = legend_title_vjust[2])
    }
    if(legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize[3]){
        guides_list$edge_color <- guide_legend(override.aes = list(edge_width = EdgeColorLegendSize),title = "Edge Group",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[3],
                                              title.position = legend_title_position[3],title.hjust = legend_title_hjust[3],title.vjust = legend_title_vjust[3])
    }
    if(legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize[4]){
        guides_list$edge_width <- guide_legend(override.aes = list(edge_width = EdgeSizeLegendSize),
                                               title = "Edge Weight",ncol = legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize[4],
                                               title.position = legend_title_position[4],title.hjust = legend_title_hjust[4],title.vjust = legend_title_vjust[4])
    }

    # 设置节点形状
    if(node_shape){
        node_group_shape <- setNames(rep(21:25, length.out = length(levels(V(graph)$node_group))), levels(V(graph)$node_group))       
    }

    
    if(plot_circular){
        
        # 环
        if(is.null(node_center)){
            # 没有设置中心节点
            # 先生成 layout_df
            layout_df <- create_layout(graph, layout = layout_type, circular = TRUE)
            layout_df$direction <- atan2(layout_df$y, layout_df$x) * 180 / pi  # 计算每个节点的法线角度
    
            scale_factor <- 2.5   # 缩放倍数
            layout_df$x <- layout_df$x * scale_factor
            layout_df$y <- layout_df$y * scale_factor
            layout_df$x2 <- layout_df$x * 0.97
            layout_df$y2 <- layout_df$y * 0.97

            edges_long <- ggraph::get_edges("long")(layout_df) %>%
                  dplyr::group_by(edge.id) %>%
                  dplyr::mutate(
                    pos = dplyr::row_number() - 1,
                    # 如果是终点 (pos == 1)，用 layout_df 的 x2
                    x = ifelse(pos == 1, layout_df$x2[node], x),
                    y = ifelse(pos == 1, layout_df$y2[node], y)
                  ) %>%
                  dplyr::ungroup()
            
            if(edge_arrow){
                p <-  ggraph(layout_df) +
                      geom_edge_loop(aes(color = edge_group, edge_width = edge_weight, direction = node1.direction), # 自环
                                       strength = edge_radian,
                                       alpha = 0.8,
                                       arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'), 
                                                           type = edge_arrow_type,
                                                           angle = edge_arrow_angle
                                                          ),
                                       end_cap = circle(2, "mm")
                                      ) +      
                      geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight), strength = edge_radian, alpha = 0.95, # 它环箭头
                                     arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'),
                                                           type = edge_arrow_type,
                                                           angle = edge_arrow_angle
                                                          ),
                                     edge_width = NA,
                                     end_cap = circle(1.5, "mm")
                                    ) + 
                      geom_edge_arc2(data = edges_long,aes(color = edge_group, edge_width = edge_weight),  # 它环曲线
                                     strength = edge_radian,
                                     alpha = 0.8,
                                     arrow = NULL,
                                     end_cap = circle(3, "mm")
                                    ) +
                      coord_fixed(clip = "off") 
            }else{
                p <-  ggraph(layout_df) +
                    geom_edge_loop(aes(color = edge_group, edge_width = edge_weight,direction = node1.direction), # 自环
                                   strength = edge_radian,alpha = 0.8) +
                    geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8) +  # 它环
                    coord_fixed(clip = "off")                  
            }

        }else{
            # 设置了中心节点
            # 所有节点
            node_names <- V(graph)$name
            
            # 拆分中心节点 & 外圈节点
            center_nodes <- intersect(node_center, node_names)  # 过滤无效节点
            outer_nodes <- setdiff(node_names, center_nodes)

            # 设置节点顺序
            if(!is.null(node_levels)){
                if(all(unique(node_levels) %in%  unique(node_names))){
                    outer_nodes <- setdiff(node_levels, center_nodes)
                }else{
                    cat("\n",magenta,"node_levels not all in  node_df.",reset,"\n")
                    stop('Please check the node_levels')
                }
            }
            
            # 坐标初始化
            layout_df <- data.frame(name = node_names, x = NA, y = NA, stringsAsFactors = FALSE)
            
            # ==== 1. 内圈中心节点 ====
            if (length(center_nodes) == 1) {
              layout_df[layout_df$name == center_nodes, c("x", "y")] <- c(0, 0)
            } else {
              theta_center <- seq(0, 2 * pi, length.out = length(center_nodes) + 1)[-1]
              r_center <- 0.4  # 内圈半径（不要太大）
              layout_df[layout_df$name %in% center_nodes, "x"] <- r_center * cos(theta_center)
              layout_df[layout_df$name %in% center_nodes, "y"] <- r_center * sin(theta_center)
            }
            
            # ==== 2. 外圈节点 ====
            theta_outer <- seq(0, 2 * pi, length.out = length(outer_nodes) + 1)[-1]
            r_outer <- 1  # 外圈半径
            idx <- match(outer_nodes, layout_df$name)  # 重新排序
            layout_df[idx, "x"] <- r_outer * cos(theta_outer)
            layout_df[idx, "y"] <- r_outer * sin(theta_outer)
            
            # 保证节点顺序匹配
            layout_df <- layout_df[match(node_names, layout_df$name), ]
            rownames(layout_df) <- NULL

            scale_factor <- 2.5   # 缩放倍数
            layout_df$x <- layout_df$x * scale_factor
            layout_df$y <- layout_df$y * scale_factor
            
            layout_df$direction <- atan2(layout_df$y, layout_df$x) * 180 / pi
            V(graph)$direction <- layout_df$direction
            
            if(edge_arrow){
                p <- ggraph(graph, layout = "manual", x = layout_df$x, y = layout_df$y) +
                        geom_edge_loop(aes(color = edge_group, edge_width = edge_weight,direction = node1.direction), # 自环
                                       span = 90,
                                       strength = edge_radian,
                                       alpha = 0.8,
                                       arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'),
                                                           type = edge_arrow_type,
                                                           angle = edge_arrow_angle
                                                          ),
                                       end_cap = circle(2, "mm")
                                      ) +
                        geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8, # 它环
                                       arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'),
                                                           type = edge_arrow_type,
                                                           angle = edge_arrow_angle
                                                          ),
                                       end_cap = circle(2, "mm")
                                      ) +
                        coord_fixed(clip = "off")
            }else{
                p <- ggraph(graph, layout = "manual", x = layout_df$x, y = layout_df$y) +
                        geom_edge_loop(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8) +  # 自环
                        geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8) +  # 它环
                        coord_fixed(clip = "off")                
            }
        }
        
    }else{
        # 非环
        layout_df <- create_layout(graph, layout = layout_type, circular = FALSE)
        layout_df$direction <- atan2(layout_df$y, layout_df$x) * 180 / pi  # 计算每个节点的法线角度

        scale_factor <- 2.5   # 缩放倍数
        layout_df$x <- layout_df$x * scale_factor
        layout_df$y <- layout_df$y * scale_factor 
        
        if(edge_arrow){

            p <- ggraph(layout_df) +
                    geom_edge_loop(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8,  # 自环
                                   arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'), 
                                                       type = edge_arrow_type,
                                                       angle = edge_arrow_angle
                                                      ),
                                   end_cap = circle(2, "mm")
                                  ) +
                    geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8,  # 它环曲线和箭头
                                   arrow = grid::arrow(length = unit(edge_arrow_length, 'mm'), 
                                                       type = edge_arrow_type,
                                                       angle = edge_arrow_angle
                                                      ),
                                   end_cap = circle(2, "mm")
                                  ) +            
                    coord_fixed(clip = "off")
        }else{
            p <- ggraph(layout_df) +
                    geom_edge_loop(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8) + # 自环
                    geom_edge_arc2(aes(color = edge_group, edge_width = edge_weight),strength = edge_radian,alpha = 0.8) + # 它环
                    coord_fixed(clip = "off")            
        }

        
    }
    # plot title
    p <- p + labs(title = title_name)
    
    # plot point       
    if(node_shape){
        p <- p + geom_node_point(aes(size = node_size,color = node_group, fill = node_group,shape = node_group)) +
                 scale_shape_manual(values = node_group_shape,guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[1])
    }else{
        p <- p + geom_node_point(aes(size = node_size,color = node_group, fill = node_group))
    }
    
    # plot theme     
     p <- p +
             scale_color_manual(values = node_group_color,
                                limits = levels(V(graph)$node_group),
                                guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[1]) +
             scale_fill_manual(values = node_group_color,
                               limits = levels(V(graph)$node_group),
                               guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[1]) +
             scale_size_continuous(range = node_size_range,
                                   guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[2],
                                   breaks = node_size_breaks) +
             scale_edge_color_manual(values = edge_group_color,
                                     limits = levels(E(graph)$edge_group),
                                     guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[3]) + 
             scale_edge_width(range = edge_width_range,
                              guide = legend_guide_NodeColor_NodeSize_EdgeColor_EdgeSize[4],
                              breaks = edge_size_breaks) +
             do.call(guides, guides_list) +
             theme_void(base_size = base_size) +
             theme(legend.position = legend_position,
                   legend.box.spacing = unit(legend_box_spacing, "mm"),
                   legend.title = element_text(face = "bold"),
                   plot.title = element_text(hjust = 0.5,vjust = 1,size = base_size * 1.2,margin = ggplot2::margin(b = title_margin_bottom)),
                   plot.margin = ggplot2::margin(b = plot_margin[1],l = plot_margin[2],t = plot_margin[3],r = plot_margin[4],unit = "mm")
                  )
    # 标签映射   
    if(plot_circular){
        # 设置非中心节点标签
        p <- p + geom_node_text(data = subset(layout_df, !(name %in% node_center)),
                    mapping = aes_text,size = label_size,show.legend = F,fontface = "bold", repel = label_repel
                   )
        # 单独设置中心节点标签
        if (!is.null(node_center)) {
            p <- p + geom_node_text(data = subset(layout_df, name %in% node_center),
                                    mapping = aes_text_center,size = label_size,show.legend = F,fontface = "bold", repel = label_repel
                                   )
        }
        
    }else{
        p <- p + geom_node_text(mapping = aes_text,size = label_size,show.legend = F,fontface = "bold",repel = label_repel)
    }    
        
    # 图例整体布局
    if(!is.null(legend_layout)){
        p <- p + theme(legend.box = legend_layout)
    }
    if(!is.null(legend_every_position)){
        p <- p + theme(legend.box.just = legend_every_position)
    }
    
    # 使用ggplotGrob 修正主图种标签会被截断的问题
    g <- ggplotGrob(p)
    g$layout$clip[g$layout$name == "panel"] <- "off"
    
    return(list(ggplot_p = p,ggplotGrob_data = g))
}
