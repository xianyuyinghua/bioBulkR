\name{network_graph_plot}
\alias{network_graph_plot}
\title{Draw a configurable node-edge network graph}
\usage{
network_graph_plot(
  edge_df,
  node_df = node_df,
  edge_from_colname = "from",
  edge_to_colname = "to",
  edge_weight_colname = "score",
  edge_group_colname = NULL,
  edge_group_levels = NULL,
  node_name_colname = "name",
  node_group_colname = NULL,
  node_size_colname = NULL,
  node_group_levels = NULL,
  node_levels = NULL,
  sources_use = NULL,
  targets_use = NULL,
  top_n = 0.1,
  plot_circular = TRUE,
  layout_type = "linear",
  ggraph_node_angle = FALSE,
  ggraph_center_node_angle = FALSE,
  label_size = 3,
  label_repel = FALSE,
  base_size = 16,
  label_radius_ratio = 1.2,
  label_align_mode = "inner",
  label_radius_center_ratio = 1,
  label_align_center_mode = "center",
  edge_directed = TRUE,
  edge_radian = 0.5,
  edge_arrow = TRUE,
  edge_arrow_length = 3,
  edge_arrow_angle = 20,
  edge_arrow_type = "closed",
  edge_width_range = c(0.2, 2),
  node_size_range = c(5, 10),
  node_shape = TRUE,
  node_center = NULL,
  node_group_color = basicR::get_colors(number = 100.1),
  edge_group_color = basicR::get_colors(number = 100.1),
  legend_position = "right",
  legend_box_spacing = 5,
  legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize = c(TRUE, TRUE, TRUE, TRUE),
  legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize = c(2, 2, 2, 2),
  legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize = c(1, 1, 1, 1),
  legend_breaknumber_NodeSize_EdgeSize = c(4, 4),
  legend_numberdigits_NodeSize_EdgeSize = c(0, 2),
  legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep("top", 4)),
  legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(0, 4)),
  legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize = c(rep(1, 4)),
  legend_layout = NULL,
  legend_every_position = NULL,
  title_name = NULL,
  title_margin_bottom = 10,
  plot_margin = c(5, 5, 5, 5),
  help = FALSE
)
}
\arguments{
\item{edge_df}{A data frame containing one row per edge.}
\item{node_df}{A node data frame. Pass \code{NULL} to derive node names,
degrees, and a common group from \code{edge_df}. The historical default is
\code{node_df}, so automatic generation requires explicitly supplying
\code{node_df = NULL}.}
\item{edge_from_colname, edge_to_colname}{Names of edge source and target
columns.}
\item{edge_weight_colname}{Optional edge-weight column name. If \code{NULL},
all weights are set to one. Because weight filtering occurs first, also set
\code{top_n = NULL} when no weight column is available.}
\item{edge_group_colname}{Optional edge-group column name, or the original
source or target column name. \code{NULL} assigns all edges to
\code{"Interaction"}.}
\item{edge_group_levels}{Optional ordering of edge groups.}
\item{node_name_colname}{Name of the node identifier column.}
\item{node_group_colname}{Optional node-group column; \code{NULL} assigns all
nodes to \code{"Node"}.}
\item{node_size_colname}{Optional numeric node-size column; \code{NULL} assigns
all nodes a size value of one.}
\item{node_group_levels}{Optional ordering of node groups.}
\item{node_levels}{Optional node order used for a manually constructed
circular layout when \code{node_center} is supplied.}
\item{sources_use, targets_use}{Optional vectors restricting source or target
nodes.}
\item{top_n}{The proportion of highest-weight edges to retain. The default
\code{0.1} retains edges at or above the 90th percentile. \code{NULL} disables
weight filtering.}
\item{plot_circular}{Logical; use a circular layout when \code{TRUE}.}
\item{layout_type}{A layout name accepted by
\code{\link[ggraph:create_layout]{ggraph::create_layout()}}.}
\item{ggraph_node_angle}{Logical; rotate outer labels with
\code{\link[ggraph:node_angle]{ggraph::node_angle()}}.}
\item{ggraph_center_node_angle}{Logical; similarly rotate center labels.}
\item{label_size, label_repel}{Label size and overlap-repelling control.}
\item{base_size}{Base theme text size.}
\item{label_radius_ratio, label_radius_center_ratio}{Radial multipliers for
outer and center labels.}
\item{label_align_mode, label_align_center_mode}{Outer and center label
alignment. Supported values are \code{"center"}, \code{"inner"}, and
\code{"outer"}.}
\item{edge_directed}{Logical; construct a directed igraph when \code{TRUE}.}
\item{edge_radian}{Curvature strength for loop and arc edges. Zero produces
straight arcs where supported; its sign controls curvature direction.}
\item{edge_arrow}{Logical; draw edge arrows.}
\item{edge_arrow_length, edge_arrow_angle, edge_arrow_type}{Arrow length in
millimetres, angle, and type. The type is normally \code{"open"} or
\code{"closed"}.}
\item{edge_width_range, node_size_range}{Two-element display ranges for edge
widths and node sizes.}
\item{node_shape}{Logical; map node groups to shapes as well as colors.}
\item{node_center}{Optional vector of nodes placed at the center of a manual
circular layout. Other nodes are placed on the outer ring.}
\item{node_group_color, edge_group_color}{Color vectors for node and edge
groups.}
\item{legend_position, legend_box_spacing}{Legend position and spacing from
the plot in millimetres.}
\item{legend_boolean_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four logical
switches for the node-color, node-size, edge-color, and edge-width legends.}
\item{legend_ncol_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four legend column
counts in the same order.}
\item{legend_keysizescale_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four legend
key-size multipliers in the same order.}
\item{legend_breaknumber_NodeSize_EdgeSize}{Two break counts for node-size and
edge-width legends.}
\item{legend_numberdigits_NodeSize_EdgeSize}{Two rounding precisions for
node-size and edge-width legend labels.}
\item{legend_title_position_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four legend
title positions.}
\item{legend_title_hjust_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four legend
title horizontal justifications.}
\item{legend_title_vjust_NodeColor_NodeSize_EdgeColor_EdgeSize}{Four legend
title vertical justifications.}
\item{legend_layout}{Optional legend-box direction, such as
\code{"horizontal"} or \code{"vertical"}.}
\item{legend_every_position}{Optional legend-box justification, such as
\code{"center"}, \code{"left"}, \code{"right"}, \code{"bottom"}, or
\code{"top"}.}
\item{title_name}{Plot title, or \code{NULL}.}
\item{title_margin_bottom}{Space below the title.}
\item{plot_margin}{Numeric bottom, left, top, and right margins in
millimetres.}
\item{help}{Logical; print the built-in parameter guide and return invisibly
without drawing when \code{TRUE}.}
}
\value{
When \code{help = FALSE}, a named list containing:
\describe{
\item{\code{ggplot_p}}{The ggraph/ggplot network object.}
\item{\code{ggplotGrob_data}}{A grob produced by
\code{\link[ggplot2:ggplotGrob]{ggplot2::ggplotGrob()}} whose panel clipping is
disabled so labels outside the panel remain visible.}
}
With \code{help = TRUE}, invisible \code{NULL}.
}
\description{
Filter and standardize edge and node tables, construct an igraph object, and
draw a circular or non-circular ggraph network with configurable labels,
arrows, node and edge aesthetics, legends, and margins.
}
\details{
Source, target, and top-weight filters are applied before columns are
standardized. \code{top_n} uses the \code{1 - top_n} quantile and retains
weights at or above that threshold. Supplied node data are reduced after
source or target filtering. If \code{node_df = NULL}, nodes and degree-like
sizes are generated from the filtered edge endpoints.

Edge columns are standardized internally to \code{from}, \code{to},
\code{edge_weight}, and \code{edge_group}. Node columns are standardized to
\code{name}, \code{node_group}, and \code{node_size}. Supplied factor-level
orders are used when they contain all observed groups; otherwise, the function
prints a message and uses observed order.

Circular plots can use a ggraph circular layout directly or, when
\code{node_center} is supplied, a manual inner/outer-ring layout.
Non-circular plots use \code{layout_type} with \code{circular = FALSE}.
Self-loops and other edges are drawn separately, with optional arrows.
}
\examples{
\dontrun{
edges <- data.frame(
  from = c("TP53", "TP53", "EGFR"),
  to = c("EGFR", "MYC", "MYC"),
  score = c(0.9, 0.8, 0.7)
)

net <- network_graph_plot(
  edge_df = edges,
  node_df = NULL,
  top_n = 1,
  plot_circular = TRUE,
  layout_type = "linear"
)

net$ggplot_p
grid::grid.draw(net$ggplotGrob_data)

# Print the function's built-in parameter guide
network_graph_plot(help = TRUE)
}
}
\seealso{
\code{\link[igraph:graph_from_data_frame]{igraph::graph_from_data_frame()}},
\code{\link[ggraph:ggraph]{ggraph::ggraph()}},
\code{\link[ggraph:create_layout]{ggraph::create_layout()}}
}
