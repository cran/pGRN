#' @import ggplot2
#' @import ggraph
#' @import visNetwork
#' @importFrom tidygraph as_tbl_graph mutate centrality_degree as.igraph groups filter activate
#' @importFrom igraph components
#' @importFrom utils combn
#' @importFrom stats cutree dist hclust quantile sd
#' 
NULL


#' Convert distance matrix to adjacency dataframe
#' 
#' Convert distance matrix to adjacency dataframe for network construction.
#' 
#' @param data distance matrix
#' @param quantile_cutoff an integer value (1-99) for quantile cutoff
#' 
#' @return adjacency dataframe (with columns "from, to, distance,direction, similarity")
#' 
#' @export
#' 
#' @examples
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#' adj_df <- matrix2adj(dtw_dist_matrix)
#' 
matrix2adj <- function(data,
                       quantile_cutoff=5){
  quantile_cutoff <- paste0(quantile_cutoff,"%")
  fm <- flat_matrix(data)
  q_v <- quantile(abs(fm$distance),seq(0,1,0.01))[[quantile_cutoff]]
  fm_s <- fm[(fm$distance < q_v) & (fm$distance > -q_v),]
  d <- fm_s$distance
  d[d>0] <- 1
  d[d<0] <- -1
  fm_s$direction <- d
  fm_s$similarity <- 1 - abs(fm_s$distance)
  row.names(fm_s) <- NULL
  return(fm_s)
}


# TODO: interactive network: visNetwork or networkD3

#' Get the list of sub-networks
#' 
#' Get sub-networks based on given adjacency data.frame input
#'  
#' @param data adjacency data.frame
#' @param centrality_degree_mod mode of centrality degree for popularity calculation
#' @param components_mod mode of sub-network extraction methods
#' @param network_min_genes minimal number of gene elements required for extracted sub-networks
#'  
#' @return list of tabl_graph objects
#'  
#' @export
#'  
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#' adj_df <- matrix2adj(dtw_dist_matrix)
#' get_networks(adj_df,network_min_genes=5)
get_networks <- function(data,
                        centrality_degree_mod="out",
                        components_mod="weak",
                        network_min_genes=10){
  
  #colnames(data) <- c("from","to","dist","direction","similarity")
  
  tf_graph_full <- as_tbl_graph(data) %>%
    mutate(Popularity = centrality_degree(mode = centrality_degree_mod))
  
  tf_network_components <- components(tf_graph_full,mode=components_mod)
  
  # tf_groups <- tf_network_components$membership
  
  # TODO: filter subnetwork components
  
  tf_groups <- tf_network_components$membership
  
  network_list <- list()
  
  new_group_id <- 0
  for(group_id in unique(tf_groups)){
    
    from <- NULL
    tf_select <- names(tf_groups[tf_groups == group_id])
    
    tf_targets_df_s <- data %>% filter(from %in% tf_select)
    
    if(nrow(tf_targets_df_s) > network_min_genes){
      
      tf_graph <- as_tbl_graph(tf_targets_df_s) %>% 
        mutate(Popularity = centrality_degree(mode = centrality_degree_mod))
      new_group_id <- new_group_id + 1
      network_list[[new_group_id]] <- tf_graph
    }
  }
  
  return(network_list)
}


#' Plot stationary network
#' 
#' Plot stationary network through ggraph
#' 
#' @param graph a tbl_graph object
#' @param ... other parameters for ggraph
#' 
#' @return ggraph
#' 
#' @export
#' 
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#' nets <- module_networks(dtw_dist_matrix,k=1,quantile_cutoff=50)
#' plot_network(nets[["module1"]])
plot_network <- function(graph, ...){
  index <- Popularity <- name <- NULL
  ggraph(graph, layout = 'kk', ...) + 
    geom_edge_fan(aes(alpha = after_stat(index)), show.legend = T) + 
    geom_node_point(aes(size = Popularity)) + 
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white',base_family = 'Helvetica')+
    geom_node_text(aes(label = name, size=Popularity), repel = TRUE, point.padding = unit(0.2, "lines"),  colour="red")
}

#' Plot interactive network
#' 
#' Plot interactive network based on igraph layout input
#' 
#' @param graph igraph layout object
#' @param save_file file name of the saved file, not save if NULL
#' 
#' @return visNetwork htmlwidget
#' 
#' @export
#' 
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#' nets <- module_networks(dtw_dist_matrix,k=1,quantile_cutoff=50)
#' plot_network_i(nets[["module1"]])
plot_network_i <- function(graph, save_file = NULL){
  #net_vis <- visIgraph(graph,...)
  edges <- data.frame(graph %>% activate(edges))
  nodes <- data.frame(graph %>% activate(nodes))
  
  nodes$id <- row.names(nodes)
  nodes$label <- nodes$name
  nodes$value <- nodes$Popularity # size
  net_vis <- visNetwork(nodes,edges)
  if(!is.null(save_file)){
    visSave(net_vis, file=save_file)
  }
  return(net_vis)
}




flat_matrix <- function(data){
  xy <- t(combn(colnames(data), 2))
  data_f <- data.frame(xy, dist=data[xy])
  colnames(data_f) <- c("from","to","distance")
  return(data_f)
}


######## Clustering ##########
get_modules <- function(data,
                        k=5){
  dist_df <- dist(data)
  clusters <- hclust(dist_df)
  group_id <- cutree(clusters,k=k)
  return(group_id)
}

#' Get module level networks
#' 
#' Given a distance matrix, calculate gene modules based on hierarchical clustering method and then get module level networks
#' 
#' @param data distance matrix
#' @param k number of gene clusters for module inference
#' @param quantile_cutoff distance cutoff based on quantile(1-99) for edge identification
#' @param centrality_degree_mod "in" or "out" for nodes popularity calculation
#' @param components_mod "weak" or "strong" for sub-network components inference
#' @param network_min_genes minial number of genes required for a network
#' 
#' @return a list networks for each module
#' 
#' @export
#' 
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#' nets <- module_networks(dtw_dist_matrix,k=1,quantile_cutoff=50)
#' plot_network(nets[["module1"]])
module_networks <- function(data,
                            k=10,
                            quantile_cutoff=10,
                            centrality_degree_mod = "out",
                            components_mod = "weak",
                            network_min_genes = 10){
  modules <- get_modules(data,k=k)
  all_genes <- names(modules)
  net_list <- list()
  n <- 0
  for(m in unique(modules)){
    m_genes <- all_genes[modules==m]
    matrix_m <- data[m_genes,m_genes]
    adj_m <- matrix2adj(matrix_m,quantile_cutoff = quantile_cutoff)
    nets_m <- get_networks(adj_m,
                           centrality_degree_mod = centrality_degree_mod,
                           components_mod = components_mod,
                           network_min_genes = network_min_genes)
    for(net in nets_m){
      n <- n+1
      net_list[[paste0("module",n)]] <- net
    }
  }
  return(net_list)
}