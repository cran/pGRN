#' pGRN: creates gene regulatory network based on single cell pseudotime information
#' 
#' Given single cell matrix and pseudotime, construct gene regulatory network (GRN)
#' 
#' @param expression_matrix expression matrix data
#' @param pseudotime_list list of pseudotime
#' @param method method for GRN construction: DTW, granger
#' @param slide_window_size	 sliding window size
#' @param slide_step_size sliding window step size
#' @param centrality_degree_mod (for DTW method) mode of centrality degree for popularity calculation
#' @param components_mod (for DTW method) mode of sub-network extraction methods (weak or strong)
#' @param network_min_genes minimal number of gene elements required for extracted sub-networks
#' @param quantile_cutoff	 an integer value (1-99) for quantile cutoff
#' @param order (for granger method) integer specifying the order of lags to include in the auxiliary regression 
#' @param cores number of cores for parallel computing
#'
#' @return a list of tabl_graph objects
#'
#' @export
#'  
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' 
#' # try DTW method
#' nets <- pGRN(expression_matrix,
#'              pseudotime_list, 
#'              method= "DTW",
#'              quantile_cutoff=50,
#'              cores=1)
#' plot_network(nets[[1]])
#' 
#' # plot the network interactively
#' plot_network_i(nets[[1]])
#' 
pGRN <- function(expression_matrix,
                 pseudotime_list,
                 method="DTW",
                 slide_window_size = 20,
                 slide_step_size = 10,
                 centrality_degree_mod="out",
                 components_mod="weak",
                 network_min_genes=10,
                 quantile_cutoff = 5,
                 order=1,
                 cores=1){
  if(method == "DTW"){
    adj_df <- run_dtw(expression_matrix, 
                      pseudotime_list,
                      slide_window_size=slide_window_size,
                      slide_step_size=slide_step_size,
                      quantile_cutoff=quantile_cutoff,
                      cores=cores)
    nets <- get_networks(adj_df,
                         components_mod = components_mod,
                         network_min_genes=network_min_genes) # a list of tbl_graph
    return(nets)
  }else if(method == "granger"){
    # time-consuming step for larger number of genes
    gt_adj_df <- run_granger_test(expression_matrix, 
                                  pseudotime_list,
                                  slide_window_size=slide_window_size,
                                  slide_step_size=slide_step_size,
                                  order=order)
    nets <- get_networks(gt_adj_df,
                         components_mod = components_mod,
                         network_min_genes=network_min_genes)
    return(nets)
  }else{
    stop("Error: method not found, current allowed method: DTW, granger")
  }
}