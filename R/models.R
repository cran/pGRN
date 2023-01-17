#' @import lmtest
#' @importFrom future plan multisession
#' 
NULL



#' Use Granger-causality Test to get gene-gene regulatory relationship
#' 
#' Based on single-cell gene expression matrix and pseudotime, calculate Granger-causality Test
#' based gene-gene regulatory relationship
#' 
#' @param data gene expression matrix (Gene * Cells)
#' @param ptime pseudotime matched with the column cells of the gene expression matrix
#' @param slide_window_size sliding window size
#' @param slide_step_size sliding window step size
#' @param pvalue_cutoff cutoff for the pvalue from transfer entropy test
#' @param order integer specifying the order of lags to include in the auxiliary regression
#' @param ... other parameters for grangertest function in lmtest
#' 
#' @return adjacency data frame
#' 
#' @export
#' 
#' @examples
#' 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' gt_adj_df <- run_granger_test(expression_matrix, pseudotime_list)
#'

run_granger_test <- function(data,
                             ptime,
                             slide_window_size = 20,
                             slide_step_size = 10,
                             pvalue_cutoff=0.01,
                             order=1,
                             ...){
  start_time <- Sys.time()
  data_tr <- data_transform(data,
                            ptime,
                            slide_window_size = slide_window_size,
                            slide_step_size = slide_step_size)
  nodes_info_df <- data.frame()
  n <- nrow(data_tr)
  genes <- rownames(data_tr)
  k <- 0
  for(i in c(1:(n-1))){
    for(j in c((i+1):n)){
      x <- as.numeric(data_tr[i,])
      y <- as.numeric(data_tr[j,])
      te_result <- grangertest(x,
                              y,
                              order=order,
                              ...)
      gene_i <- genes[i]
      gene_j <- genes[j]
      p_1 <- te_result$`Pr(>F)`[2]
      
      te_result2 <- grangertest(y,
                               x,
                               order=order,
                               ...)
      p_2 <- te_result2$`Pr(>F)`[2]
      if(p_1 < pvalue_cutoff){
        nodes_info <- c(gene_i,gene_j,p_1)
        k <- k + 1
        if(k == 1){
          nodes_info_df <- nodes_info
        }else{
          nodes_info_df <- rbind(nodes_info_df, nodes_info)
        }
      }else if(p_2 < pvalue_cutoff){
        nodes_info <- c(gene_j, gene_i, p_2)
        k <- k + 1
        if(k == 1){
          nodes_info_df <- nodes_info
        }else{
          nodes_info_df <- rbind(nodes_info_df, nodes_info)
        }
      }
    }
  }
  colnames(nodes_info_df) <- c("from","to","pvalue")
  
  nodes_info_df2 <- data.frame(from=nodes_info_df[,"from"],
                               to=nodes_info_df[,"to"],
                               pvalue=as.numeric(nodes_info_df[,"pvalue"]))
  end_time <- Sys.time()
  message("Time relapsed for granger test: ", (end_time -start_time), " seconds!")
  return(nodes_info_df2)
}