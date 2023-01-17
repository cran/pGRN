#' @import dtw
#' @importFrom proxy pr_DB
#' @import doParallel
#' @import foreach
#' @import bigmemory
#' @importFrom parallel detectCores makePSOCKcluster stopCluster
#'
NULL

#' Sliding Window Average
#'
#' Get sliding windows average values for given vector/list
#'
#' @param data list of expression
#' @param window sliding window size
#' @param step sliding window step size
#'
#' @return list/vector of sliding windows with average expression value
#'
#' @export
#'
#' @examples
#' slideWindows(c(1:1000),window=200,step=100)
#' slideWindows(c(1:1000),window=100,step=50)
#'
slideWindows <- function(data,
                         window=2,
                         step=1){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window-1)])
  }
  return(result)
}

#' Pseudotime based Expression Data Transformation
#'
#' Based on single-cell pseudotime information, get the sliding window average expression,
#' and then standard normlize the expression for each gene
#'
#' @param data expression matrix data
#' @param pseudotime list of pseudotime
#' @param slide_window_size sliding window size
#' @param slide_step_size sliding window step size
#'
#' @return Transformed new matrix
#'
#' @export
#'
#' @examples
#' data <- matrix(1,100,1000)
#' ptime <- seq(1:1000)
#' data_transform(data,
#'                ptime,
#'                slide_window_size=100,
#'                slide_step_size=50)
#'
data_transform <- function(data,
                           pseudotime,
                           slide_window_size=100,
                           slide_step_size=50){

  data_ordered <- data[,order(pseudotime,decreasing = F)] # warning message need to be corrected

  cell_num <- ncol(data_ordered)

  data_new <- apply(data_ordered,1,function(x){
    exp <- as.numeric(x)
    exp_sw <- slideWindows(exp,slide_window_size,slide_step_size)
    exp_sd <- sd(exp_sw)
    if(exp_sd == 0){
      exp_n <- (exp_sw-mean(exp_sw))/1
    }else(
      exp_n <- (exp_sw-mean(exp_sw))/exp_sd
    )

  })
  return(t(data_new))
}

#' Bidirectional DTW Distance
#'
#' Get bidirectional DTW distance.
#'
#' @param x list of x input
#' @param y list of y input
#'
#' @export
#'
#' @return numeric
#'
#' @examples
#' get_dtw_dist_bidirectional(c(1:1000),c(1:1000))
#'
get_dtw_dist_bidirectional <- function(x,y){
  dist1 <- dtw(x,y)$normalizedDistance
  dist2 <- dtw(x,-y)$normalizedDistance
  if(dist1 <= dist2){
    return(dist1)
  }else{
    return(-dist2)
  }
}


#' DTW distance matrix for all genes
#'
#' Get DTW distance matrix for all genes using pseudotime based sliding window transfromation,
#' parallel computing allowed.
#'
#' @param data gene expression matrix (Gene * Cells)
#' @param ptime pseudotime matched with the column cells of the gene expression matrix
#' @param slide_window_size sliding window size
#' @param slide_step_size sliding window step size
#' @param cores number of cores for parallel computing
#'
#' @export
#' 
#' @return bidirectional DTW distance matrix
#'
#' @examples
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix,
#'                                     pseudotime_list,
#'                                     cores=1)
#'
get_dtw_dist_mat <- function(data,
                                ptime,
                                slide_window_size = 50,
                                slide_step_size = 25,
                                cores=2){
  
  # fix for foreach dopar in windows, within foreach functions outside this function will not be found
  get_dtw_dist_bidirectional <- function(x,y){
    dist1 <- dtw::dtw(x,y)$normalizedDistance
    dist2 <- dtw::dtw(x,-y)$normalizedDistance
    if(dist1 <= dist2){
      return(dist1)
    }else{
      return(-dist2)
    }
  }
  
  data_tr <- data_transform(data,
                            ptime,
                            slide_window_size = slide_window_size,
                            slide_step_size = slide_step_size)
  
  # init distance matrix with NA default value
  n <- nrow(data_tr)
  dtw_dist_mat <- matrix(data=NA,nrow=n,ncol=n)
  
  # detect number of available cores
  numCores <- detectCores(logical = TRUE)
  if(cores > numCores){
    cores <- numCores
  }
  
  cl <- makePSOCKcluster(cores)
  
  registerDoParallel(cl)
  
  # get dtw distance for each gene pair
  i <- NULL
  results <- foreach (i=c(1:n), .combine="c") %dopar% {
    ds <- c()
    for(j in c(i:n)){
      d <- get_dtw_dist_bidirectional(data_tr[i,],data_tr[j,])
      ds <- c(ds, d)
    }
    ds
  }
  stopCluster(cl)
  
  # reformat to matrix
  idx <- 0
  for(i in c(1:n)){
    for(j in c(i:n)){
      idx <- idx + 1
      d <- results[idx]
      dtw_dist_mat[i,j] <- d
      dtw_dist_mat[j,i] <- -d
    }
  }
  
  # add row and column names
  rownames(dtw_dist_mat) <- rownames(data_tr)
  colnames(dtw_dist_mat) <- rownames(data_tr)
  
  return(dtw_dist_mat)
}



#' Get network adjacency dataframe based on DTW method
#'
#' Use DTW to calcuate gene-gene distance based on their expression and pseudotime
#'
#' @param expression_matrix expression matrix data
#' @param pseudotime_list list of pseudotime
#' @param slide_window_size	 sliding window size
#' @param slide_step_size sliding window step size
#' @param quantile_cutoff	 an integer value (1-99) for quantile cutoff
#'
#' @param cores number of cores for parallel computing
#' 
#' @return adjacency dataframe (with columns "from, to, distance,direction, similarity")
#'
#' @export
#'  
#' @examples 
#' example_data <- pGRNDB
#' expression_matrix <- example_data[["expression"]]
#' pseudotime_list <- example_data[["ptime"]]$PseudoTime
#' adj_df <- run_dtw(expression_matrix,
#'                   pseudotime_list,
#'                   quantile_cutoff=50,
#'                   cores=1)
#'                   
run_dtw <- function(expression_matrix,
                    pseudotime_list,
                    slide_window_size = 50,
                    slide_step_size = 25,
                    quantile_cutoff = 5,
                    cores=1){
  dtw_dist_matrix <- get_dtw_dist_mat(expression_matrix, 
                                      pseudotime_list,
                                      slide_window_size=slide_window_size,
                                      slide_step_size=slide_step_size,
                                      cores=cores) # time-consuming step
  adj_df <- matrix2adj(dtw_dist_matrix, quantile_cutoff=quantile_cutoff)
  return(adj_df)
}

