## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----loading data-------------------------------------------------------------
# loading pGRN package
library(pGRN)

# loading pre-build data
example_data <- pGRNDB
names(example_data)
expression_matrix <- example_data[["expression"]]
pseudotime_list <- example_data[["ptime"]]$PseudoTime

## ----call pGRN----------------------------------------------------------------
# try DTW method
nets_dtw <- pGRN(expression_matrix,pseudotime_list, method= "DTW",quantile_cutoff=50)

# try granger method
nets_gg <- pGRN(expression_matrix,pseudotime_list, method= "granger")


## ----plot network-------------------------------------------------------------
# plot network stationarily
plot_network(nets_dtw[[1]])
plot_network(nets_gg[[1]])

# plot network interactively
#plot_network_i(nets_dtw[[1]])
#plot_network_i(nets_gg[[1]])

