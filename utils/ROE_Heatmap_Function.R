library(cowplot)
library(RColorBrewer)
library(reshape2)

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        if(m2[i,j] == 0){
          div.result[i,j] <- 0
        }else{
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}
ROIE <- function(crosstab, filter = NULL){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  if(is.null(filter)){filter = 10}
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  rowsum.matrix[rowsum.matrix <= filter] <- 0
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  colsum.matrix[colsum.matrix <= filter] <- 0
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  roie <- roie[rowSums(roie)!=0,]
  return(roie)
}

# ###绘图
# library(cowplot)
# library(RColorBrewer)
# library(reshape2)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:7]) # brewer.pal(9, "YlOrRd")[1:6]
