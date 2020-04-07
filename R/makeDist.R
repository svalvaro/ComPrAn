#' Make disstance matrix
#'
#' This function calculates distance matrix for a data frame, column by column
#' requires uncenteredCor function to work
#'
#' @param df data frame, contains columns:
#'              `Protein Group Accessions` character
#'              `Protein Descriptions` character
#'               isLabel character ('TRUE'/'FALSE')
#'               columns 1 to n, numeric, n is the total number of fractions/slices, each of this columns
#'                           contains `Precursor Area` values in a given fraction(columns) for a protein(rows)
#' @param centered centered: logical,if TRUE return dist matrix based on centered Pearson correlation (uses R cor() function, fast)
#'                 ,if FALSE return dist matrix based on uncentered Pearson correlation (uses custom uncenteredCor() function, slow)
#'
#' @importFrom stats as.dist cor
#'
#' @return matrix
#' @export
makeDist <- function(df,centered = FALSE){
  if (centered){
    answer <- (1 - cor(df))
  } else {
    nc = ncol(df)
    answer = matrix(ncol=nc,nrow=nc)
    for (i in 1:nc){
      for (j in 1:nc){
        answer[i,j] = 1-(uncenteredCor(df[,i],df[,j]))
      }
    }
    rownames(answer) <- colnames(df)
    colnames(answer) <- colnames(df)
  }
  return(as.dist(answer))
}

