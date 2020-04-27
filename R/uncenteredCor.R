#' Perform uncentered correlation
#'
#' @param xx numeric vector
#' @param yy numeric vector
#'
#' @return vector
uncenteredCor <- function (xx,yy) {
    #Calculate uncentered Pearson correlation
    sd.x1 <- sqrt(sum((xx-0)^2)/(length(xx)-1))
    sd.y1 <- sqrt(sum((yy-0)^2)/(length(yy)-1))
    
    corVal <- sum(((xx - 0)/sd.x1)*((yy - 0)/sd.y1))/(length(xx)-1)
    return(corVal)
}