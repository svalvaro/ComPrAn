#' Convert extractRepPeps output to a Matrix
#'
#' Convert the dataframe as output from extractRepPeps() to matrix-like table
#' return normalized or raw values of Precursor Area, by default return normalized values
#'
#' @param .data a dataframe
#' @param applyNormalization logical apply normalization or not
#'
#' @return a matrix
#' @export
normalizeTable<- function(.data, applyNormalization = T){
  #convert data frame from extractRepPeps function to matrix-like table
  #return normalized or raw values of `Precursor Area`, by default return normalized values
  # .data <- protTableCombined
  #.data <- protTableLabelled

  labelSuffix <- unique(as.character(.data$isLabel))

  .data$LabelFraction <- paste(.data$Fraction, .data$isLabel, sep = '_')

  .data %>%
    select(-c(isLabel,Fraction)) %>%
    spread(LabelFraction,`Precursor Area`) -> .data

  #remove NA with 0 to be able to find maximum
  .data[is.na(.data)] <- 0

  #add column containing maximum value for each row
  # .data[, "max"] <- apply(.data[, 3:length(.data)], 1, max)

  #put back NAs
  #.data[.data==0] <- NA

  #convert to matrix
  # tempMatrix <- data.matrix(.data[3:(length(.data)-1)], rownames.force = T)
  # tempMatrix <- as.matrix(.data[3:(length(.data)-1)])
  tempMatrix <- as.matrix(.data[-c(1,2)])
  # dimnames(tempMatrix)
  # rownames(tempMatrix) <- .data$`Protein Group Accessions`
  # rownames(tempMatrix)
  # maxValues <- .data$max

  #divide each value by max if we want normalized values
  if (applyNormalization){
    # tempMatrix <- tempMatrix/maxValues

    tempMatrix <- t(apply(tempMatrix, 1, function(x) {x/max(x)}))
  }

  #reorder columns so they are from 1 to n
  if (length(labelSuffix) == 1){
    tempMatrix <- tempMatrix[,paste0(seq(1:ncol(tempMatrix)),"_",labelSuffix)]
  } else {
    tempMatrix <- tempMatrix[,paste0(1:(ncol(tempMatrix)/2),"_",rep(c(labelSuffix[1],labelSuffix[2]),each=ncol(tempMatrix)/2))]
  }


  #convert back to data frame
  tempMatrix <- as.data.frame(tempMatrix)

  #for unknown reasons ';' in rownames is sometimes automaticaly switched to '.', convert it back in such cases
  # rownames(tempMatrix) <- gsub("\\.", ";", rownames(tempMatrix))
  .data <- cbind(.data[,1:2], tempMatrix)

  #tempMatrix$`Protein Descriptions` <- .data$`Protein Descriptions`
  #.data <- tempMatrix

  return(.data)
}
