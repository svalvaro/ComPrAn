#' Convert extractRepPeps output to a Matrix
#'
#' Convert the dataframe as output from extractRepPeps() to matrix-like table
#' return normalized or raw values of Precursor Area, by default return
#' normalized values
#'
#' @param .data a dataframe
#' @param applyNormalization logical apply normalization or not
#'
#' @return a matrix
normalizeTable<- function(.data, applyNormalization = TRUE){
    labelSuffix <- unique(as.character(.data$isLabel))
    
    .data$LabelFraction <- paste(.data$Fraction, .data$isLabel, sep = '_')
    
    .data %>%
        select(-c(isLabel,Fraction)) %>%
        spread(LabelFraction,`Precursor Area`) -> .data
    
    #remove NA with 0 to be able to find maximum
    .data[is.na(.data)] <- 0
    
    #convert to matrix
    tempMatrix <- as.matrix(.data[-c(1,2)])

    if (applyNormalization){
        tempMatrix <- t(apply(tempMatrix, 1, function(x) {x/max(x)}))
    }
    
    #reorder columns so they are from 1 to n#->no reason to reorder columns here
    # if (length(labelSuffix) == 1){
    #     tempMatrix <- tempMatrix[,paste0(seq_len(ncol(tempMatrix)),
    #                                         "_",labelSuffix)]
    # } else {
    #     tempMatrix <- tempMatrix[,paste0(seq_len(ncol(tempMatrix)/2),
    #                                         "_",
    #                                     rep(c(labelSuffix[1],labelSuffix[2]),
    #                                         each=ncol(tempMatrix)/2))]
    # }
    #convert back to data frame
    tempMatrix <- as.data.frame(tempMatrix)
    #for unknown reasons ';' in rownames is sometimes automaticaly switched 
    #to '.', convert it back in such cases
    # rownames(tempMatrix) <- gsub("\\.", ";", rownames(tempMatrix))
    .data <- cbind(.data[,c(1,2)], tempMatrix)
    
    return(.data)
}