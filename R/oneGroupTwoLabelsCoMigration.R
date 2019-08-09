#' Compare a Single Group of Proteins Between Two Label States
#'
#' This function creates a ?scatter plot? for a subset of proteins in dataFrame specified in groupData
#' Intended use of the function - using scenario A data, compare shape of the migration profile
#' for a SINGLE GROUP of proteins BETWEEN the two LABEL STATES
#'
#' @param dataFrame dataFrame: data frame, data frame of normalised values for proteins from SCENARIO A,
#'            contains columns:
#'           `Protein Group Accessions` character
#'           `Protein Descriptions` character
#'            Fraction integer
#'            isLabel character ('TRUE'/'FALSE' values)
#'            `Precursor Area` double
#'            scenario character
#' @param max_frac numeric, total number of fractions
#' @param groupData character vector, contins list of Protein Group Accessions that belong to the group we want to plot
#' @param groupName character, name that should be used for the group specified in groupData
#' @param meanLine logical, specifies whether to plot a mean line for all values in the group
#' @param medianLine logical, specifies whether to plot a median line for all values in the group
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param labelled character, label to be used for isLabel == TRUE
#' @param unlabelled character, label to be used for isLabel == FALSE
#' @param jitterPoints numeric
#' @param pointSize numeric, size of the point in the plot
#' @param grid logical, specifies presence/absence of gridline in the plot
#' @param titleAlign character, one of the 'left', 'center'/'centre', 'right', specifies alignment of the title in plot
#' @param alphaValue numeric, transparency of the point, values 0 to 1
#'
#' @return plot
#' @export
oneGroupTwoLabelsCoMigration <- function(dataFrame, max_frac, groupData = NULL, groupName = 'group1', meanLine = FALSE, medianLine = FALSE,
                                         ylabel = 'Relative Protein Abundance', xlabel = 'Fraction',
                                         legendLabel = 'Condition', labelled = 'Labeled', unlabelled = 'Unlabeled',
                                         jitterPoints = 0.3, pointSize = 2.5, grid = FALSE,
                                         titleAlign = 'left', alphaValue = 0.5){

  if(is.null(groupData)) {
    stop('Please provide a list of proteins you would like to plot')
  }

  #filter only scenario A values:
  dataFrame <- dataFrame[dataFrame$scenario == "A",]
  dataFrame %>%
    select(-scenario) ->dataFrame


  dataFrame %>%
    filter(`Protein Group Accessions` %in% groupData) %>%
    filter(!is.na(`Precursor Area`))-> dataFrame

  dataFrame %>%
    group_by(Fraction, isLabel) %>%
    mutate (meanValue = mean(`Precursor Area`, na.rm = T)) %>%
    mutate (medianValue = median(`Precursor Area`, na.rm = T)) %>%
    ungroup() -> dataFrame

  p <- ggplot(dataFrame, aes(x = Fraction, y = `Precursor Area`, col = isLabel)) +
    geom_point(position = position_jitter(jitterPoints),alpha = alphaValue, size = pointSize) +
    scale_color_manual(legendLabel, values = col_vector_proteins,
                       labels = c("TRUE" = labelled,
                                  "FALSE" = unlabelled)) +
    scale_fill_manual(legendLabel, values = col_vector_proteins,
                      labels = c("TRUE" = labelled,
                                 "FALSE" = unlabelled)) +
    ylab(ylabel) +
    xlab(xlabel) +
    scale_x_continuous(breaks=1:max_frac,minor_breaks = NULL)+
    scale_y_continuous(breaks=seq(0,1,0.2))+
    labs(title = groupName)

  #define linetype_vector
  linetype_vector <- c('twodash', 'solid')
  names(linetype_vector) <- c('mean','median')

  #add mean line
  if(meanLine) {  ## add line that is a mean of all protein values
    p <- p + geom_line(aes(y=meanValue, colour = isLabel, linetype = 'mean'), size = 1, na.rm = T) +
      scale_linetype_manual('Line type', values = linetype_vector)
  }

  #add median line
  if (medianLine) { ##  add line that is a median of all protein values
    p <- p + geom_line(aes(y = medianValue, colour = isLabel, linetype = 'median'),size=1, na.rm = T)+
      scale_linetype_manual('Line type', values = linetype_vector)
  }

  #add grid
  if(grid){
    p<- p +theme_minimal() +
      theme(panel.grid.minor = element_blank())
  } else {
    p<- p +theme_classic()
  }

  #title alignment settings
  if (titleAlign == 'left'){
    adjust <- 0
  } else if ((titleAlign == 'centre')|(titleAlign=='center')) {
    adjust <- 0.5
  } else if(titleAlign == 'right'){
    adjust <- 1
  }

  #adjust position of title
  p <- p + theme(plot.title = element_text(hjust = adjust))

  return(p)
}
