#' Compare a Single Group of Proteins Between Two Label States
#'
#' This function creates a ?scatter plot? for a subset of proteins in dataFrame
#' specified in groupData. Intended use of the function - using scenario A data,
#' compare shape of the migration profile for a SINGLE GROUP of proteins BETWEEN
#' the two LABEL STATES.
#'
#' @param dataFrame dataFrame: data frame, data frame of normalised values for
#'  proteins from SCENARIO A,
#'            contains columns:
#'           `Protein Group Accessions` character
#'           `Protein Descriptions` character
#'            Fraction integer
#'            isLabel character ('TRUE'/'FALSE' values)
#'            `Precursor Area` double
#'            scenario character
#' @param max_frac numeric, total number of fractions
#' @param groupData character vector, contins list of Protein Group Accessions 
#' that belong to the group we want to plot
#' @param groupName character, name that should be used for the group specified 
#' in groupData
#' @param meanLine logical, specifies whether to plot a mean line for all values
#'  in the group
#' @param medianLine logical, specifies whether to plot a median line for all 
#' values in the group
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param labelled character, label to be used for isLabel == TRUE
#' @param unlabelled character, label to be used for isLabel == FALSE
#' @param jitterPoints numeric
#' @param pointSize numeric, size of the point in the plot
#' @param grid logical, specifies presence/absence of gridline in the plot
#' @param titleAlign character, one of the 'left', 'center'/'centre', 'right',
#'  specifies alignment of the title in plot
#' @param alphaValue numeric, transparency of the point, values 0 to 1
#' @param controlSample character, either labelled or unlabelled, this setting
#' will adjust plot coloring based on which sample is a control
#'
#' @importFrom stats median
#'
#' @return plot
#' @export
#' 
#' @examples 
#' 
#' ##Use example normalised proteins file
#' inputFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
#' #read file in and change structure of table to required format
#' forAnalysis <- protImportForAnalysis(inputFile)
#' ##example plot:
#' groupDV <- c("Q16540","P52815","P09001","Q13405","Q9H2W6")
#' groupName <- 'group1' 
#' max_frac <- 23 
#' oneGroupTwoLabelsCoMigration(forAnalysis, max_frac, groupDV,groupName)
#' 
oneGroupTwoLabelsCoMigration <- function(
    dataFrame, max_frac, groupData = NULL, groupName = 'group1', 
    meanLine = FALSE,medianLine = FALSE,ylabel = 'Relative Protein Abundance',
    xlabel = 'Fraction',legendLabel = 'Condition', labelled = 'Labeled',
    unlabelled = 'Unlabeled',jitterPoints = 0.3, pointSize = 2.5,
    grid = FALSE,titleAlign = 'left', alphaValue = 0.5,controlSample = ""){
    if (controlSample == "labelled"|controlSample == "labeled"){
        col_vector_proteins <- c("TRUE" = "#ff9d2e", "FALSE" = "#07b58a")
    }else if(controlSample == "unlabelled"|controlSample == "unlabeled"){
        col_vector_proteins <- c("FALSE" = "#ff9d2e", "TRUE" = "#07b58a")}
    if(is.null(groupData)) {
        stop('Please provide a list of proteins you would like to plot')}
    dataFrame <- dataFrame[dataFrame$scenario == "A",]#filter only scenario A
    dataFrame %>%
        select(-scenario) %>%
        filter(`Protein Group Accessions` %in% groupData) %>%
        filter(!is.na(`Precursor Area`)) %>% group_by(Fraction, isLabel) %>%
        mutate (meanValue = mean(`Precursor Area`, na.rm = TRUE)) %>%
        mutate (medianValue = median(`Precursor Area`, na.rm = TRUE)) %>%
        ungroup() -> dataFrame
    linetype_vector <- c('twodash', 'solid')     #define linetype_vector
    names(linetype_vector) <- c('mean','median')
    p <- ggplot(dataFrame, aes(x = Fraction, y = `Precursor Area`, 
                                col = isLabel)) +
        geom_point(position = position_jitter(jitterPoints),alpha = alphaValue,
                    size = pointSize) +
        scale_color_manual(legendLabel, values = col_vector_proteins,
                            labels = c("TRUE" = labelled,"FALSE" = unlabelled))+
        scale_fill_manual(legendLabel, values = col_vector_proteins,
                            labels = c("TRUE" =labelled, "FALSE" =unlabelled))+
        ylab(ylabel) + xlab(xlabel) +
        scale_x_continuous(breaks=seq_len(max_frac),minor_breaks = NULL)+
        scale_y_continuous(breaks=seq(0,1,0.2))+ labs(title = groupName)+
        scale_linetype_manual('Line type', values = linetype_vector)
    if(meanLine) {  ## add line that is a mean of all protein values
        p <- p + geom_line(aes(y=meanValue, colour = isLabel, linetype = 'mean')
                            ,size = 1, na.rm = TRUE)}
    if (medianLine) { ##  add line that is a median of all protein values
        p <- p + geom_line(
            aes(y = medianValue, colour = isLabel, linetype = 'median'),
            size=1, na.rm = TRUE)}
    if(grid){p<- p +theme_minimal() +     #add grid
            theme(panel.grid.minor = element_blank())
    } else {p<- p +theme_classic()}
    if (titleAlign == 'left'){adjust <- 0     #title alignment settings
    } else if ((titleAlign == 'centre')|(titleAlign=='center')) {adjust <- 0.5
    } else if(titleAlign == 'right'){adjust <- 1}
    p <- p + theme(plot.title = element_text(hjust = adjust))
    return(p)
}