#' Create Line Plots
#'
#' This function creates a line plot for a proteins in dataFrame specified
#' by protein
#'
#' @param dataFrame data frame, contains columns: `Protein Group Accessions` 
#' character; `Protein Descriptions` character;bFraction integer;
#' isLabel character ("TRUE"/"FALSE" values);`Precursor Area` double;
#' scenario character
#' @param protein character the protein of interest
#' @param max_frac integer total number of fractions
#' @param grid logical specifies presence/absence of gridline in the plot
#' @param titleLabel character, if it is 'all' or 'GN', it specifies whether to 
#' show full label or just the gene name, if any other character is used, the
#' value of titleLabel will be used as plot title
#' @param titleAlign character one of the 'left', 'center'/'centre', 'right',
#'  specifies alignment of the title in plot
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param labelled character label to be used for isLabel == TRUE
#' @param unlabelled character label to be used for isLabel == FALSE
#' @param controlSample character, either labelled or unlabelled, this setting
#' will adjust plot coloring based on which sample is a control
#' @param textSize numeric, size of text in the plot 
#' @param axisTextSize numeric, size of axis labels in the plot
#' 
#' @importFrom stringr str_extract
#'
#' @return a plot
#' @export
#' 
#' @examples
#' 
#' ##Use example normalised proteins file
#' inputFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
#' #read file in and change structure of table to required format
#' forAnalysis <- protImportForAnalysis(inputFile)
#' ##example plot:
#' protein <- "P52815"
#' max_frac <- 23
#' proteinPlot(forAnalysis, protein, max_frac)
#' 
proteinPlot <- function(dataFrame, protein, max_frac, grid = TRUE, 
                        titleLabel = 'all', titleAlign = 'left',
                        ylabel = 'Relative Protein Abundance', 
                        xlabel = 'Fraction',legendLabel = 'Condition',
                        labelled = "Labeled", unlabelled = "Unlabeled",
                        controlSample = "",textSize = 12, axisTextSize = 8){
    if (controlSample == "labelled"|controlSample == "labeled"){
        col_vector_proteins <- c("TRUE" = "#ff9d2e", "FALSE" = "#07b58a")
    }else if(controlSample == "unlabelled"|controlSample == "unlabeled"){
        col_vector_proteins <- c("FALSE" = "#ff9d2e", "TRUE" = "#07b58a")}
    dataFrame %>%
        filter(scenario == "B",`Protein Group Accessions`== protein)-> dataFrame
    description <- dataFrame$`Protein Descriptions`[1]
    p <- ggplot(dataFrame,aes(Fraction, `Precursor Area`, colour = isLabel)) +
        geom_line(na.rm = TRUE) +
        geom_point(na.rm = TRUE) +
        scale_x_continuous(breaks = seq_len(max_frac), limits = c(0,max_frac))+
        scale_color_manual( legendLabel, values = col_vector_proteins,
                            labels = c("TRUE" = labelled,"FALSE" = unlabelled))+
        ylab(ylabel) + xlab(xlabel)
    if(grid){p<- p +theme_minimal() +  #add grid
        theme(panel.grid.minor = element_blank())
    } else {p<- p +theme_classic()}
    if (titleAlign == 'left'){adjust <- 0    #title alignment settings
    } else if ((titleAlign == 'centre')|(titleAlign=='center')) {adjust <- 0.5
    } else if(titleAlign == 'right'){adjust <- 1}
    #add title to plot according to arguments
    description <- dataFrame$`Protein Descriptions`[1]
    if (titleLabel == 'all'){
        p <- p + labs(title = str_remove(
            str_extract(description, "^.* OS"), " OS"),
            subtitle = str_extract(
                str_remove(description, " PE=.*$"), "OS=.*$"),
            caption = paste('UniProt ID:', protein, sep ='')) +
            theme(plot.title = element_text(hjust = adjust))
    } else if (titleLabel == 'GN') {
        p <- p + labs(title = str_remove(
            str_extract(description, "GN=[:alnum:]*"), "GN="),
            caption = paste('UniProt ID:', protein, sep ='')) +
            theme(plot.title = element_text(hjust = adjust))
    } else {p <- p + labs(  title = titleLabel,
                            caption = paste('UniProt ID:', protein, sep ='')) +
        theme(plot.title = element_text(hjust = adjust))}
    if(str_detect(description,'\\|')){
        p <- p + labs(title = protein, subtitle = 'Multiple proteins group')  }
    p <- p + theme(text = element_text(size = textSize),
                   axis.text=element_text(size = axisTextSize))
    return(p)
}