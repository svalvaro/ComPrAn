#' Create scatter plot
#' 
#' This function creates a plot of all peptides that belong to a single protein 
#' 
#' @param .listDF list, list containing data frames of peptides for each
#'  protein indexed by `Protein Group Accessions`
#' @param protein character, `Protein Group Accession` to show in the plot
#' @param max_frac numeric, total number of fractions
#' @param meanLine logical, specifies whether to plot a mean line 
#' @param repPepLine logical, specifies whether to plot a representative 
#' peptide line
#' @param separateLabStates logical, specifies whether label states will be 
#' separated into facets
#' @param titleLabel character, what to call the plot
#' @param titleAlign character, one of the 'left', 'center'/'centre', 'right', 
#' specifies alignment of the title in plot
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param grid logical, specifies presence/absence of gridline in the plot
#' @param labelled character, label to be used for isLabel == TRUE 
#' @param unlabelled character, label to be used for isLabel == FALSE
#' @param controlSample character, either labelled or unlabelled, this setting
#' will adjust coloring based on which sample is a control
#' @param textSize numeric, size of text in the plot 
#' @param axisTextSize numeric, size of axis labels in the plot
#' 
#' @importFrom tibble rowid_to_column
#'
#' @return plot
#' @export
#' 
#' @examples 
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extData", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
#' peptides <- cleanData(peptides, fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#' ## Pick representative peptide for each protein for both scenarios
#' peptide_index <- pickPeptide(peptides)
#' 
#' ##create a plot showing all peptides of selected protein
#' protein <- "P52815"
#' max_frac <- 23
#' #default plot
#' allPeptidesPlot(peptide_index,protein, max_frac = max_frac)
#' #other plot version
#' allPeptidesPlot(peptide_index,protein, max_frac = max_frac,
#' repPepLine = TRUE, grid = FALSE, titleAlign = "center")
#' #other plot version
#' allPeptidesPlot(peptide_index,protein, max_frac = max_frac,
#' repPepLine = TRUE, meanLine = TRUE, separateLabStates =TRUE,
#' titleLabel = "GN")
allPeptidesPlot <- function(.listDF, protein, max_frac, meanLine = FALSE, 
    repPepLine = FALSE,separateLabStates = FALSE,grid = TRUE,titleLabel = 'all',
    titleAlign = 'left',ylabel = 'Precursor Area', xlabel = 'Fraction',
    legendLabel = 'Condition', labelled = "Labeled", unlabelled = "Unlabeled",
    controlSample = "",textSize = 12, axisTextSize = 8){
    #if controlSample specified, adjust colouring so control is always same 
    if (controlSample == "labelled"|controlSample == "labeled"){
        col_vector_peptides <- c("TRUE" = "#ffc125", "FALSE" = "#a020f0")
    }else if(controlSample == "unlabelled"|controlSample == "unlabeled"){
        col_vector_peptides <- c("FALSE" = "#ffc125", "TRUE" = "#a020f0")
    }
    dataFrame <- .listDF[[protein]]
    #Next lines: edit data frame to neccessary format before plotting
    description <- dataFrame$`Protein Descriptions`[1]
    #add mean and repPepValue columns to DF, keep only necessary rows
    data.frame(Value = NA, Fraction = seq_len(max_frac)) %>% 
        spread(Fraction, "Value") -> padding
    dataFrame %>% 
        select(Fraction, `Precursor Area`, isLabel, repPepA) %>%
        rowid_to_column() %>% spread(Fraction, `Precursor Area`) %>%  
        merge(padding, all.x = TRUE)%>%
        gather(Fraction, `Precursor Area`, -c(isLabel,repPepA,rowid))%>% 
        group_by(Fraction, isLabel) %>% 
        mutate (meanValue = mean(`Precursor Area`, na.rm = TRUE)) %>% 
        mutate (repPepValue = ifelse(all(is.na(`Precursor Area`)),
                                        NA, mean(`Precursor Area`[repPepA], 
                                        na.rm=TRUE))) %>% 
        ungroup() -> dataFrame
    dataFrame$Fraction <- as.numeric(as.character(dataFrame$Fraction))
    if (meanLine | repPepLine){alphaValue <- 0.20    
    } else {alphaValue <- 1}  #transparent dots in case of drawing a line
    linetype_vector <- c('twodash', 'solid')     #define linetype_vector 
    names(linetype_vector) <- c('mean','representative peptide')
    p <- ggplot(dataFrame,aes(Fraction, `Precursor Area`, colour = isLabel)) +
        geom_point(na.rm = TRUE, alpha = alphaValue) +
        scale_y_log10() +
        scale_x_continuous(breaks = seq_len(max_frac), limits = c(0,max_frac))+
        scale_colour_manual(legendLabel, values = col_vector_peptides,
                            labels = c("TRUE" = labelled,
                                        "FALSE" = unlabelled)) +
        ylab (ylabel) + xlab (xlabel) +
        scale_linetype_manual('Line type', values = linetype_vector)
    if(grid){p<- p +theme_minimal() +     #add grid 
            theme(panel.grid.minor = element_blank())
    } else {p<- p +theme_classic()}
    if (titleAlign == 'left'){adjust <- 0     #title alignment settings
    } else if ((titleAlign == 'centre')|(titleAlign=='center')) {adjust <- 0.5
    } else if(titleAlign == 'right'){adjust <- 1}
    if (titleLabel == 'all'){
        p <- p+labs(title=str_remove(str_extract(description, "^.* OS")," OS"),
            subtitle = str_extract(str_remove(description," PE=.*$"), "OS=.*$"),
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
    if(meanLine) {  ## add line that is a mean of all peptide values
        p <- p + geom_line(aes(y=meanValue, colour = isLabel,linetype ='mean'),
                            size = 1, na.rm = TRUE)
    }
    if (repPepLine) { ## add rep pep line fir scenario A
        p <- p + geom_line(aes(y = repPepValue, colour = isLabel,
                linetype = 'representative peptide'),size=1, na.rm = TRUE)+
            geom_point(aes(y=repPepValue,colour=isLabel),na.rm=TRUE,alpha = 1)
    }
    #edit title in case of a multiprotein group
    if(str_detect(description,'\\|')){
        p <- p + labs(title = protein,
                        subtitle = 'Multiple proteins group')+
            theme(plot.title = element_text(hjust = adjust))
    }
    #make facets if one wishes to separate label states
    if(separateLabStates){
        p <- p + facet_wrap(isLabel ~ ., nrow =2,scales = "free_x")+
            theme(strip.text = element_blank())
    }
    p <- p + theme(text = element_text(size = textSize), #adjust text size
                    axis.text=element_text(size = axisTextSize))
    return (p)
}