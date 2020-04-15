#' Create scatter plot
#' This function creates a scatter plot of all peptides that belong to a single proteins 
#' 
#' @param .listDF list, list containing data frames of peptides for each protein
#'       indexed by `Protein Group Accessions`
#' @param protein character, `Protein Group Accession` to show in the plot
#' @param max_frac numeric, total number of fractions
#' @param meanLine logical, specifies whether to plot a mean line 
#' @param repPepLine logical, specifies whether to plot a representative peptide line
#' @param separateLabStates logical, specifies whether label states will be separated into facets
#' @param titleLabel character, what to call the plot
#' @param titleAlign character, one of the 'left', 'center'/'centre', 'right', specifies alignment of the title in plot
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param grid logical, specifies presence/absence of gridline in the plot
#' @param labelled character, label to be used for isLabel == TRUE 
#' @param unlabelled character, label to be used for isLabel == FALSE
#' 
#' @importFrom tibble rowid_to_column
#'
#' @return plot
#' @export
allPeptidesPlot <- function(.listDF, protein, max_frac, 
                            meanLine = FALSE, repPepLine = FALSE, 
                            separateLabStates = FALSE,
                            grid = TRUE, 
                            titleLabel = 'all', titleAlign = 'left',
                            ylabel = 'Precursor Area', xlabel = 'Fraction',
                            legendLabel = 'Condition',
                            labelled = "Labeled", unlabelled = "Unlabeled"){
    
    dataFrame <- .listDF[[protein]]
    #Next 18 lines: edit data frame to neccessary format before plotting
    description <- dataFrame$`Protein Descriptions`[1]
    
    #add columns conatining mean and repPepValue to the DF, keep only neccessary rows
    data.frame(Value = NA,
               Fraction = 1:max_frac) %>% 
        spread(Fraction, "Value") -> padding
    
    dataFrame %>% 
        select(Fraction, `Precursor Area`, isLabel, repPepA) %>%
        rowid_to_column() %>% 
        spread(Fraction, `Precursor Area`) %>%  
        merge(padding, all.x = T)%>%
        gather(Fraction, `Precursor Area`, -c(isLabel,repPepA,rowid))%>% 
        group_by(Fraction, isLabel) %>% 
        mutate (meanValue = mean(`Precursor Area`, na.rm = T)) %>% 
        mutate (repPepValue = ifelse(all(is.na(`Precursor Area`)),
                                     NA,
                                     mean(`Precursor Area`[repPepA], na.rm=T))) %>% 
        ungroup() -> dataFrame
    
    dataFrame$Fraction <- as.numeric(as.character(dataFrame$Fraction))
    
    #transparent dots in case of drawing a line
    if (meanLine | repPepLine){
        alphaValue <- 0.20
    } else {
        alphaValue <- 1
    }
    
    #create a basic plot
    p <- ggplot(dataFrame,aes(Fraction, `Precursor Area`, colour = isLabel)) +
        geom_point(na.rm = T, alpha = alphaValue) +
        scale_y_log10() +
        scale_x_continuous(breaks = 1:max_frac, limits = c(0,max_frac))+
        scale_colour_manual(legendLabel, values = col_vector_peptides,
                            labels = c("TRUE" = labelled,
                                       "FALSE" = unlabelled)) +
        ylab (ylabel) +
        xlab (xlabel)
    
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
    
    #add title to plot according to arguments
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
                      # subtitle = '',
                      caption = paste('UniProt ID:', protein, sep ='')) +
            theme(plot.title = element_text(hjust = adjust))
    }
    
    #define linetype_vector 
    linetype_vector <- c('twodash', 'solid')
    names(linetype_vector) <- c('mean','representative peptide')
    
    #add mean line 
    if(meanLine) {  ## add line that is a mean of all peptide values
        p <- p + geom_line(aes(y=meanValue, colour = isLabel,linetype = 'mean'), size = 1, na.rm = T) +
            scale_linetype_manual('Line type', values = linetype_vector)
        
    }
    
    #add rep pep line
    ## add line that goes throug values of chosen representative peptide
    if (repPepLine) { 
        p <- p + geom_line(
            aes(y = repPepValue, colour = isLabel,
                linetype = 'representative peptide'),
            size=1, na.rm = T)+
            geom_point(
                aes(y = repPepValue, colour = isLabel),na.rm = T, alpha = 1)+
            scale_linetype_manual('Line type', values = linetype_vector)
    }
    
    #edit title in case of a multiprotein group
    #this hides information but there doesn't seems to be an easy way 
    #to deal with these cases in a nice way at this point
    if(str_detect(description,'\\|')){
        p <- p + labs(title = protein,
                      subtitle = 'Multiple proteins group')+
            theme(plot.title = element_text(hjust = adjust))
    }
    
    #make facets if one wishes to separate label states
    if(separateLabStates){
        p <- p + 
            facet_wrap(isLabel ~ ., nrow =2,scales = "free_x")+
            theme(strip.text = element_blank())
    }
    
    return (p)
}