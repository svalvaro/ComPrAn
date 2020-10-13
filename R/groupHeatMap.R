#' Make heatmap
#'
#' This function creates a heatmap for a subset of proteins in dataFrame
#' specified in groupData, heatmap is divided into facets according to isLabel
#'
#' @param dataFrame data frame, contains columns:
#'           `Protein Group Accessions` character
#'           `Protein Descriptions` character
#'            Fraction integer
#'            isLabel character ('TRUE'/'FALSE' values)
#'            `Precursor Area` double
#'            scenario character
#' @param groupData data frame, mandatory column: 
#' `Protein Group Accessions` character - this column is used for filtering
#'            optional columns: any other column of type character that should
#'             be used for renaming
#' @param groupName character, name that should be used for the group specified
#'  in groupData
#' @param titleAlign character, one of the 'left', 'center'/'centre', 'right', 
#' specifies alignment of the title in plot
#' @param newNamesCol character, if groupData contains column for re-naming and 
#' you want to use it, specify the column name in here
#' @param colNumber numeric, values of 1 or 2, specifies whether facets will be 
#' shown side-by-side or above each other
#' @param ylabel character
#' @param xlabel character
#' @param legendLabel character
#' @param legendPosition character, one of "right" or "bottom"
#' @param grid logical, specifies presence/absence of gridline in the plot
#' @param labelled character, label to be used for isLabel == TRUE
#' @param unlabelled character, label to be used for isLabel == FALSE
#' @param orderColumn character, if groupData contains column for re-ordering
#' and you want to use it, specify the column name in here
#'
#' @import ggplot2
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
#' groupDfn <- system.file("extdata", "exampleGroup.txt", package = "ComPrAn")
#' groupName <- 'group1'
#' groupData <- data.table::fread(groupDfn)
#' groupHeatMap(forAnalysis[forAnalysis$scenario == "B",], groupData, groupName)
#' 
groupHeatMap <- function(dataFrame, groupData, groupName, titleAlign = "left", 
    newNamesCol = NULL, colNumber = 2, ylabel = "Protein", xlabel = "Fraction",
    legendLabel = "Relative Protein Abundance", legendPosition = "right",
    grid = TRUE,labelled = "labeled",unlabelled="unlabeled",orderColumn = NULL){
    #join DF and group data - proteins present in group but absent in the data 
    #will be shown as empty
    groupData %>% 
        select(`Protein Group Accessions`, newNamesCol,orderColumn) -> groupData
    right_join(dataFrame, groupData) -> dataFrame
    if(sum(is.na(dataFrame$isLabel))>0){
        dataFrame[is.na(dataFrame$isLabel),]$isLabel <- FALSE}
    if(!is.null(newNamesCol)){ycolumn <- newNamesCol #rename proteins
    } else {ycolumn <- 'Protein Group Accessions'}
    if(is.null(orderColumn)){     #draw basic plot
        p <- ggplot(dataFrame, aes(x = Fraction,
                                y = get(ycolumn),fill = `Precursor Area`)) + 
        geom_tile(na.rm = TRUE)  +
        facet_wrap(isLabel ~ ., ncol = colNumber, 
            labeller = labeller(isLabel=c("TRUE"=labelled,"FALSE"=unlabelled)))+
        labs(title = groupName) + ylab(ylabel) + xlab(xlabel) +
        scale_fill_gradient(legendLabel,low='#cacde8',high ='#0019bf',
            na.value="grey60",breaks = seq(0,1,0.25), limits = c(0,1)) +
        coord_cartesian(expand = 0)
    }else{
        protsOrder <- dataFrame[,orderColumn]
        dataFrame[[ycolumn]]<-fct_reorder(dataFrame[[ycolumn]],desc(protsOrder))
        p <-  ggplot(dataFrame, aes(x = Fraction,y = get(ycolumn),
                                    fill = `Precursor Area`)) + 
        geom_tile(na.rm = TRUE)+
        facet_wrap(isLabel ~ .,ncol=colNumber, 
            labeller=labeller(isLabel=c("TRUE"=labelled,"FALSE"=unlabelled)))+
        labs(title = groupName) + ylab(ylabel) +  xlab(xlabel) +
        scale_fill_gradient(legendLabel, low = '#cacde8',high = '#0019bf', 
            na.value="grey60", breaks = seq(0,1,0.25), limits = c(0,1)) +
        coord_cartesian(expand = 0)}
    if(grid){p<- p +theme_minimal() +     #add grid
            theme(panel.grid.minor = element_blank())
    } else {p<- p +theme_classic()}
    if (titleAlign == 'left'){adjust <- 0     #title alignment settings
    } else if ((titleAlign == 'centre')|(titleAlign=='center')) {adjust <- 0.5
    } else if(titleAlign == 'right'){adjust <- 1}
    p <- p + theme(plot.title = element_text(hjust = adjust),
                    legend.position = legendPosition)
    if(legendPosition == "bottom"){
        p <- p+guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,
                    draw.ulim = FALSE, draw.llim = FALSE,barwidth = 12))
    }else{p<-p+guides(fill=guide_colourbar(title.position="top",title.hjust=0.5,
                draw.ulim = FALSE, draw.llim = FALSE,barwidth = 2))}
    return(p)
}