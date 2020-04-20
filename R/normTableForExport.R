#' Convert Normalized Dataframe to Export format
#'
#' This is a convenient function for plotting
#'
#' @param labTab a dataframe
#' @param unlabTab a dataframe
#' @param comboTab a dataframe
#'
#' @importFrom stringr str_remove
#' @importFrom tidyr gather separate
#'
#' @return a dataframe
normTableForExport <- function(labTab, unlabTab, comboTab) {

  # Clean up labeled
  names(labTab) <- str_remove(names(labTab), "_.*$")
  max_frac <- max(suppressWarnings(as.numeric(names(labTab))), na.rm = TRUE)
  labTab$scenario <- "A"
  labTab$label <- TRUE
  labTab <- labTab[c("Protein Group Accessions", "Protein Descriptions", "scenario", "label", 1:max_frac)]

  # Clean up unlabeled
  names(unlabTab) <- str_remove(names(unlabTab), "_.*$")
  unlabTab$scenario <- "A"
  unlabTab$label <- FALSE
  unlabTab <- unlabTab[c("Protein Group Accessions", "Protein Descriptions", "scenario", "label", 1:max_frac)]

  # Clean up combined
  comboTab %>%
    gather("key", "value", -c(`Protein Group Accessions`, `Protein Descriptions`)) %>%
    separate("key", c("Fraction", "label"), "_") %>%
    spread(Fraction, "value") %>%
    mutate(scenario = "B",
           label = as.logical(label)) -> comboTab
  comboTab <- comboTab[c("Protein Group Accessions", "Protein Descriptions", "scenario", "label", 1:max_frac)]

  return(labTab %>%
           bind_rows(unlabTab) %>%
           bind_rows(comboTab))

}
