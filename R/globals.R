utils::globalVariables(
    c("Fraction", "PSM Ambiguity", "Confidence Level",
        "Protein Group Accessions", "# Protein Groups", "Precursor Area", 
        "Modifications", "Sequence", "Protein Descriptions", "isLabel",
        "cluster", "peptide_index", "LabelFraction", "label", 
        "col_vector_proteins", 
        "meanValue", "median", "medianValue", "scenario", "UniqueCombinedID_A", 
        "UniqueCombinedID_B", "maxArea", "maxAreaPeptide", "maxN", "n_2", 
        "repPepA", "repPepB", "group", "protein", "col_vector_peptides", 
        "repPepValue", "rowid"))

# Declare colors for plots
col_vector_peptides <- c("TRUE" = "#ffc125", "FALSE" = "#a020f0")
col_vector_proteins <- c("TRUE" = "#ff9d2e", "FALSE" = "#07b58a")