###################
# sidebar.R
#
###################
sidebar <- dashboardSidebar(
  sidebarMenu(id = "menu0",
              menuItem("Introduction", tabName = "intro", icon = icon("hand-point-right")),
              menuItem("Import", tabName = "import", icon = icon("file-upload")),
              menuItem("Part1: Peptide-to-protein", tabName = "ptp", #icon = icon("tools"),
                       menuSubItem("Summary", tabName = "summary", icon = icon("eye")),
                       menuSubItem("Filter and Select", tabName = "filter", icon = icon("filter")),
                       menuSubItem("Rep Peptides", tabName = "bylabelstate", icon = icon("cogs")),
                       menuSubItem("Normalize", tabName = "normalize", icon = icon("ruler-horizontal"))

                       ),
              menuItem("Part 2: Protein workflow", tabName = "pw", #icon = icon("tools"),
                       menuSubItem("Normalized Proteins", tabName = "proteinNormViz", icon = icon("chart-bar")),
                       menuSubItem("Heatmaps", tabName = "heatMaps", icon = icon("chart-bar")),
                       menuSubItem("Co-migration plots", tabName = "coMigration", icon = icon("chart-bar")),
                       menuSubItem("Cluster", tabName = "cluster", icon = icon("shapes"))
              ),
              # menuItem("Visualize II", tabName = "dataviz", icon = icon("chart-bar")),
              menuItem("Q & A", tabName ="QuestionsAndAnswers", icon = icon("question")),
              menuItem("Feedback", tabName="feedback", icon = icon("comment", lib="glyphicon"))
  ),
  helpText("Developed by Rick Scavetta and Petra Palenikova as part of the R ComPrAn package",
           style="padding-left:1em; padding-right:1em;position:absolute; bottom:1em; ")
)
