# Import tab content
import_tab = tabItem(tabName = "import",
                     fixedRow(
                       align = "center",
                       box(
                         title = "QIIME2 -- MOTHUR -- DADA2 -- Kraken2/Bracken -- Centrifuge -- MetaPhlAn3 -- HUMAnN3",
                         width = 12,
                         collapsible = TRUE,
                         collapsed = TRUE,
                         icon = icon("folder"),
                         status = "primary",
                         tags$style(HTML("
                         .box.box-primary>.box-header {
                         color:#404040; 
                         background:#EBF5FB}
                         .box.box-primary{
                         border-top-color:#2E86C1;
                         background:#EBF5FB}
                         
                         .box.box-info>.box-header {
                         color:#404040; 
                         background:#FFFFFF}
                         .box.box-info{
                         border-top-color:#2E86C1;
                         background:#FFFFFF}")),
                         h5(style="font-family:Helvetica; font-weight: bold; color:#404040;",
                            HTML("Upload files, remove undesired patterns from the sample names and 
                         establish in which order the levels of the categorical variables will be displayed in the plots.")),
                         hr(),
                         box(title = "MOTHUR, DADA2 or Kraken2/Bracken files",
                             align = "left",
                             width = 6,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             icon = icon("folder"),
                             status = "info",
                             uiOutput('input_files1'),
                             hr(),
                             actionButton('reset1', 'Reset',
                                          style="color: #fff; background-color: #009ACD; border-color: #507786",
                                          icon = icon("redo-alt"))
                         ),
                         box(title = "MetaPhlAn3 file",
                             align = "left",
                             width = 6,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             icon = icon("folder"),
                             status = "info",
                             uiOutput('input_files2'),
                             hr(),
                             actionButton('reset2', 'Reset',
                                          style="color: #fff; background-color: #009ACD; border-color: #507786",
                                          icon = icon("redo-alt"))
                         ),
                         box(title = "HUMAnN3 file",
                             align = "left",
                             width = 6,
                             collapsible = TRUE,
                             collapsed = TRUE,
                             icon = icon("folder"),
                             status = "info",
                             uiOutput('input_files3'),
                             hr(),
                             actionButton('reset3', 'Reset',
                                          style="color: #fff; background-color: #009ACD; border-color: #507786",
                                          icon = icon("redo-alt"))

                         ),
                         fixedRow(align = "left",
                           column(width = 12,
                                  div(class="col-md-4", 
                                      fileInput(inputId = "meta",
                                                label = HTML("Metadata (Optional)"),
                                                accept = c(".xls", ".xlsx", ".xlsm", ".csv", ".tsv", ".txt")
                                      ),
                                      h6(style="font-family:Helvetica; 
                                      margin-top:-20px; font-weight: bold; color: red;", 
                                         "Metadata's 1st column MUST contain the sample names.")
                                  ),
                                  div(class="col-md-4",
                                      fileInput(inputId = "tree",
                                                label = HTML("Phylogenetic tree (Optional)"),
                                                accept = c(".qza", ".tree")
                                      )
                                  ),
                                  div(class="col-md-4", 
                                      textInput(inputId = "sep", 
                                                label = HTML("Field separator for 'txt' files (regex format)"), 
                                                value = "\t")
                                  )
                           )
                         ),
                         HTML("<b> Sample names </b>"),
                         h5(style="font-family:Helvetica; font-weight: bold; color:red;", "MUST match metadata sample names"),
                         uiOutput("samples"),
                         br(),
                         uiOutput('str2rm'),
                         fluidRow(
                           column(width = 12,
                                  uiOutput('variables')
                           ),
                           column(width = 12,
                                  uiOutput('levels')
                           )
                         ),
                         column(width = 12,
                                div(class="col-md-4", NULL),
                                div(class="col-md-4", textInput(inputId = "name", label = "Name of the phyloseq object", value = "data_phylo")),
                                div(class="col-md-4", NULL)
                                ),
                         actionButton("create", "Create phyloseq object!",
                                      style="color: #fff; background-color: #7ab733; border-color: #527b22",
                                      icon = icon("pencil-alt")),
                         br(),
                         column(width = 12,
                                div(class="col-md-4", NULL),
                                div(class="col-md-4", uiOutput('ui.physeqnames')),
                                div(class="col-md-4", NULL)
                                )
                       ),
                       box(title = "CheckM",
                           width = 12,
                           collapsible = TRUE,
                           collapsed = TRUE,
                           icon = icon("folder"),
                           tags$style(HTML("
                         .box.box-danger>.box-header {
                         color:#404040; 
                         background:#fdf7f7}
                         .box.box-danger{
                         border-top-color:#cc0000;
                         background:#fdf7f7}
                         
                         .box.box-navy>.box-header {
                         color:#404040; 
                         background:#FFFFFF}
                         .box.box-navy{
                         border-top-color:#cc0000;
                         background:#FFFFFF}")),
                           status = "danger"
                       )
                     )
)