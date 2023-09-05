# Library QC tab content
libqc_tab = tabItem(tabName = "libqc",
                    fluidRow(
                      tabBox(
                        title = "Library QC", width = 12,
                        tabPanel("Histogram",
                                 helpText("Histogram of library sizes."),
                                 sidebarLayout(
                                   uiOutput("inputs_histo"),
                                   mainPanel(uiOutput("ui.histo"))
                                 )
                        ),
                        tabPanel("Rarefaction curves",
                                 helpText("Plot of the number of OTUs/ASVs/species against the number of reads for each sample."),
                                 uiOutput("ui.rare"),
                                 hr(),
                                 fixedRow(
                                   column(12, 
                                          uiOutput("create_rare"),
                                          uiOutput("down_rare"),
                                          )
                                 ),
                                 hr(),
                                 fixedRow(
                                   column(3, 
                                          uiOutput("inputs1_rare")
                                   ),
                                   column(3,
                                          uiOutput('inputs2_rare')
                                   ),
                                   column(6,
                                          uiOutput("legend_rare"),
                                          
                                   )
                                 )
                        ),
                        tabPanel("Depth of domains", 
                                 helpText("Stacked bar plot of the number of reads belonging to each domain for each sample."),
                                 uiOutput("ui.depth"),
                                 hr(),
                                 fixedRow(
                                   column(12, 
                                          uiOutput("create_depth"),
                                          uiOutput("down_depth"),
                                   )
                                 ),
                                 hr(),
                                 fixedRow(
                                   column(4, 
                                          uiOutput("inputs_depth")
                                   ),
                                   column(4,
                                          uiOutput("top_var_depth"),
                                          uiOutput("top_col_depth")
                                   ),
                                   column(4,
                                          uiOutput("bot_var_depth"),
                                          uiOutput("bot_col_depth")
                                   )
                                 )
                        )
                      )
                    )
                    )