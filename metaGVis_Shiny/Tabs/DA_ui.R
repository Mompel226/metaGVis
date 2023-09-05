# Beta diversity tab content
DA_tab = tabItem(tabName = "DA",
                 fluidRow(
                   tabBox(
                     title = "Differential abundance", width = 12,
                     tabPanel("Differential abundance", 
                              helpText("."),
                              fluidRow(
                                column(12, 
                                       uiOutput("ui.DA")
                                )
                              ),
                              hr(),
                              fluidRow(
                                column(12, 
                                       uiOutput("create_DA"),
                                       uiOutput("down_DA")
                                )
                              ),
                              hr(),
                              fluidRow(
                                column(3,
                                       uiOutput("inputs1_DA"),
                                       uiOutput("inputs2_DA")
                                       ),
                                uiOutput("inputs3_DA"),
                              )
                     )
                   )
                 )
)