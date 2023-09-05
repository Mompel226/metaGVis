# Taxonomic profile tab content
profile_tab = tabItem(tabName = "profile",
                      fluidRow(
                        tabBox(
                          title = "Taxonomic Profile", width = 12,
                          tabPanel("Taxonomic Profile", 
                                   helpText("."),
                                   uiOutput("ui.profile"),
                                   hr(),
                                   fixedRow(
                                     column(12, 
                                            uiOutput("create_profile"),
                                            uiOutput("down_profile"),
                                     )
                                   ),
                                   hr(),
                                   fixedRow(
                                     uiOutput("inputs1_profile"),
                                     uiOutput("inputs2_profile")
                                   )
                          )
                        )
                      )
)