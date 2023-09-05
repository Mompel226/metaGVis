# Functional profile tab content
function_tab = tabItem(tabName = "function",
                      fluidRow(
                        tabBox(
                          title = "Functional Profile", width = 12,
                          tabPanel("Functional Profile", 
                                   helpText("."),
                                   uiOutput("ui.function"),
                                   uiOutput("ui.legend"),
                                   hr(),
                                   fixedRow(
                                     column(12, 
                                            uiOutput("create_function"),
                                            uiOutput("down_function"),
                                     ),
                                     column(12, 
                                            div(class="col-md-3", NULL),
                                            uiOutput("down_legend")
                                     )
                                   ),
                                   hr(),
                                   fixedRow(
                                     uiOutput("inputs1_function"),
                                     uiOutput("inputs2_function")
                                   )
                          )
                        )
                      )
)