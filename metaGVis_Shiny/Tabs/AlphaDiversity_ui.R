# Alpha diversity tab content
alpha_tab = tabItem(tabName = "alpha",
                    fluidRow(
                      tabBox(
                        title = "Alpha Diversity", width = 12,
                        tabPanel("Alpha Diversity Measures", 
                                 helpText("Point plot showing the diversity of species (species richness) within a single sample.
                                          Samples can be further grouped on the horizontal axis to better estimate the diversity range of a particular microbiome (use boxplot).
                                          For meaningful results use unfiltered, non-normalized count data."),
                                 uiOutput("ui.alpha"),
                                 hr(),
                                 uiOutput('stat.tabs_alpha'),
                                 hr(),
                                 fixedRow(
                                   column(12, 
                                          uiOutput("create_alpha"),
                                          uiOutput("down_alpha"),
                                   )
                                 ),
                                 hr(),
                                 fixedRow(
                                   uiOutput("inputs1_alpha"),
                                   uiOutput("inputs2_alpha")
                                 )
                        ),
                        tabPanel("Whittaker", 
                                 helpText("Boxplot-based rank abundance curve (~ Whittaker plot). 
                                 Chart used to visualize species richness and species evenness across the different group of samples.
                                          The boxplot helps at visualising the relative abundance variance for each species."),
                                 uiOutput("ui.whit"),
                                 hr(),
                                 fixedRow(
                                   column(12, 
                                          uiOutput("create_whit"),
                                          uiOutput("down_whit"),
                                   )
                                 ),
                                 hr(),
                                 fixedRow(
                                   uiOutput("inputs1_whit"),
                                   uiOutput("inputs2_whit")
                                 )
                        )
                      )
                    )
)