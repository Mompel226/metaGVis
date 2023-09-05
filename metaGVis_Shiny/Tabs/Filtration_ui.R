# Filtration tab content
filt_tab = tabItem(tabName = "filt",
                   fluidRow(
                     tabBox(
                       title = "Filtration", width = 12,
                       tabPanel("Filtration",
                                useShinyjs(), 
                                extendShinyjs(text = "shinyjs.button = function() {window.scrollTo(0,document.body.scrollHeight);}", functions = "button"),
                                uiOutput("inputs_filt")
                       ),
                       tabPanel("Prevalence", 
                                helpText("Plot total ASV/OTU abundance vs the fraction of samples in which an ASV/OTU is observed."),
                                hr(),
                                fluidRow(
                                  column(width = 12,
                                         uiOutput('physeq_prev'),
                                         uiOutput('rank_create_prev')
                                         ),
                                  uiOutput("ui.prev"),
                                  hr(),
                                  column(width = 12,
                                         uiOutput("down_prev")
                                         )
                                )
                                )
                       )
                     )
                   )