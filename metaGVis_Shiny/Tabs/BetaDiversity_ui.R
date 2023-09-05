# Beta diversity tab content
beta_tab = tabItem(tabName = "beta",
                   fluidRow(
                     tabBox(
                       title = "Beta diversity", width = 12,
                       tabPanel("Ordination",
                                hr(),
                                tabsetPanel(
                                  tags$style(HTML("
    .tabbable > .nav > li > a[data-value='2D ordination plot'] {background-color: #4f3f04;   color:white}
    .tabbable > .nav > li > a[data-value='Stress/Shepard plot'] {background-color: #F9AD6A;  color:black}
    .tabbable > .nav > li > a[data-value='3D ordination plot'] {background-color: #F6D55C;   color:black}
    .tabbable > .nav > li > a[data-value='Statistics'] {background-color: #D46C4E; color:black}
  ")),
                                  tabPanel("2D ordination plot",
                                           wellPanel(style = "background-color:#e2dcd1;border-color:#4f3f04;",
                                                     helpText(""),
                                                     uiOutput("ui.beta2d_1"),
                                                     hr(),
                                                     fluidRow(
                                                       column(12, 
                                                              uiOutput("create_beta2d"),
                                                              uiOutput("down_beta2d_1")
                                                       )
                                                     ),
                                                     hr(),
                                                     fluidRow(
                                                       uiOutput("inputs1_beta2d"),
                                                       uiOutput("inputs2_beta2d")
                                                     )
                                                     
                                           )
                                  ),
                                  tabPanel("Stress/Shepard plot",
                                           helpText("In ordination analysis, stress and Shepard plots are used to assess 
                                                    the quality of the ordination plot and the accuracy of the distance matrix used to generate it.
                                                    If the stress value is high or the Shepard plot shows a poor fit between the observed and predicted 
                                                    distances, it suggests that the ordination plot may not accurately represent the relationships 
                                                    between samples in the original data."),
                                           wellPanel(style = "background-color:#fff1e5;border-color:#F9AD6A;",
                                                     helpText("."),
                                                     uiOutput("ui.beta2d_2"),
                                                     hr(),
                                                     fluidRow(
                                                       column(12,
                                                              uiOutput("down_beta2d_2")
                                                       )
                                                     )
                                           )
                                  ),
                                  tabPanel("3D ordination plot", 
                                           wellPanel(style = "background-color:#fff8e4;border-color:#F6D55C;",
                                                     helpText("."),
                                                     uiOutput("ui.beta3d"),
                                                     hr(),
                                                     fluidRow(
                                                       column(12,
                                                              uiOutput("down_beta3d"),
                                                       )
                                                     )
                                           )
                                  ),
                                  tabPanel("Statistics",
                                           wellPanel(style = "background-color:#fce6e0;border-color:#D46C4E;",
                                                     hr(),
                                                     HTML(paste(
                                                       "In ordination analysis, homogeneity of dispersion and PERMANOVA statistical tests are used to assess the significance of the differences between groups of samples.
                                                       Use PERMANOVA to detect differences in the locations (centroids) of different groups.",
                                                       "When using the PERMANOVA test, it specifically tests the null hypothesis: 'the centroids of the groups are equal for all groups.'",
                                                       "If the null hypothesis is rejected, then we can assume that the microbial communities of the different groups differ.",
                                                       "However, the method may confound location and dispersion effects: significant differences may be caused by different within-group variation (dispersion) instead of different mean values of the groups.",
                                                       "PERMDISP is a common test completed in conjunction with PERMANOVA and tests the null hypothesis of 'no difference in dispersion between groups.'",
                                                       "This test can identify if it is the dispersion of the group data from the centroids that is driving the significance (of the PERMANOVA test) or if it is the centroids of the group data themselves.",
                                                       "However, for both the PERMANOVA and PERMDISP tests, it is ideal to have equal sample sizes.",
                                                       "Unbalanced experimental designs can either increase rejection rates or the test can become more conservative.",
                                                       sep="<br/>")),
                                                     hr(),
                                                     fluidRow(
                                                       column(12, 
                                                              div(class="col-md-6", HTML(paste("GLOBAL HOMOGENEITY OF DISPERSION test:",
                                                                                               "If p > 0.05 then the groups have the same dispersion.",
                                                                                               sep="<br/>")),
                                                                  hr(),
                                                                  verbatimTextOutput("homo.g")),
                                                              div(class="col-md-6", HTML(paste("GLOBAL PERMANOVA test:",
                                                                                               "If p < 0.05 then the groups have different centroids.",
                                                                                               sep="<br/>")),
                                                                  hr(),
                                                                  verbatimTextOutput("perma.g"))
                                                       )
                                                     ),
                                                     hr(),
                                                     fluidRow(
                                                       column(12, 
                                                              div(class="col-md-6", HTML(paste("PAIRWISE HOMOGENEITY OF DISPERSION test:",
                                                                                               "If adjusted p > 0.05 then the groups have the same dispersion.",
                                                                                               sep="<br>")),
                                                                  hr(),
                                                                  verbatimTextOutput("homo.p")),
                                                              div(class="col-md-6", HTML(paste("PAIRWISE PERMANOVA test:",
                                                                                               "If adjusted p < 0.05 then the groups have different centroids.",
                                                                                               sep="<br>")),
                                                                  hr(),
                                                                  verbatimTextOutput("perma.p"))
                                                       )
                                                     )
                                           )
                                  )
                                )
                       ),
                       tabPanel("Venn/Euler diagrams", 
                                helpText("Sum of each species' total relative abundance, which in turn was computed by averaging the species' mean relative abundance found in each group of samples"),
                                fluidRow(
                                  column(12, 
                                         div(class="col-md-6", uiOutput("ui.venn")),
                                         div(class="col-md-6", uiOutput("ui.venn.df"))
                                  )
                                ),
                                hr(),
                                fluidRow(
                                  column(12, 
                                         uiOutput("create_venn"),
                                         uiOutput("down_venn")
                                  ),
                                  column(12, 
                                         uiOutput("update_venn"),
                                  )
                                ),
                                hr(),
                                fluidRow(
                                  uiOutput("inputs1_venn"),
                                  uiOutput("inputs2_venn")
                                )
                       )
                     )
                   )
)

