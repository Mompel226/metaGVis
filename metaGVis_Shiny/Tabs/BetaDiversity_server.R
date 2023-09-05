# Beta diversity tab content - SERVER
################################################################################
# Plot
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create,
             { 
############## 2D Ordination plot ############### 
               
               # Display create plot button
               uiCreate("beta2d")
               
               # Reactive inputs 1
               output$inputs1_beta2d <- renderUI({
                 column(3,
                        selectInput(inputId = "data_beta2d", label = "Phyloseq object", choices = n.datalist(), width = "100%"),
                        textInput(inputId = "title_beta2d", label = "Title", value = NULL, width = "100%"),
                        checkboxInput(inputId = "names_beta2d", "Display names", value = FALSE, width = "100%"),
                        textInput(inputId = "str2rm_beta2d", label = "Patterns to remove from names", value = NULL,
                                  placeholder = "Comma-separated patterns", width = "100%"),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-6", numericInput(inputId = "pointsize_beta2d", label = "Point size", value = 3)),
                                 div(class="col-md-6", numericInput(inputId = "namesize_beta2d", label = "Name size", value = 2.5)),
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-6",                         
                                     numericInput(inputId = "axes.size_beta2d", label = "Axes size", value = 14)),
                                 div(class="col-md-6", uitheme("theme_beta2d"))
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-6",                         
                                     selectInput(inputId = "position_beta2d", label = "Legend",
                                                 choices = c("top", "bottom", "left", "right", "none"),
                                                 selected = "bottom")),
                                 div(class="col-md-6", 
                                     numericInput(inputId = "leg.size_beta2d", label = "Legend size", value = 14))
                          )
                        )
                 )
               })

               # To allow variables to change when data is changed
               observeEvent(input$data_beta2d,
                            { # Get selected phyloseq object
                              ps_beta2d <<- reactive({ datalist[[input$data_beta2d]] })
                              
                              # Get categorical variables
                              var_beta2d <<- getVar(input$data_beta2d)
                              
                              # Reactive inputs 2
                              output$inputs2_beta2d <- renderUI({
                                # Display only calculable distances
                                if(!is.null(ps_beta2d()@phy_tree)){
                                  distances <- c("euclidean", "bray", "unifrac", "wunifrac", "philr")
                                } else {
                                  distances <- c("euclidean", "bray")
                                }
                                # Display only normalisation and transformation  methods available per distance
                                if("philr" %in% distances){
                                  transformations <- c("None")
                                  normalisations <- c("None")
                                } else {
                                  transformations <- c("None", "Hellinger", "logChord_0plus1", "logChord_1p", "CLR_0plus1", "CLR_min", "CLR_runif")
                                  normalisations <- c("None", "TSS", "Rarefy", "UQ", "CSS", "DESeq", "TMM")
                                }
                                # UIs
                                list(
                                  column(3,
                                         # Select normalisation for ordination
                                         selectInput(inputId = "norm_beta2d", 
                                                     label = "Normalisation", 
                                                     choices = normalisations,
                                                     selected = "TSS", width = "100%"),
                                         # Select transformation for ordination
                                         selectInput(inputId = "transf_beta2d", 
                                                     label = "Transformation",
                                                     choices = transformations,
                                                     selected = "CLR_runif", width = "100%"),
                                         # Select distance for ordination
                                         selectInput(inputId = "dist_beta2d", 
                                                     label = "Distance", 
                                                     choices = distances,
                                                     selected = "euclidean", width = "100%"),
                                         # Select method for ordination
                                         selectInput(inputId = "method_beta2d", 
                                                     label = "Method", 
                                                     choices = c("PCA", "PCoA", "NMDS", "CCA", "RDA", "CAP"),
                                                     selected = "PCA", width = "100%"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-6", 
                                                      checkboxInput(inputId = "sqrt.dist", "Sqrt dist (adonis2)", value = FALSE, width = "100%"),),
                                                  div(class="col-md-6", 
                                                      checkboxInput(inputId = "add", "Add (adonis2)", value = FALSE, width = "100%"))
                                           )
                                         ),
                                         # Select multiple comparisons adjustment method for P-values
                                         selectInput(inputId = "adj_beta2d", 
                                                     label = "Correction method", 
                                                     choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                                 "fdr", "none"),
                                                     selected = "BH", width = "100%")
                                  ),
                                  column(6,
                                         textInput(inputId = "formula_beta2d", 
                                                   label = "RHS formula with spaces (for adonis2 and contrained CCA-RDA-CAP)",
                                                   width = "100%"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-4", selectInput(inputId = "col_var_beta2d", label = "Color variable", 
                                                                                    choices = c(var_beta2d()), width = "100%")),
                                                  div(class="col-md-4", selectInput(inputId = "shape_var_beta2d", label = "Shape variable", 
                                                                                     choices = c("None", var_beta2d()),
                                                                                     selected = "None", width = "100%")),
                                                  div(class="col-md-4", style = "margin-top: 23px",
                                                      checkboxInput(inputId = "ellipse_beta2d", "Display ellipses", 
                                                                    value = TRUE, width = "100%"))
                                           )
                                         ),
                                         hr(),
                                         column(6,
                                                ColOptions("beta2d", "Color"),
                                                div(class="col-md-12", uiOutput("col_beta2d")), 
                                         ),
                                         column(6,
                                                uiOutput("fill_options"),
                                                div(class="col-md-12", uiOutput("fill_beta2d")), 
                                         )
                                  )
                                )
                              })
                            })
               
               # Display color selection for grouping/color variable levels
               observeEvent(list(input$colset_beta2d, input$custom_col_beta2d),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_beta2d)){
                                uiColLev("col_beta2d", "col_var_beta2d", title = "Custom colors for points")
                                # Create color vector for color variable levels
                                col_beta2d <<- ColVec("col_beta2d", "col_var_beta2d")
                              } else {
                                output$col_beta2d <- renderUI({ NULL })
                                col_beta2d <<- reactive({ 
                                  col <- RColorBrewer::brewer.pal(nrow(unique(ps_beta2d()@sam_data[,input$col_var_beta2d])), 
                                                                  input$colset_beta2d) 
                                  lev.length <- length(levels(get_variable(ps_beta2d(), input$col_var_beta2d)))
                                  alpha(col[1:lev.length], input$alpha_beta2d)
                                })
                              }
                            })
               # Display fill options if ellipses are chosen to be displayed
               observeEvent(input$ellipse_beta2d,
                            { if(isTRUE(input$ellipse_beta2d)){
                                output$fill_options <- renderUI({ ColOptions("fill_beta2d", "Fill") })
                              } else {
                                output$fill_options <- renderUI({ NULL })
                              }
                            })
               # Display fill selection for grouping/color variable levels
               observeEvent(list(input$colset_fill_beta2d, input$custom_col_fill_beta2d),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_fill_beta2d) & isTRUE(input$ellipse_beta2d)){
                                uiColLev("fill_beta2d", "col_var_beta2d", title = "Custom fills for ellipses")
                                # Create color vector for color variable levels
                                fill_beta2d <<- ColVec("fill_beta2d", "col_var_beta2d")
                              } else {
                                output$fill_beta2d <- renderUI({ NULL })
                                fill_beta2d <<- reactive({ 
                                  col <- RColorBrewer::brewer.pal(nrow(unique(ps_beta2d()@sam_data[,input$col_var_beta2d])), 
                                                                  input$colset_fill_beta2d) 
                                  lev.length <- length(levels(get_variable(ps_beta2d(), input$col_var_beta2d)))
                                  alpha(col[1:lev.length], input$alpha_fill_beta2d)
                                })
                              }
                            })
               
               # Use the color vector provided
               observeEvent(req(input$hex_col_beta2d),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_col_beta2d, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_beta2d(), input$col_var_beta2d)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                col_beta2d <<- reactive({ col.vect })
                              } else {
                                col_beta2d <<- ColVec("col_beta2d", "col_var_beta2d")
                              }
                            })
               # Use the fill vector provided
               observeEvent(req(input$hex_fill_beta2d),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_fill_beta2d, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_beta2d(), input$col_var_beta2d)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                fill_beta2d <<- reactive({ col.vect })
                              } else {
                                fill_beta2d <<- ColVec("fill_beta2d", "col_var_beta2d")
                              }
                            })
               
               # Create a vector with the strings to remove from the sample names
               str2rm_beta2d <- reactive({
                 ifelse(input$str2rm_beta2d!="", str2rm <- unlist(strsplit(input$str2rm_beta2d, ",")), str2rm <- "")
                 return(str2rm)
               })
               
               # Display download button when plot is created
               observeEvent(input$create_beta2d, {
                 coef <- input$dimension[2]/input$dimension[1]
                 output$down_beta2d_1 <- renderUI({
                   dim_and_down("beta2d_1", 
                                round(2*coef*(input$dimension[2]/42),0), 
                                round(coef*(input$dimension[2]/65),0), ggplot = TRUE)
                 })
                 output$down_beta2d_2 <- renderUI({
                   dim_and_down("beta2d_2", 
                                round(2*coef*(input$dimension[2]/42),0), 
                                round(coef*(input$dimension[2]/65),0), ggplot = FALSE)
                 })
               })
               
               # 2D Ordination plot
               beta2d <- eventReactive(input$create_beta2d, {
                 # Check that the formula components are correct
                 if(input$formula_beta2d != ""){
                   formula_components <- unlist(strsplit(input$formula_beta2d, " "))
                   var_and_symbols <- c(var_beta2d(), "+", "-", ":", "*", "|")
                   if(any(is.na(match(formula_components, var_and_symbols)))){
                     m <- c()
                   } else {
                     m <- match(formula_components, var_and_symbols)
                   }
                   if(length(formula_components) != length(m)){
                     showNotification2("Formula is incorrect.<br/>
                                      Possible reasons:<br/>
                                      - Typos in variable names<br/>
                                      - Symbols not accepted <br/>
                                      (only '+', '-', ':', '*', and '|')<br/>
                                      - Spaces not included", 
                                       type = "error", 
                                       duration = 8)
                     return(NULL)
                     stop()
                   }
                 }
                 show_modal_spinner(
                   spin = "fading-circle",
                   color = "purple",
                   text = HTML("Generating 2D plot...")
                 )
                 try(plot <- ordination2D(physeq = datalist[[input$data_beta2d]],
                                          norm = input$norm_beta2d,
                                          transf = input$transf_beta2d,
                                          method = input$method_beta2d,
                                          distance = input$dist_beta2d,
                                          formula = input$formula_beta2d,
                                          color_var= input$col_var_beta2d,
                                          shape_var = switch((input$shape_var_beta2d=="None") + 1, input$shape_var_beta2d, NULL), 
                                          colors = col_beta2d(), 
                                          display.names = input$names_beta2d, 
                                          display.ellipses = input$ellipse_beta2d, 
                                          group.ellipses = input$col_var_beta2d, 
                                          ellipse.fill = fill_beta2d(),
                                          title = input$title_beta2d,
                                          str2rm = str2rm_beta2d(),
                                          point.size = input$pointsize_beta2d,
                                          names.size = input$namesize_beta2d,
                                          p.adj.method = input$adj_beta2d,
                                          sqrt.dist = input$sqrt.dist,
                                          add = input$add),
                     silent = T
                     )
                 try(plot[[1]] <- plot[[1]] + 
                       shiny_phyloseq_ggtheme_list[[input$theme_beta2d]] +
                       theme(axis.title = element_text(size = input$axes.size_beta2d + 2),
                             axis.text = element_text(size = input$axes.size_beta2d),
                             legend.position = input$position_beta2d,
                             legend.title = element_text(size = input$leg.size_beta2d + 2),
                             legend.text = element_text(size = input$leg.size_beta2d),
                             legend.key.size = unit(input$axes.size_beta2d * 0.06, 'cm')),
                     silent = T
                     )
                 remove_modal_spinner()
                 return(plot)
               })
               # Display ordination plot
               observeEvent(beta2d(), {
                 display(beta2d()[[1]], "beta2d_1", "Ordination2D", ggplot = TRUE, shiny.height = "100%")
               })
               # Display stress/shepard plot
               observeEvent(beta2d(), {
                 display(cowplot::ggdraw(beta2d()[[2]]), "beta2d_2", "Stress2D", ggplot = FALSE, shiny.height = "100%")
               })
               # Statistic report
               output$homo.g <- renderText({
                  paste0(capture.output(beta2d()[[3]][[1]]), collapse = "\n")
               })
               output$perma.g <- renderText({
                 paste0(capture.output(beta2d()[[3]][[3]]), collapse = "\n")
               })
               output$homo.p <- renderText({
                 paste0(capture.output(beta2d()[[3]][[2]]), collapse = "\n")
               })
               output$perma.p <- renderText({
                 paste0("","","","","","","",capture.output(beta2d()[[3]][[4]]), collapse = "\n")
               })
               
############## 3D Ordination plot ############### 
               
               beta3d <- eventReactive(beta2d(), {
                 show_modal_spinner(
                   spin = "fading-circle",
                   color = "purple",
                   text = HTML("Generating 3D plot...")
                 )
                 input_data <- normalisation(datalist[[input$data_beta2d]], input$norm_beta2d, input$col_var_beta2d)
                 input_data <- transformation(input_data, input$transf_beta2d)
                 try(plot <- metaGVis::MDSin3D(physeq = input_data, 
                                   explanatory_variable = input$col_var_beta2d, 
                                   distance = input$dist_beta2d, 
                                   method = input$method_beta2d, 
                                   ordination.object = beta2d()[[4]],
                                   colors = col_beta2d(),
                                   display.ellipses = input$ellipse_beta2d,
                                   constrained = switch((input$formula_beta2d=="") + 1, TRUE, FALSE), 
                                   formula = input$formula_beta2d))
                 remove_modal_spinner()
                 return(plot)
               })
               
               # Display ordination plot
               observeEvent(beta3d(), {
                 output$beta3d <- plotly::renderPlotly({ beta3d() })
                 output$ui.beta3d <- renderUI({
                   plotly::plotlyOutput(outputId = "beta3d", height = "100%")
                 })
                 output$down_beta3d <- renderUI({
                   downloadButton("down_3D_html", "Download html!")
                 })
                 # Display download button when plot is created
                 output$down_3D_html <- downloadHandler(
                   filename = function(file) { "Ordination3D.html"},
                   content = function(file) {
                     # export plotly html widget as a temp file to download.
                     saveWidget(as_widget(beta3d()), file, selfcontained = TRUE)
                   }
                 )
                 })

############## Venn/Euler diagram ############### 
               
               # Display create plot button
               uiCreate("venn")
               
               # Reactive inputs 1
               output$inputs1_venn <- renderUI({
                 column(4,
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-12",
                                   selectInput(inputId = "data_venn", label = "Phyloseq object", choices = n.datalist()),
                                   uiOutput("title_venn")
                                 )
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-6", numericInput(inputId = "numbersize_venn", label = "Number size", value = 1.3, step = 0.1)),
                                 div(class="col-md-6", numericInput(inputId = "labelsize_venn", label = "Label size", value = 1.2, step = 0.1))
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-6", colourpicker::colourInput(inputId = "numbercolor_venn", label = "Number color", value = "#000000"),),
                                 div(class="col-md-6",
                                     colourpicker::colourInput(inputId = "labelcolor_venn", label = "Label color", value = "#000000"),
                                     div(style = "margin-top: -17px"),
                                     checkboxInput(inputId = "labelcolors_venn", label = "Same as fill", value = FALSE))
                          )
                        ),
                        fluidRow(
                          column(width = 12,
                                 HTML("<br/> <b> Legend options: </b>"),
                                 hr(),
                                 div(class="col-md-4", numericInput(inputId = "legtxtsize_venn", label = "Text size", value = 10, step = 1)),
                                 div(class="col-md-4", numericInput(inputId = "legheight_venn", label = "Height", value = 1, step = 0.5)),
                                 div(class="col-md-4", numericInput(inputId = "legwidth_venn", label = "Width", value = 1, step = 0.5))
                          )
                        )
                 )
               })
               
               # To allow variables to change when data is changed
               observeEvent(input$data_venn,
                            { # Get selected phyloseq object
                              ps_venn <<- reactive({ datalist[[input$data_venn]] })
                              
                              # Get categorical variables with 7 or less levels
                              var_venn <<- reactive({
                                  sam.data <- datalist[[input$data_venn]]@sam_data
                                  if (!is.null(sam.data)){
                                    if (any(sapply(sam.data, is.factor))){
                                      variables <- colnames(sam.data[,sapply(sam.data, is.factor) & sapply(sam.data, nlevels) <= 7 ])
                                      return(variables)
                                    } else {
                                      return("Select categorical variables")
                                    }
                                  }
                                  return("Upload Metadata file")
                                })
                              
                              # Reactive inputs 2
                              output$inputs2_venn <- renderUI({
                                list(
                                  column(4,
                                         fluidRow(
                                           column(width = 12,
                                                  div(class="col-md-4", radioButtons(inputId = "type_venn", label = "Plot type",
                                                                                     choices = c("Venn", "Euler"), selected = "Venn")),
                                                  uiOutput("pack_venn"),
                                                  uiOutput("pattern_venn")
                                           )
                                         ),
                                         uiOutput("heatmapcol_venn"),
                                         fluidRow(
                                           column(width = 12,
                                                  div(class="col-md-6", selectInput(inputId = "label_venn", label = "Data shown (max 2)", 
                                                                                    choices = c("counts", "percent", "both", "abundance"), 
                                                                                    multiple = TRUE, selected = "counts")),
                                                  div(class="col-md-6", selectInput(inputId = "var_venn", label = "Grouping variable", 
                                                                                    choices = c(var_venn())))
                                           )
                                         ),
                                         uiOutput("relative_venn"),
                                         ColOptions("venn"),
                                         fluidRow(
                                           column(width = 12,
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "edgecolor_venn", label = "Edge color", value = "#000000"),
                                                      div(style = "margin-top: -17px"),
                                                      checkboxInput(inputId = "edgecolors_venn", label = "Same as fill", value = FALSE)),
                                                  div(class="col-md-4", selectInput(inputId = "edgelty_venn", label = "Edge type", 
                                                                                    choices = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                                                                                    selected = "solid")),
                                                  div(class="col-md-4", numericInput(inputId = "edgesize_venn", label = "Edge size", value = 0.7, min = 0.1, max = 5, step = 0.1))
                                           )
                                         ),
                                         fluidRow(
                                           column(width = 12,
                                                  div(class="col-md-12", selectInput(inputId = "align_venn", label = "Alignment of plots", 
                                                                                    choices = c("horizontal", "vertical"), 
                                                                                    selected = "horizontal",
                                                                                    width = "100%"))
                                                  )
                                         )
                                  ),
                                  column(4,
                                         uiOutput("intersect"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-12", uiOutput("col_venn")), 
                                           )
                                         )
                                  )
                                )
                              })
                            })
               
               # Display dynamic UIs
               observeEvent(req(input$type_venn), 
                            {# Only display "ggVennDiagram" when Venn is selected
                              if(input$var_venn == "Select categorical variables"){
                                type <- c("ggVennDiagram", "eulerr")
                              } else if (req(input$type_venn) == "Euler"){
                                type <- c("eulerr")
                              } else if (length(levels(get_variable(ps_venn(), input$var_venn))) <= 5){
                                type <- c("ggVennDiagram", "eulerr")
                              } else {
                                type <- c("ggVennDiagram")
                              }
                              output$pack_venn <- renderUI({ 
                                div(class="col-md-4", radioButtons(inputId = "pack_venn", label = "Package", 
                                                                   choices = type))
                              })
                            })
               observeEvent(req(input$pack_venn), 
                            {# Display heatmap option when "ggVennDiagram" is selected
                              output$pattern_venn <- renderUI({
                                if (input$pack_venn == "ggVennDiagram"){
                                  patt <- c("normal", "heatmap")
                                } else {
                                  patt <- c("normal")
                                }
                                div(class="col-md-4", radioButtons(inputId = "pattern_venn", label = "Color pattern", choices = patt))
                              })
                            })
               observeEvent(input$pattern_venn, 
                            {# Display heatmap gradient color options
                              output$heatmapcol_venn <- renderUI({
                                if(input$pattern_venn == "heatmap"){
                                  fluidRow(
                                    column(div(HTML("<b> Colors for heatmap gradient </b>")), width = 12,
                                           div(class="col-md-6", colourpicker::colourInput(inputId = "col1_venn", label = "", value = "#F8C9BE")),
                                           div(class="col-md-6", colourpicker::colourInput(inputId = "col2_venn", label = "", value = "#FF0000"))
                                    )
                                  )
                                }
                              })
                            })
               observeEvent(input$var_venn,
                            {# Display available levels for the intersections
                              output$intersect <- renderUI({
                                if (input$var_venn != "Select categorical variables"){
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-12",
                                               selectInput(inputId = "inter_venn", label = "Intersections to display in table", multiple = TRUE,
                                                           choices = levels(get_variable(datalist[[input$data_venn]], input$var_venn)),
                                                           selected = levels(get_variable(datalist[[input$data_venn]], input$var_venn)),
                                                           width = "100%")
                                           ),
                                           div(class="col-md-12",
                                               selectInput(inputId = "rank_venn", label = "Rank to agglomerate data", 
                                                           choices = as.list(c("None", rank_names(ps_venn(), errorIfNULL=FALSE))),
                                                           selected = "None"))
                                    )
                                  )
                                }
                              })
                            })
               observeEvent(input$label_venn,
                            {# Display relative abundance options when label "abundance" is selected
                              output$relative_venn <- renderUI({ 
                                if(any(input$label_venn == "abundance")){
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", numericInput(inputId = "fraction_venn", label = "Filter fraction", value = 0, min = 0, max = 1, step = 0.05)),
                                           div(class="col-md-6", style = "margin-top: 21px",
                                               checkboxInput(inputId = "relative_venn", label = "Relative abundance", value = TRUE))
                                    )
                                  )
                                } else {
                                  NULL
                                }
                              })
                              # Display two titles if more than two labels are selected
                              output$title_venn <- renderUI({
                                if(length(input$label_venn) > 1){
                                  list(
                                    textInput(inputId = "title1_venn", label = "Title 1", value = "default", width="100%"),
                                    textInput(inputId = "title2_venn", label = "Title 2", value = "default", width="100%")
                                  )
                                } else {
                                  textInput(inputId = "title1_venn", label = "Title", value = "default", width="100%")
                                }
                              })
                            })
               
               # Display color selection for Color variable levels
               observeEvent(list(input$custom_col_venn, input$colset_venn),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_venn)){
                                uiColLev("col_venn", "var_venn")
                                # Create color vector for grouping/color variable levels
                                col_venn <<- ColVec("col_venn", "var_venn")
                              } else {
                                output$col_venn <- renderUI({ NULL })
                                col_venn <<- reactive({ 
                                  RColorBrewer::brewer.pal(nrow(unique(ps_venn()@sam_data[,input$var_venn])), input$colset_venn) })
                              }
                            })
               
               # Use the color vector provided
               observeEvent(req(input$hex_col_venn),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_col_venn, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_venn(), input$var_venn)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                col_venn <<- reactive({ col.vect })
                              } else {
                                col_venn <<- ColVec("col_venn", "var_venn")
                              }
                            })
               
               # Venn/Euler plot
               venn <- eventReactive(input$create_venn, {
                 show_modal_spinner(
                   spin = "fading-circle",
                   color = "purple",
                   text = HTML("Generating plot...")
                 )
                 output$update_venn <- renderUI({
                   div(class='col-md-3',
                       actionButton("update_venn", "Update table",
                                    style="color: #fff; background-color: #de9849; border-color: #4a2700; margin-top: 25px;",
                                    icon = icon("redo-alt")))
                 })
                 # Make sure the number of colors matches the number of levels
                 lev <- levels(get_variable(ps_venn(), input$var_venn))
                 col <- col_venn()[1:length(lev)]
                 
                 try(plot <- multi_venn(physeq = datalist[[input$data_venn]], 
                                        group = input$var_venn, 
                                        euler = switch((input$type_venn!="Euler") + 1, TRUE, FALSE), 
                                        relative = switch((is.null(input$relative_venn)) + 1, input$relative_venn, FALSE), 
                                        package = input$pack_venn, 
                                        fraction = switch((is.null(input$fraction_venn)) + 1, input$fraction_venn, 0), 
                                        label = input$label_venn,
                                        show2plots = switch((length(input$label_venn) == 1) + 1, TRUE, FALSE), 
                                        heatmap = switch((input$pattern_venn == "normal") + 1, TRUE, FALSE),
                                        label_color = switch((isTRUE(input$labelcolors_venn)) + 1, input$labelcolor_venn, col), 
                                        label_size = input$labelsize_venn, 
                                        number_color = input$numbercolor_venn, 
                                        number_size = input$numbersize_venn, 
                                        fill_color = col,
                                        fill_alpha = input$alpha_venn,
                                        edge_color = switch((isTRUE(input$edgecolors_venn)) + 1, input$edgecolor_venn, col), 
                                        edge_size = input$edgesize_venn, 
                                        edge_lty = input$edgelty_venn,
                                        gradient2colors = switch((is.null(input$col1_venn)) + 1, c(input$col1_venn, input$col2_venn), c("#F8C9BE", "#FF0000")), 
                                        plotly = FALSE,
                                        plot.align = input$align_venn,
                                        title1 = input$title1_venn, 
                                        title2 = switch((is.null(input$title2_venn)) + 1, input$title2_venn, "default"),
                                        leg.txt.size = input$legtxtsize_venn,
                                        leg.key.height = input$legheight_venn,
                                        leg.key.width = input$legwidth_venn),
                     silent = T
                 )
                 remove_modal_spinner()
                 return(plot)
               })
               
               observeEvent(venn(), {
                 # Display download button when plot is created
                 output$down_venn <- renderUI({
                   ifelse(length(input$label_venn) > 1, mult <- 2, mult <- 1)
                   coef <- input$dimension[2]/input$dimension[1]
                   dim_and_down("venn", 
                                round(mult*coef*(input$dimension[2]/42),0), 
                                round(coef*(input$dimension[2]/50),0), ggplot = TRUE)
                 })
                 # Display venn plot
                 display(venn(), "venn", input$type_venn, ggplot = TRUE, shiny.height = "100%")
               })
               
               # Venn/Euler data frame
               observeEvent(list(venn(), input$update_venn), {
                 if (length(input$label_venn) == 1){
                   output$table_venn <- DT::renderDataTable({
                     intersect.df(physeq = datalist[[input$data_venn]],
                                  group = input$var_venn,
                                  levels = input$inter_venn,
                                  rank = input$rank_venn,
                                  relative = switch((is.null(input$relative_venn)) + 1, input$relative_venn, FALSE),
                                  weighted = ifelse(input$label_venn == "abundance", TRUE, FALSE),
                                  fraction = switch((is.null(input$fraction_venn)) + 1, input$fraction_venn, 0))
                   },
                   options = list(
                     autoWidth = TRUE,
                     scrollX = TRUE,
                     scrollY = "400px"
                   ), rownames= FALSE)
                   output$ui.venn.df <- renderUI({ DT::dataTableOutput("table_venn") })
                 } else {
                   output$ui.venn.df <- renderUI({ NULL })
                 }
               })
             }, once = TRUE)