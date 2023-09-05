# # Functional Profile - SERVER
################################################################################
# Functional Profile tab
# SERVER code for the reactive UI # Activates when phyloseq object is created
################################################################################
observeEvent(input$create,
             { # Display create plot button
               uiCreate("function")
               
               # Reactive inputs 1
               output$inputs1_function <- renderUI({
                   column(3,
                          selectInput(inputId = "data_function", label = "Phyloseq object", choices = n.datalist()),
                          div(HTML("<b> Distribution of sub-plots </b>")),
                          fluidRow(
                            column(width = 12,
                                   div(class="col-md-6", numericInput(inputId = "x_function", label = "x", value = NULL)),
                                   div(class="col-md-6", numericInput(inputId = "y_function", label = "y", value = NULL))
                            )
                          ),
                          div(HTML("<b> Labels size </b>")),
                          fluidRow(
                            column(width = 12,
                                   div(class="col-md-6", numericInput(inputId = "labelsize_function", label = "Levels", value = 1.8)),
                                   div(class="col-md-6", numericInput(inputId = "numbersize_function", label = "Numbers", value = 1))
                            )
                          )
                   )
                 })
               
               # To allow variables to change when data is changed
               observeEvent(input$data_function,
                            { # Get categorical variables
                              var_function <- getVar(input$data_function)
                              
                              ps <- datalist[[input$data_function]]
                              taxa.ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
                              function.ranks <- c("Superclass1", "Superclass2", "pathway")
                              
                              # Reactive inputs 2
                              output$inputs2_function <- renderUI({
                                list(
                                  # Display categorical variables, taxonomic and functional ranks available
                                  column(3,
                                         selectInput(inputId = "var_function", label = "Grouping variable",
                                                     choices = var_function(),
                                                     selected = NULL),
                                         selectInput("rank_function", "Function rank", 
                                                     choices = as.list(rank_names(ps, errorIfNULL=FALSE)[rank_names(ps, errorIfNULL=FALSE) %in% function.ranks])),
                                         selectInput("taxa_function", "Taxonomic rank", 
                                                     choices = as.list(rank_names(ps, errorIfNULL=FALSE)[rank_names(ps, errorIfNULL=FALSE) %in% taxa.ranks]))
                                  ),
                                  # Display extra options
                                  column(3,
                                         fluidRow(
                                           checkboxInput(inputId = "unclass_function", "Show unclassified functions", value = FALSE),
                                           checkboxInput(inputId = "unclass_taxa_function", "Show unclassified taxa", value = TRUE),
                                           column(div(HTML("<b> Color palette seed for functions </b>")), width = 12,
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col4_function", label = "", value = "#ff0000")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col5_function", label = "", value = "#00ff00")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col6_function", label = "", value = "#0000ff"))
                                           ),
                                           column(div(HTML("<b> Color palette seed for taxa </b>")), width = 12,
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col1_function", label = "", value = "#ff0000")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col2_function", label = "", value = "#00ff00")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col3_function", label = "", value = "#0000ff"))
                                           )
                                         )
                                  )
                                )
                              })
                            })
               
               # Display download button when plot is created
               observeEvent(input$create_function, 
                            { output$down_function <- renderUI({
                              dim_and_down("function", 16, 7, ggplot = FALSE)
                              })
                            output$down_legend <- renderUI({
                              dim_and_down("legend", 16, 7, ggplot = FALSE)
                            })
                            })

               # Generate main functional profile plot
               mainplot <- eventReactive(input$create_function,
                                         { show_modal_spinner(
                                           spin = "fading-circle",
                                           color = "purple",
                                           text = HTML("Generating plot... <br/>
                                           The smaller the rank the longer it takes. <br/>
                                                       Please be really patient.")
                                         )
                                           # Establish the matrix dimensions for the plot
                                           mat <- reactive({
                                             lev <- length(levels(get_variable(datalist[[input$data_function]], input$var_function)))
                                             if (is.na(input$x_function)){
                                               y <- if((lev %% 2) == 0) { round(lev/2) } else { round(lev/2) + 1 }
                                               x <- round(lev/2)
                                               mat <- c(x,y)
                                             } else {
                                               mat <- c(input$x_function, input$y_function)
                                             }
                                             return(mat)
                                           })
                                           try(plot <- FunctionalProfile(physeq = datalist[[input$data_function]],
                                                                     taxa.rank = input$taxa_function,
                                                                     func.rank = input$rank_function,
                                                                     variable = input$var_function,
                                                                     show.unclass.taxa = input$unclass_taxa_function, 
                                                                     show.unclass.func = input$unclass_function, 
                                                                     taxa.color.palette = c(input$col1_function, input$col2_function, input$col3_function),
                                                                     function.color.palette = c(input$col3_function, input$col4_function, input$col5_function),
                                                                     plot.mat = mat(),
                                                                     labelsize = input$labelsize_function,
                                                                     numbersize = input$numbersize_function,
                                                                     return.legend = FALSE)
                                               )
                                           remove_modal_spinner()
                                           return(plot)
                                         })
               # Generate legend of functional profile plot
               legend <- eventReactive(input$create_function,
                                         { show_modal_spinner(
                                           spin = "fading-circle",
                                           color = "purple",
                                           text = HTML("Generating the legend... <br/>
                                                       Please be patient.")
                                         )
                                           plot <- try(FunctionalProfile(physeq = datalist[[input$data_function]],
                                                                     taxa.rank = input$taxa_function,
                                                                     func.rank = input$rank_function,
                                                                     variable = input$var_function,
                                                                     show.unclass.taxa = input$unclass_taxa_function, 
                                                                     show.unclass.func = input$unclass_function,
                                                                     taxa.color.palette = c(input$col1_function, input$col2_function, input$col3_function),
                                                                     function.color.palette = c(input$col3_function, input$col4_function, input$col5_function),
                                                                     plot.mat = NULL,
                                                                     labelsize = input$labelsize_function,
                                                                     numbersize = input$numbersize_function,
                                                                     return.legend = TRUE)
                                                       )
                                           remove_modal_spinner()
                                           return(plot)
                                         })
               
               # Render Main plot
               observeEvent(mainplot(), {
                 display(mainplot(), "function", "FunctionalProfile", ggplot = FALSE, shiny.height = "100%")
               })
               # Render Legend
               observeEvent(legend(), {
                 display(legend(), "legend", "Legend", ggplot = FALSE, shiny.height = "100%", data.suffix = "function")
               })
             }, once = TRUE)