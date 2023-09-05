# Library QC tab content - SERVER
################################################################################
# To plot histogram
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create,{
  # Select phylo object used for histogram
  output$inputs_histo <- renderUI({
    sidebarPanel(
      selectInput(inputId = "data_histo", label = "Phyloseq object", choices =  n.datalist()),
      textInput(inputId = "title_histo", label = "Title", value = 'Histogram of library sizes'),
      textInput(inputId = "xlab_histo", label = "X-axis title", value = "Library size (reads)"),
      fluidRow(column(width = 12,
                      div(class="col-md-4", 
                          colourpicker::colourInput(inputId = "col_histo", label = "Colour", value = '#0030A080', allowTransparent = TRUE)),
                      div(class="col-md-4", 
                          numericInput(inputId = "bins_histo", "Bins", value = 10)),
                      div(class='col-md-4', 
                          uitheme("theme_histo"))
      ),
      column(width = 12,
             div(class="col-md-4", 
                 numericInput("width_histo", "Plot width", value = 1, min = 1, max = 100, step = 1)),
             div(class="col-md-4", 
                 numericInput("height_histo", "Plot height", value = 1, min = 1, max = 100, step = 1)),
             div(class='col-md-4', 
                 graphicTypeUI("downtype_histo"))
      ),
      column(width = 12, align="center", 
             div(class='col-md-12', 
                 downloadButton("download_histo", 'Download plot!', style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
      )
      )
    )
  })
  # Histogram
  histo <- reactive({
    plot = NULL
    try(plot <- ggplot(data.frame(x=colSums(datalist[[input$data_histo]]@otu_table)), aes(x = x)) + 
          geom_histogram(bins = input$bins_histo, fill = input$col_histo) + 
          ggtitle(input$title_histo) +
          xlab(input$xlab_histo) + 
          ylab("Frequency") +
          shiny_phyloseq_ggtheme_list[[input$theme_histo]] +
          theme(axis.title.y = element_text(margin = margin(0,10,0,0)),
                axis.title.x = element_text(margin = margin(10,0,0,0))),
        silent = TRUE)
    return(plot)
  })
  # Display histogram plot
  observeEvent(histo(), 
               { display(histo(), "histo", "Histogram", ggplot = TRUE)
               })
  ################################################################################
  # To plot rarecurves
  ################################################################################
  # Reactive UI - reactive parameters depending on the phylo object
  # ONLY when phyloseq object is created
  # Display create plot button
  uiCreate("rare")
  # Display download button when plot is created
  observeEvent(input$create_rare, {
    output$down_rare <- renderUI({
      dim_and_down("rare", 16, 7, ggplot = FALSE)
    })
  })
  # Display basic UI inputs for the rarecurves plot
  output$inputs1_rare <- renderUI({
    list(
      selectInput(inputId = "data_rare", label = "Phyloseq object", choices = n.datalist()),
      textInput(inputId = "title_rare", label = "Title", value = "CHANGE TITLE"),
      textInput(inputId = "y_lab_rare", label = "Y-axis title", value = "Number of different ASVs"),
      checkboxInput(inputId = "labels", label = "Display sample names", value = TRUE),
      textInput(inputId = "str2rm_rare", label = "Patterns to remove from sample names", value = NULL,
                placeholder = "Comma-separated patterns")
    )
  })
  # Display only categorical variables in color and shape input boxes
  observeEvent(input$data_rare,
               { # Get selected phyloseq object
                 ps_rare <<- reactive({ datalist[[input$data_rare]] })
                 # Get defined categorical variables
                 var_rare <- getVar(input$data_rare)
                 # Display step and categorical variable for COLOR and LINE-TYPE
                 output$inputs2_rare <- renderUI({
                   list(
                     sliderInput(inputId = "step", label = "Step size", min = 100, max = 10000, value = 1000),
                     selectInput(inputId = "rank_rare", label = "Taxonomic rank", 
                                 choices = as.list(rank_names(ps_rare(), errorIfNULL=FALSE))),
                     selectInput(inputId = "col_var_rare", label = "Color variable", choices = c("None", var_rare())),
                     uiOutput("col_rare_options"),
                     fluidRow(
                       column(12,
                              div(class="col-md-12", uiOutput("col_rare")), 
                       )
                     ),
                     selectInput(inputId = "line_var_rare", label = "Line type variable", choices = c("None", var_rare()))
                   )
                 })
               })
  # Display legend parameters when COLOR or LINE-TYPE variables are selected
  output$legend_rare <- renderUI({
    if (req(input$col_var_rare) != "None" | req(input$line_var_rare) != "None"){
      list(
        # Legend title and background colors
        column(width = 12,
               div(class="col-md-6",
                   colourpicker::colourInput(inputId = "titlecol", label = "Legend title colour", value = "#00008B")),
               div(class="col-md-6",
                   colourpicker::colourInput(inputId = "back", label = "Legend background colour", value = "#D3D3D3", allowTransparent = TRUE))
               ),
        # Text size, position, width and heigth of legend boxes 
        column(width = 12,
               div(class="col-md-6", 
                   numericInput(inputId = "text_size", "Legend text size", value = 1.4)),
               div(class="col-md-6", 
                   numericInput("x", "X-axis start of the legend boxes", value = 0.8, min = 0.01, max = 1)),
        ),
        column(width = 12,
               div(class="col-md-6", 
                   numericInput("y1", "Y-axis start legend box 1", value = 0.88, min = 0.1, max = 1)),
               div(class="col-md-6", 
                   numericInput("y2", "Y-axis start legend box 1", value = 0.42, min = 0.1, max = 1)),
        ),
        column(width = 12,
               div(class="col-md-6", 
                   numericInput("width1", "Legend box 1 width", value = 0.45, min = 0.1, max = 5)),
               div(class="col-md-6", 
                   numericInput("width2", "Legend box 2 width", value = 0.45, min = 0.1, max = 5)),
        ),
        column(width = 12,
               div(class="col-md-6", 
                   numericInput("height1", "Legend box 1 height", value = 1.5, min = 0.1, max = 10)),
               div(class="col-md-6", 
                   numericInput("height2", "Legend box 2 height", value = 1.5, min = 0.1, max = 10)),
        )
      )
    }
  })
  
  # Create reactive values to avoid errors
  col_rare <<- reactive({ NULL })
  # Display color options for COLOR variable levels
  observeEvent(input$col_var_rare,
               { # Refresh col_var_rare()
                 col_rare <<- reactive({ NULL })
                 if(input$col_var_rare != "None"){
                   output$col_rare_options <- renderUI({ ColOptions("col_rare", label = "Color") })
                 } else {
                   output$col_rare_options <- renderUI({ NULL })
                   output$col_rare <- renderUI({ NULL })
                 }
               })
  # Display custom colors for COLOR variable levels
  observeEvent(list(req(input$colset_col_rare), input$custom_col_col_rare),
               {# Display color and fill selection for grouping/color variable levels
                 if(isTRUE(input$custom_col_col_rare)){
                   uiColLev("col_rare", "col_var_rare")
                   col_rare <<- ColVec("col_rare", "col_var_rare")
                 } else {
                   output$col_rare <- renderUI({ NULL })
                   col_rare <<- reactive({ 
                     print(input$colset_col_rare)
                     col <- RColorBrewer::brewer.pal(nrow(unique(ps_rare()@sam_data[,input$col_var_rare])), 
                                                     input$colset_col_rare)
                     lev.length <- length(levels(get_variable(ps_rare(), input$col_var_rare)))
                     alpha(col[1:lev.length], input$alpha_col_rare)
                   })
                 }
               })
  # Use the COLOR vector provided
  observeEvent(req(input$hex_col_rare),
               {# Make sure text input are HEX codes
                 col.vect <- strsplit(input$hex_col_rare, ",")[[1]]
                 lev.length <- length(levels(get_variable(ps_rare(), input$col_var_rare)))
                 if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                   col_rare <<- reactive({ col.vect })
                 } else {
                   col_rare <<- ColVec("col_rare", "col_var_rare")
                 }
               })
  # Create a vector with the strings to remove from the sample names
  str2rm_rare <- reactive({
    ifelse(input$str2rm_rare!="", str2rm <- unlist(strsplit(input$str2rm_rare, ",")), str2rm <- "")
    return(str2rm)
  })
  # Rarecurves plot
  vegan.data <- eventReactive(input$create_rare, { 
    show_modal_spinner(
      spin = "fading-circle",
      color = "purple",
      text = HTML("Generating plot... <br/>
      The smaller the step size the longer it takes. . <br/>
      Please be patient.")
    )
    try(vegan.data <- rare.curve(physeq = datalist[[input$data_rare]],
                           step = input$step, 
                           rank = input$rank_rare,
                           return.vegandata = TRUE),
        silent = TRUE)
    remove_modal_spinner()
    return(vegan.data)
  })
  rarecurve <- reactive({
    if (req(input$create_rare) > 0){
      print(col_rare())
      print(col_rare())
      try(plot <- rare.curve(physeq = datalist[[input$data_rare]],
                             vegandata = vegan.data(),
                             rank = input$rank_rare,
                             title = input$title_rare,
                             y_lab = input$y_lab_rare,
                             col_var = switch((input$col_var_rare=="None") + 1, input$col_var_rare, NULL), 
                             colors = switch((is.null(col_rare())) + 1, col_rare(), NULL),
                             line_var = switch((input$line_var_rare=="None") + 1, input$line_var_rare, NULL), 
                             linetypes = NULL,
                             labels = input$labels, 
                             legend.box1.width = input$width1,
                             legend.box2.width = input$width2,
                             legend.box1.height = input$height1,
                             legend.box2.height = input$height2,
                             legend.text.size = input$text_size,
                             legend.titlecol = input$titlecol,
                             legend.back = input$back,
                             legend.box.col = c(input$x, input$y1),
                             legend.box.line = c(input$x, input$y2),
                             str2rm = str2rm_rare()),
          silent = TRUE)
      return(plot)
    } else {
      return(NULL)
    }
  })
  
  # Display rarefaction curves plot
  observeEvent(rarecurve(), {
    display(rarecurve(), "rare", "RareCurves", ggplot = FALSE, shiny.height = "100%")
  })
  ################################################################################
  # To plot sample depth
  ################################################################################
  # Reactive UI - reactive parameters depending on the phylo object
  # ONLY when phyloseq object is created
  # Display create plot button
  uiCreate("depth")
  # Display download button when plot is created
  observeEvent(input$create_depth, 
               { output$down_depth <- renderUI({
                 dim_and_down("depth", 16, 7, ggplot = TRUE)})
               })

  # Display basic UI inputs for the rarecurves plot
  output$inputs_depth <- renderUI({
    list(
      selectInput(inputId = "data_depth", label = "Phyloseq object", choices = n.datalist()),
      textInput(inputId = "title_depth", label = "Title", value = "CHANGE TITLE"),
      textInput(inputId = "str2rm_depth", label = "Patterns to remove from sample names", value = NULL,
                placeholder = "Comma-separated patterns")
    )
  })
  observeEvent(input$data_depth,
               { # Get defined categorical variables
                 var_depth <- getVar(input$data_depth)
                 # Select categorical variable for TOP strip
                 output$top_var_depth <- renderUI({
                   selectInput(inputId = "top_var_depth", label = "Top strip variable", 
                               choices = c("None", var_depth()),
                               selected = "None")
                 })
                 # Select categorical variable for BOTTOM strip
                 output$bot_var_depth <- renderUI({
                   selectInput(inputId = "bot_var_depth", label = "Bottom strip variable", 
                               choices = c("None", var_depth()),
                               selected = "None")
                 })
               })
  # Display color selection for TOP variable levels
  observeEvent(input$top_var_depth, 
               { uiColLev("top_col_depth", "top_var_depth")
               }) 
  # Display color selection for BOTTOM variable levels
  observeEvent(input$bot_var_depth, 
               { uiColLev("bot_col_depth", "bot_var_depth")
               }) 
  # Create color vector for TOP variable levels
  topcols_depth <- ColVec("top_col_depth", "top_var_depth")
  # Create color vector for BOTTOM variable levels
  botcols_depth <- ColVec("bot_col_depth", "bot_var_depth")
  
  # Create a vector with the strings to remove from the sample names
  str2rm_depth <- reactive({
    ifelse(input$str2rm_depth!="", str2rm <- unlist(strsplit(input$str2rm_depth, ",")), str2rm <- "")
    return(str2rm)
  })
  
  # Sample depth plot
  depth <- eventReactive(input$create_depth, {
    show_modal_spinner(
      spin = "fading-circle",
      color = "purple",
      text = HTML("Generating plot... <br/>
      Please be patient.")
    )
    try(plot <- Taxonomic.analyis(data.count.plot = datalist[[input$data_depth]],
                                  abundance.plot = "blank", 
                                  count.plot = TRUE, 
                                  title = input$title_depth,
                                  Top.var = switch((input$top_var_depth=="None") + 1, input$top_var_depth, NULL),
                                  Top.strip.fill = switch((is.null(topcols_depth())) + 1, topcols_depth(), "blank"), 
                                  Bot.var = switch((input$bot_var_depth=="None") + 1, input$bot_var_depth, NULL),
                                  Bot.strip.fill = switch((is.null(botcols_depth())) + 1, botcols_depth(), "blank"),
                                  str2rm = str2rm_depth(), 
                                  return.plot = TRUE),
        silent = TRUE
    )
    remove_modal_spinner()
    return(plot)
  })
  # Display depth plot
  observeEvent(depth(), {
    display(depth(), "depth", "LibraryDepth", ggplot = TRUE, shiny.height = "100%")
  })
}, once = TRUE)