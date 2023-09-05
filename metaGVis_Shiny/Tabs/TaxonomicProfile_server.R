# Taxonomic profile tab content - SERVER
################################################################################
# Plot
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create,
             { # Display create plot button
               uiCreate("profile")
               
               # Reactive inputs 1
               output$inputs1_profile <- renderUI({
                 column(3,
                        selectInput(inputId = "data_profile", label = "Phyloseq object", choices = n.datalist()),
                        textInput(inputId = "title_profile", label = "Title", value = "CHANGE TITLE"),
                        textInput(inputId = "str2rm_profile", label = "Patterns to remove from sample names", value = NULL,
                                  placeholder = "Comma-separated patterns"),
                        div(HTML("<b> Labels size </b>")),
                        br(),
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-4", numericInput(inputId = "taxasize_profile", label = "Taxa", value = 11)),
                                 div(class="col-md-4", numericInput(inputId = "samplesize_profile", label = "Samples", value = 8)),
                                 div(class="col-md-4", numericInput(inputId = "stripsize_profile", label = "Strip", value = 10))
                          ),
                          column(width = 12,
                                 div(class="col-md-6", uiOutput("bs_size")),
                                 div(class="col-md-6", uiOutput("bs_nudge"))
                          )
                        )
                 )
               })
               
               # To allow variables to change when data is changed
               observeEvent(input$data_profile,
                            { # Get categorical variables
                              var_profile <<- getVar(input$data_profile)

                              # Get selected phyloseq object
                              ps_profile <<- reactive({ datalist[[input$data_profile]] })
                              
                              # Reactive inputs 2
                              output$inputs2_profile <- renderUI({
                                list(
                                  column(3,
                                         radioButtons(inputId = "abund_plot", "Type of plot", 
                                                      choices = c("barplot", "pointplot"), selected = "barplot", inline = TRUE),
                                         selectInput(inputId = "rank_profile", label = "Rank", 
                                                     choices = as.list(rank_names(ps_profile(), errorIfNULL=FALSE))),
                                         uiOutput("taxa_profile"),
                                         fluidRow(
                                           column(div(HTML("<b> Color palette seed for taxa </b>")), width = 12,
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col1_profile", label = "", value = "#ff0000")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col2_profile", label = "", value = "#00ff00")),
                                                  div(class="col-md-4", colourpicker::colourInput(inputId = "col3_profile", label = "", value = "#0000ff"))
                                           )
                                         ),
                                         checkboxInput(inputId = "unident_profile", "Show upper rank for unclassified", value = FALSE),
                                         uiOutput("format_profile")
                                  ),
                                  column(3,
                                         selectInput(inputId = "merge.by_profile", label = "Merge by", 
                                                     choices = c("None", var_profile()),
                                                     selected = "None"),
                                         checkboxInput(inputId = "count_plot", "Add library depth plot", value = FALSE),
                                         uiOutput("count_inputs"),
                                         checkboxInput(inputId = "dendo_plot", "Add dendogram plot", value = FALSE),
                                         uiOutput("dendo_inputs")
                                  ),
                                  column(3,
                                         selectInput(inputId = "top_var_profile", label = "Strip/Tile variable", 
                                                     choices = c("None", var_profile()),
                                                     selected = "None"),
                                         uiOutput("TopColOptions"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-12", uiOutput("top_col_profile")), 
                                           )
                                         ),
                                         uiOutput("bot_var_profile"),
                                         uiOutput("BotColOptions"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-12", uiOutput("bot_col_profile")), 
                                           )
                                         ),
                                         uiOutput("icons_var_profile"),
                                         uiOutput("icons_files"),
                                         uiOutput("IconsColOptions"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-12", uiOutput("icons_col_profile")), 
                                           )
                                         ),
                                  )
                                )
                              })
                            })
             
               # Display number of taxa per selected rank:               
               output$taxa_profile <- renderUI({
                 n.taxa <- length(unique(ps_profile()@tax_table[,input$rank_profile]))
                 # Define input
                 sliderInput(inputId = "ntaxa_profile", label = "Nº of taxa to display", 
                             min = 0, max = n.taxa, value = ifelse(n.taxa < 15, n.taxa, 15), step = 1)
               })
               
               # Show format of renamed unclassified taxa
               output$format_profile <- renderUI({
                 if ( isTRUE(input$unident_profile) ){
                   list(
                     # Show format of renamed unclassified taxa
                     radioButtons(inputId = "format_profile", "Format", choices = c("short", "long"), selected = "short"),
                     # Show the number of ASVs regrouped under each taxa
                     checkboxInput(inputId = "asv_profile", "Show ASVs regrouped under each taxa", value = TRUE) 
                   )
                 }
               })
               
               # Count plot inputs
               output$count_inputs <- renderUI({
                 if( isTRUE(input$count_plot) ){
                   list(
                     # Select phylo object used for the count plot 
                     selectInput(inputId = "data_count_profile", label = "Phyloseq obj. for lib. depth plot", choices = n.datalist()),
                     # Make count data relative
                     checkboxInput(inputId = "count_rel", "Relative abundance", value = FALSE)
                   )
                 }
               })
               
               # Dendogram inputs
               output$dendo_inputs <- renderUI({
                 if( isTRUE(input$dendo_plot) ){
                   list(
                     # Select transformation for dendogram
                     selectInput(inputId = "transf_profile", label = "Transformation",
                                 choices = c("Hellinger", "logChord_0plus1", "logChord_1p", "CLR_0plus1", "CLR_min", "CLR_runif"),
                                 selected = "CLR_runif"),
                     # Select distance for dendogram
                     selectInput(inputId = "dist_profile", 
                                 label = "Distance", 
                                 choices = c("Bray", "wUniFrac", "UniFrac", "Euclidean"),
                                 selected = "Euclidean"),
                     # Select clustering method for dendogram
                     selectInput(inputId = "clust_profile", label = "Clustering method", 
                                 choices = c("single", "complete", "average", "ward.D"),
                                 selected = "ward.D"),
                     # Select number of clusters for dendogram
                     numericInput(inputId = "nclust_profile", label = "Nº of clusters", value = 2),
                     # Select number of bootstraps for dendogram
                     numericInput(inputId = "nboot_profile", label = "Nº of bootstrap", value = 500),
                     # Select if the clusteing is performed using agglomerated data
                     checkboxInput(inputId = "aglom_clust", "Perform clustering using agglomerated taxa instead of ASVs/OTUs", value = FALSE),
                     # Apply y-axis nudge to bootstrap values
                     numericInput(inputId = "bs_nudge", label = "Y-nudge bs label", value = 0),
                     # Size of bootstrap label
                     numericInput(inputId = "bs_size", label = "Bootstrap", value = 2.5),
                     # Add icons to dendogram
                     checkboxInput(inputId = "dendo_icons", "Add sample icons to dendogram", value = FALSE)
                   )
                 }
               })
               
               # Select variable used to assign icons
               output$icons_var_profile <- renderUI({
                 if( isTRUE(input$dendo_plot) & isTRUE(input$dendo_icons) ){
                   selectInput(inputId = "icons_var_profile", label = "Icons variable",
                               choices = var_profile(),
                               selected = NULL) 
                 } else {
                   NULL
                 }
               })
               
               # Write PATH used to assign icons
               output$icons_files <- renderUI({
                 if( isTRUE(input$dendo_plot) & isTRUE(input$dendo_icons) ){
                   list(fileInput(inputId = "icons_files",
                                  label = "Icons files",
                                  accept = c(".png"),
                                  multiple = TRUE),
                        checkboxInput(inputId = "icons_with_color", "Colored icons", value = FALSE)
                   )
               }})

               # Create reactive values to avoid errors
               top_col_profile <<- reactive({ NULL })
               bot_col_profile <<- reactive({ NULL })
               icons_col_profile <<- reactive({ NULL })
               
               # Display color options for TOP variable levels
               observeEvent(input$top_var_profile,
                            { # Reset all color vectors
                              top_col_profile <<- reactive({ NULL })
                              bot_col_profile <<- reactive({ NULL })
                              icons_col_profile <<- reactive({ NULL })
                              if(input$top_var_profile != "None"){
                                output$TopColOptions <- renderUI({ ColOptions("top_profile", label = "Color") })
                              } else {
                                output$TopColOptions <- renderUI({ NULL })
                                output$top_col_profile <- renderUI({ NULL })
                              }
                            })
               # Display custom colors for TOP variable levels
               observeEvent(list(req(input$colset_top_profile), input$custom_col_top_profile, input$top_var_profile),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_top_profile)){
                                uiColLev("top_col_profile", "top_var_profile")
                                top_col_profile <<- ColVec("top_col_profile", "top_var_profile")
                              } else {
                                output$top_col_profile <- renderUI({ NULL })
                                top_col_profile <<- reactive({ 
                                  col <- RColorBrewer::brewer.pal(nrow(unique(ps_profile()@sam_data[,input$top_var_profile])), 
                                                                  input$colset_top_profile)
                                  lev.length <- length(levels(get_variable(ps_profile(), input$top_var_profile)))
                                  alpha(col[1:lev.length], input$alpha_top_profile)
                                })
                              }
                            })
               
               # Display BOTTOM strip variable
               output$bot_var_profile <- renderUI({
                 if(input$top_var_profile != "None" & !isTRUE(input$dendo_icons)){
                   selectInput(inputId = "bot_var_profile", label = "Bottom strip variable", 
                               choices = c("None", var_profile()[!var_profile() %in% input$top_var_profile]),
                               selected = "None")
                 } else {
                   output$BotColOptions <- renderUI({ NULL })
                   output$bot_col_profile <- renderUI({ NULL })
                   NULL
                 }
               })
               # Display color options for BOTTOM variable levels
               observeEvent(input$bot_var_profile,
                            { # Refresh bot_col_profile()
                              bot_col_profile <<- reactive({ NULL })
                              if(input$bot_var_profile != "None"){
                                output$BotColOptions <- renderUI({ ColOptions("bot_profile", label = "Color") })
                              } else {
                                output$BotColOptions <- renderUI({ NULL })
                                output$bot_col_profile <- renderUI({ NULL })
                              }
                            })
               # Display custom colors for BOTTOM variable levels
               observeEvent(list(req(input$colset_bot_profile), input$custom_col_bot_profile, input$bot_var_profile),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_bot_profile)){
                                uiColLev("bot_col_profile", "bot_var_profile")
                                bot_col_profile <<- ColVec("bot_col_profile", "bot_var_profile")
                              } else {
                                output$bot_col_profile <- renderUI({ NULL })
                                bot_col_profile <<- reactive({ 
                                  col <- RColorBrewer::brewer.pal(nrow(unique(ps_profile()@sam_data[,input$bot_var_profile])), 
                                                                  input$colset_bot_profile)
                                  lev.length <- length(levels(get_variable(ps_profile(), input$bot_var_profile)))
                                  alpha(col[1:lev.length], input$alpha_bot_profile)
                                })
                              }
                            })
               # Display color options for ICON variable levels
               observeEvent(list(input$icons_with_color, input$dendo_icons),
                            { # Reset the icon color vector
                              icons_col_profile <<- reactive({ NULL })
                              if(isTRUE(input$dendo_icons) & isTRUE(input$icons_with_color)){
                                output$IconsColOptions <- renderUI({ ColOptions("icons_profile", label = "Color") })
                              } else {
                                output$IconsColOptions <- renderUI({ NULL })
                                output$icons_col_profile <- renderUI({ NULL })
                              }
                            })
               # Display custom colors for ICONS variable levels
               observeEvent(list(req(input$colset_icons_profile), input$custom_col_icons_profile, 
                                 input$icons_var_profile,input$icons_with_color, input$dendo_icons),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_icons_profile)){
                                uiColLev("icons_col_profile", "icons_var_profile")
                                icons_col_profile <<- ColVec("icons_col_profile", "icons_var_profile")
                              } else if (isTRUE(input$icons_with_color)){
                                output$icons_col_profile <- renderUI({ NULL })
                                icons_col_profile <<- reactive({ 
                                  col <- RColorBrewer::brewer.pal(nrow(unique(ps_profile()@sam_data[,input$icons_var_profile])), 
                                                                  input$colset_icons_profile)
                                  lev.length <- length(levels(get_variable(ps_profile(), input$icons_var_profile)))
                                  alpha(col[1:lev.length], input$alpha_icons_profile)
                                })
                              } else {
                                bot_col_profile <<- reactive({ NULL })
                                output$icons_col_profile <- renderUI({ NULL })
                                icons_col_profile <<- reactive({ NULL })
                              }
                            })
               
               # Use the TOP color vector provided
               observeEvent(req(input$hex_top_col_profile),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_top_col_profile, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_profile(), input$top_var_profile)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                top_col_profile <<- reactive({ col.vect })
                              } else {
                                top_col_profile <<- ColVec("top_col_profile", "top_var_profile")
                              }
                            })
               # Use the BOTTOM color vector provided
               observeEvent(req(input$hex_bot_col_profile),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_bot_col_profile, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_profile(), input$bot_var_profile)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                bot_col_profile <<- reactive({ col.vect })
                              } else {
                                bot_col_profile <<- ColVec("bot_col_profile", "bot_var_profile")
                              }
                            })
               # Use the ICONS color vector provided
               observeEvent(req(input$hex_icons_col_profile),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_icons_col_profile, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_profile(), input$icons_var_profile)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                icons_col_profile <<- reactive({ col.vect })
                              } else {
                                icons_col_profile <<- ColVec("icons_col_profile", "icons_var_profile")
                              }
                            })
               
               # Create a vector with the strings to remove from the sample names
               str2rm_profile <- reactive({
                 ifelse(input$str2rm_profile!="", str2rm <- unlist(strsplit(input$str2rm_profile, ",")), str2rm <- "")
                 return(str2rm)
               })
               
               # Display download button when plot is created
               observeEvent(input$create_profile, {
                 output$down_profile <- renderUI({
                   dim_and_down("profile", 16, 7, ggplot = TRUE)
                 })
               })
               
               # Taxonomic profile plot
               profile <- eventReactive(input$create_profile, {
                 show_modal_spinner(
                   spin = "fading-circle",
                   color = "purple",
                   text = HTML("Generating plot... <br/>
                   Unfiltered data takes considerably longer. <br/>
                   The smaller the rank the longer it takes. <br/>
                   The bigger the bootstrap value the longer it takes. <br/>
                               Please be really patient.")
                 )
       
                 # Create a copy of files with original names in tmp folder
                 if (isTRUE(input$dendo_icons) & !is.null(names(input$icons_files))){
                   icons_lev <- levels(get_variable(datalist[[input$data_profile]], input$icons_var_profile))
                   file_names <- tools::file_path_sans_ext(input$icons_files$name)
                   if (length(icons_lev) == length(icons_lev[icons_lev %in% file_names])){
                     for (i in seq_along(icons_lev)){
                       file <- eval(parse(text = paste0("input$icons_files$datapath[", i, "]")))
                       name <- eval(parse(text = paste0("input$icons_files$name[", i, "]")))
                       icons.folder <- dirname(file)
                       file.copy(file, paste0(icons.folder, "/", name), overwrite = TRUE)
                     }
                     dendo_icons <- TRUE
                   } else {
                     showNotification("Icon names do not match the 'icon variable' level names", type = "error")
                     dendo_icons <- FALSE
                   }
                 } else {
                   # Make sure dendo_icons object exists
                   dendo_icons <- FALSE
                 }
                 
                 # Make sure tile.col object exists
                 tile.col <- "blank" # Reset object
                 ifelse(!is.null(top_col_profile()), tile.col <- top_col_profile(), tile.col <- "blank")
                 
                 # Make sure icons.col object exists
                 icons.col <- "blank" # Reset object
                 ifelse(!is.null(icons_col_profile()), icons.col <- icons_col_profile(), icons.col <- "blank")
                 
                 
                 # Make sure bot_var_profile object exists
                 bot_var_profile <- "None" # Reset object
                 ifelse(is.null(input$bot_var_profile), bot_var_profile <- "None", 
                        bot_var_profile <- input$bot_var_profile)
                 # Make sure bot_var_profile is equal to "None" when icons are displayed
                 ifelse(isTRUE(input$dendo_icons), bot_var_profile <- "None", "")

                 try(plot <- Taxonomic.analyis(physeq = datalist[[input$data_profile]],
                                               merge.by = input$merge.by_profile,
                                               title = input$title_profile,
                                               abundance.plot = input$abund_plot,
                                               rank = input$rank_profile,
                                               ntaxa = input$ntaxa_profile,
                                               show.unidentified = input$unident_profile,
                                               name.format = input$format_profile,
                                               ASV.count = input$asv_profile,
                                               count.plot = input$count_plot,
                                               data.count.plot = datalist[[input$data_count_profile]],
                                               count.plot.rel = input$count_rel,
                                               dendogram = input$dendo_plot,
                                               transf = input$transf_profile,
                                               distance = input$dist_profile,
                                               clust.method = input$clust_profile,
                                               n.clusters = input$nclust_profile,
                                               nboot = input$nboot_profile,
                                               aglom.clust = input$aglom_clust,
                                               nudge.bs.label = switch((input$bs_nudge != 0 + 1), NULL, input$bs_nudge),
                                               bs.label.size = input$bs_size,
                                               icons.dendogram = dendo_icons,
                                               icons.variable = switch(("icons_var_profile" %in% names(input) + 1), NULL, input$icons_var_profile),
                                               icons.folder = icons.folder,
                                               taxa.txt.size = input$taxasize_profile,
                                               sample.txt.size = input$samplesize_profile,
                                               strip.text.size = input$stripsize_profile,
                                               col.palette = c(input$col1_profile, input$col2_profile, input$col3_profile),
                                               Top.var = switch((input$top_var_profile=="None") + 1, input$top_var_profile, NULL),
                                               Top.strip.fill = switch((is.null(top_col_profile())) + 1, top_col_profile(), "blank"), 
                                               Bot.var = switch((bot_var_profile=="None") + 1, bot_var_profile, NULL),
                                               Bot.strip.fill = switch((is.null(bot_col_profile())) + 1, bot_col_profile(), "blank"),
                                               tile.col = tile.col,
                                               icons.col = icons.col,
                                               str2rm = str2rm_profile(), 
                                               return.plot = TRUE),
                     silent = F
                 )
                 remove_modal_spinner()
                 return(plot)
               })
               # Display profile plot
               observeEvent(profile(), {
                 display(profile(), "profile", "TaxonomicProfile", ggplot = TRUE, shiny.height = "100%")
               })
             }, once = TRUE)
               
                 
               