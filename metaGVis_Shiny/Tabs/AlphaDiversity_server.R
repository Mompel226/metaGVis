# AlphaDiversity tab content - SERVER
################################################################################
# To plot Alpha Diversity mesasures
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(list(input$create,input$filter),
             { # Display create plot button
               uiCreate("alpha")
               
               # Reactive inputs 1
               output$inputs1_alpha <- renderUI({
                   column(3,
                          selectInput(inputId = "data_alpha", label = "Phyloseq object", choices = n.datalist()),
                          textInput(inputId = "title_alpha", label = "Title", value = "CHANGE TITLE"),
                          uitheme("theme_alpha"),
                          fluidRow(
                            column(width = 12,
                                   div(class="col-md-6", numericInput(inputId = "samplesize_alpha", label = "Samples size", value = 16)),
                                   div(class="col-md-6", numericInput(inputId = "tipsize_alpha", label = "Tip length", value = 0.01))
                            )
                          )
                   )
                 })
               # To allow variables to change when data is changed
               observeEvent(input$data_alpha,
                            { # Get selected phyloseq object
                              ps_alpha <<- reactive({ datalist[[input$data_alpha]] })
                              
                              # Get categorical variables
                              var_alpha <- getVar(input$data_alpha)
                              
                              # Reactive inputs 2
                              output$inputs2_alpha <- renderUI({
                                list(
                                  column(3,
                                         selectInput(inputId = "measure_alpha", 
                                                     label = "Alpha-diversity measure", 
                                                     choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"),
                                                     multiple = TRUE,
                                                     selected = c("Observed", "Shannon", "Simpson")),
                                         checkboxInput(inputId = "boxplot", "Boxplot", value = TRUE),
                                         checkboxInput(inputId = "pair_alpha", "Pairwise comparison", value = FALSE),
                                         uiOutput("pair_method_alpha")
                                         
                                  ),
                                  column(3,
                                         selectInput(inputId = "var_alpha", label = "Group variable", 
                                                     choices = c(var_alpha()),
                                                     selected = var_alpha()),
                                         selectInput(inputId = "col_var_alpha", label = "Color variable", 
                                                     choices = c("None", var_alpha()),
                                                     selected = "None"),
                                         selectInput(inputId = "shape_var_alpha", label = "Shape variable", 
                                                     choices = c("None", var_alpha()),
                                                     selected = "None"),
                                         selectInput(inputId = "legend_alpha", label = "Displayed legend", 
                                                     choices = c("colors", "labels", "both"),
                                                     selected = "colors")
                                  ),
                                  column(3,
                                         ColOptions("alpha", "Color"),
                                         fluidRow(
                                           column(12,
                                                  div(class="col-md-12", uiOutput("col_alpha")), 
                                           )
                                         )
                                  )
                                )
                              })
                            })
               
               # Display pairwise comparison methods when pair_alpha is TRUE
               output$pair_method_alpha <- renderUI({
                 if( isTRUE(input$pair_alpha) ){
                   list(
                     fluidRow(
                       column(12,
                              div(class="col-md-12", 
                                  selectInput(inputId = "pair_method_alpha", label = "Mean comparison method", 
                                              choices = list('Parametric' = list("t.test"), 'Non-parametric' = list("wilcox.test")),
                                              selected = "wilcox.test"),
                                  div(style = "margin-top: -17px"),
                                  checkboxInput(inputId = "paired_alpha", label = "Paired", value = FALSE))
                       )
                     ),
                     fluidRow(
                       column(12,
                              selectInput(inputId = "p.adjust_alpha", label = "Method for adjusting p values", 
                                          choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                                          selected = "BH"),
                              checkboxInput(inputId = "hide_alpha", "Hide ns symbol", value = FALSE),
                              checkboxInput(inputId = "rm_alpha", "Remove ns brackets", value = TRUE)
                       )
                     )
                   )
                 }
               })
            
               # Display color selection for Color variable levels
               observeEvent(list(input$colset_alpha, input$custom_col_alpha),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_alpha)){
                                uiColLev("col_alpha", "col_var_alpha")
                                # Create color vector for color variable levels
                                col_DA <<- ColVec("col_alpha", "col_var_alpha")
                              } else {
                                output$col_alpha <- renderUI({ NULL })
                                col_alpha <<- reactive({ 
                                  RColorBrewer::brewer.pal(nrow(unique(ps_alpha()@sam_data[,input$col_var_alpha])), 
                                                           input$colset_alpha) })
                              }
                            })
               # Use the color vector provided
               observeEvent(req(input$hex_col_alpha),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_col_alpha, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_alpha(), input$col_var_alpha)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                col_alpha <<- reactive({ col.vect })
                              } else {
                                col_alpha <<- ColVec("col_alpha", "col_var_alpha")
                              }
                            })
               
               # Display download button when plot is created
               observeEvent(input$create_alpha, {
                 output$down_alpha <- renderUI({
                   dim_and_down("alpha", 16, 7, ggplot = TRUE)
                 })
                 
                 # Alpha diversity measures plot
                 alpha <- eventReactive(input$create_alpha, {
                   show_modal_spinner(
                     spin = "fading-circle",
                     color = "purple",
                     text = HTML("Generating plot... <br/>
                               Please be patient.")
                   )
                   
                   if (isTRUE(input$pair_alpha)){
                     m <- sample_data(ps_alpha())
                     m <- m %>% 
                       group_by_at(input$var_alpha) %>%
                       summarise(Nrows = length((!!sym(input$var_alpha))))
                     if (isTRUE(input$paired_alpha) & length(unique(m$Nrows)) != 1){
                       text = "Levels of the group variable must have the same number of samples when paired option is selected."
                       p <- ggplot() + 
                         annotate("text", x = 4, y = 25, size=8, color="red", label = text) + 
                         theme_void()
                       remove_modal_spinner()
                       plot <- list(p, NULL)
                       return(plot)
                     }
                   }
                   
                   try(plot <- Alpha_Diversity_plot(physeq = datalist[[input$data_alpha]], 
                                                    color_variable = switch((input$col_var_alpha != "None") + 1, NULL, input$col_var_alpha), 
                                                    shape_variable = switch((input$shape_var_alpha != "None") + 1, NULL, input$shape_var_alpha), 
                                                    group_variable = switch((input$var_alpha != "None") + 1, NULL, input$var_alpha),
                                                    boxplot = input$boxplot,
                                                    measures = input$measure_alpha, 
                                                    colors = col_alpha(), 
                                                    pair_compare = input$pair_alpha, 
                                                    hide.ns = input$hide_alpha, 
                                                    remove.ns = input$rm_alpha, 
                                                    pair_compare_method = input$pair_method_alpha,
                                                    paired = input$paired_alpha,
                                                    p.adjust = input$p.adjust_alpha,
                                                    title = input$title_alpha,
                                                    size.text.samples = input$samplesize_alpha,
                                                    tip.length = input$tipsize_alpha,
                                                    ggtheme = shiny_phyloseq_ggtheme_list[[input$theme_alpha]],
                                                    shown.legend = input$legend_alpha,
                                                    return.p.values = input$pair_alpha),
                       silent = TRUE
                   )
                   remove_modal_spinner()
                   return(plot)
                 })
                 
                 # Display Alpha diversity measures plot
                 observeEvent(alpha()[[1]], {
                   display(alpha()[[1]], "alpha", "AlphaDiversity", ggplot = TRUE, shiny.height = "100%")
                 })
                 # Display statistical results as tabs
                 observe({
                   if(isTRUE(input$pair_alpha) & !is.null(alpha()[[2]])){
                     # Display stats if calculated for the specified Grouping Variable
                     lev <- levels(get_variable(ps_alpha(), input$var_alpha))
                     g1 <- unlist(alpha()[[2]][[1]]["group1"])
                     g2 <- unlist(alpha()[[2]][[1]]["group2"])
                     g <- unique(c(g1,g2))
                     if (all(lev %in% g)){
                       lapply(seq_len(length(alpha()[[2]])), function(i) {
                         local({
                           measure <- input$measure_alpha[i]
                           richness <- estimate_richness(datalist[[input$data_alpha]], 
                                                         measures = measure)
                           out <- shapiro.test(richness[,1]) # Normality test
                           tablename <- paste0("Stats_", measure)
                           output[[tablename]] <- renderTable({
                             as.data.frame(alpha()[[2]][[i]])
                           })
                           normname <- paste0("Norm_", measure)
                           output[[normname]] <- renderTable({ 
                             data.frame("W" = out$statistic, "p-value" = out$p.value)
                           })
                         })
                       })
                     } else {
                       lapply(seq_len(length(alpha()[[2]])), function(i) {
                         local({
                           measure <- input$measure_alpha[i]
                           tablename <- paste0("Stats_", measure)
                           output[[tablename]] <- renderTable({ NULL })
                           normname <- paste0("Norm_", measure)
                           output[[normname]] <- renderTable({ NULL })
                         })
                       })
                     }
                   }
                 })
                 output$stat.tabs_alpha <- renderUI({
                   if(isTRUE(input$pair_alpha) & !is.null(alpha()[[2]])){
                     measure <- input$measure_alpha
                     myTabs = lapply(measure, function(x){
                       N.test <- paste0('Norm_', x)
                       P.stats <- paste0('Stats_', x)
                       tabPanel(P.stats,
                                helpText(HTML("Normality test performed with the shapiro.test() function from the stats R package. <br/> 
                                                If the test is significant (p < 0.05), the distribution is non-normal. <br/>")),
                                tableOutput(N.test),
                                br(),
                                helpText(HTML("Statistical analysis performed with the compare_means() function from the ggpubr R package. <br/> ")),
                                tableOutput(P.stats)
                       )
                     })
                     do.call(tabsetPanel, myTabs)
                   }
                 })
               })
               
################################################################################
# To plot Boxplot-based rank abundance curve (~ Whittaker plot)
################################################################################  
               
               uiCreate("whit")
               
               # Reactive inputs 1
               output$inputs1_whit <- renderUI({
                 column(4,
                        selectInput(inputId = "data_whit", label = "Phyloseq object", choices = n.datalist()),
                        textInput(inputId = "title_whit", label = "Title", value = "CHANGE TITLE"),
                        uitheme("theme_whit")
                 )
               })
               
               # To allow variables to change when data is changed
               observeEvent(input$data_whit,
                            { # Get categorical variables
                              var_whit <- getVar(input$data_whit)
                              
                              # Reactive inputs 2
                              output$inputs2_whit <- renderUI({
                                list(
                                  column(4,
                                         selectInput(inputId = "var_whit", label = "Group variable", 
                                                     choices = c(var_whit()),
                                                     selected = var_whit())
                                  ),
                                  column(4,
                                         numericInput(inputId = "colum_whit", label = "Nº of columns", value = 3)
                                  )
                                )
                              })
                            })
                 
               # Whittaker plot
               whit <- eventReactive(input$create_whit, {
                 if (input$var_whit != "Upload Metadata file"){
                   show_modal_spinner(
                     spin = "fading-circle",
                     color = "purple",
                     text = HTML("Generating plot... <br/>
                               Please be patient.")
                   )
                   df <- as.data.frame(datalist[[input$data_whit]]@otu_table)
                   df[df > 0] <- 1
                   group.max.x = max(colSums(df)) + (max(colSums(df)) * 0.1)
                   try(plot <- whittaker(physeq = datalist[[input$data_whit]],
                                         group = input$var_whit,
                                         max.x = group.max.x,
                                         ncol = input$colum_whit,
                                         title = input$title_whit,
                                         ggtheme = shiny_phyloseq_ggtheme_list[[input$theme_whit]]),
                       silent = F
                   )
                   remove_modal_spinner()
                   return(plot)
                 }
               })
               
               observeEvent(whit(), {
                 # Display download button when plot is created
                 output$down_whit <- renderUI({
                   coef <- input$dimension[2]/input$dimension[1]
                   dim_and_down("whit", 
                                round(2*coef*(input$dimension[2]/42),0), 
                                round(coef*(input$dimension[2]/50),0), ggplot = TRUE)
                 })
                 # Display whittaker plot
                 display(whit(), "whit", "Whittaker", ggplot = TRUE, shiny.height = "100%")
               })
             }, once = TRUE)