# Filtration - SERVER
################################################################################
# Filtration tab: to filter data 
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create, {
  output$inputs_filt <- renderUI({ 
    js$button()
    list(
      fluidRow(
        column(width = 6,
               h4(style="font-family:Helvetica;", HTML("<br/> <b> Select phyloseq object </b>")),
               helpText(HTML("<br/> Only unfiltered phyloseq objects are available.")),
               selectInput(inputId = "data_filt", label = "", choices = n.datalist()[!grepl("filt_", n.datalist())])
        ),
        column(width = 6, style = 'border-left: 1px solid',
               h4(style="font-family:Helvetica;", HTML("<br/> <b> Default filtering </b>")),
               helpText(HTML("Remove chloroplasts, mitochondria and phylum taxa with NA annotation from the data.
               <br/> Optionally, remove undesired domains before any further filtration steps.")),
               div(class="col-md-2", radioButtons("default_filt", "", choices = c("Yes", "No"), selected = "Yes")),
               div(class="col-md-6", uiOutput("filt_domain"))                                           
        )
      ),
      fluidRow(
        column(width = 6,
               h4(style="font-family:Helvetica;", HTML("<br/> <b> Failed samples and low prevalence taxa filtering </b>")),
               div(class="col-md-3", numericInput(inputId = "min_depth", label = "Min. depth", value = 1000)),
               div(class="col-md-3", numericInput(inputId = "min_prev", label = "Min. prevalence", value = 2)),
               div(class="col-md-4", numericInput(inputId = "min_abund", label = "Min. abundance (CPM)", value = 10)),
               div(class="col-md-2", radioButtons(inputId = "abund_type", "", choices = c("Mean", "Total"), selected = "Mean"))
        ),
        column(width = 6, style = 'border-left: 1px solid',
               h4(style="font-family:Helvetica;", HTML("<br/> <b> Sample subsetting </b>")),
               div(class="col-md-5", selectInput(inputId = "filtvar", 
                                                 label = "Variable used to filter samples", 
                                                 choices = c("None", outvar()))),
               
               uiOutput('filtlev'),
        )
      ),
      h4(style="font-family:Helvetica;", HTML("<br/> <b> Taxa subsetting (sequencially done)</b>")),
      tags$style(HTML('#buttons {display: flex; align-items: center; justify-content: center;}')),
      fluidRow(
        column(width = 12, id="buttons",
               div(class="col-md-4", uiOutput('selected_rank')),
               div(class="col-md-4", uiOutput('selected_taxa')),
               div(class="col-md-4", uiOutput('keep_or_remove'))
        )
      ),
      actionButton("add_taxa_filter", "+ Add"),
      tags$head(tags$style(HTML('#add_taxa_filter{font-weight: bold}'))),
      actionButton("remove_taxa_filter", "- Remove"),
      tags$head(tags$style(HTML('#remove_taxa_filter{font-weight: bold}'))),
      br(),
      br(),
      fluidRow(
        column(width = 12,
               div(class="col-md-2", actionButton("filter", "Perform Filtration",
                                                  style="color: #fff; background-color: #7ab733; border-color: #527b22; margin-top: 25px;",
                                                  icon = icon("filter"))),
               div(class="col-md-2", textInput(inputId = "filtsuffix", label = "Suffix (optional)", value = ""))
        )
      ),
      br(),
      uiOutput('ui.physeqnames_filt'),
      uiOutput('ui.taxa_count')
    )
  })

  # Selected phyloseq object
  ps <- reactive({ datalist[[input$data_filt]] })
  
  # Display domains to remove before any further filtration steps
  output$filt_domain <- renderUI({
    selectInput("filt_domain", "Undesired Domains", 
                get_taxa_unique(ps(), colnames(ps()@tax_table)[1]), 
                selected = NULL, multiple = TRUE)
  })
  
  # Display filtering variable levels
  output$filtlev <- renderUI({
    ifelse(input$filtvar != "None", lev_choices <- unique(META()[,input$filtvar]), lev_choices <- "")
    list(
      div(class="col-md-3", selectInput(inputId = "keeplev", 
                                        label = paste0("Levels to retain"), 
                                        choices = lev_choices, 
                                        multiple = TRUE, 
                                        selected = NULL)),
      div(class="col-md-4", radioButtons("filtorder", label = "", 
                                         choices = c("Before next filtrations", "After next filtrations"), 
                                         selected = "Before next filtrations"))
    )
  })
  
  # Track the number of taxa filters added
  taxa_filt_counter <- reactiveValues(n = 1)
  observeEvent(input$add_taxa_filter, 
               { if (taxa_filt_counter$n >= 1) { taxa_filt_counter$n <- taxa_filt_counter$n + 1 }
               })
  observeEvent(input$remove_taxa_filter, 
               { if (taxa_filt_counter$n > 1) { taxa_filt_counter$n <- taxa_filt_counter$n - 1 }
               })
  
  # Display taxa ranks
  output$selected_rank <- renderUI({
    # Get value of "+ Add more" button, which represents number of times pressed
    inputsToShow <- taxa_filt_counter$n
    # Create rank list
    rank_list = list("None")
    # Populate the rank list
    lapply(1:inputsToShow, function(i){
      # Define unique input id
      selected_rankId <- paste0("selected_rank", i)
      selected_rank <- "None"
      if (selected_rankId %in% names(input)){
        selected_rank <- input[[selected_rankId]]
        # Define list of choices
        if( "data_filt" %in% names(input) ){
          rank_list <- c(rank_list, as.list(rank_names(ps(), errorIfNULL=FALSE)))
          rank_list <- c(rank_list, list(Rownames="Rownames"))
        } else {
          rank_list = list("None", "Upload data")
        }
      }
      # Define new input
      selectInput(selected_rankId, "Rank", rank_list, selected_rank)
    })
  })
  
  # Display available taxa per selected rank           
  output$selected_taxa <- renderUI({
    # Get value of "+ Add more" button, which represents number of times pressed
    inputsToShow <- taxa_filt_counter$n
    # Create taxa list
    taxa_list = list("")
    # Populate the taxa list
    lapply(1:inputsToShow,function(i){
      # Get selected rank (& avoid warning message)
      selected_rank <- "None"
      if (paste0("selected_rank", i) %in% names(input)){
        selected_rank <- input[[paste0("selected_rank", i)]]
      }
      # Define list of choices
      if (selected_rank != "None") {
        if(selected_rank == "Rownames"){
          taxa_list <- c(taxa_list, as.list(taxa_names(ps())))
        } else {
          taxa_list <- c(taxa_list, as.list(get_taxa_unique(ps(), selected_rank)))
        }
      }
      # Define unique input id and check for pre-existing selected taxa 
      selected_taxaId <- paste0("selected_taxa", i)
      newTaxa <- NULL
      if (selected_taxaId %in% names(input)){
        newTaxa <- input[[selected_taxaId]]
      }
      # Define new input
      selectInput(selected_taxaId, "Taxa", taxa_list, selected = newTaxa, multiple = TRUE)
    })
    })
  
  # Display keep or remove options
  output$keep_or_remove <- renderUI({
    # Get value of "+ Add more" button, which represents number of times pressed
    inputsToShow <- taxa_filt_counter$n
    # Initialize list of inputs and taxa list
    inputTagList <- tagList()
    # Populate the taxa list
    lapply(1:inputsToShow,function(i){
      # Define unique input id and label
      keep_or_removeId <- paste0("keep_or_remove", i)
      newChoice <- "Keep"
      if (keep_or_removeId %in% names(input)){
        newChoice <- input[[keep_or_removeId]]
      }
      # Define new input
      newInput <- radioButtons(keep_or_removeId, "", choices = c("Keep", "Remove"), selected = newChoice)
      # Append new input to list of existing inputs
      inputTagList <<- tagAppendChild(inputTagList, newInput)
    })
    # Return updated list of inputs
    inputTagList 
  })
  
  # Create list with filtered phyloseq objects 
  observeEvent(input$filter, 
               { # Scroll down the page once the plot is created
                 js$button()
                 # Filter taxa to retain
                 filt.ps <- orig.ps <- ps()
                 # Filter undesired domains
                 if(!is.null(input$filt_domain)){
                   oldMA <- as(tax_table(filt.ps), "matrix")
                   oldDF <- data.frame(oldMA)
                   newDF <- subset(oldDF, !eval(parse(text = colnames(ps()@tax_table)[1])) %in% input$filt_domain)
                   newMA <- as(newDF, "matrix")
                   tax_table(filt.ps) <- tax_table(newMA)
                 }
                 # Default filtering
                 if( input$default_filt == "Yes"){
                   filt.ps %<>% subset_taxa((Order !="Chloroplast") | is.na(Order))
                   filt.ps %<>% subset_taxa((Family != "Mitochondria") | is.na(Family))
                   filt.ps %<>% subset_taxa(!is.na(Phylum))
                 }
                 # Filter samples before any further filtrations
                 if(input$filtvar != "None" & !is.null(input$keeplev)){
                   if(input$filtorder == "Before next filtrations"){
                     sam.filt.ps <- data.frame(filt.ps@sam_data)
                     samples.kept <- rownames(sam.filt.ps[sam.filt.ps[,input$filtvar] %in% input$keeplev,])
                     filt.ps <- prune_samples(samples.kept, filt.ps)
                   }  
                 }
                 # Remove samples with less than 'min.depth' reads
                 filt.ps = prune_samples(sample_sums(filt.ps) >= input$min_depth, filt.ps)
                 # Filter taxa (species/ASVs/OTUs) with less than 'min.abund' reads
                 prevdf <- calc.prev(filt.ps)
                 filt.t <- (prevdf$Prevalence >= input$min_prev & prevdf[,input$abund_type] >= input$min_abund)
                 taxatokeep <- rownames(prevdf)[filt.t]
                 filt.ps <- prune_taxa(taxatokeep, filt.ps)
                 # Custom taxa filtration
                 if( input$selected_rank1 != "None" ){
                   keepTaxa = NULL
                   lapply(1:taxa_filt_counter$n, function(i){
                     rank <- input[[paste0("selected_rank", i)]]
                     taxa <- input[[paste0("selected_taxa", i)]]
                     keep_or_remove <- input[[paste0("keep_or_remove", i)]]
                     if (rank == "Rownames" & keep_or_remove == "Keep"){
                       keepTaxa = taxa
                     } else if (rank == "Rownames" & keep_or_remove == "Remove"){
                       keepTaxa = !taxa_names(filt.ps) %in% taxa
                     } else if (keep_or_remove == "Keep"){
                       TT = as(tax_table(filt.ps), "matrix")
                       keepTaxa = TT[, rank] %in% taxa 
                     } else if (keep_or_remove == "Remove"){
                       TT = as(tax_table(filt.ps), "matrix")
                       keepTaxa = !TT[, rank] %in% taxa
                     }
                     filt.ps <<- prune_taxa(keepTaxa, filt.ps)
                   })
                 }
                 # Filter samples after all previous filtrations
                 if(input$filtvar != "None" & !is.null(input$keeplev)){
                   if(input$filtorder == "After next filtrations"){
                     sam.filt.ps <- data.frame(filt.ps@sam_data)
                     samples.kept <- rownames(sam.filt.ps[sam.filt.ps[,input$filtvar] %in% input$keeplev,])
                     filt.ps <- prune_samples(samples.kept, filt.ps)
                   }  
                 }
                 filt.ps <- prune_taxa(taxa_sums(filt.ps) != 0, filt.ps) # Remove taxa with 0 total abundance
                 if (inherits(filt.ps, "phyloseq")){
                   filt.ps <- list(filt.ps)
                   new.name <- paste0("filt_", input$data_filt, input$filtsuffix)
                   if (new.name %in% names(datalist)){
                     rpt <- gsub(paste0("(", new.name,").*"), "\\1", names(datalist))
                     id <- table(rpt)[new.name]
                     names(filt.ps) <- paste0(new.name, ".", id)
                   } else {
                     names(filt.ps) <- new.name
                   }
                   datalist <<- c(filt.ps, datalist)
                   # To display created phyloseq objects
                   output$physeqnames_filt <- renderText({ 
                     filt.names <- names(datalist)[grepl("filt_", names(datalist))]
                     sapply(filt.names, function(x){ paste(x, "\n") })
                   })
                   output$ui.physeqnames_filt <- renderUI({
                     list(
                       HTML("<br/> <b> Filtered objects: </b>"),
                       verbatimTextOutput("physeqnames_filt")
                     )
                   })
                 }
                 # Display features per sample and taxa count histograms (unfiltered and filtered)
                 output$taxa_count <- renderPlot({
                   df <- as.data.frame(orig.ps@otu_table)
                   df[df > 0] <- 1
                   group.max.y = max(colSums(df)) + max(colSums(df)) * 0.1
                   Before.fil <- filtering_QC(orig.ps, "Before filtration", max.y = group.max.y)
                   After.fil <- filtering_QC(filt.ps[[1]], "After filtration", max.y = group.max.y)
                   Final.plot <- cowplot::plot_grid(Before.fil, After.fil, nrow = 1, ncol = 2, align = "hv", rel_widths = c(1, 1))
                   Final.plot
                 })
                 output$ui.taxa_count <- renderUI({ plotOutput(outputId = "taxa_count", height = "500px") })
               })
}, once = TRUE)
               
################################################################################
# Prevalence tab: prevalence plots 
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create, 
             { # Select phylo objects used for the prevalence plots
               output$physeq_prev <- renderUI({ 
                 # Keep previous choices
                 choice.unfilt <- choice.filt <- "None"
                 if ("data_prev.o" %in% names(input)){
                   choice.unfilt <- input$data_prev.o
                 }
                 if ("data_prev.f" %in% names(input)){
                   choice.filt <- input$data_prev.f
                 }
                 list(
                   # Unfiltered phylo objects
                   div(class="col-md-4", selectInput(inputId = "data_prev.o", 
                                                     label = "Unfiltered phyloseq object", 
                                                     choices = c("None", n.datalist()[!grepl("filt_", n.datalist())]),
                                                     selected = choice.unfilt)
                   ),
                   # Filtered phylo objects
                   div(class="col-md-4", selectInput(inputId = "data_prev.f", 
                                                     label = "Filtered phyloseq object", 
                                                     choices = c("None", n.datalist()[grepl("filt_", n.datalist())]),
                                                     selected = choice.filt)
                   )
                 )
               })
               
               # Get taxonomic ranks of unfiltered or filtered phyloseq object
               output$rank_create_prev <- renderUI({
                 if (req(input$data_prev.o) != "None" | req(input$data_prev.f) != "None"){
                   if (input$data_prev.o != "None"){
                     cols <- colnames(datalist[[input$data_prev.o]]@tax_table)
                   } else {
                     cols <- colnames(datalist[[input$data_prev.f]]@tax_table)
                   }
                   # Keep previous choice
                   choice.rank <- "None"
                   if ("rank_prev" %in% names(input)){
                     choice.rank <- input$rank_prev
                   }
                   list(
                     div(class="col-md-4", selectInput(inputId = "rank_prev", 
                                                       label = "Taxonomic rank", 
                                                       choices = c("None", cols),
                                                       selected = choice.rank)
                     )
                   )
                 }
               })
               
               # Prevalence plot
               prev <- reactive({
                 if (req(input$data_prev.o) != "None" | req(input$data_prev.f) != "None"){
                   if(req(input$rank_prev) != "None"){
                     show_modal_spinner(
                       spin = "fading-circle",
                       color = "purple",
                       text = "Generating plot..."
                     )
                     # Display download button when plot is created
                     output$down_prev <- renderUI({
                       dim_and_down("prev", 16, 16, ggplot = TRUE)
                     })
                     # Create plots
                     plot.list <- list(NULL)
                     if (input$data_prev.o != "None"){
                       ps <- datalist[[input$data_prev.o]]
                       prevdf <- calc.prev(ps)
                       # Original data prevalence plot
                       df.P = subset(prevdf, eval(as.symbol(input$rank_prev)) %in% get_taxa_unique(ps, input$rank_prev))
                       plot.list[[1]] <- ggplot(df.P, aes(Total, Prevalence / nsamples(ps), color = eval(as.symbol(input$rank_prev)))) +
                         facet_wrap(~ eval(as.symbol(input$rank_prev))) +
                         geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
                         geom_point(size = 2, alpha = 0.7) +
                         scale_x_log10() +
                         labs(title = paste0(input$rank_prev, " prevalence", 
                                             ifelse(grepl("filt_", input$data_prev.o), " after filtration.", " before filtration.")),
                              subtitle = ifelse((input$default_filt == "Yes" & grepl("filt_", input$data_prev.o)), 
                                                "No chloroplast, mitochondria and NA phyla.", ""),
                              x = "Total Abundance (CPM)",
                              y = "Prevalence [fraction] (the dotted line indicates a 0.05 prevalence)") +
                         theme(legend.position = "none")
                     }
                     # Filtered data prevalence plot
                     if (input$data_prev.f != "None"){
                       filt.ps <- datalist[[input$data_prev.f]]
                       filt.prevdf <- calc.prev(filt.ps)
                       df.P.filt = subset(filt.prevdf, eval(as.symbol(input$rank_prev)) %in% get_taxa_unique(filt.ps, input$rank_prev))
                       plot.list[[2]] <- ggplot(df.P.filt, aes(Total, Prevalence / nsamples(filt.ps), color = eval(as.symbol(input$rank_prev)))) +
                         facet_wrap(~ eval(as.symbol(input$rank_prev))) +
                         geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
                         geom_point(size = 2, alpha = 0.7) +
                         scale_x_log10() +
                         labs(title = paste0(input$rank_prev, " prevalence", 
                                             ifelse(grepl("filt_", input$data_prev.f), " after filtration.", " before filtration.")),
                              subtitle = ifelse((input$default_filt == "Yes" & grepl("filt_", input$data_prev.f)), 
                                                "No chloroplast, mitochondria and NA phyla.", ""),
                              x = "Total Abundance (CPM)",
                              y = "Prevalence [fraction] (the dotted line indicates a 0.05 prevalence)") +
                         theme(legend.position = "none")
                     }
                     plot <- ggpubr::ggarrange(plotlist = plot.list, nrow = length(plot.list), ncol = 1, heights = c(1, 0.85))
                     remove_modal_spinner()
                     return(plot)
                   } else { 
                     # Remove download button
                     output$down_prev <- renderUI({ NULL })
                     return("None") }
                 } else { 
                   # Remove download button
                   output$down_prev <- renderUI({ NULL })
                   return("None") }
               })
               
               # Display prevalence plot
               observeEvent(prev(), {
                 display(prev(), "prev", "Prevalence", ggplot = TRUE, shiny.height = "100%")
               })
             }, once = TRUE)