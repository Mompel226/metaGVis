# Differential abundance tab content - SERVER
################################################################################
# Plot
################################################################################
# Reactive UI - reactive parameters depending on the phylo object
# ONLY when phyloseq object is created
observeEvent(input$create,
             { 
               # Display create plot button
               uiCreate("DA")
               # Reactive inputs 1
               output$inputs1_DA <- renderUI({
                        fluidRow(
                          column(width = 12,
                                 div(class="col-md-12",
                                   selectInput(inputId = "data_DA", label = "Phyloseq object", 
                                               choices = n.datalist(), width="100%"),
                                   selectInput(inputId = "type_DA", label = "Method used", 
                                               choices = c("LEfSe", "ANCOM-BC"), selected = "ANCOM-BC",
                                               width="100%")
                                 )
                          )
                        )
                 })
               
               # To allow variables to change when data is changed
               observeEvent(req(input$data_DA),
                            { # Get selected phyloseq object
                              ps_DA <<- reactive({ datalist[[input$data_DA]] })
                              
                              # Get categorical variables
                              var_DA <- getVar(input$data_DA)
                              
                              # Reactive inputs 2
                              output$inputs2_DA <- renderUI({
                                list(
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", selectInput(inputId = "var_DA", label = "Group", 
                                                                             choices = var_DA())),
                                           div(class="col-md-6", selectInput(inputId = "subvar_DA", label = "Subgroup", 
                                                                             choices = c("", var_DA()), selected = ""))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-12", 
                                               selectInput(inputId = "rank_DA", label = "Rank", 
                                                           choices = as.list(rank_names(ps_DA(), errorIfNULL=FALSE)),
                                                           width = "100%"))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", numericInput(inputId = "tit.size_DA", label = "Title size", value = 12, step = 1)),
                                           div(class="col-md-6", numericInput(inputId = "atit.size_DA", label = "Axis titles size", value = 12, step = 1))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", numericInput(inputId = "x.size_DA", label = "X axis size", value = 12, step = 1)),
                                           div(class="col-md-6", numericInput(inputId = "y.size_DA", label = "Y axis size", value = 12, step = 1))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", 
                                               selectInput(inputId = "lpos_DA", label = "Legend position", 
                                                           choices = c("bottom", "right"),
                                                           selected = "bottom")),
                                           div(class="col-md-6", 
                                               numericInput(inputId = "lsize_DA", label = "Legend expansion", value = 0.1, step = 0.1))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-6", numericInput(inputId = "ltit.size_DA", label = "Legend title size", value = 14, step = 1)),
                                           div(class="col-md-6", numericInput(inputId = "ltxt.size_DA", label = "Legend text size", value = 12, step = 1))
                                    )
                                  ),
                                  fluidRow(
                                    column(width = 12,
                                           div(class="col-md-12", uitheme("theme_DA", width = "100%"))
                                    )
                                  )
                                )
                              })
                            })
               
               # Reactive inputs 3
               output$inputs3_DA <- renderUI({
                 list(
                   column(6,
                          fluidRow(
                            column(12,
                                   div(class="col-md-6", 
                                       numericInput(inputId = "ntaxa_DA", 
                                                    label = "Maximum number of taxa to show per group", 
                                                    value = 20, step = 1, min = 1)),
                                   div(class="col-md-6", uiOutput("plot_DA"))
                            )
                          ),
                          uiOutput("ancomopt_DA"),
                          uiOutput("lefseopt_DA")
                   ),
                   column(3,
                          ColOptions("DA"),
                          fluidRow(
                            column(12,
                                   div(class="col-md-12", uiOutput("col_DA")), 
                            )
                          )
                   )
                 )
               })
               
               # Display dynamic UIs
               observeEvent(input$type_DA, 
                            {# Display options available for ANCOM-BC analysis
                              if(input$type_DA == "ANCOM-BC"){
                                output$ancomopt_DA <- renderUI({
                                  list(
                                    div(HTML("<b> ANCOM-BC options </b>")),
                                    hr(),
                                    fluidRow(
                                      column(12,
                                             div(class="col-md-6", 
                                                 radioButtons(inputId = "focus_DA",
                                                              label = "Focus on taxa differentialy abundant against several groups",
                                                              choices = c("Yes", "No"), 
                                                              selected = "Yes", inline = T)),
                                             div(class="col-md-6",
                                                 radioButtons(inputId = "showlefse_DA", 
                                                              label = "Prioritize taxa matching LEfSe results (shown in bold)",
                                                              choices = c("Yes", "No"),
                                                              selected = "Yes", inline = T))
                                      ),
                                      column(12,
                                             div(class="col-md-6", NULL),
                                             div(class="col-md-6",
                                                 radioButtons(inputId = "showlefse_only", 
                                                              label = "Show only LEfSe matches",
                                                              choices = c("Yes", "No"),
                                                              selected = "No", inline = T))
                                      )
                                    )
                                  )
                                })
                              } else {
                                output$ancomopt_DA <- renderUI({ NULL })
                              }
                            })
               observeEvent(list(req(input$type_DA, req(input$showlefse_DA))), 
                            {# Display options available for LEfSe analysis
                              if(input$type_DA == "LEfSe" | input$showlefse_DA == "Yes"){
                                output$lefseopt_DA <- renderUI({
                                  list(
                                    div(HTML("<b> LEfSe options </b>")),
                                    hr(),
                                    fluidRow(
                                      column(12,
                                             div(class="col-md-6", 
                                                 numericInput(inputId = "LDA_DA", 
                                                              label = "Threshold on the logarithmic LDA score", 
                                                              value = 2, step = 0.5, min = 0)),
                                             div(class="col-md-6", 
                                                 numericInput(inputId = "boot_DA", 
                                                              label = "Number of bootstrap iteration for LDA", 
                                                              value = 30, step = 1, min = 1))
                                      )
                                    ),
                                    fluidRow(
                                      column(12,
                                             div(class="col-md-6", 
                                                 numericInput(inputId = "KW_DA", 
                                                              label = "p value cutoff of Kruskal-Wallis test among groups", 
                                                              value = 0.01, step = 0.01, max = 0.1, min = 0.001)),
                                             div(class="col-md-6", 
                                                 numericInput(inputId = "W_DA", 
                                                              label = "p value cutoff of pairwise Wilcoxon test between subgroups", 
                                                              value = 0.01, step = 0.01, max = 0.1, min = 0.001))
                                             )
                                      ),
                                      fluidRow(
                                        column(12,
                                               div(class="col-md-6", 
                                                   radioButtons(inputId = "same_DA",
                                                                label = "Perform Wilcoxon test only among the subgroups with the same name", 
                                                                choices = c("Yes", "No"), selected = "No")),
                                               div(class="col-md-6", 
                                                   numericInput(inputId = "min_DA", 
                                                                label = "Min nº of samples per subgroup to perform Wilcoxon test", 
                                                                value = 10, step = 1, min = 1))
                                               )
                                        ),
                                    fluidRow(
                                      column(12,
                                             div(class="col-md-12", 
                                                 radioButtons(inputId = "multi_DA",
                                                          label = "Setting for a group with more than 2 levels",
                                                          choices = c("One-against one (more strict)", "One-against all (less strict)"), 
                                                          selected = "One-against one (more strict)", width="100%"))
                                      )
                                    )
                                  )
                                })
                              } else {
                                output$lefseopt_DA <- renderUI({ NULL })
                              }
                            })
               
               # Display plots available per type of analysis
               output$plot_DA <- renderUI({ 
                 if(is.null(ps_DA()@phy_tree)){
                   lefse.plot <- c("LEfSe barplot")
                 } else {
                   lefse.plot <- c("LEfSe barplot", "LEfSe cladogram")
                 }
                 if (input$type_DA == "ANCOM-BC"){
                   choices = c(lefse.plot, "All bars", "Significant bars", "Positive bars")
                 } else {
                   choices = lefse.plot
                 }
                 selectInput(inputId = "plot_DA", label = "Type of plot", 
                             choices = choices,
                             selected = "Positive bars")
               })
               
               # Create color vector depending on the color option selected
               observeEvent(list(input$colset_DA, input$custom_col_DA),
                            {# Display color and fill selection for grouping/color variable levels
                              if(isTRUE(input$custom_col_DA)){
                                uiColLev("col_DA", "var_DA")
                                col_DA <<- ColVec("col_DA", "var_DA")
                              } else {
                                output$col_DA <- renderUI({ NULL })
                                col_DA <<- reactive({ 
                                  RColorBrewer::brewer.pal(nrow(unique(ps_DA()@sam_data[,input$var_DA])), input$colset_DA) })
                              }
                            })
               # Use the color vector provided
               observeEvent(req(input$hex_col_DA),
                            {# Make sure text input are HEX codes
                              col.vect <- strsplit(input$hex_col_DA, ",")[[1]]
                              lev.length <- length(levels(get_variable(ps_DA(), input$var_DA)))
                              if (length(col.vect) == lev.length & all(grepl("^#[[:xdigit:]]{6}$", col.vect))){
                                col_DA <<- reactive({ col.vect })
                              } else {
                                col_DA <<- ColVec("col_DA", "var_DA")
                              }
                            })
               
               # Named col vector to display all samples in legend
               namedcol_DA <- reactive({ 
                 list(input$colset_DA, input$custom_col_DA, input$hex_col_DA)
                 a.col <- alpha(col_DA(), input$alpha_DA)
                 lev <- levels(get_variable(ps_DA(), input$var_DA))
                 names(a.col) <- lev
                 a.col <- a.col[1:length(lev)]
                 return(a.col)
                 })
               
               # Differential abundance plot
               DA <- eventReactive(input$create_DA, {
                 show_modal_spinner(
                   spin = "fading-circle",
                   color = "purple",
                   text = HTML("Performing LEfSe analysis...<br/>
                   This analysis takes a long time. <br/>
                    Please be patient.")
                 )
                 # LEfSe analysis
                 group.levels <<- levels(get_variable(ps_DA(), input$var_DA))
                 agg.data <- metaGVis::unclassified_ranks(ps_DA(), input$rank_DA, format = "short", agglomerate = TRUE, 
                                                          count = TRUE, output = "physeq")
                 agg.data@sam_data[,input$var_DA] <- factor(get_variable(ps_DA(), input$var_DA), levels = group.levels)
                 colnames(agg.data@tax_table)[1] <- "Kingdom"
                 lefse.out <- NULL
                 try(
                   lefse.out <- microbiomeMarker::run_lefse(agg.data, 
                                                            group = input$var_DA,
                                                            subgroup = switch((input$subvar_DA=="") + 1, input$subvar_DA, NULL),
                                                            norm = "CPM", 
                                                            taxa_rank = switch((input$rank_DA=="Domain") + 1, input$rank_DA, "Kingdom"), 
                                                            wilcoxon_cutoff = input$W_DA, 
                                                            kw_cutoff = input$KW_DA, 
                                                            multigrp_strat = switch((input$multi_DA!="One-against one (more strict)") + 1, TRUE, FALSE), 
                                                            lda_cutoff = input$LDA_DA, 
                                                            bootstrap_n = input$boot_DA,
                                                            sample_min = input$min_DA,
                                                            only_same_subgrp = switch((input$same_DA=="No") + 1, TRUE, FALSE)),
                   silent = TRUE
                 )
                 remove_modal_spinner()
                 if(!is.null(lefse.out)){
                   if(!is.null(lefse.out@marker_table)){
                     lefse.out <- data.frame(lefse.out@marker_table)
                     lefse.out%<>%transmute(taxon = feature, direction = enrich_group, lda = ef_lda, p.value = padj)
                     
                     # Order levels
                     lefse.out$direction <- factor(lefse.out$direction, levels = levels(get_variable(ps_DA(), input$var_DA)))
                     lefse.out <- with(lefse.out, lefse.out[rev(order(lda)),])
                     lefse.out <- with(lefse.out, lefse.out[order(direction),])
                     lefse.out$taxon <- factor(lefse.out$taxon, levels = rev(lefse.out$taxon))
                   } else {
                     lefse.out <- NULL
                   }
                 } else {
                   lefse.out <- NULL
                 }
                 if (input$type_DA == "ANCOM-BC"){
                   show_modal_spinner(
                     spin = "fading-circle",
                     color = "purple",
                     text = HTML("Performing ANCOM-BC analysis...<br/>
                   This analysis takes a long time. <br/>
                    Please be patient.")
                   )
                   # ANCOM-BC analysis
                   ANCOM.rank <<- input$rank_DA
                   ANCOM.data.list <- list()
                   try({
                     for (i in group.levels){
                         
                         print(paste0("Performing ANCOM-BC using ", i, " as reference."))
                         
                         # Re-arranging levels order to set new reference (first value is the reference)
                         new.levels.order <- i
                         temp.levels.order <- group.levels[group.levels != i]
                         new.levels.order <- c(new.levels.order, temp.levels.order)
                         agg.data@sam_data[,input$var_DA] <- factor(get_variable(ps_DA(), input$var_DA), 
                                                                    levels = new.levels.order)
                         
                         # ANCOM-BC analysis (no structural zeros are detected, otherwise all taxa are significant)
                         out <- ANCOMBC::ancombc(agg.data, formula = input$var_DA, p_adj_method = "holm", 
                                                 zero_cut = 0.9, lib_cut = 0, group = input$var_DA, 
                                                 struc_zero = FALSE, neg_lb = FALSE, tol = 1e-05, max_iter = 100, 
                                                 conserve = TRUE, alpha = 0.05, global = FALSE)
                         print(paste0("Completed"))
                         # Extract the relevant information from the ANCOM-BC output and merge it into a new data frame (dat.fig)
                         # Structural zeros are not detected (for our analysis) but the code to extract the information is written.
                         # W 
                         W <- out$res$W
                         colnames(W) <- sub(input$var_DA, "", colnames(W))
                         rank.names <- as.data.frame(unclassified_ranks(agg.data, input$rank_DA, output = "taxnames"))
                         rank.names$ID <- rownames(rank.names)
                         sub.rank.names <- rank.names[rank.names$ID %in% rownames(W),]
                         sub.rank.names <- sub.rank.names[-2]
                         colnames(sub.rank.names) <- "Taxa"
                         W <- merge(W, sub.rank.names, by="row.names", all=TRUE)
                         rownames(W) <- W[[1]]
                         W <- W[-1]
                         W <- reshape2::melt(W, value.name = "W", id.vars = "Taxa")
                         colnames(W) <- c("Taxa", "Group", "W")
                         # q_value
                         q_val <- out$res$q_val
                         colnames(q_val) <- sub(input$var_DA, "", colnames(q_val))
                         q_val <- merge(q_val, sub.rank.names, by="row.names", all=TRUE)
                         rownames(q_val) <- q_val[[1]]
                         q_val <- q_val[-1]
                         q_val <- reshape2::melt(q_val, value.name = "q_val", id.vars = "Taxa")
                         colnames(q_val) <- c("Taxa", "Group", "q_val")
                         # se
                         se <- out$res$se
                         colnames(se) <- sub(input$var_DA, "", colnames(se))
                         se <- merge(se, sub.rank.names, by="row.names", all=TRUE)
                         rownames(se) <- se[[1]]
                         se <- se[-1]
                         se <- reshape2::melt(se, value.name = "se", id.vars = "Taxa")
                         colnames(se) <- c("Taxa", "Group", "se")
                         # Merge data.frames
                         dat.fig=cbind(W, se, q_val)
                         dat.fig=dat.fig[,-(4:5)]
                         dat.fig=dat.fig[,-(5:6)]
                         alpha.adj=0.05/nrow(out$res$se)
                         critic.val=qnorm(1-alpha.adj/2)
                         dat.fig=dat.fig%>%transmute(Taxa,
                                                     Group,
                                                     W = -W,
                                                     se,
                                                     q_val,
                                                     ci.low=W-critic.val*se, 
                                                     ci.up=W+critic.val*se, 
                                                     star=ifelse(q_val<.001, "***", 
                                                                 ifelse(q_val<.01, "**",
                                                                        ifelse(q_val<.05, "*", ""))))
                         # List data
                         dat.fig$Reference <- rep(paste0(i), nrow(dat.fig))
                         ANCOM.data.list[[i]] <- dat.fig
                     }
                     },
                     silent = TRUE)
                   remove_modal_spinner()
                   if (is_empty(ANCOM.data.list)){
                     ANCOM.out <- NULL
                   } else {
                     print("Cleaning ANCOM-BC results...")
                     ANCOM.out <- bind_rows(ANCOM.data.list)
                     # Remove taxa not present (ANCOM fails to do so and creates false data)
                     agg.data_melted <- reshape2::melt((otu_table(agg.data)))
                     agg.data_melted <- merge(agg.data_melted, sample_data(agg.data), by.x = "Var2", by.y = "row.names")
                     agg.data_agg <- stats::aggregate(as.formula(paste("value ~ Var1 +", input$var_DA)),
                                                      data = agg.data_melted,
                                                      function(x) (sum(x > 0)/length(x) >= 0) * mean(x))
                     agg.data_mat <- reshape2::dcast(as.formula(paste("Var1 ~ ", input$var_DA)), 
                                                     data = agg.data_agg, 
                                                     value.var = "value")
                     rownames(agg.data_mat) <- agg.data_mat[, 1]
                     agg.data_mat <- agg.data_mat[, -1]
                     df <- merge(agg.data_mat, tax_table(agg.data), by=0)
                     for (g in group.levels){
                       for (i in df[,rank]){
                         row.data = df[df[,rank] == i,]
                         if (row.data[,g] == 0){
                           RowsToRm <- rownames(ANCOM.out[ANCOM.out$Taxa == i & ANCOM.out$Reference == g,])
                         } else {
                           RowsToRm <- ""
                         }
                         ANCOM.out %<>% dplyr::filter(!row.names(ANCOM.out) %in% RowsToRm)
                       }
                     }
                     print("Done")
                   }
                 } else {
                   ANCOM.out <- NULL
                 }
                 results <- list(lefse.out, ANCOM.out)
                 return(results)
               })
               
               # Reactive plots
               DA.plot <- reactive({
                 if (req(input$create_DA) > 0 & !is.null(DA()[[1]])){
                   list(input$colset_DA, input$custom_col_DA, input$hex_col_DA)
                   if(input$plot_DA == "LEfSe barplot"){
                     # Create LEfSe plots
                     show_modal_spinner(
                       spin = "fading-circle",
                       color = "purple",
                       text = HTML("Generating plot...")
                     )
                     data <- DA()[[1]] %>% group_by(direction) %>% slice_max(order_by = lda, n = input$ntaxa_DA)
                     p<-ggplot(subset(data,lda >= input$LDA_DA),
                               aes(x=taxon, y=lda, group=direction)) + 
                       geom_bar(aes(fill=direction), stat="identity", width=0.7, position="dodge")+
                       geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5)+
                       labs(title = "LEfSe scores", x=NULL, y="LDA score (log10)") + coord_flip()+
                       scale_fill_manual(name=paste0("Group:"), values = alpha(col_DA(), input$alpha_DA)) +
                       shiny_phyloseq_ggtheme_list[[input$theme_DA]] +
                       theme(plot.title = element_text(hjust = 0.5, size = input$tit.size_DA),
                             axis.text.x = element_text(size = input$x.size_DA),
                             axis.text.y = element_text(size = input$y.size_DA),
                             axis.title = element_text(size = input$atit.size_DA),
                             legend.title = element_text(size = input$ltit.size_DA), 
                             legend.text = element_text(size = input$ltxt.size_DA),
                             panel.grid.major = element_line(),
                             panel.grid.minor = element_blank(),
                             strip.background = element_rect(fill="white"))
                     remove_modal_spinner()
                     return(p)
                   } else if(!is.null(DA()[[2]])){
                     # Create ANCOM-BC plots
                     show_modal_spinner(
                       spin = "fading-circle",
                       color = "purple",
                       text = HTML("Generating plot...")
                     )
                     Plots <- list()
                     for (l in group.levels){
                       n1 <- filter(DA()[[2]], Reference == l) # Filter individual reference
                       n1$Taxa <- gsub("Unidentified ", "NA", n1$Taxa) # Remove "Unidentified ..." prefix
                       n1.pos <- n1%>%filter(W > 0,) # Filter taxa with positive W values
                       # Aggregate W values of each taxa
                       n1.agg <- aggregate(n1.pos$W, by=list(Category=n1.pos$Taxa), FUN=mean) 
                       n2 <- n1%>%filter(Taxa %in% unique(n1.agg$Category))
                       n3 <- n2%>%filter(W > 0,) # Select positive W values (for plots)
                       n3 <- n3%>%filter(star != "") # Select taxa that are significatively DE
                       n3 <- with(n3, n3[rev(order(W)),]) # Order W values from bigger to smaller
                       # Give priority to taxa that are differentially abundant against several groups
                       if (input$focus_DA == "Yes" & input$showlefse_DA == "No"){
                         maxn <- length(levels(n3$Group))-1
                         new.order <- vector()
                         for(k in maxn:1){
                           x <- n3 %>% count(Taxa) %>% filter(n == k)
                           new.order <- c(new.order, as.vector(x$Taxa))
                         }
                         n3 <- n3%>%filter(Taxa %in% new.order[1:input$ntaxa_DA])
                       }
                       # Create a vector with the taxa also identified by LEfSe
                       if (input$showlefse_DA == "Yes" & !is.null(DA()[[1]])){
                         Lefse.taxa <- as.vector(as.matrix(DA()[[1]]%>%
                                                             filter(direction == l)%>%
                                                             transmute(Lefse = taxon)))
                         Lefse.taxa <- gsub("Unidentified ", "NA", Lefse.taxa)
                         Lefse.taxa <- n3%>%filter(Taxa %in% Lefse.taxa)
                         Lefse.taxa <- as.vector(unique(Lefse.taxa$Taxa))
                         # Keep taxa names also found in LEfSe
                         not.in.Lefse <- n3%>%filter(!Taxa %in% Lefse.taxa)
                         not.in.Lefse$Taxa <- factor(not.in.Lefse$Taxa, levels = unique(not.in.Lefse$Taxa))
                         ifelse(length(Lefse.taxa) > input$ntaxa_DA, minus <- input$ntaxa_DA, minus <- length(Lefse.taxa))
                         if (input$focus_DA == "Yes" & length(Lefse.taxa) < input$ntaxa_DA){
                           ancom.taxa <- input$ntaxa_DA - length(Lefse.taxa)
                           maxn <- length(levels(n3$Group))-1
                           new.order <- vector()
                           for(k in c(maxn:1)){
                             x <- not.in.Lefse %>% count(Taxa) %>% filter(n == k)
                             new.order <- c(new.order, as.vector(x$Taxa))
                           }
                           Lefse.Ancom.taxa <- c(Lefse.taxa, new.order[1:ancom.taxa])
                         } else {
                           Lefse.Ancom.taxa <- c(Lefse.taxa, levels(not.in.Lefse$Taxa)[1:(input$ntaxa_DA-minus)])
                         }
                         n3 <- n3%>%filter(Taxa %in% Lefse.Ancom.taxa)
                       }
                       # Remove "Unidentified ..." prefix to make names shorter
                       if (input$showlefse_DA == "Yes" & !is.null(DA()[[1]])){
                         Lefse.taxa <- gsub("Unidentified ", "NA", Lefse.taxa)
                       }
                       # Limit the number of taxa shown
                       n2$Taxa <- factor(n2$Taxa, levels = unique(n2$Taxa)) # Transform taxa names to levels
                       if (length(levels(n2$Taxa)) > input$ntaxa_DA){
                         n2 <- n2%>%filter(Taxa %in% levels(n2$Taxa)[1:input$ntaxa_DA])
                         n2$Taxa <- factor(n2$Taxa, levels = unique(n2$Taxa)) # Re-transform taxa names to levels
                       }
                       n3$Taxa <- factor(n3$Taxa, levels = unique(n3$Taxa)) # Transform taxa names to levels
                       if (length(levels(n3$Taxa)) > input$ntaxa_DA){
                         n3 <- n3%>%filter(Taxa %in% levels(n3$Taxa)[1:input$ntaxa_DA])
                         n3$Taxa <- factor(n3$Taxa, levels = unique(n3$Taxa)) # Re-transform taxa names to levels
                       }
                       # Reverse Group factors for a better plotting
                       n2$Group <- factor(n2$Group, levels = rev(group.levels))
                       n3$Group <- factor(n3$Group, levels = rev(group.levels))
                       
                       # Display only LEfSe taxa
                       if (input$showlefse_only == "Yes"){
                         filt.taxa <- Lefse.taxa[1:input$ntaxa_DA]
                         n2 <- n2%>%filter(Taxa %in% filt.taxa)
                         n2$Taxa <- factor(n2$Taxa, levels = filt.taxa) # Re-transform taxa names to levels
                         n3 <- n3%>%filter(Taxa %in% filt.taxa)
                         n3$Taxa <- factor(n3$Taxa, levels = filt.taxa) # Re-transform taxa names to levels
                       }
                       
                       # Bold names for taxa names also found in LEfSe
                       if (input$showlefse_DA == "Yes" & !is.null(DA()[[1]])){
                         bold.labels <- ifelse(rev(levels(n3$Taxa)) %in% Lefse.taxa, yes = "bold", no = "plain")
                       } else {
                         bold.labels <- NULL
                       }
                       # Plot with bars for each group and stars showing significance
                       if(input$plot_DA == "All bars"){
                         data <- n2
                         width <- 1
                       } else if(input$plot_DA == "Significant bars"){
                         data <- n2%>%filter(star != "")
                         width <- 0.7
                       } else if(input$plot_DA == "Positive bars"){
                         data <- n3
                         width <- 0.7
                       }
                       if(all(names(namedcol_DA()) %in% group.levels)){
                         p<-ggplot(data, 
                                   aes(x=Taxa, y=W, ymin=ci.low, ymax=ci.up, group=Group)) + 
                           geom_bar(aes(fill=Group), stat="identity", width=width, position="dodge") +
                           geom_errorbar(width=0.2, size=0.25, position=position_dodge(width=width)) +
                           geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
                           labs(title = paste0("Ref: ", l, " - Rank: ", input$rank_DA), 
                                x = NULL, y = "Log Fold Change") + 
                           coord_flip() +
                           scale_fill_manual(name="Group:", values = namedcol_DA()) +
                           scale_x_discrete(limits = rev(levels(data$Taxa))) +
                           shiny_phyloseq_ggtheme_list[[input$theme_DA]] +
                           theme(plot.title = element_text(hjust = 0.5, size = input$tit.size_DA),
                                 axis.text.x = element_text(size = input$x.size_DA),
                                 axis.text.y = ggtext::element_markdown(size = input$y.size_DA, face = bold.labels),
                                 axis.title = element_text(size = input$atit.size_DA),
                                 legend.position=input$lpos_DA,
                                 legend.title = element_text(size = input$ltit.size_DA), 
                                 legend.text = element_text(size = input$ltxt.size_DA),
                                 panel.grid.major = element_line(),
                                 panel.grid.minor = element_blank(),
                                 strip.background = element_rect(fill="white"))
                         
                         legend <- cowplot::get_legend(p)
                         p <- p + theme(legend.position = "none")
                         if(input$plot_DA == "All bars"){
                           p <- p +
                             geom_text(aes(y=ci.up + max(W)*0.07, label=star), vjust=.75, size = 3, 
                                       color="black", position=position_dodge(width = width))
                         } else {
                           p <- p +
                             geom_text(aes(y=ci.up + sign(W)*1.5, label=star), vjust=.7, size = 5, 
                                       color="black", position=position_dodge(width = width))
                         }
                         Plots[[l]] <- p
                       }
                     }

                     if(all(names(namedcol_DA()) %in% group.levels) & ANCOM.rank == input$rank_DA){
                       p <- ggpubr::ggarrange(plotlist = Plots, align = "hv")
                       if (input$lpos_DA == "bottom"){
                         p <- ggpubr::ggarrange(p, legend, ncol = 1, nrow = 2, heights = c(1,input$lsize_DA))
                       } else {
                         p <- ggpubr::ggarrange(p, legend, ncol = 2, nrow = 1, widths = c(1,input$lsize_DA))
                       }
                     } else {
                       p <- ggplot() + theme_void()
                     }
                     remove_modal_spinner()
                     return(p)
                   }
                 }
               })
            
               observeEvent(DA.plot(), {
                 # Display download button when plot is created
                 output$down_DA <- renderUI({
                   coef <- input$dimension[2]/input$dimension[1]
                   dim_and_down("DA", 
                                round(2*coef*(input$dimension[2]/42),0), 
                                round(coef*(input$dimension[2]/50),0), ggplot = TRUE)
                 })
                 # Display DA plot
                 display(DA.plot(), "DA", input$type_DA, ggplot = TRUE, shiny.height = "100%")
               })
             }, once = TRUE)