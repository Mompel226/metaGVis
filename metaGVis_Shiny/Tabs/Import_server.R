# Import tab content - SERVER
# Create reactive values for the file paths to be able to reset them
path <- reactiveValues(asv = NULL, taxa = NULL, biom = NULL, uncla = NULL, mphlan = NULL, humann = NULL)
observe({ path$asv <- input$asv$datapath[1] })
observe({ path$taxa <- input$taxa$datapath[1] })
observe({ path$biom <- input$biom$datapath[1] })
observe({ path$uncla <- input$uncla$datapath[1] })
observe({ path$mphlan <- input$mphlan$datapath[1] })
observe({ path$humann <- input$humann$datapath[1] })

# MOTHUR -- DADA2 -- Kraken2/Bracken files (having a reactive UI allows to reset the files!!)
output$input_files1 <- renderUI({
  input$reset1
  path$asv <- path$taxa <- path$biom <- NULL
  list(
    fileInput(inputId = "asv",
              label = HTML("Abundance table (Required)"),
              accept = c(".qza", ".csv", ".tsv", ".txt")),
    fileInput(inputId = "taxa",
              label = HTML("Taxonomic table (Required)"),
              accept = c(".qza", ".csv", ".tsv", ".txt")),
    h5(style="font-family:Helvetica; font-weight: bold; color:#404040;", "OR"),
    fileInput(inputId = "biom",
              label = HTML("BIOM file (Required)"),
              accept = c(".biom")),
    hr(),
    fileInput(inputId = "uncla",
              label = HTML("Unclassified file (Optional for Bracken)"),
              accept = c(".tsv", ".txt")),
    h6(style="margin-top:-20px; font-family:Helvetica; color:grey;", 
       "The Bracken report removes the 'Unclassified' taxon. 
                           To add this information extract it from the Kraken report 
                           and create a txt file with the the sample names in the first column 
                           and the corresponding 'Unclassified' values in the second column.")
  )
})
# MetaPhlAn file (having a reactive UI allows to reset the file!!)
output$input_files2 <- renderUI({
  input$reset2
  path$mphlan <- NULL
  fileInput(inputId = "mphlan",
            label = HTML("Merged tables file (Required)"),
            accept = c(".tsv", ".txt"))
})
# HUMAnN3 file (having a reactive UI allows to reset the file!!)
output$input_files3 <- renderUI({
  input$reset3
  path$humann <- NULL
  fileInput(inputId = "humann",
            label = HTML("Merged tables file (Required)"),
            accept = c(".tsv", ".txt"))
})
# To see sample names in the Import data tabs
outsamples <- eventReactive(list(path$asv, path$taxa, path$biom, path$mphlan, path$humann, input$str2rm),
                            { if(!is.null(path$biom)){
                                OTU <- import_biom(path$biom, parseFunction=parse_taxonomy_greengenes, parallel=TRUE, version=1.0)
                              } else if (!is.null(path$asv)){
                                ifelse(tools::file_ext(path$asv) == "csv", {OTU <- read.csv(path$asv, sep = ",")},
                                       ifelse(tools::file_ext(path$asv) == "tsv", {OTU <- read.csv(path$asv, sep = "\t")},
                                              ifelse(tools::file_ext(path$asv) == "txt", {OTU <- read.csv(path$asv, sep = input$sep)},
                                                     ifelse(tools::file_ext(path$asv) == "qza", {OTU <- qiime2R::read_qza(path$asv)}, ""))))
                              } else if (!is.null(path$mphlan)){
                                mphlan <- read.csv(path$mphlan, sep = "\t", strip.white = T, stringsAsFactors = F, skip = 1, row.names = 1)
                                mphlan <- mphlan[,-1]
                                colnames(mphlan) %<>% sub('_metaphlan3', "", .)
                              } else if (!is.null(path$humann)){
                                humann <- file2meco::humann2meco(abund_table = path$humann, db = "MetaCyc")
                                humann <- file2meco::meco2phyloseq(humann)
                                for (i in c('_Abundance-RPKs', '_Abundance', '_Coverage')){
                                  colnames(humann@otu_table) %<>% sub(i, "", .)
                                }
                              }
                              ifelse(!is.null(path$biom), samples <- colnames(otu_table(OTU)),
                                     ifelse(!is.null(path$asv), samples <- colnames(OTU)[-1], 
                                            ifelse(!is.null(path$mphlan), samples <- colnames(mphlan), 
                                                   ifelse(!is.null(path$humann), samples <- colnames(otu_table(humann)), ""))))
                              
                              ifelse(is.null(path$biom) & is.null(path$asv) & is.null(path$mphlan) & is.null(path$humann), 
                                     samples <- "To see the sample names upload the abundance table, 
         the BIOM file, the metaPhlAn3 file or the HUMAnN3 file", "")
                              if (!is.null(input$str2rm)){
                                for (i in unlist(strsplit(input$str2rm, ","))){
                                  samples <- sub(i, "", samples)
                                }
                              }
                              samples <- as.list(samples)
                              return(samples)
                            })
  
output$samples <- renderUI({ outsamples() })

# Display textInput that allows to change sample names
output$str2rm <- renderUI({
  if (!is.null(path$asv) | !is.null(path$biom) | !is.null(path$mphlan) | !is.null(path$humann)){
    textInput(inputId = "str2rm", label = "Patterns to remove from sample names (OPTIONAL)", value = NULL,
              placeholder = "Comma-separated patterns (see outcome above)")
  }
})

# Get metadata for internal use
META <- reactive({
  if(!is.null(input$meta)){
    path <- as.character(input$meta$datapath[1])
    ifelse(tools::file_ext(path) %in% c("xls", "xlsx", "xlsm"), {META  <- readxl::read_excel(path)},
           ifelse(tools::file_ext(path) == "csv", {META  <- read.csv(path, sep = ",")},
                  ifelse(tools::file_ext(path) == "tsv", {META  <- read.csv(path, sep = "\t")},
                         ifelse(tools::file_ext(path) == "txt", {META  <- read.csv(path, sep = input$sep)},""))))
    return(META)
  }})

# Display metadata variables
outvar <- reactive({
  if(!is.null(input$meta)){ coln <- colnames(META()) } else { coln <- "Upload metadata file" }
  return(coln)
})
observeEvent(META(),
             { output$variables <- renderUI({
               div(class="col-md-12", selectInput(inputId = "var", label = "Select categorical variables (REQUIRED)", 
                                                  choices = outvar(), multiple = TRUE))
             })
             })
# Display levels for each categorical variable selected
observeEvent(input$var,
             { output$levels <- renderUI({
               list(
                 h5(style="font-family:Helvetica; font-weight: bold; color:red;", "Order ALL levels for each selected variable (REQUIRED)"),
                 lapply(input$var, function(i) {
                   div(class="col-md-6", selectInput(inputId = i, label = paste0("Levels of ", i), 
                                                     choices = unique(META()[,i]), multiple = TRUE))
                 })
               )
             })
             })

# To create phyloseq object
datalist <- list()
observeEvent(input$create, 
             { # To associate levels order to variables
               show_modal_spinner(
                 spin = "fading-circle",
                 color = "purple",
                 text = HTML("Generating phyloseq object... <br/>
                             Please be patient.")
               )
               factor.lev <- list()
               for (i in input$var){ 
                 l <- list(input[[i]])
                 names(l) <- i
                 factor.lev <- c(factor.lev, l)
               }
               ifelse(input$str2rm!="", str2rm <- unlist(strsplit(input$str2rm, ",")), str2rm <- "")
               new.physeq <- try(Files2Phyloseq(BIOM.file = path$biom,
                                                OTU.file = path$asv,
                                                TAXA.file = path$taxa,
                                                META.file = input$meta$datapath[1],
                                                TREE.file = input$tree$datapath[1],
                                                unclass.file = path$uncla,
                                                mphlan.file = path$mphlan,
                                                humann.file = path$humann,
                                                sep = ifelse(!is.null(input$sep), input$sep, "\t"),
                                                str2rm = str2rm,
                                                factor.lev = factor.lev,
                                                filt.var = NULL,
                                                lev.kept = NULL))
               
               if (inherits(new.physeq, "phyloseq")){
                 new.physeq <- list(new.physeq)
                 if (input$name %in% names(datalist)){
                   rpt <- gsub(paste0("(", input$name,").*"), "\\1", names(datalist))
                   id <- table(rpt)[input$name]
                   names(new.physeq) <- paste0(input$name, ".", id)
                 } else {
                   names(new.physeq) <- input$name
                 }
                 datalist <<- c(new.physeq, datalist)
                 # To display created phyloseq objects as a list
                 output$physeqnames <- renderText({ 
                   unfilt.names <- names(datalist)[!grepl("filt_", names(datalist))]
                   sapply(unfilt.names, function(x) {paste(x,"\n")}) 
                   })
                 output$ui.physeqnames <- renderUI({
                   list(
                     HTML("<br/> <b> Created objects: </b>"),
                     verbatimTextOutput("physeqnames")
                   )
                 })
               }
               remove_modal_spinner()
             })
# Update phyloseq objects displayed in select input ui 
n.datalist <- eventReactive(input$create | input$filter, { names(datalist) })