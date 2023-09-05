# # Globally defined lists and functions
################################################################################
# Lines of code re-used from the package shiny-phyloseq
# Original code: https://github.com/joey711/shiny-phyloseq/blob/master/ui.R

# List of ggplot themes 
shiny_phyloseq_ggtheme_list <- list(
  bl_wh = theme_bw(),
  thin = theme_linedraw(),
  light = theme_light(),
  minimal = theme_minimal(),
  classic = theme_classic(),
  gray = theme_gray()
)
# To specify the format of downloaded files
graphicTypeUI <- function(inputId, label="Format", choices=graphicFormats, selected="png"){
  selectInput(inputId, label, choices, selected, multiple = FALSE, selectize = TRUE)
}

# UI function to define ggplot2 themes. Reused in many plots.
uitheme <- function(id, default="bl_wh", width = NULL){
  selectInput(id, "Theme",
              choices = names(shiny_phyloseq_ggtheme_list),
              selected = default, width = width
  )
}
################################################################################
# Original code
################################################################################
# Display create plot button
uiCreate <- function(suffix){
  output[[paste0("create_", suffix)]] <- renderUI({
    div(class='col-md-3', style = "margin-top: 25px;", 
        actionButton(paste0("create_", suffix), "Create/update plot!",
                     style="color: #fff; background-color: #7ab733; border-color: #527b22",
                     icon = icon("pencil-alt")))
  })
}

# Get defined categorical variables
getVar <- function(data){
  reactive({
    sam.data <- datalist[[data]]@sam_data
    if (!is.null(sam.data)){
      if (any(sapply(sam.data, is.factor))){
        variables <- colnames(sam.data[,sapply(sam.data, is.factor)])
        return(variables)
      } else {
        return("Select categorical variables")
      }
    }
    return("Upload Metadata file")
  })
}

# Display color options (alpha, set, custom)
ColOptions <- function(Id, label = "Fill"){
  list(
    fluidRow(
      column(width = 12,
             div(class="col-md-6", numericInput(inputId = paste0("alpha_",Id), label = paste0(label, " alpha"), 
                                                value = 0.8, min = 0, max = 1, step = 0.1)),
             div(class="col-md-6", selectInput(inputId = paste0("colset_",Id), label = paste0(label, " set"),
                                               choices = c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"),
                                               selected = "Set2"),
                 div(style = "margin-top: -17px"),
                 checkboxInput(inputId = paste0("custom_col_", Id), label = "Customise", value = FALSE)
             )
      )
    )
  )
}
# Get a color picker UI for each level of a defined categorical variable 
uiColLev <- function(Id, variable, title = "Custom colors"){
  if (input[[variable]] != "None"){
    output[[Id]] <- renderUI({
      suffix <- tail(unlist(strsplit(Id, "_")), n=1)
      lev <- levels(get_variable(datalist[[input[[paste0("data_", suffix)]]]], input[[variable]]))
      list(
        div(HTML(paste0("<b>", title, "</b>"))),
        hr(),
        fluidRow(
          column(width = 12,
                 div(class="col-md-12",
                     textInput(inputId = paste0("hex_",Id), label = "HEX color codes", 
                               value = "Comma-separated",
                               width = "100%"))
          )
        ),
        lapply(seq_along(lev), function(i) {
          div(class="col-md-6", colourpicker::colourInput(inputId = paste0(Id, i),
                                                          label = paste0(lev[i]))
          )
        })
      )
    })
  } else {
    output[[Id]] <- renderUI({ NULL })
  }
}

# Create a vector with the chosen colors using the color picker UIs
ColVec <- function(Id, variable){
  reactive({
    if(input[[variable]] != "None"){
      suffix <- tail(unlist(strsplit(Id, "_")), n=1)
      lev <- levels(get_variable(datalist[[input[[paste0("data_", suffix)]]]], input[[variable]]))
      cols <- vector(length = length(lev))
      for (i in seq_along(lev)){
        cols[i] <- input[[paste0(Id, i)]]
      }
      return(cols)
    } else {
      return(NULL)
    }
  })
}

# UI to download a plot
dim_and_down <- function(suffix, width = 0, height = 0, ggplot = TRUE){
  if (isTRUE(ggplot)) {
    list(
      div(class='col-md-3', style = "margin-top: 25px;", downloadButton(paste0("download_", suffix), 'Download plot!', 
                                                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
      div(class="col-md-2", numericInput(paste0("width_", suffix), "Plot width", value = width, min = 0, max = 100, step = 1)),
      div(class="col-md-2", numericInput(paste0("height_", suffix), "Plot height", value = height, min = 0, max = 100, step = 1)),
      div(class='col-md-2', graphicTypeUI(paste0("downtype_", suffix)))
    )
  } else {
    list(
      div(class='col-md-3', style = "margin-top: 25px;", downloadButton(paste0("download_", suffix), 'Download plot!', 
                                                                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
      div(class="col-md-2", numericInput(paste0("width_", suffix), "Plot width", value = width, min = 1, max = 100, step = 1)),
      div(class="col-md-2", numericInput(paste0("height_", suffix), "Plot height", value = height, min = 1, max = 100, step = 1))
    )
  }
}

# Server code to display and dowload a plot
display <- function(plot, suffix, file.name, ggplot = TRUE, shiny.height = "500px", data.suffix){
  if (missing(data.suffix)){
    data.suffix <- suffix
  }
  output[[suffix]] <- renderPlot({ 
    plot },
    width=function(){ 72*input[[paste0("width_", suffix)]] },
    height=function(){ 72*input[[paste0("height_", suffix)]] }
  )
  # Conditional UI plot 
  output[[paste0("ui.", suffix)]] <- renderUI({
    if (plot[1] == "None"){
      NULL
    } else {
      plotOutput(outputId = suffix, height = shiny.height)
    }
  })
  output[[paste0("download_", suffix)]] <- downloadHandler(
    filename = function(){
      if (isTRUE(ggplot)){
        paste0(file.name, "_", input[[paste0("data_", data.suffix)]], ".", input[[paste0("downtype_", suffix)]])
      } else {
        paste0(file.name, "_", input[[paste0("data_", data.suffix)]], ".png")
      }
    },
    content = function(file){
      if (isTRUE(ggplot)){
        ggsave2(filename = file,
                plot = plot,
                device = input[[paste0("downtype_", suffix)]],
                width = ifelse(input[[paste0("width_", suffix)]] == 1, 18, input[[paste0("width_", suffix)]]), 
                height = ifelse(input[[paste0("height_", suffix)]] == 1, 7, input[[paste0("height_", suffix)]]), 
                dpi = 600L, 
                units = "in")
      } else {
        png(file, 
            width = ifelse(input[[paste0("width_", suffix)]] == 1, 18, input[[paste0("width_", suffix)]]),
            height = ifelse(input[[paste0("height_", suffix)]] == 1, 7, input[[paste0("height_", suffix)]]), 
            units = "in", 
            res = 600)
        print(plot)
        dev.off()
      }
    }
  )
}
# Server code to allow newlines in notifications
showNotification2 <- function (ui, action = NULL, duration = 5, closeButton = TRUE, 
                               id = NULL, type = c("default", "message", "warning", "error"), 
                               session = shiny:::getDefaultReactiveDomain()) {
  if (is.null(id)) 
    id <- shiny:::createUniqueId(8)
  res <- shiny:::processDeps(HTML(ui), session)
  actionRes <- shiny:::processDeps(action, session)
  session$sendNotification("show", list(html = res$html, action = actionRes$html, 
                                        deps = c(res$deps, actionRes$deps), duration = if (!is.null(duration)) duration * 
                                          1000, closeButton = closeButton, id = id, type = match.arg(type)))
  id
}
# Server code to normalise data
normalisation <- function(physeq, norm, group_var){
  if (norm == "Rarefy"){
    ps <- phyloseq::rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                                      rngseed = 711, replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
  } else if (norm == "TSS"){
    ps <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x))
  } else if (norm == "UQ"){
    # Since every gene contains at least one zero, DESeq cannot compute log geometric means
    # To overcome this issue add a pseudocount to all genes
    # HOWEVER:  If you have very atypical samples that really do have zeros in every gene 
    # (really common in microbial ecology amplicon datasets), 
    # don't just blindly add pseducounts to force the data through a workflow that is not designed for your kind of data.
    ps <- physeq
    ps@otu_table <- ps@otu_table + 0.01
    df <- as.data.frame(ps@otu_table)
    m <- as.matrix(df)
    group = get_variable(physeq, group_var)
    dge <- edgeR::DGEList(counts = m, group = group, genes = tax_table(ps), remove.zeros = TRUE)
    y <- edgeR::calcNormFactors(dge, method = "upperquartile")
    CPM <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
    ps@otu_table <- otu_table(CPM, taxa_are_rows = TRUE)
  } else if (norm == "CSS"){
    ps <- metagMisc::phyloseq_transform_css(physeq, norm = TRUE, log = FALSE)
  } else if (norm == "DESeq"){
    ps <- physeq
    ps@otu_table <- ps@otu_table + 1
    DESeq.object = phyloseq::phyloseq_to_deseq2(ps, as.formula(paste0("~", group_var)))
    dds <- estimateSizeFactors(DESeq.object)
    normalized_counts <- counts(dds, normalized=TRUE)
    ps@otu_table <- otu_table(normalized_counts, taxa_are_rows = TRUE)
  } else if (norm == "TMM"){
    group = get_variable(physeq, group_var)
    df <- as.data.frame(physeq@otu_table)
    m <- as.matrix(df)
    dge <- edgeR::DGEList(counts = m, group = group, genes = tax_table(physeq))
    y <- edgeR::calcNormFactors(dge, method = "TMM")
    CPM <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
    ps <- physeq
    ps@otu_table <- otu_table(CPM, taxa_are_rows = TRUE)
  } else {
    ps <- physeq
  }
  return(ps)
}
# Server code to transform data
transformation <- function(physeq, transf){
  if (transf == "CLR_min"){
    t <- metaGVis::clr_min_percent(physeq)
  } else if (transf == "CLR_0plus1"){
    t <- metaGVis::clr_0_plus_1(physeq)
  } else if (transf == "CLR_runif"){
    t <- metaGVis::clr_runif(physeq)
  } else if (transf == "Hellinger"){
    t <- metaGVis::hellinger(physeq)
  } else if (transf == "logChord_0plus1"){
    t <- metaGVis::logchord_0_plus_1(physeq)
  } else if (transf == "logChord_1p"){
    t <- metaGVis::logchord_1p(physeq)
  } else {
    t <- physeq
  }
  return(t)
}
