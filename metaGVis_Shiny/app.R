# GENERAL UI
setwd("~/Dropbox/PhD/R/Packages/metaGVis_Shiny")
################################################################################
# Analysis related libraries
################################################################################
library(metaGVis)
library(magrittr)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(data.table)
library(file2meco)
################################################################################
# Shiny related libraries
################################################################################
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyFiles)
library(colourpicker)
library(plotly)
library(htmlwidgets) 
################################################################################
# UI tabs
source("ggsave.R", local = TRUE)
source("Tabs/Import_ui.R", local = TRUE)
source("Tabs/Filtration_ui.R", local = TRUE)
source("Tabs/LibraryQC_ui.R", local = TRUE)
source("Tabs/AlphaDiversity_ui.R", local = TRUE)
source("Tabs/TaxonomicProfile_ui.R", local = TRUE)
source("Tabs/BetaDiversity_ui.R", local = TRUE)
source("Tabs/DA_ui.R", local = TRUE)
source("Tabs/FunctionalProfile_ui.R", local = TRUE)
################################################################################
# Define the overall user-interface "ui"
################################################################################
ui = dashboardPage(
  scrollToTop = TRUE,
  skin = "blue",
  options = list(sidebarExpandOnHover = TRUE),
  header = dashboardHeader(title = "metaGVis"),
  sidebar = dashboardSidebar(minified = TRUE, 
                             collapsed = FALSE,
                             # Code to get screen pixel dimensions
                             dashboardSidebar(
                               tags$head(tags$script('
                                var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                            ')),
                               sidebarMenu(style = "position:fixed;",
                                 menuItem("Import data", tabName = "import", icon = icon("file-import", lib = "font-awesome")),
                                 menuItem("Filtration", tabName = "filt", icon = icon("filter", lib = "font-awesome")),
                                 menuItem("Library QC", tabName = "libqc", icon = icon("chart-bar", lib = "font-awesome")),
                                 menuItem("Alpha diversity", tabName = "alpha", icon = icon("adn", lib = "font-awesome")),
                                 menuItem("Beta diversity", tabName = "beta", icon = icon("project-diagram", lib = "font-awesome")),
                                 menuItem("Taxonomic profile", tabName = "profile", icon = icon("chart-pie", lib = "font-awesome")),
                                 menuItem("Differential abundance", tabName = "DA", icon = icon("map-signs", lib = "font-awesome")),
                                 menuItem("Functional profile", tabName = "function", icon = icon("hubspot", lib = "font-awesome"))
                               )
                             )),
  body = dashboardBody(
    tabItems(
      # Import tabItem
      import_tab,
      # Filtration tabItem
      filt_tab,
      # Library QC tabItem
      libqc_tab,
      # Alpha diversity tabItem
      alpha_tab,
      # Beta diversity tabItem
      beta_tab,
      # Taxonomic profile tabItem
      profile_tab,
      # Differential abundance tabItem
      DA_tab,
      # Functional profile tabItem
      function_tab
      )
    )
  )
################################################################################
# Define the overall "server"
################################################################################
server = function(input, output, session) {
  
  # Globally defined functions 
  source("server_functions.R", local = TRUE)
  # Globally defined plot functions 
  source("plot_functions.R", local = TRUE)
  # Import data
  source("Tabs/Import_server.R", local = TRUE)
  # Filtration
  source("Tabs/Filtration_server.R", local = TRUE)
  # Library QC
  source("Tabs/LibraryQC_server.R", local = TRUE)
  # Alpha Diversity
  source("Tabs/AlphaDiversity_server.R", local = TRUE)
  # Beta Diversity
  source("Tabs/BetaDiversity_server.R", local = TRUE)
  # Taxonomic Profile
  source("Tabs/TaxonomicProfile_server.R", local = TRUE)
  # Differential abundance
  source("Tabs/DA_server.R", local = TRUE)
  # Functional Profile
  source("Tabs/FunctionalProfile_server.R", local = TRUE)
  
}

# Run app
shinyApp(ui, server)
