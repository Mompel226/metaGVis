############ GIFS
## https://cran.r-project.org/web/packages/shinybusy/shinybusy.pdf
library(shiny)
library(shinybusy)
ui <- fluidPage(
  # Use this function somewhere in UI
  use_busy_gif(
    src = "https://jeroen.github.io/images/banana.gif",
    height = 70, width = 70
  ),
  actionButton("play", "Play GIF"),
  actionButton("stop", "Stop GIF")
)
server <- function(input, output, session) {
  observeEvent(input$play, {
    play_gif()
  })
  observeEvent(input$stop, {
    stop_gif()
  })
}
if (interactive()) {
  shinyApp(ui, server)
}

################# POP-UP WINDOWS

rm(list = ls())
library(shiny)
library(DT)

ui <- fluidPage(
  
  # Inlcude the line below in ui.R so you can send messages
  tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
  titlePanel("Pop-up Alerts"),
  sidebarPanel(
    sliderInput("my_slider", "Range Slider:", min = 0, max = 150, value = 40, step=1),
    dateInput('my_daterange',label = '',value = Sys.Date()),
    actionButton("run","Execute")),
  mainPanel(DT::dataTableOutput('tbl'))
)

server <- function(input, output, session) {
  
  # Alert below will trigger if the slider is over 100
  observe({
    if(input$my_slider >= 100)
    {
      my_slider_check_test <- "Your slider value is above 100 - no data will be displayed"
      js_string <- 'alert("SOMETHING");'
      js_string <- sub("SOMETHING",my_slider_check_test,js_string)
      session$sendCustomMessage(type='jsCode', list(value = js_string))
    }
  })
  
  
  # Alert below about dates will notify you if you selected today
  observe({
    if (is.null(input$run) || input$run == 0){return()}
    isolate({
      input$run
      if(input$my_daterange == Sys.Date())
      {
        my_date_check_test <- "Today Selected"
        js_string <- 'alert("SOMETHING");'
        js_string <- sub("SOMETHING",my_date_check_test,js_string)
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
      # Alert will also trigger and will notify about the dates
      if(input$my_daterange == Sys.Date())
      {
        my_date_check_test <- paste0("You selected: ",input$my_daterange)
        js_string <- 'alert("SOMETHING");'
        js_string <- sub("SOMETHING",my_date_check_test,js_string)
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
      
    })
  })
  
  my_data <- reactive({
    
    if(input$run==0){return()}
    isolate({
      input$run
      if(input$my_slider >= 100)
      {
        # Alert below will trigger if you adjusted the date but slider is still 100
        my_slider_check_test <- "Slider is still over 100"
        js_string <- 'alert("SOMETHING");'
        js_string <- sub("SOMETHING",my_slider_check_test,js_string)
        session$sendCustomMessage(type='jsCode', list(value = js_string))
      }
      if(input$my_slider < 100)
      {
        iris[1:input$my_slider,]
      }
    })  
  })
  output$tbl = DT::renderDataTable(my_data(), options = list(lengthChange = FALSE))
}

shinyApp(ui = ui, server = server)

############## CONDITIONAL PANELS

if (interactive()) {
  ui <- fluidPage(
    sidebarPanel(
      selectInput("plotType", "Plot Type",
                  c(Scatter = "scatter", Histogram = "hist")
      ),
      # Only show this panel if the plot type is a histogram
      conditionalPanel(
        condition = "input.plotType == 'hist'",
        selectInput(
          "breaks", "Breaks",
          c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
        ),
        # Only show this panel if Custom is selected
        conditionalPanel(
          condition = "input.breaks == 'custom'",
          sliderInput("breakCount", "Break Count", min = 1, max = 50, value = 10)
        )
      )
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
  
  server <- function(input, output) {
    x <- rnorm(100)
    y <- rnorm(100)
    
    output$plot <- renderPlot({
      if (input$plotType == "scatter") {
        plot(x, y)
      } else {
        breaks <- input$breaks
        if (breaks == "custom") {
          breaks <- input$breakCount
        }
        
        hist(x, breaks = breaks)
      }
    })
  }
  
  shinyApp(ui, server)
}