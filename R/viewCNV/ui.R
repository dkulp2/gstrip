#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  singleton(
    tags$head(tags$script(src = "message-handler.js"))
  ),
  
  # Application title
  titlePanel("Copy Number Variants"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      tags$table(tags$td(selectInput("saved.site","Saved Sites",list(),width='350px')),
                 tags$td(actionButton("saved.prev","Prev")),
                 tags$td(actionButton("saved.next","Next"))),
      actionButton("reset","Reset"),actionButton("delete","Delete"),
      actionButton("save.button","Save Site"),
      tags$table(tags$td(selectInput("candidate","GenomeStrip CNVs", list(), width='350px')),
                 tags$td(actionButton("candidate.prev","Prev")),
                 tags$td(actionButton("candidate.next","Next"))),
      hr(),
      tags$table(tags$td(textInput("seg.chr","Chromosome",NA)),
                 tags$td(numericInput("seg.start","Start Position", NA)),
                 tags$td(numericInput("seg.end", "End Position", NA)),
                 tags$td(textOutput("seg.size")),
                 tags$td(actionButton("zoom.out", "Zoom Out"))),
      selectizeInput("seg.sample", "Samples", choices=list(), multiple=TRUE),
      textareaInput("note","Note (first line is title)", "", rows = 3),
      h3("Show"),
      tags$table(tags$td(checkboxInput("show_cnv","CNV Extents", value=TRUE)),
                 tags$td(checkboxInput("show_winGeno","Windowed CN Genotypes", value=TRUE)),
                 tags$td(checkboxInput("show_frag","Fragment Profiles", value=TRUE)),
                 tags$td(checkboxInput("show_readPairs","Read Pair Counts", value=FALSE))),
      conditionalPanel("input.show_winGeno",
                       sliderInput("winGeno_nudge_L","Nudge Windowed CN Genotypes (inches)", min=0, max=3, value=0.75, step=0.25)),
      conditionalPanel("input.show_frag",
                       sliderInput("win.size", "Number of Bins in Profile Window:", min = 1, max = 200, value = 20, step=5)),
      numericInput("pad","Left/Right Padding (bases)", min=0, max=500000, value=5000, step=1000),
      conditionalPanel("input.show_cnv",
                       sliderInput("min.cnv.len", "Minimum Length of Displayed CNV Extents", min=0, value=500, max=5000, step=500)),
      sliderInput("paired.reads","Minimum Number of Paired Reads", min=0, max=50, value=0, step=1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("allThree", height="750px"),
      plotOutput("expPlot",height="200px")
      # conditionalPanel("input.show_cnv",plotOutput("cnvPlot",height=300)),
      # conditionalPanel("input.show_winGeno",plotOutput("winGenoPlot",height=250)),
      # conditionalPanel("input.show_frag",plotOutput("fragPlot",height=300))
    )
  )
))
