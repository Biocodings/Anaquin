library(shiny)

shinyUI(fluidPage(
  titlePanel("Anaquin - Uploading Sequence Files!"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose File',
                accept=c('text/csv', 
								 'text/comma-separated-values,text/plain', 
								 '.csv')),
      tags$hr()
    ),
    mainPanel(
          img(src="http://www.anaquin.org/wp-content/uploads/2015/10/logo_white12.svg")
    )
  )
))
