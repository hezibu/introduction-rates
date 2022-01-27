library(shiny)
library(shinydashboard)
library(ggplot2)

# THIS APP DOES NOT WORK ITS OWN, I NEED TO UNDERSTAND HOW SOURCE() WORKS WITH SHINY


ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(sliderInput("sliderb0", "b0", min= -2, max= 3   , step=0.1  , value= -0.5),
                   sliderInput("sliderb1", "b1", min= -0.2, max= 0.2   , step=0.001  , value=  0.01),
                   sliderInput("sliderg0", "g0", min= -2, max= 2   , step=0.1  , value=   -1),
                   sliderInput("sliderg1", "g1", min=  0, max= 0.1   , step=0.0001  , value=    0),
                   sliderInput("sliderg2", "g2", min=  0, max= 0.01, step=0.000001, value=    0)),
  dashboardBody( 
    fluidRow(column(6,plotOutput('waveplot')))
  ))

server <- function(input, output, session) { 
  
  source("02_global_functions.R", local = TRUE)
  
  output$waveplot <- renderPlot({
    x <- seq_along(150)
    params <- set_params(c(input$sliderb0,
                           input$sliderb1,
                           input$sliderg0,
                           input$sliderg1,
                           input$sliderg2))
    set.seed(1234)
    y <- sim(N = length(x), params = params, constant = 0)
    df <- data.frame(x,y)
    ggplot()+
      #geom_line(aes(x = x, y = cumsum(test)), size = 1.1, color = "grey")+
      geom_line(data = df, aes(x=x,y=cumsum(y)), size = 1.3)+ 
      scale_x_continuous()+
      #scale_y_continuous(limits = c(0, max(cumsum(test))))+
      theme_classic()
  })
}



shinyApp(ui, server)
