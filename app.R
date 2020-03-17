#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)

baseparam = function() {
  param = list()
  param$aI = 0.5
  param$d = 1/6
  param$m = 0.03
  param$r = 1/6
  
  param$tStart = 30
  param$tEnd = 99
  param$eff = 0.3
  
  param$tMax = 100
  
  return(param)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    # Add google analytics tracking:
    includeHTML(("googleanalytics.html")),
    # Make rules widers:
    tags$style(HTML("hr {border-top: 1px solid #444444;}"))
  ),
  # Application title
  titlePanel("Corona infection simulator"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3('Quarantine:')
      ,
      sliderInput("tStart",
                  "Day quarantine starts:",
                  min = 0,
                  max = 100,
                  value = 30),
      sliderInput("tQuarantine",
                  "Length of quarantine:",
                  min = 0,
                  max = 99,
                  value = 14),
      sliderInput("eff",
                  "Transmission reduction during quarantine:",
                  min = 0,
                  max = 1,
                  value = 0.3)
      ,
      h3('Simulation parameters')
      ,
      sliderInput("aI",
                  "Transmission rate (per day):",
                  min = 0,
                  max = 1,
                  value = 0.5),
      sliderInput("d_recip",
                  "Time from infection to symptions (days):",
                  min = 0,
                  max = 12,
                  value = 6),
      sliderInput("m",
                  "Mortality (% of diseased dying):",
                  min = 0,
                  max = 20,
                  value = 3),
      sliderInput("r_recip",
                  "Time to recovery  (days):",
                  min = 0,
                  max = 12,
                  value = 6)
      #),
      # sidebarPanel(
      
    ),
    # Show plots
    mainPanel(
      tabsetPanel(
        tabPanel('About',
                 uiOutput("about")),
        tabPanel('Simulation results',
                 plotOutput("plotEpidemic")),
        selected='Simulation results'
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  sim = eventReactive({
    input$aI
    input$d_recip
    input$m
    input$r_recip
    input$tStart
    input$tQuarantine
    input$eff
  },
  {
    # Set all parameters
    param = list()
    param$a = input$aI
    param$d = 1/input$d_recip
    param$m = input$m/100
    param$r = 1/input$r_recip
    param$tMax = 100
    param$tStart = input$tStart
    param$tEnd = min(param$tMax-1, input$tQuarantine+param$tStart)
    param$eff = input$eff
    
    # Simulate
    return( runCorona(param) )   
  })
  
  output$plotEpidemic <- renderPlot({
    plotCorona( sim(), input$tStart, min(99, input$tQuarantine+input$tStart) )
  }, height=800)
  
  output$about <- renderText("The simulator is only for illustration and should not be used for decision support <br>
See https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
for in-depth analysis <br>
                             <br>
Code available on https://github.com/Kenhasteandersen/Corona <br><br>
                             Made by Ken H Andersen, kha@aqua.dtu.dk<br>")
}

runCorona <- function(param) {
  
  derivatives <- function(t,y,param) {
    S = y[1]
    I = y[2]
    D = y[3]
    M = y[4]
    R = y[5]
    
    transmission = param$a*(I+D)*S
    diseased = param$d*I
    death = param$m*param$r*D
    recovery = param$r*D
    
    dSdt = -transmission
    dIdt = transmission - diseased
    dDdt = diseased - death - recovery
    dMdt = death
    dRdt = recovery
    
    dydt = list(c(dSdt, dIdt, dDdt, dMdt, dRdt))
  }
  
  y0 = c(1, 1/5e6, 0,0,0)
  # Run before quarantine:
  out = as.data.frame( ode(y0, seq(1, param$tStart, by=1), derivatives, parms = param) )
  # Run during quarantine:
  if (param$tEnd>param$tStart) {
    param$a = param$a * param$eff # Reduce transmission
    n = dim(out)[1]
    y0 = as.numeric( out[ n, 2:6] )
    out = rbind( out[ 1:(n-1), ],
                 as.data.frame( ode(y0, seq(param$tStart, param$tEnd, by=1), derivatives, parms = param)))
    # Run after end of quarantine:
    param$a = param$a / param$eff # Reduce transmission
  }
  n = dim(out)[1]
  y0 = as.numeric( out[ n, 2:6] )
  out = rbind( out[ 1:(n-1), ],
               as.data.frame( ode(y0, seq(param$tEnd, param$tMax, by=1), derivatives, parms = param)))
  
  names(out) = c("time", "S", "I", "D", "M", "R")
  return(out)
}

plotCorona = function(out, tStart, tEnd) {
  col = c("black","orange","red","blue","green")
  par(mfcol=c(2,1))
  #
  # Linear plot
  #
  plot(out$time, out$S, type="l", lwd=3, ylim=c(0,1), 
       xlab="Time (days)", ylab="Fraction of population", col=col[1])
  #
  # Quarantine patch
  #
  if (tEnd > tStart) {
    polygon( c(tStart, tEnd, tEnd, tStart), c(0,0,1,1), col=grey(0.8), border=NA ) 
    text( tStart + 0.5*(tEnd-tStart), 0.95, labels="Quarantine")
  }
  #lines( tStart*c(1,1), c(0,1), lty=3 )
  #lines( tEnd*c(1,1), c(0,1), lty=3 )
  lines(out$time, out$S, lwd=3)
  lines(out$time, out$I, col=col[2], lwd=3)
  lines(out$time, out$D, col=col[3], lwd=3)
  lines(out$time, out$M, col=col[4], lwd=3)
  lines(out$time, out$R, col=col[5], lwd=3)
  
  legend( x="right",
          legend=c("Susceptible (healthy)", "Infected (no symptoms)", "Diseased (symptoms)", "Dead", "Recovered"),
          col=col, lty=rep(1,5), lwd=3, bty="n")
  
  ix = length(out$time)
  text(x=out$time[ix], y=out$M[ix]+0.03, 
       adj=1, col="blue", 
       labels=paste( format(out$M[ix], digits=2), "% dead"))
  #
  # Semilog plot
  #
  #plot(out$time, out$S, type="l", lwd=3, ylim=c(1e-3,1), log="y",
  #     xlab="Time (days)", ylab="Fraction of population", col=col[1])
  #lines(out$time, out$I, col=col[2], lwd=3)
  #lines(out$time, out$D, col=col[3], lwd=3)
  #lines(out$time, out$M, col=col[4], lwd=3)
  #lines(out$time, out$R, col=col[5], lwd=3)
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

