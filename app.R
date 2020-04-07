#
# This is a Shiny web application. 
#
library(shiny)
library(deSolve)

baseparam = function() {
  param = list()
  #
  # Basic SIR parameters
  #
  param$d = 1/6 # rate of developing symptoms
  param$m = 0.01 # The risk of dying if the disease is contracted
  param$r = 1/6 # recovery rate
  param$aI = 2.5/(1/param$r+1/param$d)
  #
  # Quarantine parameters
  #
  param$tStart = 30
  param$tEnd = 45
  param$eff = 0.3
  param$tEnd2 = 66
  param$eff2 = 0.15
  #
  # ICU (per person)
  #
  param$ICU = 0.0001 # US; https://www.medpagetoday.com/hospitalbasedmedicine/generalhospitalpractice/84845
  param$ICUAdmissionFraction = 0.0025  # https://www.medpagetoday.com/hospitalbasedmedicine/generalhospitalpractice/84845
  #
  # Simulation time
  #
  param$tMax = 200
  
  return(param)
}
#
# Define UI for application
#
ui <- fluidPage(
  tags$head(
    # Add google analytics tracking:
    includeHTML(("googleanalytics.html")),
    # Make rules widers:
    tags$style(HTML("hr {border-top: 1px solid #444444;}"))
  ),
  # Application title
  titlePanel("Corona quarantine simulator"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3('Quarantine:')
      ,
      sliderInput("tStart",
                  "Day quarantine starts:",
                  min = 2,
                  max = 100,
                  value = 41),
      sliderInput("tQuarantine",
                  "Length of quarantine:",
                  min = 0,
                  max = 99,
                  value = 21),
      sliderInput("eff",
                  "Transmission reduction during quarantine (% reduction):",
                  min = 0,
                  max = 100,
                  value = 50),
      sliderInput("tQuarantine2",
                  "Length of reduced quarantine:",
                  min = 0,
                  max = 99,
                  value = 40),
      sliderInput("eff2",
                  "Transmission reduction during reduced quarantine (% reduction):",
                  min = 0,
                  max = 100,
                  value = 25)
      ,
      hr(),
      checkboxInput("bParamICU",
                    label="Intensive care capacity: ",
                    value=FALSE)
      ,
      conditionalPanel(
        condition = "input.bParamICU==true"
        ,
        h3('Intensive care capacity')
        ,
        sliderInput("ICU",
                    "No. of intensive care units per million",
                    min=0,
                    max=1000,
                    value=300),
        sliderInput("ICUAdmissionFraction",
                    "% needing intensive care",
                    min=0,
                    max=1.5,
                    value=0.25,
                    step=0.05)
        ,
        hr(),
        checkboxInput("bMoreParam",
                      label="Show infection parameters",
                      value=FALSE)
        ,
        conditionalPanel(
          condition = "input.bMoreParam==true"
          ,
          h3('Epidemic parameters')
          ,
          
          sliderInput("aI",
                      "Transmission rate (per day):",
                      min = 0,
                      max = 1,
                      value = 0.2),
          sliderInput("d_recip",
                      "Time from infection to symptoms (days):",
                      min = 0,
                      max = 12,
                      value = 6),
          sliderInput("m",
                      "Mortality (% of diseased dying):",
                      min = 0,
                      max = 5,
                      value = 1),
          sliderInput("r_recip",
                      "Time to recovery  (days):",
                      min = 0,
                      max = 12,
                      value = 6)
        )
      )
      #),
      # sidebarPanel(
      
    ),
    # Show plots
    mainPanel(
      tabsetPanel(
        tabPanel('About',
                 uiOutput("about")),
        tabPanel('Simulation results',
                 plotOutput("plotEpidemic"),
                 plotOutput("plotR")),
        selected='Simulation results'
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  sim = eventReactive({
    input$tStart
    input$tQuarantine
    input$eff
    input$tQuarantine2
    input$eff2
    input$ICU
    input$ICUAdmissionFraction
    input$aI
    input$d_recip
    input$m
    input$r_recip
    
  },
  {
    # Set all parameters
    param = baseparam()
    
    param$tStart = input$tStart
    param$tEnd = min(param$tMax-1, input$tQuarantine+param$tStart)
    param$eff = 1-input$eff/100
    
    param$tEnd2 = min(param$tMax-1, param$tEnd + input$tQuarantine2)
    param$eff2 = 1-input$eff2/100
    
    param$ICU = input$ICU*1e-6
    param$ICUAdmissionFraction = input$ICUAdmissionFraction/100
    
    param$a = input$aI
    param$d = 1/input$d_recip
    param$m = input$m/100
    param$r = 1/input$r_recip
    
    # Simulate
    return( list(sim=runCorona(param), param=param) )   
  })
  
  output$plotEpidemic <- renderPlot({
    plotCorona( sim()$sim, sim()$param )
  }, height=400)
  
  output$plotR <- renderPlot({
    plotR( sim()$sim, sim()$param )
  }, height=250)
  
  output$about <- renderUI({tagList(p("The simulator is only for illustration and should not be used for decision support"),
                                    p("The simulator is based on a standard SIR (Susceptible-Infected-Recovered) epidemics model, though with an added distinction between
                                    those being infection (without symptoms) and those with symptoms."),
                                    p("Parameters are set based on rough estimates from wikipedia etc. Please let me know if any are wrong"),
                                    "See ", a("here", href="https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf") ,
                                    " for in-depth analysis. ",
                                    ("Code available on "),
                                    a("github", href="https://github.com/Kenhasteandersen/Corona"),
                                    ". Made by ",
                                    a("Ken H Andersen", href="http://ken.haste.dk"))
  })
  
  #"The simulator is only for illustration and should not be used for decision support <br><br>
  #The simulator is based on a standard SIR (Susceptible-Infected-Recovered) epidemics model, though with an added distinction between
  #those being infection (without symptoms) and those with symptoms.<br><br>
  #Parameters are set based on rough estimates from wikipedia etc. Please let me know if any are wrong<br><br>
  #See https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
  #for in-depth analysis <br>
  #                             <br>
  #Code available on https://github.com/Kenhasteandersen/Corona <br><br>
  #                             Made by Ken H Andersen, kha@aqua.dtu.dk, http://ken.haste.dk.")
  
}

runCorona <- function(param) {
  
  derivatives <- function(t,y,param) {
    S = y[1]
    I = y[2]
    D = y[3]
    M = y[4]
    R = y[5]
    #Dintensive = y[6]
    
    transmission = param$a*(I+D)*S
    diseased = param$d*I
    #intensive_care = param$ICUAdmissionFraction*param$d*I
    
    #D_missing_ICU = max(0, D*param$ICUAdmissionFraction-param$ICU)
    
    death_normal = param$m*param$r*D
    #death_intensive = param$m*param$r*Dintensive
    #if (Dintensive > param$ICU) 
    #  death_intensive = death_intensive + (Dintensive - param$ICU)
    recovery = param$r*D 
    #recovery_intensive = 0.5*param$r*Dintensive
    
    dSdt = -transmission
    dIdt = transmission - diseased
    dDdt = diseased - death_normal - recovery
    dMdt = death_normal
    dRdt = recovery 
    #dDintensive_dt = intensive_care - death_intensive - recovery_intensive
    
    dydt = list(c(dSdt, dIdt, dDdt, dMdt, dRdt))
  }
  
  y0 = c(1, 500/1e6, 0,0,0)
  # Run before quarantine:
  out = as.data.frame( ode(y0, seq(1, param$tStart, by=1), derivatives, parms = param) )
  # Run during quarantine:
  if (param$tEnd>param$tStart) {
    param$a = param$a * param$eff # Reduce transmission
    n = dim(out)[1]
    y0 = as.numeric( out[ n, 2:6] )
    out = rbind( out[ 1:(n-1), ],
                 as.data.frame( ode(y0, seq(param$tStart, param$tEnd, by=1), derivatives, parms = param)))
    # Run during reduced quarantine:
    param$a = param$a/param$eff * param$eff2 # Reduce transmission
    n = dim(out)[1]
    y0 = as.numeric( out[ n, 2:6] )
    print(param$eff2)
    print(param$eff)
    out = rbind( out[ 1:(n-1), ],
                 as.data.frame( ode(y0, seq(param$tEnd, param$tEnd2, by=1), derivatives, parms = param)))
    
    # Run after end of quarantine:
    param$a = param$a / param$eff2 # Reduce transmission
  }
  n = dim(out)[1]
  y0 = as.numeric( out[ n, 2:6] )
  out = rbind( out[ 1:(n-1), ],
               as.data.frame( ode(y0, seq(param$tEnd2, param$tMax, by=1), derivatives, parms = param)))
  
  names(out) = c("time", "S", "I", "D", "M", "R")
  #
  # Calc R
  #
  a = rep(param$a, length(out$time))

  ixQuarantine = out$time>param$tStart & out$time<param$tEnd
  a[ixQuarantine] = a[ixQuarantine]*param$eff
  
  ixQuarantine2 = out$time>=param$tEnd & out$time<param$tEnd2
  a[ixQuarantine2] = a[ixQuarantine2]*param$eff2
  
  out$RR = a * out$S * (1/param$r+1/param$d)
  print(param$a)
  print(1/param$r+1/param$d)
  return(out)
}

plotCorona = function(out, param) {
  col = c("blue","blue","red","black","green")
  #par(mfcol=c(2,1))
  #
  # Main plot
  #
  # plot(out$time, out$S, type="l", lwd=3, ylim=c(0,1), 
  #      xlab="Time (days)", ylab="Fraction of population", col=col[1])
  # #
  # # Quarantine patch
  # #
  # if (param$tEnd > param$tStart) {
  #   polygon( c(param$tStart, param$tEnd, param$tEnd, param$tStart), c(0,0,1,1), col=grey(0.8), border=NA ) 
  #   text( param$tStart + 0.5*(param$tEnd-param$tStart), 0.95, labels="Quarantine")
  # }
  # #
  # # SIR output:
  # #
  # lines(out$time, out$S, lwd=3)
  # lines(out$time, out$I, col=col[2], lwd=3)
  # patchBelow(out$time, 0, out$D, col=rgb(1,0,0,0.5))
  # #?polygonlines(out$time, out$D, col=col[3], lwd=3)
  # lines(out$time, out$M, col=col[4], lwd=3)
  # lines(out$time, out$R, col=col[5], lwd=3)
  # 
  # legend( x="right",
  #         legend=c("Susceptible (healthy)", "Infected (no symptoms)", "Diseased (symptoms)", "Dead", "Recovered"),
  #         col=col, lty=rep(1,5), lwd=3, bty="n")
  # 
  # ix = length(out$time)
  # text(x=out$time[ix], y=out$M[ix]+0.03, 
  #      adj=1, col="blue", 
  #      labels=paste( format(out$M[ix]*100, digits=2), "% dead"))
  #
  # Main plot #2
  #
  plot(out$time, out$S, type="n", lwd=3, ylim=c(0,0.4), 
       xlab="Time (days)", ylab="Fraction of population", col=col[1])
  #
  # Quarantine patch
  #
  if (param$tEnd > param$tStart) {
    polygon( c(param$tStart, param$tEnd, param$tEnd, param$tStart), c(0,0,1,1), col=grey(0.8), border=NA ) 
    text( param$tStart + 0.5*(param$tEnd-param$tStart), 0.39, labels="Quarantine")
  }
  
  if (param$tEnd > param$tStart) {
    polygon( c(param$tEnd, param$tEnd2, param$tEnd2, param$tEnd), 
             c(0,0,1,1), col=grey(0.9), border=NA ) 
    text( param$tEnd + 0.5*(param$tEnd2-param$tEnd), 0.39, 
          labels="(reduced)")
  }
  
  #
  # SIR output:
  #
  #lines(out$time, out$S, lwd=3)
  patchBelow(out$time, 0, out$D, col=rgb(1,0,0,0.5))
  patchBelow(out$time, param$ICU/param$ICUAdmissionFraction, pmax(param$ICU/param$ICUAdmissionFraction, out$D), col=rgb(1,0,0,1))
  lines(out$time, out$I, col=col[2], lwd=3)
  lines(range(out$time), param$ICU/param$ICUAdmissionFraction*c(1,1), col="red")
  #?polygonlines(out$time, out$D, col=col[3], lwd=3)
  lines(out$time, out$M, col=col[4], lwd=3)
  #lines(out$time, 100*out$Dintensive, col="yellow")
  text(x=0, y=param$ICU/param$ICUAdmissionFraction, label="Limit of intensive\ncare capacity", adj=c(-0.01,-0.2), col="red")
  #lines(out$time, out$R, col=col[5], lwd=3)
  # lines(range(out$time), param$ICU*c(1,1))
  # text(x=0, y=param$ICU, label="No. of intensive care units", adj=c(-0.01,-0.2))
  #
  legend( x="topright",
          legend=c("","Infected (no symptoms)", "Diseased (symptoms)", "Dead"),
          col=col[1:4], lty=rep(1,4), lwd=c(0,3,3,3), bty="n")
  
  ix = length(out$time)
  text(x=out$time[ix], y=out$M[ix]+0.015, 
       adj=1, col=col[4], 
       labels=paste( format(out$M[ix]*100, digits=2), "% dead"))
  #
  # Patients in ICUs:
  #
  # plot(out$time, out$D*param$ICUAdmissionFraction, type="nS", lwd=3, ylim=c(0,1e-3), 
  #      xlab="Time (days)", ylab="Fraction of population", col=col[1],
  #      main="Patients in need of intensive care units")
  # # Quarantine patch
  # if (param$tEnd > param$tStart) {
  #   polygon( c(param$tStart, param$tEnd, param$tEnd, param$tStart), c(0,0,1,1), col=grey(0.8), border=NA ) 
  #   text( param$tStart + 0.5*(param$tEnd-param$tStart), 0.95, labels="Quarantine")
  # }
  # patchBelow(out$time, 0, out$D*param$ICUAdmissionFraction, col=rgb(1,0,0,0.5))
  # patchBelow(out$time, param$ICU, pmax(param$ICU, out$D*param$ICUAdmissionFraction), col=rgb(1,0,0,1))
  # lines(range(out$time), param$ICU*c(1,1))
  # text(x=0, y=param$ICU, label="No. of intensive care units", adj=c(-0.01,-0.2))
  #
  # Semilog plot
  #
  # plot(out$time, out$S, type="l", lwd=3, ylim=c(1e-4,1), log="y",
  #      xlab="Time (days)", ylab="Fraction of population", col=col[1])
  # #
  # # Quarantine patch
  # #
  # if (param$tEnd > param$tStart) {
  #   polygon( c(param$tStart, param$tEnd, param$tEnd, param$tStart), c(1e-5,1e-5,1,1), col=grey(0.8), border=NA ) 
  #   #text( param$tStart + 0.5*(param$tEnd-param$tStart), 0.95, labels="Quarantine")
  # }  
  # lines(out$time, out$S, lwd=3)
  # lines(out$time, out$I, col=col[2], lwd=3)
  # #lines(out$time, out$D, col=col[3], lwd=1)
  # patchBelow(out$time, 1e-5, out$D, col=rgb(1,0,0,0.5))
  # patchBelow(out$time, 1e-5, out$D*param$ICUAdmissionFraction, col=rgb(1,0,0,1))
  # #lines(out$time, out$D*param$ICUAdmissionFraction, col=col[3], lwd=3)
  # lines(out$time, out$M, col=col[4], lwd=3)
  # lines(out$time, out$R, col=col[5], lwd=3)
  # #
  # # ICUs
  # #
  # lines(range(out$time), param$ICU*c(1,1))
}

plotR = function(out, param) {
  plot(out$time, out$RR, type="l", lwd=3, col="blue",
      xlab="Time (days)", ylab="Epidemic reproductive number, R",
      ylim=c(0, 3))
  #
  # Quarantine patch
  #
  if (param$tEnd > param$tStart) {
    polygon( c(param$tStart, param$tEnd, param$tEnd, param$tStart), 
             c(0,0,5,5), col=grey(0.8), border=NA ) 
  }
  if (param$tEnd > param$tStart) {
    polygon( c(param$tEnd, param$tEnd2, param$tEnd2, param$tEnd), 
             c(0,0,5,5), col=grey(0.9), border=NA ) 
  }

  lines(out$time, out$RR, type="l", lwd=3, col="blue")  
  lines(out$time, 1+0*out$time, lty=3)
  
  text(param$tMax, 1.3, adj=c(1,0), labels="Evolving epidemic", col="red")
  text(param$tMax, 0.8, adj=c(1,1), labels="Controlled epidemic")
}

patchBelow <- function(x,y0,y,col) {
  ix = seq(length(x),1,by = -1)
  polygon(c(x, x[ix]), c(y*0+y0, y[ix]), col=col, border=NA)
}

# Run the application 
shinyApp(ui = ui, server = server)

