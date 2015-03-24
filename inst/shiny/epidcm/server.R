##
## Server File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
## Run online: http://statnet.shinyapps.io/epidcm/
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output) {

  ## Main reactive functions
  param <- reactive({
    param.dcm(inf.prob = input$inf.prob,
              act.rate = input$act.rate,
              rec.rate = input$rec.rate,
              b.rate = input$b.rate,
              ds.rate = input$ds.rate,
              di.rate = input$di.rate,
              dr.rate = input$dr.rate)
  })
  init <- reactive({
    if (input$modtype == "SIR") {
      init.dcm(s.num = input$s.num,
               i.num = input$i.num,
               r.num = input$r.num)
    } else {
      init.dcm(s.num = input$s.num,
               i.num = input$i.num)
    }
  })
  control <- reactive({
    control.dcm(type = input$modtype,
                nsteps = input$nsteps,
                dt = input$dt,
                verbose = FALSE)
  })
  mod <- reactive({
    input$runMod
    isolate(dcm(param(), init(), control()))
  })
  sleg <- reactive({
    l <- ifelse(input$showleg == TRUE, "full", "n")
    return(l)
  })


  ## Main Plot tab
  output$MainPlot <- renderPlot({
    par(mar=c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
    if (input$compsel == "Compartment Prevalence") {
      plot(mod(),
           leg.cex = 1.1,
           alpha = input$alpha,
           lwd = 3.5,
           leg = sleg(),
           main = "")
    }
    if (input$compsel == "Compartment Size") {
      plot(mod(),
           popfrac = FALSE,
           leg.cex = 1.1,
           alpha = input$alpha,
           lwd = 3.5,
           leg = sleg(),
           main = "")
    }
    if (input$compsel == "Disease Incidence") {
      plot(mod(),
           y = "si.flow",
           popfrac = FALSE,
           leg.cex = 1.1,
           alpha = input$alpha,
           lwd = 3.5,
           leg = sleg(),
           main = "")
    }
  })
  output$dlMainPlot <- downloadHandler(
    filename = "MainPlot.pdf",
    content = function(file) {
      pdf(file = file, height = 6, width = 10)
      par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
      if (input$compsel == "Compartment Prevalence") {
        plot(mod(),
             leg.cex = 1.1,
             alpha = input$alpha,
             lwd = 3.5,
             leg = sleg(),
             main = "")
      }
      if (input$compsel == "Compartment Size") {
        plot(mod(),
             popfrac = FALSE,
             leg.cex = 1.1,
             alpha = input$alpha,
             lwd = 3.5,
             leg = sleg(),
             main = "")
      }
      if (input$compsel == "Disease Incidence") {
        plot(mod(),
             y = "si.flow",
             popfrac = FALSE,
             leg.cex = 1.1,
             alpha = input$alpha,
             lwd = 3.5,
             leg = sleg(),
             main = "")
      }
      dev.off()
    }
  )

  ## Summary and Compartment plot tab
  # Outfrom from summary
  output$outSummary <- renderPrint({
    if (is.na(input$summTs)) {
      summat <- 1
    } else {
      summat <- input$summTs
    }
    summary(mod(),
            at = summat,
            digits = input$summDig)
  })

  # comp_plot
  output$CompPlot <- renderPlot({
    if (is.na(input$summTs)) {
      summat <- 1
    } else {
      summat <- input$summTs
    }
    comp_plot(mod(),
              at = summat,
              digits = input$summDig)
  })

  # Download for comp_plot
  output$dlCompPlot <- downloadHandler(
    filename = "CompPlot.pdf",
    content = function(file) {
      pdf(file = file, height = 6, width = 10)
      comp_plot(mod(),
                at = input$summTs,
                digits = input$summDig)
      dev.off()
    }
  )


  ## Data tab
  output$outData <- renderDataTable({
    as.data.frame(mod())
  }, options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))

  output$dlData <- downloadHandler(
    filename = "ModelData.csv",
    content = function(file) {
      write.csv(as.data.frame(mod()), file)
    }
  )

})
