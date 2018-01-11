##
## Server File for epiicm Shiny Application
##
## Run local: epiweb(class = "icm")
## Run online: http://statnet.shinyapps.io/epiicm/
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output) {

  ## Main reactive functions
  param <- reactive({
    param.icm(inf.prob = input$inf.prob,
              act.rate = input$act.rate,
              rec.rate = input$rec.rate,
              b.rate = input$b.rate,
              ds.rate = input$ds.rate,
              di.rate = input$di.rate,
              dr.rate = input$dr.rate)
  })
  init <- reactive({
    if (input$modtype == "SIR") {
      init.icm(s.num = as.numeric(input$s.num),
               i.num = as.numeric(input$i.num),
               r.num = as.numeric(input$r.num))
    } else {
      init.icm(s.num = as.numeric(input$s.num),
               i.num = as.numeric(input$i.num))
    }
  })
  control <- reactive({
    control.icm(type = input$modtype,
                nsteps = input$nsteps,
                nsims = input$nsims,
                verbose = FALSE)
  })
  mod <- reactive({
    input$runMod
    isolate(icm(param(), init(), control()))
  })
  showqnts <- reactive({
    ifelse(input$qntsrng == 0, FALSE, input$qntsrng)
  })

  ## Main Plot tab
  output$MainPlot <- renderPlot({
    par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
    if (input$compsel == "Compartment Prevalence") {
      plot(mod(),
           popfrac = TRUE,
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           legend = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
    if (input$compsel == "Compartment Size") {
      plot(mod(),
           popfrac = FALSE,
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           legend = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
    if (input$compsel == "Disease Incidence") {
      plot(mod(),
           y = "si.flow",
           popfrac = FALSE,
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           legend = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
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
             mean.line = input$showmean,
             sim.lines = input$showsims,
             qnts = showqnts(),
             legend = input$showleg,
             leg.cex = 1.1,
             lwd = 3.5)
      }
      if (input$compsel == "Compartment Size") {
        plot(mod(),
             mean.line = input$showmean,
             sim.lines = input$showsims,
             qnts = showqnts(),
             legend = input$showleg,
             popfrac = FALSE,
             leg.cex = 1.1,
             lwd = 3.5)
      }
      if (input$compsel == "Disease Incidence") {
        plot(mod(),
             y = "si.flow",
             popfrac = FALSE,
             mean.line = input$showmean,
             sim.lines = input$showsims,
             qnts = showqnts(),
             legend = input$showleg,
             leg.cex = 1.1,
             lwd = 3.5)
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

  # Comp_plot
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
  output$simnoControl <- renderUI({
    input$runMod
    maxsims <- isolate(input$nsims)
    sliderInput(inputId = "datasim",
                label = strong("Simulation Number"),
                min = 1,
                max = maxsims,
                value = 1,
                step = 1)
  })
  output$outData <- renderDataTable({
    if (input$datasel == "Means") {
      round(as.data.frame(mod()), input$tabdig)
    } else if (input$datasel == "Standard Deviations") {
      round(as.data.frame(mod(), out = "sd"), input$tabdig)
    } else if (input$datasel == "Simulations") {
      as.data.frame(mod(), out = "vals", sim = max(1, input$datasim))
    }
  }, options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10))
  output$dlData <- downloadHandler(
    filename = "ModelData.csv",
    content = function(file) {
      if (input$datasel == "Means") {
        write.csv(as.data.frame(mod()), file)
      } else if (input$datasel == "Standard Deviations") {
        write.csv(as.data.frame(mod(), out = "sd"), file)
      } else if (input$datasel == "Simulations") {
        write.csv(as.data.frame(mod(), out = "vals", sim = input$datasim), file)
      }

    }
  )

})
