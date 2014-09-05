
library(EpiModel)
library(shiny)

shinyServer(function(input, output, session) {

  ## Main reactive functions
  param <- reactive({
    vital <- ifelse(input$b.rate > 0 |
                    input$ds.rate > 0 |
                    input$di.rate > 0 |
                    input$dr.rate > 0, TRUE, FALSE)
    l <- list(inf.prob = input$inf.prob,
              act.rate = input$act.rate,
              rec.rate = input$rec.rate,
              b.rate = input$b.rate,
              ds.rate = input$ds.rate,
              di.rate = input$di.rate,
              dr.rate = input$dr.rate,
              vital = vital,
              groups = 1)
    class(l) <- "param.dcm"
    return(l)
  })
  init <- reactive({
    if (input$modtype == "SIR") {
      l <- list(s.num = input$s.num,
                i.num = input$i.num,
                r.num = input$r.num)
    } else {
      l <- list(s.num = input$s.num,
                i.num = input$i.num)
    }
    class(l) <- "init.dcm"
    return(l)
  })
  control <- reactive({
    l <- list(type = input$modtype,
              dt = seq(1, input$nsteps, input$dt),
              odemethod = "rk4",
              print.mod = FALSE,
              verbose = FALSE)
    class(l) <- "control.dcm"
    return(l)
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
    summary(mod(),
            at = input$summTs,
            digits = input$summDig)
  })

  # comp_plot
  output$CompPlot <- renderPlot({
    comp_plot(mod(),
              at = input$summTs,
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
  }, options = list(aLengthMenu = c(10, 25, 50, 100), iDisplayLength = 10))

  output$dlData <- downloadHandler(
    filename = "ModelData.csv",
    content = function(file) {
      write.csv(as.data.frame(mod()), file)
    }
  )


})
