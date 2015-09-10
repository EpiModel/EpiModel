##
## Server File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output, session) {

  #Main reactive objects
  net <- reactive({
    network.initialize(n = input$num, directed = input$directed)
  })
  coef.diss <- reactive({
    dissolution_coefs(dissolution = as.formula(input$dissolution),
                      duration = input$dur)
  })
  target.stats <- reactive({
    eval(parse(text = paste0("c(", input$form.targets, ")")))
  })
  fit <- reactive({
    if(input$runMod == 0){return()}
    isolate({
      fit.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(fit.progress$close())

      fit.progress$set(value = 0.5, message = "Fitting model")
      netest(net(), formation = as.formula(input$formation),
           target.stats = target.stats(),
           coef.diss = coef.diss(),
           verbose = FALSE)
    })
  })
  dxsim <- reactive({
    if(input$runMod == 0){return()}
    input$runDx
    isolate({
      dx.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(dx.progress$close())

      dx.progress$set(value = 0.5, message = "Diagnosing fit")
      netdx(fit(), nsims = input$dx.nsims, nsteps = input$dx.nsteps,
            keep.tedgelist = FALSE, verbose = FALSE)
    })
  })

  param <- reactive({
    param.net(inf.prob = input$inf.prob, act.rate = input$act.rate,
              rec.rate = input$rec.rate, b.rate = input$b.rate,
              ds.rate = input$ds.rate, di.rate = input$di.rate,
              dr.rate = input$dr.rate)
  })
  init <- reactive({
    init.net(i.num = input$i.num, r.num = input$r.num)
  })
  control <- reactive({
    control.net(type = input$modtype, nsims = input$epi.nsims,
                nsteps = input$epi.nsteps, verbose.int = 0)
  })
  episim <- reactive({
    if(input$runEpi == 0){return()}
    isolate({
      epi.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(epi.progress$close())

      epi.progress$set(value = 0.5, message = "Simulating Epidemic")
      netsim(fit(), param = param(), init = init(), control = control())
    })
  })
  showqnts <- reactive({
    ifelse(input$qntsrng == 0, FALSE, input$qntsrng)
  })

  #Output objects
  output$dxplot <- renderPlot({
    par(mar = c(5, 4, 2, 2))
    if(!is.null(dxsim())){
      plot(dxsim(), type = input$dxtype)
    }
  })
  output$modeldx <- renderPrint({
    if(!is.null(dxsim())){
      dxsim()
    }
  })
  output$epiplot <- renderPlot({
    if(input$runEpi == 0){return()}
    par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
    if (input$compsel == "Compartment Prevalence") {
      plot(episim(),
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           leg = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
    if (input$compsel == "Compartment Size") {
      plot(episim(),
           popfrac = FALSE,
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           leg = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
    if (input$compsel == "Disease Incidence") {
      plot(episim(),
           y = "si.flow",
           popfrac = FALSE,
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           leg = input$showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
  })
  outputOptions(output, "epiplot", suspendWhenHidden = FALSE)
  output$sumtimeui <- renderUI({
    numericInput("sumtime", label = "Time Step",
                 value = 1, min = 1, max = input$epi.nsteps)
  })
  output$episum <- renderPrint({
    if(is.null(input$sumtime)) {return()}
    summary(episim(), at = input$sumtime)
  })


})