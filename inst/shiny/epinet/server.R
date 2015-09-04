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
  fit <- reactive({
    if(input$runMod == 0){return()}
    isolate({
      fit.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(fit.progress$close())

      fit.progress$set(value = 0.5, message = "Fitting model")
      netest(net(), formation = as.formula(input$formation),
           target.stats = input$form.targets,
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
    param.net(inf.prob = input$infprob, act.rate = input$actrate,
              rec.rate = input$recrate)
  })
  init <- reactive({
    init.net(i.num = input$inum, r.num = input$rnum)
  })
  control <- reactive({
    control.net(type = input$modtype, nsims = input$epi.nsims,
                nsteps = input$epi.nsteps, verbose.int = 0)
  })
  episim <- reactive({
    if(input$runEpi == 0){return()}
    isolate(
      netsim(fit(), param = param(), init = init(), control = control())
    )
  })

  #Output objects
  output$modelsum <- renderPrint({
    if(!is.null(fit())){
      summary(fit())
    }
  })
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
  output$episum <- renderPrint({
    episim()
  })
  output$epiplot <- renderPlot({
    plot(episim())
  })

})