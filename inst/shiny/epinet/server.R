##
## Server File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output) {

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
    isolate(
      netest(net(), formation = as.formula(input$formation),
             target.stats = input$form.targets,
             coef.diss = coef.diss())
    )
  })
  dxsim <- reactive({
    if(input$runMod == 0){return()}
    input$runDx
    isolate(
      netdx(fit(), nsims = input$dx.nsims, nsteps = input$dx.nsteps,
            keep.tedgelist = FALSE)
    )
  })

  #Output objects
  output$modelsum <- renderPrint({
    if(!is.null(fit())){
      summary(fit())
    }
  })
  output$dxplot <- renderPlot({
    if(!is.null(dxsim())){
      plot(dxsim(), type = input$dxtype)
    }
  })

})