##
## Server File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyServer(function(input, output, session) {

# Reactive Objects --------------------------------------------------------

  ##Network Diagnostics
  net <- reactive({
    network.initialize(n = input$num, directed = FALSE)
  })

  lowestconc <- reactive({
    n <- input$num
    e <- input$edge.target
    s <- 0
    conc <- 0
    if(e > n/2){
      if(e <= n-1){
        # can put all edges on one node
        conc <- 1
      } else {
        # e >= n
        # one node has been filled
        # first extra edge creates two concurrent nodes
        conc <- 3
        leftovers <- e - n
        if(leftovers > 0){
          moreconc <- floor((-1 + sqrt(1 + 8*leftovers))/2)
        }
        conc <- conc + moreconc
      }
    }
    if(e > .3 * n){
      # add buffer 10%
      conc <- conc + (.1 * n)
      if(e > .5 * n){
        #add buffer 10% for md over 1
        conc <- conc + (.1 * n)
      }
    }

    conc/n
  })

  #link duration slider and numeric input
  observeEvent(input$meandur, {
    updateNumericInput(session, "dur",
                       label = "Edge Durations",
                       value = input$meandur)
  })
  observeEvent(input$dur, {
    updateSliderInput(session, "meandur",
                      label = "Mean Partnership Duration",
                      value = input$dur, min = 1, max = 100)
  })

  #link mean degree slider and edge target numeric input
  observeEvent(input$meandeg, {
    updateNumericInput(session, "edge.target",
                       label = "Target: edges",
                       value = input$num * input$meandeg / 2,
                       step = 0.1)
  })
  observeEvent(input$num, {
    updateNumericInput(session, "edge.target",
                       label = "Target: edges",
                       value = input$num * input$meandeg / 2,
                       step = 0.01)
  })
  observeEvent(input$edge.target, {
    updateNumericInput(session, "meandeg",
                       label = "Mean Degree",
                       value = input$edge.target * 2 / input$num,
                       min = 0.1,
                       max = 1.5,
                       step = 0.01)
  })

  #link concurrency dropdown and formation formula and nwstats to track
  observeEvent(input$conc, {
    form <- ifelse(input$conc == "Concurrency not included in model",
                   yes = "~edges",
                   no = "~edges + concurrent")
    updateSelectInput(session, "formation",
                      label = "Formation Formula",
                      choices = c("~edges", "~edges + concurrent"),
                      selected = form)
  })
  observeEvent(input$formation, {
    conc <- ifelse(input$formation == "~edges",
                   yes = "Concurrency not included in model",
                   no = "Target % concurrency")
    updateSelectInput(session, "conc",
                      label = "Concurrency Rule",
                      choices = c("Concurrency not included in model",
                                  "Target % concurrency"),
                      selected = conc)
  })
  observeEvent(input$formation, {
    if(input$formation == "~edges"){
      track <- "edges"
    } else {
      track <- c("edges", "concurrent")
    }
    updateSelectInput(session, "nwstats",
                      label = "Network Stats to Track",
                      choices = c("edges",
                                  "concurrent",
                                  "isolates",
                                  "mean degree" = "meandeg"),
                      selected = track)
  })

  #link percent concurrent with concurrent target
  observeEvent(input$percConc, {
    updateNumericInput(session, "conc.target",
                       label = "Target: concurrent",
                       value = input$percConc * input$num / 100)
  })
  observeEvent(input$num, {
    updateNumericInput(session, "conc.target",
                       label = "Target: concurrent",
                       value = input$percConc * input$num / 100)
  })
  observeEvent(input$conc.target, {
    updateSliderInput(session, "percConc",
                      label = "Percent of nodes with concurrent partners",
                      value = input$conc.target / input$num * 100,
                      min = lowestconc()*100,
                      max = 50,
                      step = 5)
  })

  coef.diss <- reactive({
    dissolution_coefs(dissolution = as.formula(input$dissolution),
                      duration = input$dur)
  })
  target.stats <- reactive({
    if(input$formation == "~edges"){
      c(input$edge.target)
    } else {
      c(input$edge.target, input$conc.target)
    }
  })
  nwstats <- reactive({
    as.formula(paste0("~", paste(input$nwstats, collapse = "+")))
  })
  fit <- reactive({
    if(input$runMod == 0){return()}
    isolate({
      fit.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(fit.progress$close())
      fit.progress$set(value = NULL, message = "Fitting model")
      fit <- tryCatch({netest(net(),
                              formation = as.formula(input$formation),
                              target.stats = target.stats(),
                              coef.diss = coef.diss(),
                              verbose = FALSE)},
                      error = function(e) e)
      validate(need(!("error" %in% class(fit)),
        message = "There is a problem with model specification, please try again."))
      fit
    })
  })
  dxsim <- reactive({
    if(input$runMod == 0){return()}
    input$runDx
    isolate({
      dx.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(dx.progress$close())
      dx.progress$set(value = NULL, message = "Diagnosing fit")
      netdx(fit(),
            nsims = input$dx.nsims,
            nsteps = input$dx.nsteps,
            nwstats.formula = nwstats(),
            keep.tedgelist = FALSE,
            verbose = FALSE)
    })
  })
  dx.showqnts <- reactive({
    ifelse(input$dx.qntsrng == 0,  FALSE, input$dx.qntsrng)
  })

  ##Epidemic Simulation
  param <- reactive({
    param.net(inf.prob = input$inf.prob,
              act.rate = input$act.rate,
              rec.rate = input$rec.rate)
  })
  init <- reactive({
    init.net(i.num = input$i.num,
             r.num = input$r.num)
  })
  control <- reactive({
    control.net(type = input$modtype,
                nsims = 1,
                nsteps = input$epi.nsteps,
                verbose = FALSE)
  })
  episim <- reactive({
    if(input$runEpi == 0){return()}
    isolate({
      epi.progress <- Progress$new(session, min = 0, max = 1)
      on.exit(epi.progress$close())
      epi.progress$set(value = 0, message = "Simulating Epidemic")
      nsims <- input$epi.nsims

      epi.progress$inc(amount = 1/nsims,
                       message = "Simulating Epidemic",
                       detail = paste0("Sim 1/", nsims))
      x <- netsim(fit(),
                  param = param(),
                  init = init(),
                  control = control())
      if(nsims > 1){
        for(i in 2:nsims){
          epi.progress$inc(amount = 1/nsims,
                           detail = paste0("Sim ", i, "/", nsims))
          y <- netsim(fit(),
                      param = param(),
                      init = init(),
                      control = control())

          x <- merge(x, y)
        }
      }
      x
    })
  })
  epi.showqnts <- reactive({
    ifelse(input$epi.qntsrng == 0, FALSE, input$epi.qntsrng)
  })

# Output objects ----------------------------------------------------------

  ## netdx page
  output$percConcSlider <- renderUI({
    sliderInput("percConc",
                "Percent of nodes with concurrent partners",
                value = 10,
                min = lowestconc()*100,
                max = 50,
                step = 5,
                post = "%")
  })

  output$dxplot <- renderPlot({
    input$runMod
    par(mar = c(5, 4, 2, 2))
    if(!is.null(dxsim())){
      stats <- isolate(input$nwstats)
      plot(dxsim(),
           stats = stats,
           type = input$dxtype,
           mean.line = input$dx.showmean,
           sim.lines = input$dx.showsims,
           qnts = dx.showqnts(),
           plots.joined = input$plots.joined,
           leg = input$dx.showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
  })
  output$dxplotDL <- downloadHandler(
    filename = "netdxplot.pdf",
    content = function(file) {
      pdf(file = file, height = 6, width = 10)
      par(mar = c(5, 4, 2, 2), mgp = c(2.1, 1, 0))
      if(!is.null(dxsim())){
        input$runMod
        stats <- isolate(input$nwstats)
        plot(dxsim(),
             stats = stats,
             type = input$dxtype,
             mean.line = input$dx.showmean,
             sim.lines = input$dx.showsims,
             qnts = dx.showqnts(),
             leg = input$dx.showleg,
             leg.cex = 1.1,
             lwd = 3.5,
             main = "")
      }
      dev.off()
  })
  output$modeldx <- renderPrint({
    if(!is.null(dxsim())){
      dxsim()
    }
  })

  ## Epi page
  output$epiplot <- renderPlot({
    if(input$runEpi == 0){return()}
    par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
    if (input$compsel == "Compartment Prevalence") {
      plot(episim(),
           mean.line = input$epi.showmean,
           sim.lines = input$epi.showsims,
           qnts = epi.showqnts(),
           leg = input$epi.showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    } else if (input$compsel == "Compartment Size") {
      plot(episim(),
           popfrac = FALSE,
           mean.line = input$epi.showmean,
           sim.lines = input$epi.showsims,
           qnts = epi.showqnts(),
           leg = input$epi.showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    } else if (input$compsel == "Disease Incidence") {
      plot(episim(),
           y = "si.flow",
           popfrac = FALSE,
           mean.line = input$epi.showmean,
           sim.lines = input$epi.showsims,
           qnts = epi.showqnts(),
           leg = input$epi.showleg,
           leg.cex = 1.1,
           lwd = 3.5,
           main = "")
    }
  })
  outputOptions(output, "epiplot", suspendWhenHidden = FALSE)
  output$epiplotDL <- downloadHandler(
    filename = "epiplot.pdf",
    content =  function(file){
      pdf(file = file, height = 6, width = 10)
      par(mar = c(3.5, 3.5, 1.2, 1), mgp = c(2.1, 1, 0))
      if (input$compsel == "Compartment Prevalence") {
        plot(episim(),
             mean.line = input$epi.showmean,
             sim.lines = input$epi.showsims,
             qnts = epi.showqnts(),
             leg = input$epi.showleg,
             leg.cex = 1.1,
             lwd = 3.5,
             main = "")
      } else if (input$compsel == "Compartment Size") {
        plot(episim(),
             popfrac = FALSE,
             mean.line = input$epi.showmean,
             sim.lines = input$epi.showsims,
             qnts = epi.showqnts(),
             leg = input$epi.showleg,
             leg.cex = 1.1,
             lwd = 3.5,
             main = "")
      } else if (input$compsel == "Disease Incidence") {
        plot(episim(),
             y = "si.flow",
             popfrac = FALSE,
             mean.line = input$epi.showmean,
             sim.lines = input$epi.showsims,
             qnts = epi.showqnts(),
             leg = input$epi.showleg,
             leg.cex = 1.1,
             lwd = 3.5,
             main = "")
      }
      dev.off()
    }
  )

  output$sumtimeui <- renderUI({
    numericInput("sumtime",
                 label = "Time Step",
                 value = 1,
                 min = 1,
                 max = input$epi.nsteps)
  })
  output$episum <- renderPrint({
    if(is.null(input$sumtime)) {return()}
    summary(episim(), at = input$sumtime)
  })

  ## nw plots page
  output$nwplot <- renderPlot({
    if(input$runEpi == 0){return()}
    simno <- ifelse(input$nwplotsim == "mean",
                    yes = "mean",
                    no = as.numeric(input$nwplotsim))
    par(mar = c(0, 0, 0, 0))
    plot(episim(),
         type = "network",
         col.status = TRUE,
         at = input$nwplottime,
         sims = simno)
  })
  output$nwplot2 <- renderPlot({
    simno <- ifelse(input$nwplotsim2 == "mean",
                    yes = "mean",
                    no = as.numeric(input$nwplotsim2))
    par(mar = c(0, 0, 0, 0))
    plot(episim(),
         type = "network",
         col.status = TRUE,
         at = input$nwplottime2,
         sims = simno)
  })

  output$nwplot1DL <- downloadHandler(
    filename = "nwplot.pdf",
    content = function(file){
      pdf(file = file, height = 6, width = 6)
      par(mar = c(0, 0, 0, 0))
      plot(episim(),
           type = "network",
           col.status = TRUE,
           at = input$nwplottime,
           sims = input$nwplotsim)
      dev.off()
    }
  )
  output$nwplot2DL <- downloadHandler(
    filename = "nwplot.pdf",
    content = function(file){
      pdf(file = file, height = 6, width = 6)
      par(mar = c(0, 0, 0, 0))
      plot(episim(),
           type = "network",
           col.status = TRUE,
           at = input$nwplottime2,
           sims = input$nwplotsim2)
      dev.off()
    }
  )

  output$plotUI <- renderUI({
    if(input$secondplot){
      out <- fluidRow(
        column(6, style = "padding-right: 0px;",
               plotOutput("nwplot")),
        column(6, style = "padding-left: 0px;",
               plotOutput("nwplot2"))
      )
    } else {
      out <- plotOutput("nwplot")
    }
    out
  })
  output$plotoptionsUI <- renderUI({
    fluidRow(
      column(6,
         selectInput("nwplotsim",
                     label = "Simulation",
                     choices = c("mean", 1:input$epi.nsims)),
         numericInput("nwplottime",
                      label = "Time Step",
                      value = 1,
                      min = 1,
                      max = input$epi.nsteps,
                      step = 1),
         downloadButton("nwplot1DL", label = "Download Plot 1")
      ),
      conditionalPanel("input.secondplot",
           column(6,
              selectInput("nwplotsim2",
                          label = "Simulation",
                          choices = c("mean", 1:input$epi.nsims)),
              numericInput("nwplottime2",
                           label = "Time Step",
                           value = 1,
                           min = 1,
                           max = input$epi.nsteps,
                           step = 1),
              downloadButton("nwplot2DL", label = "Download Plot 2")
                    )
      )
    )
  })

  ## Data page
  output$simnoControl <- renderUI({
    input$runEpi
    maxsims <- isolate(input$epi.nsims)
    numericInput(inputId = "datasim",
                label = strong("Simulation Number"),
                value = 1,
                min = 1,
                max = maxsims)
  })
  output$outData <- renderDataTable({
    if (input$datasel == "Means") {
      round(as.data.frame(episim()), input$tabdig)
    } else if (input$datasel == "Standard Deviations") {
      round(as.data.frame(episim(), out = "sd"), input$tabdig)
    } else if (input$datasel == "Simulations") {
      as.data.frame(episim(), out = "vals", sim = max(1, input$datasim))
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
