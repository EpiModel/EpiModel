

# Server ------------------------------------------------------------------

server <- function(input, output) {

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
      init.icm(s.num = input$s.num,
               i.num = input$i.num,
               r.num = input$r.num)
    } else {
      init.icm(s.num = input$s.num,
               i.num = input$i.num)
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
           mean.line = input$showmean,
           sim.lines = input$showsims,
           qnts = showqnts(),
           leg = input$showleg,
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
           leg = input$showleg,
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
           leg = input$showleg,
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
             leg = input$showleg,
             leg.cex = 1.1,
             lwd = 3.5)
      }
      if (input$compsel == "Compartment Size") {
        plot(mod(),
             mean.line = input$showmean,
             sim.lines = input$showsims,
             qnts = showqnts(),
             leg = input$showleg,
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
             leg = input$showleg,
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
      as.data.frame(mod())
    } else if (input$datasel == "Standard Deviations") {
      as.data.frame(mod(), out = "sd")
    } else if (input$datasel == "Simulations") {
      as.data.frame(mod(), out = "vals", sim = max(1, input$datasim))
    }
  }, options = list(aLengthMenu = c(10, 25, 50, 100), iDisplayLength = 10))
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


# UI ----------------------------------------------------------------------

ui <- shinyUI(pageWithSidebar(

  # Header
  headerPanel("EpiModel: Stochastic Individual Contact Models"),

  # Sidebar
  sidebarPanel(

    h4("Instructions"),
    helpText("Click Run Model after changing model parameters",
             "or conditions."),
    actionButton(inputId = "runMod", "Run Model"),
    br(), br(),

    h4("Model Type"),
    selectInput(inputId = "modtype",
                label = "",
                choices = c("SI", "SIR", "SIS")),
    br(),

    h4("Initial Conditions"),
    numericInput(inputId = "s.num",
                 label = "Number Susceptible",
                 value = 1000,
                 min = 0),
    numericInput(inputId = "i.num",
                 label = "Number Infected",
                 value = 1,
                 min = 0),
    conditionalPanel("input.modtype == 'SIR'",
                     numericInput(inputId = "r.num",
                                  label = "Number Recovered",
                                  value = 0,
                                  min = 0)
    ),
    p(),

    h4("Time"),
    numericInput(inputId = "nsteps",
                 label = "Timesteps",
                 value = 500,
                 min = 1),
    numericInput(inputId = "nsims",
                 label = "Simulations",
                 value = 5,
                 min = 1),
    br(),

    h4("Parameters"),
    numericInput(inputId = "inf.prob",
                 label="Transmission Probability per Act",
                 min = 0,
                 max = 1,
                 value = 0.1),
    br(),
    numericInput(inputId = "act.rate",
                 label = "Act Rate",
                 min = 0,
                 value = 0.5),
    br(),
    conditionalPanel("input.modtype != 'SI'",
                     numericInput(inputId = "rec.rate",
                                  label = "Recovery Rate",
                                  min = 0,
                                  value = 0), br()
    ),
    numericInput(inputId = "b.rate",
                 label = "Birth Rate",
                 min = 0,
                 value = 0.0),
    br(),
    numericInput(inputId = "ds.rate",
                 label = "Death Rate (Sus.)",
                 min = 0,
                 value = 0.0),
    br(),
    numericInput(inputId = "di.rate",
                 label = "Death Rate (Inf.)",
                 min = 0,
                 value = 0.0),
    br(),
    conditionalPanel("input.modtype == 'SIR'",
                     numericInput(inputId = "dr.rate",
                                  label = "Death Rate (Rec.)",
                                  min = 0,
                                  value = 0.0)
    )

  ), #End sidebarPanel

  # Main panel
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               fluidRow(
                 h4("Plot of Model Results"),
                 plotOutput(outputId = "MainPlot"),
                 br()
               ),
               wellPanel(
                 fluidRow(
                   column(5,
                          selectInput(inputId = "compsel",
                                      label = strong("Plot Selection"),
                                      choices = c("Compartment Prevalence",
                                                  "Compartment Size",
                                                  "Disease Incidence")),
                          p()
                   )
                 ),
                 fluidRow(
                   h5("Display Options"),
                   column(3,
                          checkboxInput(inputId = "showmean",
                                        label = "Mean Line",
                                        value = TRUE)
                   ),
                   column(3,
                          checkboxInput(inputId = "showsims",
                                        label = "Sim Lines",
                                        value = TRUE)
                   ),
                   column(3,
                          checkboxInput(inputId = "showleg",
                                        label = "Legend",
                                        value = TRUE)
                   )
                 ),
                 fluidRow(
                   p(),
                   column(6,
                          sliderInput(inputId = "qntsrng",
                                      label = "Quantile Band",
                                      min = 0,
                                      max = 1,
                                      value = 0.5,
                                      step = 0.01)
                   )
                 ),
                 fluidRow(
                   br(),
                   h5("Download PDF"),
                   column(5,
                          downloadButton(outputId = "dlMainPlot",
                                         label = "Download")
                   )
                 )
               )

      ), # End tabPanel Plot


      tabPanel("Summary",
               h4("Time-Specific Model Summary"),
               helpText("Select the time step of interest and the number of
                        significant digits to output in the table and compartment
                        plot below."),
               p(),
               fluidRow(
                 column(5,
                        numericInput(inputId = "summTs",
                                     label = strong("Time Step"),
                                     value = 1,
                                     min = 1,
                                     max = 500)
                 ),
                 column(5,
                        numericInput(inputId = "summDig",
                                     label = strong("Significant Digits"),
                                     value = 3,
                                     min = 0,
                                     max = 8)
                 ),
                 p()
               ),
               fluidRow(
                 verbatimTextOutput(outputId = "outSummary"),
                 br()
               ),
               fluidRow(
                 h4("Compartment Plot"),
                 plotOutput(outputId = "CompPlot"),
                 downloadButton(outputId = "dlCompPlot",
                                label = "Download Plot"),
                 br(),br()
               )
               ), # End tabPanel Summary


      tabPanel("Data",
               h4("Model Data"),
               helpText("Select output as the time-specific means or standard
                        deviations across simulations, or individual simulation
                        values (if the last, also input the desired simulation
                        number."),
               p(),
               wellPanel(
                 fluidRow(
                   column(5,
                          selectInput(inputId = "datasel",
                                      label = strong("Data Selection"),
                                      choices = c("Means",
                                                  "Standard Deviations",
                                                  "Simulations"))
                   ),
                   column(4,
                          conditionalPanel("input.datasel == 'Simulations'",
                                           #                      numericInput(inputId = "datasim",
                                           #                      label = strong("Simulation Number"),
                                           #                      value = 1,
                                           #                      min = 1)
                                           uiOutput("simnoControl")
                          )
                   )
                 )
               ),
               fluidRow(
                 dataTableOutput("outData"),
                 br(),br()
               ),
               fluidRow(
                 downloadButton(outputId = "dlData",
                                label = "Download Data")
               )
                        ), # End tabPanel Data


      tabPanel("About",
               p("This application solves and plots a stochastic individual contact epidemic models
                 (ICMs), which are intended to serve as microsimulation analogs to deterministic
                 compartmental models (DCMs). The model simulations are driven by the",
                 a("EpiModel", href = "http://cran.r-project.org/web/packages/EpiModel/index.html"),
                 "package in R."),
               p("Models here are limited to basic one-group homogenous mixing models with
                 a limited set of parameters, initial conditions, and control settings. More
                 complex models are available in the command-line version of EpiModel. For
                 further details, including more background on the mathematics and theory behind
                 these ICMs, please consult the documentation, tutorials, and workshop materials
                 at the main", a("EpiModel website.", href = "http://statnet.github.io/EpiModel")),
               p("This web application, built with",
                 a("Shiny", href="http://shiny.rstudio.com/"), "may be lauched via an R session with
                 EpiModel and Shiny installed (see the epiweb function), or directly on any web
                 browser (no R needed)", a("here.", href = "http://statnet.shinyapps.io/epiicm/")),
               br(),
               strong("Authors"), p("Samuel M. Jenness, Department of Epidemiology,
                                    University of Washington"),
               p("Steven M. Goodreau, Department of Anthropology,
                 University of Washington"),
               p("Martina Morris, Departments of Statistics and Sociology,
                 University of Washington")
               ) # End tabPanel About

               )
      )
))


# Run app -----------------------------------------------------------------

shinyApp(ui = ui, server = server)
