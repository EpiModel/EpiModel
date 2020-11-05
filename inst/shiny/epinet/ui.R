##
## UI File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyUI(
navbarPage(title = NULL, windowTitle = "EpiModel: Network Models",
  tabPanel("About",

      column(6, offset = 1,
           h2("Stochastic Network Models with EpiModel",
              style = "color: #445555;"),
           p(a("EpiModel", href = "http://www.epimodel.org/",
               target = "_blank"),
             "is an R package that provides tools for simulating and
             analyzing mathematical models of infectious disease.
             Details about the package, including the epidemic model classes
             supported by the software can be found at the link above."),
           p("This web-based application allows for simple modeling of epidemics
              over dynamic contact networks. These stochastic network models are
             based on the statistical framework of",
             a("temporal exponential random graph models.",
               href = "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3891677/",
             target = "_blank"), "This web application is built with",
             a("Shiny,", href = "http://shiny.rstudio.com/", target = "_blank"),
             "and may be lauched via an R session with EpiModel and Shiny
             installed (see the", code("epiweb"), "function), or directly on any
             web browser (no R needed)",
             a("here.", href = "https://statnet.shinyapps.io/epinet",
               target = "_blank")),
           p("To get started, create a statistical network model in the Model
            Estimation page using one of the two model specification methods.
            This page fits a temporal ERGM using the", code("netest"),
            "function and runs diagnostics on the fitted model with the",
            code("netdx"), "function. After the model is properly specified,
            simulate an epidemic on the network using the Epidemic Simulation
            page. This runs the", code("netsim"), "function in EpiModel, and the
            epidemic parameters are described in detail in the help pages there.
            Model output may be plotted to show the epidemic time series or
            static network plots, as well as viewing numerical data
            summaries."),
           p("The author of this application is Emily Beylerian, Software
            Developer at the University of Washington Centers for Studies in
            Demoraphy and Ecology. The authors of the larger EpiModel project
            are Samuel Jenness at Emory University, and Steven Goodreau and
            Martina Morris at the University of Washington. Development of this
            software is supported by the following grants from the National
            Institutes of Health: R01HD68395 (NICHD), T32HD007543 (NICHD), and
             R24HD042828 (NICHD).")
             )
           ),
  tabPanel("Network Model Estimation",
           tagList(
             tags$head(
               tags$link(rel = "stylesheet", type = "text/css",
                         href = "style.css")
             )
           ),
           fluidRow(
             column(4,
                    br(),
                    wellPanel(
                      h4(tags$u("Network Model Estimation")),

                      actionButton("runMod", "Fit Model & Run Diagnostics",
                                   style = "margin-bottom: 10px;"),
                      fluidRow(
                        column(7, numericInput("num",
                                               label = "Number of Nodes",
                                               value = 100,
                                               min = 0))
                      ),

                      hr(),
                      h4("Specification: Method 1", style = "margin-bottom:0;"),
                      helpText("Summary Stat Targets", style = "margin-top:0;"),

                      sliderInput("meandeg",
                                  label = "Mean Degree",
                                  value = 0.5,
                                  min = 0.1,
                                  max = 1.5,
                                  step = 0.01),
                      sliderInput("meandur",
                                  label = "Mean Partnership Duration",
                                  value = 50,
                                  min = 1,
                                  max = 100),
                      selectInput("conc",
                                  label = "Concurrency Rule",
                                  choices =
                                    c("Concurrency not included in model",
                                      "Target % concurrency")),
                      conditionalPanel("input.conc == 'Target % concurrency'",
                                       helpText("Note: A mean degree greater
                                                than one always implies some
                                                level of concurrency. The model
                                                will not be run if concurrency
                                                is too low for the chosen mean
                                                degree."),
                                       uiOutput("percConcSlider")
                                       ),
                      hr(),
                      h4("Specification: Method 2", style = "margin-bottom:0;"),
                      helpText("Model and NW Stat Targets",
                               style = "margin-top:0;"),
                      fluidRow(
                        column(7, selectInput("formation",
                                              label = "Formation Formula",
                                              choices =
                                                c("~edges",
                                                  "~edges + concurrent"))),
                        column(5, numericInput("edge.target",
                                               label = "Target: edges",
                                               value = 25,
                                               step = 0.1),
                               conditionalPanel("input.formation ==
                                                '~edges + concurrent'",
                                  numericInput("conc.target",
                                               label = "Target: concurrent",
                                               value = 10)
                                                )
                               )),

                      fluidRow(
                        column(7, selectInput("dissolution",
                                              label = "Dissolution Formula",
                                              choices = c("~offset(edges)"))),
                        column(5, numericInput("dur",
                                               label = "Edge Durations",
                                               value = 50)))
                    )), #end sidebar
             #main panel
             column(8,
                   br(),
                   fluidRow(
                     column(3, numericInput("dx.nsims",
                                            label = "Simulations",
                                            value = 5, min = 1),
                            actionButton("runDx",
                                         label = "Re-Run Diagnostics")),

                     column(3, numericInput("dx.nsteps",
                                            label = "Time Steps per Sim",
                                            value = 500, min = 1)),
                     column(3, selectInput("nwstats",
                                           label = "Network Stats to Track",
                                           multiple = TRUE,
                                           choices = c("edges",
                                                       "concurrent",
                                                       "isolates",
                                                       "mean degree" =
                                                         "meandeg"),
                                           selected = "edges"))
                   ),
                   plotOutput("dxplot", height = "600px"),
                   wellPanel(
                     h4("Plot Options"),
                     fluidRow(
                       column(4,
                          selectInput("dxtype",
                               label = "Plot Type",
                               choices = c("formation", "dissolution",
                                           "duration"))),
                       column(5,
                              sliderInput(inputId = "dx.qntsrng",
                                          label = "Quantile Band",
                                          min = 0,
                                          max = 1,
                                          value = 0.5,
                                          step = 0.01))
                       ),
                       fluidRow(
                         column(3,
                                checkboxInput(inputId = "plots.joined",
                                              label = "Join Plots",
                                              value = TRUE)),
                         column(3,
                                checkboxInput(inputId = "dx.showmean",
                                              label = "Mean Line",
                                              value = TRUE)),
                         column(3,
                                checkboxInput(inputId = "dx.showsims",
                                              label = "Sim Lines",
                                              value = FALSE)),
                         column(3,
                                checkboxInput(inputId = "dx.showleg",
                                              label = "Legend",
                                              value = FALSE))),
                     fluidRow(
                       column(3,
                            downloadButton("dxplotDL", label = "Download Plot"))
                     )
                   ),

                   verbatimTextOutput("modeldx"))
           )
          ),
  tabPanel("Epidemic Simulation",
           fluidRow(
             column(4,
               br(),
               wellPanel(
                 h4(tags$u("Epidemic Simulation")),
                 actionButton("runEpi", label = "Simulate Epidemic",
                              style = "margin-bottom: 10px"),
                 helpText("Click the button above after changing model",
                          "parameters or conditions."),
                 selectInput("modtype",
                             label = "Disease Type",
                             choices = c("SI", "SIR", "SIS")),

                 h4("Initial Conditions", style = "margin-top: 25px;"),
                 numericInput(inputId = "i.num",
                              label = "Number Infected",
                              value = 1,
                              min = 0),
                 conditionalPanel("input.modtype == 'SIR'",
                                  numericInput(inputId = "r.num",
                                               label = "Number Recovered",
                                               value = 0,
                                               min = 0)),

                 h4("Time and Simulations", style = "margin-top: 25px;"),
                 numericInput("epi.nsims",
                              label = "Simulations",
                              value = 5, min = 1),
                 numericInput("epi.nsteps",
                              label = "Time Steps per Sim",
                              value = 500, min = 0),

                 h4("Parameters", style = "margin-top: 25px;"),
                 numericInput("inf.prob",
                              label = "Transmission Probability per Act",
                              min = 0,
                              max = 1,
                              value = 0.1,
                              step = 0.01),
                 numericInput(inputId = "act.rate",
                              label = "Act Rate",
                              min = 0,
                              value = 0.5,
                              step = 0.01),
                 conditionalPanel("input.modtype != 'SI'",
                                  numericInput(inputId = "rec.rate",
                                               label = "Recovery Rate",
                                               min = 0,
                                               value = 0,
                                               step = 0.01))
#                  numericInput(inputId = "a.rate",
#                               label = "Arrival Rate",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  numericInput(inputId = "ds.rate",
#                               label = "Departure Rate (Sus.)",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  numericInput(inputId = "di.rate",
#                               label = "Departure Rate (Inf.)",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  conditionalPanel("input.modtype == 'SIR'",
#                                   numericInput(inputId = "dr.rate",
#                                                label =
#                                                   "Departure Rate (Rec.)",
#                                                min = 0,
#                                                value = 0.0,
#                                                step = 0.005))
               )), #end sidebar
             #main panel
             column(8,
                tabsetPanel(
                  tabPanel("Time Series Plots",
                     plotOutput("epiplot", height = "600px"),
                     wellPanel(
                       h4("Plot Options"),
                       fluidRow(
                         column(5,
                                selectInput(inputId = "compsel",
                                            label = strong("Plot Type"),
                                            choices =
                                              c("Compartment Prevalence",
                                                "Compartment Size",
                                                "Disease Incidence"))),
                         column(5,
                                sliderInput(inputId = "epi.qntsrng",
                                            label = "Quantile Band",
                                            min = 0,
                                            max = 1,
                                            value = 0.5,
                                            step = 0.01))
                         ),
                       fluidRow(
                         column(3,
                                checkboxInput(inputId = "epi.showmean",
                                              label = "Mean Line",
                                              value = TRUE)),
                         column(3,
                                checkboxInput(inputId = "epi.showsims",
                                              label = "Sim Lines",
                                              value = FALSE)),
                         column(3,
                                checkboxInput(inputId = "epi.showleg",
                                              label = "Legend",
                                              value = TRUE))),
                       fluidRow(
                         downloadButton("epiplotDL", "Download Plot")
                         )
                       )
                           ),
                  tabPanel("Network Plots",
                           uiOutput("plotUI"),
                           br(),
                           wellPanel(
                             h4("Plot Options"),
                             helpText("Plotting the mean network shows the plot
                                      of the simulation that is closest to
                                      the mean prevalence at each time step."),
                             checkboxInput("secondplot",
                                           label = "Plot two time steps",
                                           value = FALSE),
                             uiOutput("plotoptionsUI")
                             )
                           ),
                  tabPanel("Data",
                           div(style = "margin: auto; width: 90%;",
                               br(),
                               helpText("Select output as the time-specific
                                        means or standard deviations across
                                        simulations, or individual simulation
                                        values (if the last, also input the
                                        desired simulation number)."),
                               fluidRow(
                                 column(3,
                                        selectInput(inputId = "datasel",
                                                    label =
                                                      strong("Data Selection"),
                                                    choices =
                                                      c("Means",
                                                        "Standard Deviations",
                                                        "Simulations"))),
                                 conditionalPanel("input.datasel ==
                                                  'Simulations'",
                                        column(3,
                                               uiOutput("simnoControl"))),
                                 column(3,
                                        numericInput(inputId = "tabdig",
                                                     label =
                                                       "Significant Digits",
                                                     min = 0,
                                                     value = 2))),
                               fluidRow(
                                dataTableOutput("outData")),
                               fluidRow(
                             downloadButton(outputId = "dlData",
                                            label = "Download Data"))
                           )
                           ),
                  tabPanel("Summary",
                       br(),
                       uiOutput("sumtimeui"),
                       verbatimTextOutput("episum")
                           )
                )

              ) #end main panel
            )
           ) #end epi page
  )
)
