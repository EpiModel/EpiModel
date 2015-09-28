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
  tabPanel("EpiModel: Network Models",

      column(6, offset = 1,
           h2("EpiModel: Network Models", style = "color: #445555;"),
           p(a("EpiModel", href = "http://www.epimodel.org/", target = "_blank"),
             "is an R package that provides tools for simulating and
             analyzing mathematical models of infectious disease."),
           p("More text..."),
           p("More text..."),
           p("This web application, built with",
           a("Shiny,", href = "http://shiny.rstudio.com/", target = "_blank"),
             "may be lauched via an R session with EpiModel and Shiny installed
             (see the epiweb function), or directly on any web browser (no R
             needed)", a("here.", href = "https://statnet.shinyapps.io/epinet",
                         target = "_blank"))
             ),
      column(4,
             img(src = "dxplot.png", class = "transp", title = "Network Diagnostics",
                 width = "250px", height = "151px"),
             img(src = "epiplot.png", class = "transp", title = "Disease Prevalence",
                 width = "250px", height = "151px"),
             img(src = "nwplot.png", class = "transp", title = "Network Plot",
                 width = "250px", height = "164px")
             )
           ),
  tabPanel("Model Estimation",
           tagList(
             tags$head(
               tags$link(rel="stylesheet", type="text/css", href="style.css")
             )
           ),
           fluidRow(
             column(4,
                    br(),
                    wellPanel(
                      h4("Model Estimation"),

                      actionButton("runMod", "Fit Model & Run Diagnostics",
                                   style = "margin-bottom: 10px;"),
                      fluidRow(
                        column(7, numericInput("num",
                                               label = "Number of Nodes",
                                               value = 100,
                                               min = 0))
                      ),

                      hr(),
                      h4("Model Specification: Method 1"),
                      sliderInput("meandeg",
                                  label = "Mean Degree",
                                  value = 0.5,
                                  min = 0.1,
                                  max = 1.5,
                                  step = 0.1),
                      sliderInput("meandur",
                                  label = "Mean Partnership Duration",
                                  value = 50,
                                  min = 1,
                                  max = 100),
                      selectInput("conc",
                                  label = "Concurrency Rule",
                                  choices = c("No concurrency specified",
                                              "Target % concurrency")),
                      conditionalPanel("input.conc == 'Target % concurrency'",
                                       sliderInput("percConc",
                                                   "Percent of nodes with concurrent partners",
                                                   value = 10,
                                                   min = 0,
                                                   max = 50,
                                                   step = 10,
                                                   post = "%")),
                      hr(),
                      h4("Model Specification: Method 2"),
                      fluidRow(
                        column(7, selectInput("formation",
                                              label = "Formation Formula",
                                              choices = c("~edges", "~edges + concurrent"))),
                        column(5, numericInput("edge.target",
                                               label = "Target: edges",
                                               value = 25,
                                               step = 0.1),
                               conditionalPanel("input.formation == '~edges + concurrent'",
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
                                                       "mean degree" = "meandeg"),
                                           selected = "edges"))
                   ),
                   plotOutput("dxplot", height = "600px"),
                   wellPanel(
                     h4("Plot Options"),
                     fluidRow(
                       column(4,
                          selectInput("dxtype",
                               label = "Plot Type",
                               choices = c("formation", "dissolution", "duration"))),
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
                 h4("Epidemic Simulation"),
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
#                  numericInput(inputId = "b.rate",
#                               label = "Birth Rate",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  numericInput(inputId = "ds.rate",
#                               label = "Death Rate (Sus.)",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  numericInput(inputId = "di.rate",
#                               label = "Death Rate (Inf.)",
#                               min = 0,
#                               value = 0.0,
#                               step = 0.005),
#                  conditionalPanel("input.modtype == 'SIR'",
#                                   numericInput(inputId = "dr.rate",
#                                                label = "Death Rate (Rec.)",
#                                                min = 0,
#                                                value = 0.0,
#                                                step = 0.005))
               )), #end sidebar
             #main panel
             column(8,
                tabsetPanel(
                  tabPanel("Plot",
                     plotOutput("epiplot", height = "600px"),
                     wellPanel(
                       h4("Plot Options"),
                       fluidRow(
                         column(5,
                                selectInput(inputId = "compsel",
                                            label = strong("Plot Type"),
                                            choices = c("Compartment Prevalence",
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
                                              value = TRUE)),
                         column(3,
                                checkboxInput(inputId = "epi.showleg",
                                              label = "Legend",
                                              value = TRUE))),
                       fluidRow(
                         downloadButton("epiplotDL", "Download Plot")
                         )
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
           ), #end epi page
  tabPanel("Network Plots",

           uiOutput("plotUI"),
           br(),
           div(style = "margin: auto; width: 60%;",
            wellPanel(
             h4("Plot Options"),
             helpText("Look at the network plot of any epidemic simulation at
                      any time step. Plotting the mean network shows the plot
                      of the simulation that is closest to the mean prevalence
                      at each time step."),
             checkboxInput("secondplot",
                           label = "Plot two time steps",
                           value = FALSE),
             uiOutput("plotoptionsUI")
           ))

           ), #end nw plots page
  tabPanel("Data",
      div(style = "margin: auto; width: 80%;",
           h4("Model Data"),
           helpText("Select output as the time-specific means or standard
                  deviations across simulations, or individual simulation
                  values (if the last, also input the desired simulation
                  number)."),
           p(),
           wellPanel(
             fluidRow(
               column(5,
                      selectInput(inputId = "datasel",
                                  label = strong("Data Selection"),
                                  choices = c("Means",
                                              "Standard Deviations",
                                              "Simulations")),
                      conditionalPanel("input.datasel == 'Simulations'",
                                       uiOutput("simnoControl"))),
               column(4, offset = 1,
                      numericInput(inputId = "tabdig",
                                   label = "Significant Digits",
                                   min = 0,
                                   value = 2)))
           ), # end wellPanel
           fluidRow(
             dataTableOutput("outData")),
           fluidRow(
             downloadButton(outputId = "dlData",
                            label = "Download Data")),
           br()
      )
           ) #end data page
  )
)