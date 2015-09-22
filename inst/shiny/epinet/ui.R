##
## UI File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyUI(
navbarPage("EpiModel: Network Models",

  tabPanel("Network Diagnostics",
           tagList(
             tags$head(
               tags$link(rel="stylesheet", type="text/css",href="style.css")
             )
           ),
           fluidRow(
             column(4,
                    br(),
                    wellPanel(

                      actionButton("runMod", "Fit Model & Run Diagnostics",
                                   style = "margin-bottom: 10px;"),
                      fluidRow(
                        column(7, numericInput("num",
                                               label = "Number of Nodes",
                                               value = 100,
                                               min = 0))
                      ),
                      fluidRow(
                        column(7, selectInput("formation",
                                              label = "Formation Formula",
                                              choices = c("~edges", "~edges + concurrent"))),
                        column(5, numericInput("edge.target",
                                               label = "Target: edges",
                                               value = 20),
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
                                               value = 90)))
                    )), #end sidebar
             #main panel
             column(8,
                   br(),
                   fluidRow(
                     column(3, numericInput("dx.nsims",
                                            label = "Simulations",
                                            value = 5, min = 1)),

                     column(3, numericInput("dx.nsteps",
                                            label = "Time Steps per Sim",
                                            value = 500, min = 1)),
                     column(3, actionButton("runDx",
                                            label = "Re-Run Diagnostics",
                                            style = "margin-top: 25px;"))
                   ),
                   plotOutput("dxplot"),
                   wellPanel(
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
                 actionButton("runEpi", label = "Simulate Epidemic",
                              style = "margin-bottom: 10px"),
                 helpText("Click the button above after changing model",
                          "parameters or conditions."),
                 selectInput("modtype",
                             label = "Disease Type",
                             choices = c("SI", "SIR", "SIS")),

                 h4("Initial Conditions", style = "margin-top: 25px"),
                 numericInput(inputId = "i.num",
                              label = "Number Infected",
                              value = 1,
                              min = 0),
                 conditionalPanel("input.modtype == 'SIR'",
                                  numericInput(inputId = "r.num",
                                               label = "Number Recovered",
                                               value = 0,
                                               min = 0)),

                 h4("Time and Simulations", style = "margin-top: 25px"),
                 numericInput("epi.nsims",
                              label = "Simulations",
                              value = 5, min = 1),
                 numericInput("epi.nsteps",
                              label = "Time Steps per Sim",
                              value = 500, min = 0),

                 h4("Parameters", style = "margin-top: 25px"),
                 numericInput("inf.prob",
                              label = "Transmission Probability per Act",
                              min = 0,
                              max = 1,
                              value = 0.1),
                 numericInput(inputId = "act.rate",
                              label = "Act Rate",
                              min = 0,
                              value = 0.5),
                 conditionalPanel("input.modtype != 'SI'",
                                  numericInput(inputId = "rec.rate",
                                               label = "Recovery Rate",
                                               min = 0,
                                               value = 0)),
                 numericInput(inputId = "b.rate",
                              label = "Birth Rate",
                              min = 0,
                              value = 0.0),
                 numericInput(inputId = "ds.rate",
                              label = "Death Rate (Sus.)",
                              min = 0,
                              value = 0.0),
                 numericInput(inputId = "di.rate",
                              label = "Death Rate (Inf.)",
                              min = 0,
                              value = 0.0),
                 conditionalPanel("input.modtype == 'SIR'",
                                  numericInput(inputId = "dr.rate",
                                               label = "Death Rate (Rec.)",
                                               min = 0,
                                               value = 0.0))
               )), #end sidebar
             #main panel
             column(8,
                tabsetPanel(
                  tabPanel("Plot",
                     plotOutput("epiplot"),
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
             checkboxInput("secondplot", label = "Plot two time steps",
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
                  number."),
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
           ), #end data page
  tabPanel("About")
  )
)