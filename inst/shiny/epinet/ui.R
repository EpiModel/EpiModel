##
## UI File for epinet Shiny Application
##
## Run local:
## Run online:
##

library(shiny)
library(EpiModel)

shinyUI(fluidPage(

  titlePanel("EpiModel: Network Models"),

  fluidRow(
    column(4,
      tabsetPanel(
        tabPanel("Initialize Network",
                 br(),
                 wellPanel(

                   actionButton("runMod", "Fit Model & Run Diagnostics",
                                style = "margin-bottom: 10px;"),
                   fluidRow(
                     column(7, numericInput("num",
                                            label = "Number of Nodes",
                                            value = 100,
                                            min = 0)),
                     column(5, checkboxInput("directed",
                                             label = "Directed?",
                                             value = FALSE))),
                   fluidRow(
                     column(7, selectInput("formation",
                                           label = "Formation Formula",
                                           choices = c("~edges"))),
                     column(5, numericInput("form.targets",
                                            label = "Target Statistics",
                                            value = 20))),

                   fluidRow(
                     column(7, selectInput("dissolution",
                                           label = "Dissolution Formula",
                                           choices = c("~offset(edges)"))),
                     column(5, numericInput("dur",
                                            label = "Edge Durations",
                                            value = 90)))
                 ),
                 verbatimTextOutput("modelsum")
        ),
        tabPanel("Epidemic Parameters",
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
                 ))
      )

    ), #End sidebar

    # Main panel
    column(8,
      tabsetPanel(
        tabPanel("Network Diagnostics",
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
                 selectInput("dxtype",
                             label = "Plot Type",
                             choices = c("formation", "dissolution", "duration")),
                 verbatimTextOutput("modeldx")),
        tabPanel("Epidemic Simulation",
                 plotOutput("epiplot"),
                 wellPanel(
                   h4("Plot Options"),
                   fluidRow(
                     column(5,
                            selectInput(inputId = "compsel",
                                        label = strong("Plot Selection"),
                                        choices = c("Compartment Prevalence",
                                                    "Compartment Size",
                                                    "Disease Incidence")))),
                   fluidRow(
                     column(3,
                            checkboxInput(inputId = "showmean",
                                          label = "Mean Line",
                                          value = TRUE)),
                     column(3,
                            checkboxInput(inputId = "showsims",
                                          label = "Sim Lines",
                                          value = TRUE)),
                     column(3,
                            checkboxInput(inputId = "showleg",
                                          label = "Legend",
                                          value = TRUE))),
                   fluidRow(
                     column(5,
                            sliderInput(inputId = "qntsrng",
                                        label = "Quantile Band",
                                        min = 0,
                                        max = 1,
                                        value = 0.5,
                                        step = 0.01)))
                 )),
        tabPanel("Epidemic Summary",
                 uiOutput("sumtimeui"),
                 verbatimTextOutput("episum"))
      )
    )
  )
))