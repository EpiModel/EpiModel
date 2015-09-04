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
                   selectInput("modtype",
                               label = "Disease model",
                               choices = c("SI", "SIR", "SIS")),
                   numericInput("infprob",
                                label = "Per act transmission probability",
                                value = 0.4, max = 1, min = 0),
                   numericInput("actrate",
                                label = "Acts per partnership per time step",
                                value = 2, min = 0),
                   numericInput("recrate",
                                label = "Recovery rate (inverse of disease duration)",
                                value = 0.01, min = 0),
                   numericInput("inum",
                                label = "Initially infected",
                                value = 10, min = 0, step = 1),
                   numericInput("rnum",
                                label = "Initially recovered",
                                value = 0, min = 0, step = 1),
                   numericInput("epi.nsims",
                                label = "Simulations",
                                value = 1, min = 1),
                   numericInput("epi.nsteps",
                                label = "Time Steps per Sim",
                                value = 500, min = 1)

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