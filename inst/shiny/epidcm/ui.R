##
## UI File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
## Run online: http://statnet.shinyapps.io/epidcm/
##

library(shiny)
library(EpiModel)

shinyUI(fluidPage(

  titlePanel("EpiModel: Deterministic Compartmental Models"),

  sidebarLayout(
    sidebarPanel(

      h3("Instructions", style = "margin-top: 0px"),
      helpText("Click Run Model after changing model parameters",
               "or conditions."),
      actionButton(inputId = "runMod", "Run Model"),

      h4("Type", style = "margin-top: 20px"),
      selectInput(inputId = "modtype",
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


      h4("Time", style = "margin-top: 25px"),
      numericInput(inputId = "nsteps",
                   label = "Time Steps",
                   value = 500,
                   min = 0),
      numericInput(inputId = "dt",
                   label = "dt (Step Size)",
                   value = 1,
                   min = 0),
      selectInput(inputId = "nimeth",
                  label = "Integration Method",
                  choices = c("rk4", "lsoda", "euler")),

      h4("Parameters", style = "margin-top: 25px"),
      numericInput(inputId = "inf.prob",
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
      numericInput(inputId = "a.rate",
                   label = "Arrival Rate",
                   min = 0,
                   value = 0.0),
      numericInput(inputId = "ds.rate",
                   label = "Departure Rate (Sus.)",
                   min = 0,
                   value = 0.0),
      numericInput(inputId = "di.rate",
                   label = "Departure Rate (Inf.)",
                   min = 0,
                   value = 0.0),
      conditionalPanel("input.modtype == 'SIR'",
                       numericInput(inputId = "dr.rate",
                                    label = "Departure Rate (Rec.)",
                                    min = 0,
                                    value = 0.0))
    ), # End sidebarPanel

    mainPanel(
      tabsetPanel(
        tabPanel("Plot",
         h4("Plot of Model Results"),
         plotOutput(outputId = "MainPlot"),
         br(),
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
             column(5,
                    sliderInput(inputId = "alpha",
                                label = "Line Transparency",
                                min = 0.1,
                                max = 1,
                                value = 0.8,
                                step = 0.1))),
           fluidRow(
             column(3,
                    checkboxInput(inputId = "showleg",
                                  label = "Plot Legend",
                                  value = TRUE))),
           fluidRow(
             column(5,
                    downloadButton(outputId = "dlMainPlot",
                                   label = "Download PDF")))
         ) # End wellPanel
        ), # End tabPanel


        tabPanel("Summary",
         h4("Time-Specific Model Summary"),
         helpText("Select the time step of interest and the number of
                  significant digits to output in the table and compartment
                  plot below."),
         fluidRow(
           column(5,
                  numericInput(inputId = "summTs",
                               label = strong("Time Step"),
                               value = 1,
                               min = 1,
                               max = 500)),
           column(5,
                  numericInput(inputId = "summDig",
                               label = strong("Significant Digits"),
                               value = 3,
                               min = 0,
                               max = 8))),
         fluidRow(
           verbatimTextOutput(outputId = "outSummary")),
         fluidRow(
           h4("Compartment Plot", style = "margin-top: 25px"),
           plotOutput(outputId = "CompPlot"),
           downloadButton(outputId = "dlCompPlot",
                          label = "Download Plot"))
        ), # end tabPanel Summary

        tabPanel("Data",
         h4("Model Data"),
         DT::DTOutput("outData"),
         fluidRow(
           column(4,
                  numericInput(inputId = "tabdig",
                               label = "Significant Digits",
                               min = 0,
                               value = 2))),
         downloadButton(outputId = "dlData",
                        label = "Download Data")
        ), # end tabPanel Data


        tabPanel("About",
         p("This application solves and plots a deterministic, compartmental
           epidemic models (DCMs). The model simulations are driven by the",
           a("EpiModel", href =
               "http://cran.r-project.org/web/packages/EpiModel/index.html"),
           "package in R.", style = "margin-top: 25px"),
         p("Models here are limited to basic one-group homogenous mixing models
           with a limited set of parameters, initial conditions, and control
           settings. More complex models are available in the command-line
           version of EpiModel. For further details, including more background
           on the mathematics and theory behind these DCMs, please consult the
           documentation, tutorials, and workshop materials
           at the main", a("EpiModel website.", href = "http://epimodel.org/")),
         p("This web application, built with",
           a("Shiny", href = "https://shiny.rstudio.com/"), "may be lauched via
           an R session with EpiModel and Shiny installed (see the epiweb
           function), or directly on any web browser (no R needed)",
           a("here.", href = "http://statnet.shinyapps.io/epidcm/")),
         br(),
         strong("Authors"), p("Samuel M. Jenness, Department of Epidemiology,
                              Emory University"),
         p("Steven M. Goodreau, Department of Anthropology,
           University of Washington"),
         p("Martina Morris, Departments of Statistics and Sociology,
           University of Washington")
        ) # End tabPanel About
      ) # end tabsetPanel
    ) # end mainPanel
  ) # end sidebarLayout
))
