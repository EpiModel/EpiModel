##
## UI File for epiicm Shiny Application
##
## Run local: epiweb(class = "icm")
## Run online: http://statnet.shinyapps.io/epiicm/
##

library(shiny)
library(EpiModel)

shinyUI(fluidPage(

  titlePanel("EpiModel: Stochastic Individual Contact Models"),

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

      h4("Time and Simulations", style = "margin-top: 25px"),
      numericInput(inputId = "nsteps",
                   label = "Time Steps",
                   value = 500,
                   min = 0),
      numericInput(inputId = "nsims",
                   label = "Simulations",
                   value = 5,
                   min = 1),

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
    ), #End sidebarPanel

    # Main panel
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
                                step = 0.01))),
           fluidRow(
             column(5,
                    downloadButton(outputId = "dlMainPlot",
                                   label = "Download PDF")))
         ) # end wellPanel
        ), # End tabPanel Plot

        tabPanel("Summary",
         h4("Time-Specific Model Summary"),
         helpText("Select the time step of interest and the number of
                  significant digits to output in the table and compartment
                  plot below."),
         fluidRow(
           column(5,
                  numericInput(inputId = "summTs", label = strong("Time Step"),
                               value = 1, min = 1, max = 500)),
           column(5,
                  numericInput(inputId = "summDig",
                               label = strong("Significant Digits"),
                               value = 3, min = 0, max = 8))),
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
                                            "Simulations"))),
             column(4,
                    conditionalPanel("input.datasel == 'Simulations'",
                                     uiOutput("simnoControl"))))
         ), # end wellPanel
         fluidRow(
           dataTableOutput("outData")),
         fluidRow(
           column(4,
                  numericInput(inputId = "tabdig",
                               label = "Significant Digits",
                               min = 0,
                               value = 2))),
         fluidRow(
           downloadButton(outputId = "dlData",
                          label = "Download Data"))
        ), # End tabPanel Data


        tabPanel("About",
         p("This application solves and plots a stochastic individual contact
           epidemic models (ICMs), which are intended to serve as
           microsimulation analogs to deterministic compartmental models (DCMs).
           The model simulations are driven by the",
           a("EpiModel", href =
               "http://cran.r-project.org/web/packages/EpiModel/index.html"),
           "package in R.", style = "margin-top: 25px"),
         p("Models here are limited to basic one-group homogenous mixing models
           with a limited set of parameters, initial conditions, and control
           settings. More complex models are available in the command-line
           version of EpiModel. For further details, including more background
           on the mathematics and theory behind these ICMs, please consult the
           documentation, tutorials, and workshop materials at the main",
           a("EpiModel website.", href = "http://epimodel.org/")),
         p("This web application, built with",
           a("Shiny", href = "http://shiny.rstudio.com/"), "may be lauched via
           an R session with EpiModel and Shiny installed (see the epiweb
           function), or directly on any web browser (no R needed)",
           a("here.", href = "http://statnet.shinyapps.io/epiicm/")),
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
