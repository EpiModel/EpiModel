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

  sidebarLayout(
    sidebarPanel(

      h3("Instructions", style = "margin-top: 0px"),
      helpText("Click Run Model after changing model parameters",
               "or conditions."),
      actionButton("runMod", "Run Model"),

      h4("Initialize Network", style = "margin-top: 25px"),
      numericInput("num",
                   label = "Number of Nodes",
                   value = 100,
                   min = 0),
      checkboxInput("directed",
                    label = "Directed?",
                    value = FALSE),

      selectInput("formation",
                  label = "Formation Formula",
                  choices = c("~edges")),
      numericInput("form.targets",
                   label = "Target Statistics",
                   value = 20),

      selectInput("dissolution",
                  label = "Dissolution Formula",
                  choices = c("~offset(edges)")),
      numericInput("dur",
                   label = "Edge Durations",
                   value = 90)

    ), #End sidebarPanel

    # Main panel
    mainPanel(
      tabsetPanel(
        tabPanel("Net Diagnostics",
                 fluidRow(
                   column(3, offset = 1,
                          numericInput("dx.nsims",
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
                 verbatimTextOutput("modelsum"))
      )
    )
  )
))