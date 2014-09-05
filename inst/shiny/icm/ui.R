shinyUI(pageWithSidebar(

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

  ),

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

      ),


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
      ),


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
      ),


      tabPanel("About",
         p("This application solves and plots a stochastic individual contact
            epidemic models. Models here are limited to one-group models, but
            two-group models are available directly through the icm function.
            The underlying modeling software for this application is the",
            a("EpiModel", href="http://cran.r-project.org/web/packages/EpiModel/index.html"),
            "package in R. The web application is built with",
            a("Shiny.", href="http://www.rstudio.com/shiny/")),
         br(),
         strong("Author"), p("Samuel M. Jenness, Department of Epidemiology,
                              University of Washington")
      )
    )
  )
))
