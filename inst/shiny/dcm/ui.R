shinyUI(pageWithSidebar(

  # Header
  headerPanel("EpiModel: Deterministic Compartmental Models"),

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
                 label = "Time Steps",
                 value = 500,
                 min = 0),
    numericInput(inputId = "dt",
                 label = "dt (Step Size)",
                 value = 1,
                 min = 0),
    br(), br(),

    h4("Parameters"),
    numericInput(inputId = "trans.rate",
                 label = "Transmission per Contact",
                 min = 0,
                 max = 1,
                 value = 0.5),
    br(),
    numericInput(inputId = "act.rate",
                 label = "Acts per Timestep",
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

  mainPanel(
    tabsetPanel(

      tabPanel("Plot",
         h4("Plot of Model Results"),
         plotOutput(outputId = "MainPlot"),
         br(),
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
             column(2,
                checkboxInput(inputId = "showleg",
                              label = "Legend",
                              value = TRUE)
             ),
             column(5,
                 sliderInput(inputId = "alpha",
                             label = "Line Transparency",
                             min = 0.1,
                             max = 1,
                             value = 0.8,
                             step = 0.1),
               br()
             )
           ),
           fluidRow(
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
         dataTableOutput("outData"),
         br(),br(),
         downloadButton(outputId = "dlData",
                        label = "Download Data")
      ),


      tabPanel("About",
         p("This application solves and plots a deterministic, compartmental
            epidemic models. Models here are limited to one-group models, but
            two-group models are available directly through the dcm function.
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

