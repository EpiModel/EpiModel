##
## UI File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
##

library(shiny)
library(bslib)
library(DT)
library(EpiModel)

# -- Theme --
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Source Sans Pro"),
  heading_font = font_google("Source Sans Pro"),
  "primary" = "#2C3E50",
  "success" = "#18BC9C",
  "info" = "#3498DB",
  "danger" = "#E74C3C"
)

# -- UI --
page_sidebar(
  title = "EpiModel: Deterministic Compartmental Models",
  theme = app_theme,

  # ===== Sidebar =====
  sidebar = sidebar(
    width = 320,

    # Run button
    actionButton("runMod", "Run Model",
                 class = "btn-primary btn-lg w-100 mb-3"),

    # Disease type and presets
    selectInput("modtype", "Disease Type",
                choices = c("SI", "SIR", "SIS")),

    selectInput("preset", "Scenario Preset",
                choices = c("Custom",
                            "Flu-like (SIR)",
                            "STI-like (SIS)",
                            "Measles-like (SIR)")),

    hr(),

    # Parameter accordion
    accordion(
      id = "params_accordion",
      open = c("Population & Time", "Epidemic Parameters"),

      # -- Population & Time --
      accordion_panel(
        "Population & Time",
        numericInput("s.num", "Number Susceptible",
                     value = 1000, min = 1, step = 100),
        numericInput("i.num", "Number Infected",
                     value = 1, min = 0, step = 1),
        conditionalPanel(
          "input.modtype == 'SIR'",
          numericInput("r.num", "Number Recovered",
                       value = 0, min = 0, step = 1)
        ),
        numericInput("nsteps", "Time Steps",
                     value = 500, min = 10, step = 50)
      ),

      # -- Epidemic Parameters --
      accordion_panel(
        "Epidemic Parameters",
        sliderInput("inf.prob", "Transmission Probability per Act",
                    min = 0, max = 1, value = 0.1, step = 0.01),
        sliderInput("act.rate", "Act Rate",
                    min = 0, max = 20, value = 0.5, step = 0.1),
        conditionalPanel(
          "input.modtype != 'SI'",
          sliderInput("rec.rate", "Recovery Rate",
                      min = 0, max = 1, value = 0.05, step = 0.01)
        )
      ),

      # -- Intervention --
      accordion_panel(
        "Intervention",
        checkboxInput("enable_intervention", "Enable Intervention",
                      value = FALSE),
        conditionalPanel(
          "input.enable_intervention",
          sliderInput("inter.eff", "Intervention Efficacy",
                      min = 0, max = 1, value = 0.5, step = 0.01),
          helpText("Proportional reduction in transmission probability."),
          numericInput("inter.start", "Intervention Start (time step)",
                       value = 100, min = 1, step = 10)
        )
      ),

      # -- Sensitivity Analysis --
      accordion_panel(
        "Sensitivity Analysis",
        checkboxInput("enable_sensitivity", "Enable Sensitivity Analysis",
                      value = FALSE),
        conditionalPanel(
          "input.enable_sensitivity",
          selectInput("sens_param", "Parameter to Vary",
                      choices = c("inf.prob", "act.rate", "rec.rate")),
          numericInput("sens_min", "Minimum Value",
                       value = 0.05, min = 0, step = 0.01),
          numericInput("sens_max", "Maximum Value",
                       value = 0.5, min = 0, step = 0.01),
          selectInput("sens_nruns", "Number of Values",
                      choices = c(3, 5, 7), selected = 5)
        )
      ),

      # -- Vital Dynamics --
      accordion_panel(
        "Vital Dynamics",
        checkboxInput("enable_vital", "Include Births & Deaths",
                      value = FALSE),
        conditionalPanel(
          "input.enable_vital",
          numericInput("a.rate", "Birth Rate (per person per time step)",
                       value = 0.001, min = 0, step = 0.001),
          numericInput("death_rate", "Death Rate (uniform, all compartments)",
                       value = 0.001, min = 0, step = 0.001)
        )
      )
    ) # end accordion
  ), # end sidebar


  # ===== Main Panel =====

  # Top row: plot + key measures
  layout_columns(
    col_widths = c(8, 4),

    # -- Main Plot Card --
    card(
      full_screen = TRUE,
      card_header(
        class = "d-flex justify-content-between align-items-center",
        "Epidemic Trajectory",
        div(
          class = "d-flex align-items-center gap-2",
          selectInput("compsel", NULL,
                      choices = c("Compartment Prevalence",
                                  "Compartment Size",
                                  "Disease Incidence"),
                      width = "220px"),
          checkboxInput("use_plotly", "Interactive", value = FALSE),
          downloadButton("dlMainPlot", "PDF", class = "btn-sm btn-outline-secondary")
        )
      ),
      card_body(
        conditionalPanel(
          "!input.use_plotly",
          plotOutput("MainPlot", height = "450px")
        ),
        conditionalPanel(
          "input.use_plotly",
          uiOutput("plotlyUI")
        )
      )
    ),

    # -- Right Column: Key Measures + Compartment Diagram --
    layout_columns(
      col_widths = 12,

      # Key Measures
      card(
        card_header("Key Measures"),
        card_body(
          uiOutput("keyMeasures")
        )
      ),

      # Compartment Diagram
      card(
        card_header("Model Diagram"),
        card_body(
          plotOutput("CompPlot", height = "220px"),
          sliderInput("summTs", "Time Step",
                      min = 1, max = 500, value = 1, step = 1,
                      width = "100%")
        )
      )
    )
  ), # end layout_columns top row

  # Bottom row: Summary / Data / About tabs
  navset_card_tab(
    id = "main_tabs",

    nav_panel(
      "Summary",
      verbatimTextOutput("outSummary")
    ),

    nav_panel(
      "Data",
      DTOutput("outData"),
      div(
        class = "mt-2",
        downloadButton("dlData", "Download CSV",
                       class = "btn-sm btn-outline-secondary")
      )
    ),

    nav_panel(
      "About",
      div(
        class = "p-3",
        h4("About This Application"),
        p("This application solves and visualizes deterministic compartmental
          epidemic models (DCMs) using the",
          tags$a("EpiModel", href = "https://www.epimodel.org/",
                 target = "_blank"),
          "package for R."),
        p("Models include SI, SIR, and SIS disease types with options for
          vital dynamics (births and deaths) and public health interventions.
          Explore how epidemic parameters affect disease spread by adjusting
          the controls in the sidebar."),
        p("For more complex models, including stochastic individual contact
          models and network models, see the command-line version of EpiModel.
          Tutorials and documentation are available at the",
          tags$a("EpiModel website.", href = "https://www.epimodel.org/",
                 target = "_blank")),
        p("Source code and bug reports:",
          tags$a("GitHub", href = "https://github.com/EpiModel/EpiModel",
                 target = "_blank")),
        p("R package on CRAN:",
          tags$a("EpiModel",
                 href = "https://cran.r-project.org/package=EpiModel",
                 target = "_blank")),
        hr(),
        p(tags$strong("Citation:"),
          "Jenness SM, Goodreau SM, Morris M (2018).",
          tags$em("EpiModel: An R Package for Mathematical Modeling of
                   Infectious Disease over Networks."),
          "Journal of Statistical Software, 84(8), 1-47.",
          tags$a("doi:10.18637/jss.v084.i08",
                 href = "https://doi.org/10.18637/jss.v084.i08",
                 target = "_blank")),
        hr(),
        p(tags$strong("Authors")),
        p("Samuel M. Jenness, Department of Epidemiology, Emory University"),
        p("Steven M. Goodreau, Department of Anthropology, University of Washington"),
        p("Martina Morris, Departments of Statistics and Sociology, University of Washington")
      )
    )
  ) # end navset_card_tab
) # end page_sidebar
