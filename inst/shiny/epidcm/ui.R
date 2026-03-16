##
## UI File for epidcm Shiny Application
##
## Run local: epiweb(class = "dcm")
##

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(EpiModel)

# -- Theme --
app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Atkinson Hyperlegible"),
  heading_font = font_google("Atkinson Hyperlegible"),
  "primary" = "#2C3E50",
  "success" = "#18BC9C",
  "info" = "#3498DB",
  "danger" = "#E74C3C"
)

# -- UI --
page_sidebar(
  title = div(
    class = "d-flex justify-content-between align-items-center w-100",
    span("EpiModel: Deterministic Compartmental Models"),
    actionButton("runMod", "Run Model",
                 class = "btn-light btn-sm",
                 icon = icon("play"))
  ),
  theme = app_theme,

  # ===== Sidebar =====
  sidebar = sidebar(
    width = 320,
    open = "always",

    # Parameter accordion
    accordion(
      id = "params_accordion",
      open = c("Disease Type or Scenario", "Population & Time",
               "Epidemic Parameters"),

      # -- Disease Type or Scenario --
      accordion_panel(
        "Disease Type or Scenario",
        selectInput("modtype", "Disease Type",
                    choices = c("SI", "SIR", "SIS")),
        selectInput("preset", "Scenario Preset",
                    choices = c("Custom",
                                "Flu-like (SIR)",
                                "STI-like (SIS)",
                                "Measles-like (SIR)"))
      ),

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
                      choices = c(3, 5, 7), selected = 5),
          helpText("These values override the fixed parameter value set above.")
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
          numericInput("ds.rate", "Death Rate",
                       value = 0.001, min = 0, step = 0.001),
          checkboxInput("diff_death_rates",
                        "Different death rates by compartment",
                        value = FALSE),
          conditionalPanel(
            "input.diff_death_rates",
            numericInput("di.rate", "Death Rate, Infected",
                         value = 0.001, min = 0, step = 0.001),
            conditionalPanel(
              "input.modtype == 'SIR'",
              numericInput("dr.rate", "Death Rate, Recovered",
                           value = 0.001, min = 0, step = 0.001)
            )
          )
        )
      )
    ) # end accordion
  ), # end sidebar


  # ===== Main Panel =====

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
                    width = "330px"),
        downloadButton("dlMainPlot", "PDF", class = "btn-sm btn-outline-secondary")
      )
    ),
    card_body(
      plotlyOutput("MainPlotly", height = "500px")
    )
  ),

  # Bottom row: Summary / Data / About tabs
  navset_card_tab(
    id = "main_tabs",
    height = "500px",

    nav_panel(
      "Summary",
      div(class = "p-3", style = "overflow-y: auto; height: 100%;",
          uiOutput("outSummary"))
    ),

    nav_panel(
      "Data",
      div(style = "overflow-y: auto; height: 100%;",
        DTOutput("outData"),
        div(
          class = "mt-2",
          downloadButton("dlData", "Download CSV",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    nav_panel(
      "Guide",
      div(
        class = "p-3",
        style = "max-width: 900px; overflow-y: auto; height: 100%;",

        # --- Overview ---
        h4("User Guide"),
        p("This application solves and visualizes",
          tags$strong("deterministic compartmental models (DCMs)"),
          "of infectious disease transmission using the",
          tags$a("EpiModel", href = "https://www.epimodel.org/",
                 target = "_blank"),
          "package for R. DCMs represent populations as aggregate compartments
          and track the flow of individuals between disease states over time.
          The dynamics are governed by ordinary differential equations (ODEs)
          solved numerically with an adaptive step-size solver (lsoda)."),

        # --- Disease Types ---
        hr(),
        h5("Disease Types"),
        p("Three compartmental structures are available, each representing
          a different natural history of infection:"),
        tags$ul(
          tags$li(tags$strong("SI (Susceptible-Infected):"),
                  "Individuals move from susceptible to infected and remain
                  infected and infectious permanently. There is no recovery.
                  This is appropriate for chronic infections like HIV (without
                  treatment) or herpes."),
          tags$li(tags$strong("SIR (Susceptible-Infected-Recovered):"),
                  "Infected individuals recover and gain lasting immunity.
                  Once recovered, they cannot be re-infected. This applies to
                  many acute viral infections like measles, influenza, and
                  SARS-CoV-2 (as a simplification)."),
          tags$li(tags$strong("SIS (Susceptible-Infected-Susceptible):"),
                  "Infected individuals recover but return to the susceptible
                  state with no lasting immunity. They can be re-infected.
                  This is typical for many bacterial STIs like gonorrhea
                  and chlamydia.")
        ),

        # --- Epidemic Parameters ---
        hr(),
        h5("Epidemic Parameters"),
        tags$ul(
          tags$li(tags$strong("Transmission Probability per Act:"),
                  "The probability that an infection is transmitted during a
                  single contact (act) between an infected and a susceptible
                  individual. Ranges from 0 to 1."),
          tags$li(tags$strong("Act Rate:"),
                  "The average number of contact events (acts) per person
                  per unit time. Higher act rates mean more opportunities for
                  transmission."),
          tags$li(tags$strong("Recovery Rate:"),
                  "(SIR and SIS only.) The per-capita rate at which infected
                  individuals recover per time step. In a closed population,
                  the average duration of infection is 1 / recovery rate. With
                  vital dynamics enabled, the average duration in the infected
                  state is 1 / (recovery rate + departure rate for infecteds).")
        ),
        p("Together, these parameters determine the",
          tags$strong("basic reproduction number"),
          HTML("(R<sub>0</sub>),"),
          "the average number of secondary infections caused by one infected
          individual in a fully susceptible population:"),
        p(class = "text-center",
          style = "font-size: 1.1rem;",
          HTML("R<sub>0</sub> = (transmission probability &times;
                act rate) / recovery rate")),
        p(HTML("When R<sub>0</sub> > 1, the epidemic will grow.
          When R<sub>0</sub> < 1, the infection will die out.")),

        # --- Scenario Presets ---
        hr(),
        h5("Scenario Presets"),
        p("Three built-in presets configure all parameters to reasonable
          values for common infectious diseases. In all presets, each time
          step represents one day:"),
        tags$ul(
          tags$li(tags$strong("Flu-like (SIR):"),
                  "Low per-act transmission probability (0.03) with a high
                  contact rate (10 acts/day), reflecting airborne spread.
                  Average duration of infection of about 7 days (recovery
                  rate = 0.15). Initial population of 10,001."),
          tags$li(tags$strong("STI-like (SIS):"),
                  "Moderate transmission probability (0.2) with a low contact
                  rate (0.5 acts/day), reflecting sexual transmission.
                  Average duration of infection of about 100 days (recovery
                  rate = 0.01). Initial population of 1,010."),
          tags$li(tags$strong("Measles-like (SIR):"),
                  "High transmission probability (0.5) and moderate contact
                  rate (3 acts/day), producing a high",
                  HTML("R<sub>0</sub>"),
                  "of 15. Average duration of infection of about 10 days
                  (recovery rate = 0.1). Initial population of 10,001.")
        ),
        p("Select", tags$em("Custom"), "to set parameters manually."),

        # --- Interventions ---
        hr(),
        h5("Interventions"),
        p("When enabled, an intervention reduces the transmission probability
          by a proportional amount starting at a specified time step. For
          example, an efficacy of 0.5 cuts the transmission probability in
          half. This can represent public health measures such as vaccination
          campaigns, mask mandates, contact tracing, or behavioral
          interventions."),
        p("The intervention effect is shown as a vertical dashed line on the
          plot. Compare epidemic trajectories with and without the intervention
          to assess its impact."),

        # --- Sensitivity Analysis ---
        hr(),
        h5("Sensitivity Analysis"),
        p("Sensitivity analysis runs the model multiple times, each with a
          different value for a selected parameter. The parameter is varied
          across evenly spaced values between a minimum and maximum that you
          specify. This reveals how sensitive the epidemic outcome is to
          uncertainty in that parameter."),
        p("The range of values specified in the sensitivity analysis panel
          overrides whatever fixed value is set for that parameter in the
          Epidemic Parameters section above."),
        p("When sensitivity analysis is active, the plot shows one trajectory
          per run (colored from blue to red), with each line labeled by its
          parameter value. The Summary tab reports statistics for the first
          run only."),

        # --- Vital Dynamics ---
        hr(),
        h5("Vital Dynamics"),
        p("By default, the population is closed (no births or deaths). Enabling
          vital dynamics adds a constant per-capita birth rate (new susceptibles
          entering the population) and per-capita death rates. By default, a
          single death rate is applied uniformly across all compartments.
          Checking", tags$em("Different death rates by compartment"),
          "allows setting separate death rates for susceptible, infected,
          and (for SIR models) recovered individuals. This is useful when
          disease-induced mortality differs from background mortality. Vital
          dynamics are important for modeling endemic equilibria over longer
          time horizons, where demographic turnover replenishes the susceptible
          pool."),

        # --- Reading the Output ---
        hr(),
        h5("Reading the Output"),
        p("The application provides three output views:"),
        tags$ul(
          tags$li(tags$strong("Plot:"),
                  "Interactive time-series plots (powered by plotly) showing
                  compartment sizes, prevalence, disease incidence, or specific
                  flow variables. Hover over any point to see exact values.
                  Use the plotly toolbar to zoom, pan, or download the plot
                  as a PNG."),
          tags$li(tags$strong("Summary:"),
                  HTML("Key epidemic metrics including R<sub>0</sub>, peak
                  infected count and timing, cumulative infections, incidence
                  rate, and (when applicable) attack rate. Additional sections
                  appear when interventions, vital dynamics, or sensitivity
                  analysis are enabled.")),
          tags$li(tags$strong("Data:"),
                  "The full simulation output as a searchable, sortable table.
                  Download the raw data as a CSV file for further analysis.")
        ),

        tags$h6(class = "fw-semibold mt-3", "Incidence Metrics"),
        p("The summary tab reports two measures of disease incidence:"),
        tags$ul(
          tags$li(tags$strong("Incidence Rate:"),
                  "Total new infections divided by total susceptible
                  person-time (reported per 1,000 person-timesteps). This
                  measure is always shown and is valid for all model types,
                  including SIS models with re-infection and models with vital
                  dynamics, because it accounts for the changing size of the
                  susceptible population over time."),
          tags$li(tags$strong("Attack Rate:"),
                  "Total new infections divided by the initial number of
                  susceptibles. This is a simple proportion representing the
                  fraction of the original susceptible population that became
                  infected. It is only displayed for SI and SIR models without
                  vital dynamics, because it requires a closed population with
                  no re-infection. In SIS models (where individuals can be
                  re-infected) or models with births and deaths (where the
                  susceptible pool changes), the attack rate is not
                  well-defined and is therefore hidden.")
        ),

        # --- Further Resources ---
        hr(),
        h5("Further Resources"),
        p("This application uses a subset of EpiModel's capabilities. The
          full R package supports stochastic individual contact models (ICMs)
          and stochastic network models built on exponential-family random
          graph models (ERGMs). For tutorials and documentation, visit the",
          tags$a("EpiModel website.", href = "https://www.epimodel.org/",
                 target = "_blank")),
        tags$ul(
          tags$li("Source code and bug reports:",
                  tags$a("GitHub",
                         href = "https://github.com/EpiModel/EpiModel",
                         target = "_blank")),
          tags$li("R package on CRAN:",
                  tags$a("EpiModel",
                         href = "https://cran.r-project.org/package=EpiModel",
                         target = "_blank")),
          tags$li("Course materials:",
                  tags$a("SISMID EpiModel Workshop",
                         href = "https://epimodel.github.io/sismid/",
                         target = "_blank"))
        ),

        # --- Citation ---
        hr(),
        h5("Citation"),
        p("Jenness SM, Goodreau SM, Morris M (2018).",
          tags$em("EpiModel: An R Package for Mathematical Modeling of
                   Infectious Disease over Networks."),
          "Journal of Statistical Software, 84(8), 1-47.",
          tags$a("doi:10.18637/jss.v084.i08",
                 href = "https://doi.org/10.18637/jss.v084.i08",
                 target = "_blank")),

        # --- Authors ---
        hr(),
        h5("Authors"),
        p("Samuel M. Jenness, Department of Epidemiology, Emory University"),
        p("Steven M. Goodreau, Department of Anthropology, University of Washington"),
        p("Martina Morris, Departments of Statistics and Sociology, University of Washington")
      )
    )
  ) # end navset_card_tab
) # end page_sidebar
