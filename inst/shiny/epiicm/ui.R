##
## UI File for epiicm Shiny Application
##
## Run local: epiweb(class = "icm")
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
    span("EpiModel: Stochastic Individual Contact Models"),
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
                     value = 500, min = 1, step = 50),
        numericInput("i.num", "Number Infected",
                     value = 1, min = 0, step = 1),
        conditionalPanel(
          "input.modtype == 'SIR'",
          numericInput("r.num", "Number Recovered",
                       value = 0, min = 0, step = 1)
        ),
        numericInput("nsteps", "Time Steps",
                     value = 300, min = 10, step = 50),
        numericInput("nsims", "Number of Simulations",
                     value = 5, min = 1, max = 50, step = 1),
        helpText("More simulations capture stochastic variability",
                 "but take longer to run.")
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
        downloadButton("dlMainPlot", "PDF",
                       class = "btn-sm btn-outline-secondary")
      )
    ),
    card_body(
      plotlyOutput("MainPlotly", height = "500px")
    )
  ),

  # Bottom row: Summary / Data / Guide tabs
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
        div(
          class = "mb-3",
          fluidRow(
            column(4,
              selectInput("datasel", "Data View",
                          choices = c("Means", "Standard Deviations",
                                      "Individual Simulations"))
            ),
            column(4,
              conditionalPanel(
                "input.datasel == 'Individual Simulations'",
                uiOutput("simnoControl")
              )
            )
          )
        ),
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
          tags$strong("stochastic individual contact models (ICMs)"),
          "of infectious disease transmission using the",
          tags$a("EpiModel", href = "https://www.epimodel.org/",
                 target = "_blank"),
          "package for R. ICMs are agent-based microsimulation analogs to
          deterministic compartmental models. Each person in the population
          is represented as a discrete agent, and disease transmission is a
          stochastic (random) process that depends on the contact patterns
          and transmission probability per contact event."),

        p("Because ICMs are stochastic, running the model multiple times
          with the same parameters produces different outcomes. The
          application runs multiple simulations and displays their mean
          trajectory along with the interquartile range (IQR) to show the
          distribution of possible outcomes."),

        # --- Disease Types ---
        hr(),
        h5("Disease Types"),
        p("Three compartmental structures are available, each representing
          a different natural history of infection:"),
        tags$ul(
          tags$li(tags$strong("SI (Susceptible-Infected):"),
                  "Individuals move from susceptible to infected and remain
                  infected permanently. There is no recovery. This is appropriate
                  for chronic infections like HIV (without treatment) or herpes."),
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
                  individuals recover per time step. The average duration of
                  infection is 1 / recovery rate.")
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

        # --- Stochastic Variation ---
        hr(),
        h5("Stochastic Variation and Simulations"),
        p("Unlike deterministic compartmental models (DCMs), ICMs are
          stochastic. Each simulation may produce a different epidemic
          trajectory even with identical parameters, because individual
          transmission events are governed by random draws. This means
          outcomes like peak timing, peak size, and final epidemic size
          will vary across simulations."),
        p("The", tags$strong("Number of Simulations"), "control sets how
          many independent runs are performed. More simulations give a
          better picture of the range of possible outcomes. The plot shows:"),
        tags$ul(
          tags$li(tags$strong("Mean line:"), "The average trajectory across
                  all simulations, shown as a solid line."),
          tags$li(tags$strong("IQR band:"), "The interquartile range (25th
                  to 75th percentile) shown as a shaded ribbon. This captures
                  the middle 50% of simulation outcomes at each time step.")
        ),
        p("With only one simulation, only the individual trajectory is shown."),

        # --- Scenario Presets ---
        hr(),
        h5("Scenario Presets"),
        p("Three built-in presets configure all parameters to reasonable
          values for common infectious diseases:"),
        tags$ul(
          tags$li(tags$strong("Flu-like (SIR):"),
                  "Low per-act transmission probability (0.03) with a high
                  contact rate (10 acts/time step), reflecting airborne spread.
                  Recovery in about 7 days. Population of 500."),
          tags$li(tags$strong("STI-like (SIS):"),
                  "Moderate transmission probability (0.2) with a low contact
                  rate (0.5 acts/time step), reflecting sexual transmission.
                  Slow recovery (rate = 0.01). Population of 500."),
          tags$li(tags$strong("Measles-like (SIR):"),
                  "High transmission probability (0.5) and moderate contact
                  rate (3 acts/time step), producing a high",
                  HTML("R<sub>0</sub>"),
                  "of 15. Rapid epidemic growth and quick resolution.
                  Population of 500.")
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
                  compartment sizes, prevalence, or disease incidence.
                  Mean lines and IQR bands summarize the stochastic variation
                  across simulations. Hover over any point to see exact values.
                  Use the plotly toolbar to zoom, pan, or download the plot
                  as a PNG."),
          tags$li(tags$strong("Summary:"),
                  HTML("Key epidemic metrics including R<sub>0</sub>, peak
                  infected count and timing, cumulative infections, incidence
                  rate, and (when applicable) attack rate, plus stochastic
                  variation statistics. Additional sections appear when
                  interventions or vital dynamics are enabled.")),
          tags$li(tags$strong("Data:"),
                  "The full simulation output as a searchable, sortable table.
                  Choose to view means across simulations, standard deviations,
                  or individual simulation values. Download the raw data as a
                  CSV file for further analysis.")
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

        # --- ICM vs DCM ---
        hr(),
        h5("ICM vs. DCM"),
        p("Individual contact models (ICMs) and deterministic compartmental
          models (DCMs) both track the flow of individuals between disease
          states, but they differ in key ways:"),
        tags$ul(
          tags$li(tags$strong("Stochasticity:"), "ICMs model each individual
                  agent and each contact event stochastically, while DCMs solve
                  a system of ordinary differential equations (ODEs)
                  deterministically."),
          tags$li(tags$strong("Population representation:"), "ICMs track
                  discrete individuals; DCMs track aggregate compartment sizes
                  that can take fractional values."),
          tags$li(tags$strong("Small populations:"), "ICMs naturally handle
                  small populations where integer effects and random extinction
                  matter. DCMs assume large, well-mixed populations."),
          tags$li(tags$strong("Computational cost:"), "ICMs are slower because
                  each agent is simulated individually. Use the Run Model
                  button to control when the simulation executes.")
        ),

        # --- Further Resources ---
        hr(),
        h5("Further Resources"),
        p("This application uses a subset of EpiModel's capabilities. The
          full R package supports deterministic compartmental models (DCMs)
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
