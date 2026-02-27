
#' @title EpiModel Web
#'
#' @description Launches interactive Shiny applications for exploring
#'              deterministic compartmental models (DCMs) and stochastic
#'              individual contact models (ICMs) of infectious disease
#'              transmission.
#'
#' @param class Model class, with options of `"dcm"` and `"icm"`.
#' @param ... Additional arguments passed to [shiny::runApp], such as
#'        `port`, `host`, or `launch.browser`.
#'
#' @details
#' `epiweb` launches a web-based graphical interface for configuring, running,
#' and analyzing epidemic models built on EpiModel. Each application provides
#' interactive controls for model parameters, real-time (or on-demand)
#' visualization with plotly, narrative model summaries, and downloadable
#' data tables. Both applications include a **Guide** tab with a comprehensive
#' user manual covering the underlying theory, parameter definitions, and
#' instructions for interpreting output.
#'
#' @section DCM Application (`class = "dcm"`):
#' The DCM app solves systems of ordinary differential equations for one-group
#' SI, SIR, and SIS models. Models run instantly and the interface is fully
#' reactive (no run button). Features include:
#' \itemize{
#'   \item **Scenario presets:** Flu-like (SIR), STI-like (SIS), and
#'     Measles-like (SIR) configurations that set all parameters to
#'     epidemiologically plausible values.
#'   \item **Interventions:** An optional mid-epidemic intervention that
#'     reduces the transmission probability by a user-specified efficacy,
#'     shown as a vertical dashed line on the plot.
#'   \item **Sensitivity analysis:** Vary any epidemic parameter across a
#'     range of values to compare trajectories side by side with a
#'     blue-to-red color ramp.
#'   \item **Vital dynamics:** Optional births and deaths to model endemic
#'     equilibria over longer time horizons.
#'   \item **Summary tab:** Displays the basic reproduction number
#'     (R\ifelse{html}{\out{<sub>0</sub>}}{0}), peak timing, cumulative
#'     infections, attack rate, and contextual interpretation.
#'   \item **Data tab:** Searchable table of model output with CSV download.
#' }
#'
#' @section ICM Application (`class = "icm"`):
#' The ICM app runs agent-based stochastic microsimulations for one-group
#' SI, SIR, and SIS models. Because ICMs are computationally more expensive,
#' a **Run Model** button controls when simulations execute. Features include:
#' \itemize{
#'   \item **Stochastic visualization:** Multiple simulations are summarized
#'     with mean trajectory lines and interquartile range (IQR) ribbons
#'     showing the middle 50\% of outcomes. Single simulations display
#'     individual trajectories.
#'   \item **Scenario presets:** The same Flu-like, STI-like, and Measles-like
#'     configurations as the DCM app, adapted for smaller populations suitable
#'     for agent-based simulation.
#'   \item **Interventions:** Same intervention system as DCM, natively
#'     supported by [param.icm].
#'   \item **Vital dynamics:** Optional births and deaths, as in DCM.
#'   \item **Summary tab:** Displays R\ifelse{html}{\out{<sub>0</sub>}}{0},
#'     epidemic timeline from the mean trajectory, and a stochastic variation
#'     section showing the range and standard deviation of peak infections and
#'     final prevalence across simulations.
#'   \item **Data tab:** View simulation output as means, standard deviations,
#'     or individual simulation values, with CSV download.
#' }
#'
#' @section Required Packages:
#' Both applications require the `shiny`, `bslib`, `DT`, and `plotly` packages.
#' Install any missing packages with
#' `install.packages(c("shiny", "bslib", "DT", "plotly"))`.
#'
#' @seealso
#' [dcm] and [param.dcm] for the DCM modeling API;
#' [icm] and [param.icm] for the ICM modeling API.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Deterministic compartmental models
#' epiweb(class = "dcm")
#'
#' ## Stochastic individual contact models
#' epiweb(class = "icm")
#' }
#'
epiweb <- function(class, ...) {
  if (class == "dcm") {
    if (!requireNamespace("shiny", quietly = TRUE)) {
      stop("Package \"shiny\" is required. Install with: ",
           "install.packages(\"shiny\")", call. = FALSE)
    }
    if (!requireNamespace("bslib", quietly = TRUE)) {
      stop("Package \"bslib\" is required for the DCM app. Install with: ",
           "install.packages(\"bslib\")", call. = FALSE)
    }
    shiny::runApp(system.file("shiny", "epidcm", package = "EpiModel"), ...)
  } else if (class == "icm") {
    if (!requireNamespace("shiny", quietly = TRUE)) {
      stop("Package \"shiny\" is required. Install with: ",
           "install.packages(\"shiny\")", call. = FALSE)
    }
    if (!requireNamespace("bslib", quietly = TRUE)) {
      stop("Package \"bslib\" is required for the ICM app. Install with: ",
           "install.packages(\"bslib\")", call. = FALSE)
    }
    shiny::runApp(system.file("shiny", "epiicm", package = "EpiModel"), ...)
  } else {
    stop("Specify class as either \"dcm\", \"icm\" ", call. = FALSE)
  }
}
