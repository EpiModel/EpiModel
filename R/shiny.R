
#' @title EpiModel Web
#'
#' @description Runs a web browser-based GUI of deterministic compartmental
#'              models, stochastic individual contact models, and basic network
#'              models.
#'
#' @param class Model class, with options of `"dcm"`, `"icm"`, and `"net"`.
#' @param ... Additional arguments passed to `shiny::runApp`.
#'
#' @details
#' `epiweb` runs a web-based GUI of one-group deterministic compartmental
#' models, stochastic individual contact models, and stochastic network models
#' with user input on model type, state sizes, and parameters. Model output may
#' be plotted, summarized, and saved as raw data using the core `EpiModel`
#' functionality for these model classes. These applications are built using
#' the `shiny` package framework.
#'
#' @references
#' RStudio. shiny: Web Application Framework for R. R package version 1.0.5.
#' 2015. <https://shiny.posit.co/>.
#'
#' @seealso [`dcm`], [`icm`], [`netsim`]
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
#'
#' ## Stochastic network models
#' epiweb(class = "net")
#' }
#'
epiweb <- function(class, ...) {
  if (class == "dcm") {
    shiny::runApp(system.file("shiny", "epidcm", package = "EpiModel"), ...)
  } else if (class == "icm") {
    shiny::runApp(system.file("shiny", "epiicm", package = "EpiModel"), ...)
  } else if (class == "net") {
    shiny::runApp(system.file("shiny", "epinet", package = "EpiModel"), ...)
  } else {
    stop("Specify class as either \"dcm\", \"icm\" or \"net\" ", call. = FALSE)
  }
}
