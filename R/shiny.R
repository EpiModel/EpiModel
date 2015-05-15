
#' @title EpiModel Web
#'
#' @description Runs a web browser-based GUI of deterministic compartmental
#'              models and stochastic individual contact models.
#'
#' @param class Model class, with options of \code{"dcm"} and \code{"icm"}.
#'
#' @details
#' \code{epiweb} runs a web-based GUI of a one-group \code{\link{dcm}} models
#' and \code{\link{icm}} models with user input on model type, state sizes, and
#' parameters. Model output may be plotted, summarized, and saved as raw data
#' using the core \code{EpiModel} functionality for these model classes. These
#' applications are built using the \code{shiny} package framework.
#'
#' These apps are also hosted online at Rstudio's shinyapps site here:
#' \itemize{
#'    \item DCM App: \samp{https://statnet.shinyapps.io/epidcm/}
#'    \item ICM App: \samp{https://statnet.shinyapps.io/epiicm/}
#' }
#'
#' @references
#' RStudio. shiny: Web Application Framework for R. R package version 0.9.1. 2014.
#' \url{http://www.rstudio.com/shiny/}
#'
#' @seealso \code{\link{dcm}}, \code{\link{icm}}
#'
#' @keywords GUI
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
epiweb <- function(class) {
  if (class == "dcm") {
    shiny::runApp(system.file("shiny", "epidcm", package = "EpiModel"))
  } else if (class == "icm") {
    shiny::runApp(system.file("shiny", "epiicm", package = "EpiModel"))
  } else {
    stop("Specify class as either \"dcm\" or \"icm\" ")
  }
}
