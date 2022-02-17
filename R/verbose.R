
#' @title Progress Print Module for Deterministic Compartmental Models
#'
#' @description This function prints progress from deterministic compartmental
#'              models simulated with \code{dcm} to the console.
#'
#' @param x If the \code{type} is "startup", then an object of class
#'        \code{control.icm}, otherwise the all main data object in \code{icm}
#'        simulations.
#' @param type Progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s Current simulation number, if type is "progress".
#'
#' @export
#' @keywords internal
#'
verbose.dcm <- function(x, type, s = 1) {

  if (type == "startup") {
    if (x$verbose == TRUE & x$nruns > 1) {
      cat("\nStarting DCM Simulation...")
    }
  }

  if (type == "progress") {
    if (x$verbose == TRUE & x$nruns > 1) {
      cat("\nRun = ", s, "/", x$nruns, sep = "")
    }
  }

}

#' @title Progress Print Module for Stochastic Individual Contact Models
#'
#' @description This function prints progress from stochastic individual contact
#'              models simulated with \code{icm} to the console.
#'
#' @param x If the \code{type} is "startup", then an object of class
#'        \code{control.icm}; otherwise \code{dat}, the main data object in
#'        \code{icm} simulations.
#' @param type Progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s Current simulation number, if type is "progress".
#' @param at Current time step, if type is "progress".
#'
#' @export
#' @keywords internal
#'
verbose.icm <- function(x, type, s = 1, at = 2) {

  if (type == "startup") {
    if (x$verbose == TRUE) {
      cat("\nStarting ICM Simulation...")
    }
  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep = "")
      }
      if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
        cat("\014")
        cat("\nEpidemic Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, sep = "")
        cat("\nTimestep: ", at, "/", x$control$nsteps, sep = "")
        status <- x$attr$status
        if (class(status) == "character") {
          status <- ifelse(status == "i", 1, 0)
        }
        cat("\nPrevalence:", sum(status, na.rm = TRUE))
        cat("\nPopulation Size:", sum(x$attr$active == 1))
        cat("\n----------------------------")
      }
    }
  }

}


#' @title Progress Print Module for Stochastic Network Models
#'
#' @description This function prints progress from stochastic network models
#'              simulated with \code{netsim} to the console.
#'
#' @param x If the \code{type} is "startup", then an object of class
#'        \code{control.net}; otherwise \code{dat}, the main data object in
#'        \code{netsim} simulations.
#' @param type Progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s Current simulation number, if type is "progress".
#' @param at Current time step, if type is "progress".
#'
#' @export
#' @keywords internal
#'
verbose.net <- function(x, type, s = 1, at = 2) {

  if (type == "startup" && x$ncores == 1) {
    if (x$verbose == TRUE) {
      cat("\nStarting Network Simulation...")
    }
  }

  if (type == "progress" && x$control$ncores == 1) {
    if (x$control$verbose == TRUE) {
      if (x$control$verbose.int == 0 && at == x$control$nsteps) {
        cat("\nSim = ", s, "/", x$control$nsims, sep = "")
      }
      if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
        cat("\014")
        cat("\nEpidemic Simulation")
        cat("\n----------------------------")
        cat("\nSimulation: ", s, "/", x$control$nsims, sep = "")
        cat("\nTimestep: ", at, "/", x$control$nsteps, sep = "")
        active <- x$attr$active
        status <- x$attr$status[which(active == 1)]
        if (class(status) == "character") {
          status <- ifelse(status == "i", 1, 0)
        }
        cat("\nPrevalence:", sum(status, na.rm = TRUE))
        cat("\nPopulation Size:", sum(active == 1))
        cat("\n----------------------------")
      }
    }
  }

}
