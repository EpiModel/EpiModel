
#' @title Stochastic Individual Contact Models
#'
#' @description Simulates stochastic individual contact epidemic models for
#'              infectious disease.
#'
#' @param param Model parameters, as an object of class \code{\link{param.icm}}.
#' @param init Initial conditions, as an object of class \code{\link{init.icm}}.
#' @param control Control settings, as an object of class
#'        \code{\link{control.icm}}.
#'
#' @details
#' Individual contact models are intended to be the stochastic microsimulation
#' analogs to deterministic compartmental models. ICMs simulate disease spread
#' on individual agents in discrete time as a function of processes with stochastic
#' variation. The stochasticity is inherent in all transition processes:
#' infection, recovery, and demographics. A detailed description of these models
#' may be found in the \href{http://statnet.github.io/tut/BasicICMs.html}{Basic
#' ICMs} tutorial.
#'
#' The \code{icm} function performs  modeling of both the base model types
#' and original models. Base model types include one-group and two-group
#' models with disease types for Susceptible-Infected (SI),
#' Susceptible-Infected-Recovered (SIR), and Susceptible-Infected-Susceptible (SIS).
#' Original models may be built by writing new process modules that either take
#' the place of existing modules (for example, disease recovery), or supplement
#' the set of existing processes with a new one contained in an original module.
#'
#' @return
#' A list of class \code{icm} with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{param}, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{control}, with additional controls added as necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model. Outputs for base models always include the
#'        size of each compartment, as well as flows in, out of, and between
#'        compartments.
#' }
#'
#' @keywords model
#'
#' @seealso Extract the model results with \code{\link{as.data.frame.icm}}.
#' Summarize the time-specific model results with \code{\link{summary.icm}}.
#' Plot the model results with \code{\link{plot.icm}}. Plot a compartment flow
#' diagram with \code{\link{comp_plot}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Example 1: SI Model
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
#' mod1 <- icm(param, init, control)
#' mod1
#' plot(mod1)
#'
#' ## Example 2: SIR Model
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
#' init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.icm(type = "SIR", nsteps = 500, nsims = 10)
#' mod2 <- icm(param, init, control)
#' mod2
#' plot(mod2)
#'
#' ## Example 3: SIS Model
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SIS", nsteps = 500, nsims = 10)
#' mod3 <- icm(param, init, control)
#' mod3
#' plot(mod3)
#'
#' ## Example 4: SI Model with Vital Dynamics (Two-Group)
#' param <- param.icm(inf.prob = 0.4,  inf.prob.g2 = 0.1,
#'                    act.rate = 0.25, balance = "g1",
#'                    a.rate = 1/100, a.rate.g2 = NA,
#'                    ds.rate = 1/100, ds.rate.g2 = 1/100,
#'                    di.rate = 1/50, di.rate.g2 = 1/50)
#' init <- init.icm(s.num = 500, i.num = 1,
#'                  s.num.g2 = 500, i.num.g2 = 0)
#' control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
#' mod4 <- icm(param, init, control)
#' mod4
#' plot(mod4)
#' }
#'
icm <- function(param, init, control) {

  crosscheck.icm(param, init, control)
  verbose.icm(control, type = "startup")

  # Simulation loop start
  for (s in 1:control$nsims) {

    ## Initialization module
    if (!is.null(control[["initialize.FUN"]])) {
      dat <- do.call(control[["initialize.FUN"]], list(param, init, control))
    }


    # Timestep loop
    for (at in 2:control$nsteps) {

      ## User Modules
      um <- control$user.mods
      if (length(um) > 0) {
        for (i in 1:length(um)) {
          dat <- do.call(control[[um[i]]], list(dat, at))
        }
      }

      ## Infection
      if (!is.null(control[["infection.FUN"]])) {
        dat <- do.call(control[["infection.FUN"]], list(dat, at))
      }


      ## Recovery
      if (!is.null(control[["recovery.FUN"]])) {
        dat <- do.call(control[["recovery.FUN"]], list(dat, at))
      }


      ## Departure Module
      if (!is.null(control[["departures.FUN"]])) {
        dat <- do.call(control[["departures.FUN"]], list(dat, at))
      }


      ## Arrival Module
      if (!is.null(control[["arrivals.FUN"]])) {
        dat <- do.call(control[["arrivals.FUN"]], list(dat, at))
      }


      ## Outputs
      if (!is.null(control[["get_prev.FUN"]])) {
        dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
      }


      ## Track progress
      verbose.icm(dat, type = "progress", s, at)
    }

    # Set output
    if (s == 1) {
      out <- saveout.icm(dat, s)
    } else {
      out <- saveout.icm(dat, s, out)
    }



  } # Simulation loop end


  class(out) <- "icm"
  return(out)
}

