
#' @title Stochastic Individual Contact Models
#'
#' @description Simulates stochastic individual contact epidemic models for
#'              infectious disease.
#'
#' @param param model parameters, as an object of class \code{\link{param.icm}}.
#' @param init initial conditions, as an object of class \code{\link{init.icm}}.
#' @param control control settings, as an object of class
#'        \code{\link{control.icm}}.
#'
#' @details
#' Individual contact models are intended to be the stochastic microsimulation
#' analogs to deterministic compartmental models. ICMs simulate disease spread
#' on individual agents in discrete time as a function of processes with stochastic
#' variation. The stochasticity is inherent in all transition processes:
#' infection, recovery, and demographics. A detailed description of these models
#' may be found in Section 3 of the
#' \href{http://statnet.org/EpiModel/vignette/Tutorial.pdf}{EpiModel Tutorial}.
#'
#' The \code{icm} function performs  modeling of both the built-in model types
#' and original models. Built-in model types include one-group and two-group
#' models with disease types for Susceptible-Infected (SI),
#' Susceptible-Infected-Recovered (SIR), and Susceptible-Infected-Susceptible (SIS).
#' Original models may be built by writing new process modules that either take
#' the place of existing modules (for example, disease recovery), or supplement
#' the set of existing processes with a new one contained in an original module.
#' New and replacement modules may be written and input into \code{icm} following
#' the steps outlined in the
#' \href{http://statnet.org/EpiModel/vignette/NewICMs.html}{Solving New ICMs
#' with EpiModel} tutorial.
#'
#' @return
#' A list of class \code{icm} with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{param}, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through,
#'        \code{control}, with additional controls added as necessary.
#'  \item \strong{epi:} a list of of data frames, one for each epidemiological
#'        output from the model. Outputs for built-in models always include the
#'        size of each compartment, as well as flows in, out, and between
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
#'                    b.rate = 1/100, b.rate.g2 = NA,
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
    all <- do.call(control[["initialize.FUN"]], list(param, init, control))


    # Timestep loop
    for (at in 2:control$nsteps) {

      ## User Modules
      bim <- grep(".FUN", names(formals(control.icm)), value = TRUE)
      um <- which(grepl(".FUN", names(control)) &
                    !(names(control) %in% bim))
      if (length(um) > 0) {
        for (i in 1:length(um)) {
          umn <- names(control)[um[i]]
          all <- do.call(control[[umn]], list(all, at))
        }
      }

      ## Infection
      all <- do.call(control[["infection.FUN"]], list(all, at))


      ## Recovery
      all <- do.call(control[["recovery.FUN"]], list(all, at))


      ## Mortality
      all <-  do.call(control[["deaths.FUN"]], list(all, at))


      ## Birth Module
      all <-  do.call(control[["births.FUN"]], list(all, at))


      ## Outputs
      all <-  do.call(control[["get_prev.FUN"]], list(all, at))


      ## Track progress
      verbose.icm(all, type = "progress", s, at)
    }

    # Set output
    if (s == 1) {
      out <- saveout.icm(all, s)
    } else {
      out <- saveout.icm(all, s, out)
    }



  } # Simulation loop end


  class(out) <- "icm"
  return(out)
}

