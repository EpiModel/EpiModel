
#' @title Deterministic Compartmental Models
#'
#' @description Solves deterministic compartmental epidemic models for
#'              infectious disease.
#'
#' @param param Model parameters, as an object of class \code{\link{param.dcm}}.
#' @param init Initial conditions, as an object of class \code{\link{init.dcm}}.
#' @param control Control settings, as an object of class
#'        \code{\link{control.dcm}}.
#'
#' @details
#' The \code{dcm} function uses the ordinary differential equation solver in
#' the \code{deSolve} package to model disease as a deterministic compartmental
#' system. The parameterization for these models follows the standard approach
#' in \code{EpiModel}, with epidemic parameters, initial conditions, and control
#' settings.
#'
#' The \code{dcm} function performs  modeling of both base model types and
#' original models with new structures. Base model types include one-group
#' and two-group models with disease types for Susceptible-Infected (SI),
#' Susceptible-Infected-Recovered (SIR), and Susceptible-Infected-Susceptible
#' (SIS). Both base and original models require the \code{param},
#' \code{init}, and \code{control} inputs.
#'
#' @return
#' A list of class \code{dcm} with the following elements:
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
#' @references
#' Soetaert K, Petzoldt T, Setzer W. Solving Differential Equations in
#' R: Package deSolve. Journal of Statistical Software. 2010; 33(9): 1-25.
#' \doi{10.18637/jss.v033.i09}.
#'
#' @keywords model
#'
#' @seealso Extract the model results with \code{\link{as.data.frame.dcm}}.
#' Summarize the time-specific model results with \code{\link{summary.dcm}}.
#' Plot the model results with \code{\link{plot.dcm}}. Plot a compartment flow
#' diagram with \code{\link{comp_plot}}.
#'
#' @export
#'
#' @examples
#' ## Example 1: SI Model (One-Group)
#' # Set parameters
#' param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.dcm(s.num = 500, i.num = 1)
#' control <- control.dcm(type = "SI", nsteps = 500)
#' mod1 <- dcm(param, init, control)
#' mod1
#' plot(mod1)
#'
#' ## Example 2: SIR Model with Vital Dynamics (One-Group)
#' param <- param.dcm(inf.prob = 0.2, act.rate = 5,
#'                    rec.rate = 1/3, a.rate = 1/90, ds.rate = 1/100,
#'                    di.rate = 1/35, dr.rate = 1/100)
#' init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
#' control <- control.dcm(type = "SIR", nsteps = 500)
#' mod2 <- dcm(param, init, control)
#' mod2
#' plot(mod2)
#'
#' ## Example 3: SIS Model with act.rate Sensitivity Parameter
#' param <- param.dcm(inf.prob = 0.2, act.rate = seq(0.1, 0.5, 0.1),
#'                    rec.rate = 1/50)
#' init <- init.dcm(s.num = 500, i.num = 1)
#' control <- control.dcm(type = "SIS", nsteps = 500)
#' mod3 <- dcm(param, init, control)
#' mod3
#' plot(mod3)
#'
#' ## Example 4: SI Model with Vital Dynamics (Two-Group)
#' param <- param.dcm(inf.prob = 0.4,  inf.prob.g2 = 0.1,
#'                    act.rate = 0.25, balance = "g1",
#'                    a.rate = 1/100, a.rate.g2 = NA,
#'                    ds.rate = 1/100, ds.rate.g2 = 1/100,
#'                    di.rate = 1/50, di.rate.g2 = 1/50)
#' init <- init.dcm(s.num = 500, i.num = 1,
#'                  s.num.g2 = 500, i.num.g2 = 0)
#' control <- control.dcm(type = "SI", nsteps = 500)
#' mod4 <- dcm(param, init, control)
#' mod4
#' plot(mod4)
#'
dcm <- function(param, init, control) {
  check.control.class("dcm", "EpiModel dcm")

  crosscheck.dcm(param, init, control)

  # Model selection ---------------------------------------------------------
  if (is.null(control$new.mod)) {

    if (control$type == "SI") {
      if (param$groups == 1) {
        if (param$vital == FALSE) {
          model <- mod_SI_1g_cl
          init <- c(init, si.flow = 0)
        } else {
          model <- mod_SI_1g_op
          init <- c(init, si.flow = 0, a.flow = 0, ds.flow = 0, di.flow = 0)
        }
      } else {
        if (param$vital == FALSE) {
          model <- mod_SI_2g_cl
          init <- c(init, si.flow = 0, si.flow.g2 = 0)
        } else {
          model <- mod_SI_2g_op
          init <- c(init,
                    si.flow = 0, a.flow = 0, ds.flow = 0, di.flow = 0,
                    si.flow.g2 = 0, a.flow.g2 = 0, ds.flow.g2 = 0,
                    di.flow.g2 = 0)
        }
      }
    }
    if (control$type == "SIR") {
      if (param$groups == 1) {
        if (param$vital == FALSE) {
          model <- mod_SIR_1g_cl
          init <- c(init, si.flow = 0, ir.flow = 0)
        } else {
          model <- mod_SIR_1g_op
          init <- c(init, si.flow = 0, ir.flow = 0, a.flow = 0,
                    ds.flow = 0, di.flow = 0, dr.flow = 0)
        }
      } else {
        if (param$vital == FALSE) {
          model <- mod_SIR_2g_cl
          init <- c(init, si.flow = 0, ir.flow = 0,
                    si.flow.g2 = 0, ir.flow.g2 = 0)
        } else {
          model <- mod_SIR_2g_op
          init <- c(init, si.flow = 0, ir.flow = 0, a.flow = 0,
                    ds.flow = 0, di.flow = 0, dr.flow = 0,
                    si.flow.g2 = 0, ir.flow.g2 = 0, a.flow.g2 = 0,
                    ds.flow.g2 = 0, di.flow.g2 = 0, dr.flow.g2 = 0)
        }
      }
    }
    if (control$type == "SIS") {
      if (param$groups == 1) {
        if (param$vital == FALSE) {
          model <- mod_SIS_1g_cl
          init <- c(init, si.flow = 0, is.flow = 0)
        } else {
          model <- mod_SIS_1g_op
          init <- c(init, si.flow = 0, is.flow = 0,
                    a.flow = 0, ds.flow = 0, di.flow = 0)
        }
      } else {
        if (param$vital == FALSE) {
          model <- mod_SIS_2g_cl
          init <- c(init, si.flow = 0, is.flow = 0,
                    si.flow.g2 = 0, is.flow.g2 = 0)
        } else {
          model <- mod_SIS_2g_op
          init <- c(init, si.flow = 0, is.flow = 0,
                    a.flow = 0, ds.flow = 0, di.flow = 0,
                    si.flow.g2 = 0, is.flow.g2 = 0,
                    a.flow.g2 = 0, ds.flow.g2 = 0, di.flow.g2 = 0)
        }
      }
    }

  } else {
    # Sub in new model
    model <- control$new.mod
  }


  # Print model option ------------------------------------------------------
  if (control$print.mod == TRUE) {
    return(print(model))
  }


  # Initial conditions ------------------------------------------------------
  t0 <- setNames(as.numeric(init), names(init))


  # Sensitivity parameters --------------------------------------------------
  if (control$sens.param == FALSE) {
    control$nruns <- 1
  } else {
    control$nruns <- max(sapply(param, length))
  }
  if (control$nruns > 1) {
    longv <- which(sapply(param, length) == max(sapply(param, length)))
    longvn <- names(longv)
    lim.p <- param[!(names(param) %in% longvn)]
  }


  # Model runs --------------------------------------------------------------
  verbose.dcm(control, type = "startup")
  for (s in 1:control$nruns) {

    ## Sensitivity parameter input
    if (control$nruns > 1) {
      all.p <- as.list(lim.p)
      for (j in seq_along(longvn)) {
        sens.p <- param[[longvn[j]]][s]
        all.p[[longvn[j]]] <- sens.p
      }
    } else {
      all.p <- param
    }

    ## Timesteps
    if (length(control$nsteps) == 1) {
      times <- seq(1, control$nsteps, control$dt)
    } else {
      times <- control$nsteps
    }


    ## Differential equation solvers
    if (control$dede == FALSE) {
      df <- data.frame(ode(y = t0,
                           times = times,
                           func = model,
                           parms = all.p,
                           method = control$odemethod))
    } else {
      df <- data.frame(dede(y = t0,
                            times = times,
                            func = model,
                            parms = all.p))
    }

    ## Recalculate Flows
    isFlow <- grep(".flow", x = names(df))
    if (length(isFlow) > 0) {
      for (j in isFlow) {
        df[, j] <- c(diff(df[, j]), NA)
      }
    }


    # Output ------------------------------------------------------------------

    ## Save to out object
    if (s == 1) {
      out <- saveout.dcm(df, s, param, control)
    } else {
      out <- saveout.dcm(df, s, param, control, out)
    }

    ## Progress
    verbose.dcm(control, type = "progress", s)

  }

  class(out) <- "dcm"
  invisible(out)
}
