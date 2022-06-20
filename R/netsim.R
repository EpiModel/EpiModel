
#' @title Stochastic Network Models
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease.
#'
#' @param x Fitted network model object, as an object of class \code{netest}.
#'        Alternatively, if restarting a previous simulation, may be an object
#'        of class \code{netsim}.
#' @param param Model parameters, as an object of class \code{param.net}.
#' @param init Initial conditions, as an object of class \code{init.net}.
#' @param control Control settings, as an object of class
#'        \code{control.net}.
#'
#' @details
#' Stochastic network models explicitly represent phenomena within and across
#' edges (pairs of nodes that remain connected) over time. This enables edges to
#' have duration, allowing for repeated transmission-related acts within the
#' same dyad, specification of edge formation and dissolution rates, control
#' over the temporal sequencing of multiple edges, and specification of
#' network-level features. A detailed description of these models, along with
#' examples, is found in the \href{http://www.epimodel.org/tut.html}{Basic
#' Network Models} tutorials.
#'
#' The \code{netsim} function performs modeling of both the base model types
#' and original models. Base model types include one-group and two-group models
#' with disease types for Susceptible-Infected (SI),
#' Susceptible-Infected-Recovered (SIR), and
#' Susceptible-Infected-Susceptible (SIS).
#'
#' Original models may be parameterized by writing new process modules that
#' either take the place of existing modules (for example, disease recovery), or
#' supplement the set of existing processes with a new one contained in a new
#' module. This functionality is documented in the
#' \href{http://www.epimodel.org/tut.html}{Extension Network Models} tutorials.
#' The list of modules within \code{netsim} available for modification is listed
#' in \code{\link{modules.net}}.
#'
#' @return
#' A list of class \code{netsim} with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{param}, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through
#'        \code{control}, with additional controls added as necessary.
#'  \item \strong{epi:} a list of data frames, one for each epidemiological
#'        output from the model. Outputs for base models always include the
#'        size of each compartment, as well as flows in, out of, and between
#'        compartments.
#'  \item \strong{stats:} a list containing two sublists, \code{nwstats} for any
#'        network statistics saved in the simulation, and \code{transmat} for
#'        the transmission matrix saved in the simulation. See
#'        \code{\link{control.net}} and the
#'        \href{http://www.epimodel.org/tut.html}{tutorials} for further
#'        details.
#'  \item \strong{network:} a list of \code{networkDynamic} objects,
#'         one for each model simulation.
#' }
#' If \code{control$raw.output == TRUE}: A list of the raw (pre-processed)
#' \code{netsim} \code{dat} objects, for use in simulation continuation.
#'
#' @references
#' Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical
#' Modeling of Infectious Disease over Networks. Journal of Statistical
#' Software. 2018; 84(8): 1-47.
#'
#' @seealso Extract the model results with \code{\link{as.data.frame.netsim}}.
#'          Summarize the time-specific model results with
#'          \code{\link{summary.netsim}}. Plot the model results with
#'          \code{\link{plot.netsim}}.
#'
#' @keywords model
#' @export
#'
#' @examples
#' \dontrun{
#' ## Example 1: SI Model without Network Feedback
#' # Network model estimation
#' nw <- network_initialize(n = 100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Epidemic model
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
#' mod1 <- netsim(est1, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod1
#' plot(mod1)
#' summary(mod1, at = 50)
#'
#' ## Example 2: SIR Model with Network Feedback
#' # Recalculate dissolution coefficient with departure rate
#' coef.diss <- dissolution_coefs(
#'   dissolution = ~ offset(edges), duration = 20,
#'   d.rate = 0.0021
#' )
#'
#' # Reestimate the model with new coefficient
#' est2 <- netest(nw, formation, target.stats, coef.diss)
#'
#' # Reset parameters to include demographic rates
#' param <- param.net(
#'   inf.prob = 0.3, inf.prob.g2 = 0.15,
#'   rec.rate = 0.02, rec.rate.g2 = 0.02,
#'   a.rate = 0.002, a.rate.g2 = NA,
#'   ds.rate = 0.001, ds.rate.g2 = 0.001,
#'   di.rate = 0.001, di.rate.g2 = 0.001,
#'   dr.rate = 0.001, dr.rate.g2 = 0.001
#' )
#' init <- init.net(
#'   i.num = 10, i.num.g2 = 10,
#'   r.num = 0, r.num.g2 = 0
#' )
#' control <- control.net(
#'   type = "SIR", nsteps = 100, nsims = 5,
#'   resimulate.network = TRUE, tergmLite = TRUE
#' )
#'
#' # Simulate the model with new network fit
#' mod2 <- netsim(est2, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod2
#' plot(mod2)
#' summary(mod2, at = 40)
#' }
#'
netsim <- function(x, param, init, control) {
  crosscheck.net(x, param, init, control)
  if (!is.null(control[["verbose.FUN"]])) {
    do.call(control[["verbose.FUN"]], list(control, type = "startup"))
  }

  if (control$nsims == 1) {
    control$ncores <- 1
  } else {
    control$ncores <- min(parallel::detectCores(), control$ncores)
  }

  s <- NULL
  if (control$ncores == 1) {
    sout <- lapply(seq_len(control$nsims), function(s) {
      netsim_run(x, param, init, control, s)
    })
  } else if (control$ncores > 1) {
    doParallel::registerDoParallel(ncores)
    sout <- foreach(s = seq_len(control$nsims)) %dopar% {
      netsim_loop(x, param, init, control, s)
    }
  }

  if (control$raw.output) {
    out <- sout
  } else {
    out <- process_out.net(sout)
  }

  return(out)
}

#' @title Internal Function Running the Network Simulation Loop
#'
#' @description This function runs the initialization and simulation loop for
#'              one simulation. Errors, warnings, and messages are pretty
#'              printed using the \code{netsim_cond_msg} function (utils.R)
#' @inheritParams initialize.net
#'
#' @inherit recovery.net return
#'
#' @keywords internal
#'
netsim_loop <- function(x, param, init, control, s) {
  ## Instantiate random parameters
  param <- generate_random_params(param, verbose = FALSE)

  dat <- withCallingHandlers(
    expr = {
      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        current_mod <- "initialize.FUN"
        at <- paste0("`Initialization Step` (", control$start, ")")
        dat <- do.call(control[[current_mod]], list(x, param, init, control, s))
      }

      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {
          current_mod <- "epimodel.internal"
          dat <- set_current_timestep(dat, at)
          # Applies updaters, if any
          dat <- input_updater(dat)

          ## Module order
          morder <- get_control(dat, "module.order", override.null.error = TRUE)
          if (is.null(morder)) {
            bi.mods <- get_control(dat, "bi.mods")
            user.mods <- get_control(dat, "user.mods")
            lim.bi.mods <- bi.mods[
              -which(bi.mods %in% c("initialize.FUN", "verbose.FUN"))
            ]
            morder <- c(user.mods, lim.bi.mods)
          }

          ## Evaluate modules
          for (i in seq_along(morder)) {
            current_mod <- morder[[i]]
            mod.FUN <- get_control(dat, current_mod)
            dat <- do.call(mod.FUN, list(dat, at))
          }

          ## Verbose module
          verbose <- !is.null(control[["verbose.FUN"]])
          if (verbose == TRUE) {
            current_mod <- "verbose.FUN"
            do.call(control[[current_mod]], list(dat, type = "progress", s, at))
          }
        }
      }

      dat
    },
    message = function(e) message(netsim_cond_msg("MESSAGE", current_mod, at)),
    warning = function(e) message(netsim_cond_msg("WARNING", current_mod, at)),
    error = function(e) message(netsim_cond_msg("ERROR", current_mod, at))
  )

  return(dat)
}

#' @export
netsim_initialize <- function(x, param, init, control, s = 1) {
  param <- generate_random_params(param, verbose = FALSE)
  dat <- control[["initialize.FUN"]](x, param, init, control, s)
  dat <- set_current_timestep(dat, 1)

  return(dat)
}

#' @export
netsim_run_nsteps <- function(dat, nsteps, s) {
  last_timestep <- get_current_timestep(dat) + nsteps
  while (get_current_timestep(dat) < last_timestep) {
    dat <- increment_timestep(dat)
    dat <- netsim_run_modules(dat, s)
  }

  return(dat)
}

#' @export
netsim_run_modules <- function(dat, s) {
  at <- get_current_timestep(dat)

  dat <- withCallingHandlers(
    expr = {
      current_mod <- "epimodel.internal"
      # Applies updaters, if any
      dat <- input_updater(dat)

      ## Module order
      morder <- get_control(dat, "module.order", override.null.error = TRUE)
      if (is.null(morder)) {
        bi.mods <- get_control(dat, "bi.mods")
        user.mods <- get_control(dat, "user.mods")
        lim.bi.mods <- bi.mods[
          -which(bi.mods %in% c("initialize.FUN", "verbose.FUN"))
        ]
        morder <- c(user.mods, lim.bi.mods)
      }

      ## Evaluate modules
      for (i in seq_along(morder)) {
        current_mod <- morder[[i]]
        mod.FUN <- get_control(dat, current_mod)
        dat <- do.call(mod.FUN, list(dat, at))
      }

      ## Verbose module
      if (!is.null(get_control(dat, "verbose.FUN"))) {
        current_mod <- "verbose.FUN"
        verbose.FUN <- get_control(dat, current_mod)
        do.call(get_control(dat, "verbose.FUN"), list(dat, type = "progress", s, at))
      }

      dat
    },
    message = function(e) message(netsim_cond_msg("MESSAGE", current_mod, at)),
    warning = function(e) message(netsim_cond_msg("WARNING", current_mod, at)),
    error = function(e) message(netsim_cond_msg("ERROR", current_mod, at))
  )

  return(dat)
}

netsim_run <- function(x, param, init, control, s = 1) {
  dat <- netsim_initialize(x, param, init, control)
  dat <- netsim_run_nsteps(dat, get_control(dat, "nsteps") - 1, s)
  return(dat)
}
