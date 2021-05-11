
#' @title Stochastic Network Models
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease.
#'
#' @param x Fitted network model object, as an object of class \code{netest}.
#'        Alternatively, if restarting a previous simulation, may be an object
#'        ofclass \code{netsim}.
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
#' and original models. Base model types include one-mode and two-group models
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
#'        \code{\link{control.net}} and the Tutorial for further details.
#'  \item \strong{network:} a list of \code{networkDynamic} objects,
#'         one for each model simulation.
#' }
#' If \code{control$raw.output == TRUE}: A list of the raw (pre-processed)
#' nestsim dat objects, for use in simulation continuation.
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
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
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
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
#'                                d.rate = 0.0021)
#'
#' # Reestimate the model with new coefficient
#' est2 <- netest(nw, formation, target.stats, coef.diss)
#'
#' # Reset parameters to include demographic rates
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15,
#'                    rec.rate = 0.02, rec.rate.g2 = 0.02,
#'                    a.rate = 0.002, a.rate.g2 = NA,
#'                    ds.rate = 0.001, ds.rate.g2 = 0.001,
#'                    di.rate = 0.001, di.rate.g2 = 0.001,
#'                    dr.rate = 0.001, dr.rate.g2 = 0.001)
#' init <- init.net(i.num = 10, i.num.g2 = 10,
#'                  r.num = 0, r.num.g2 = 0)
#' control <- control.net(type = "SIR", nsteps = 100, nsims = 5,
#'                        resimulate.network = TRUE, tergmLite = TRUE)
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

  nsims <- control$nsims
  ncores <- ifelse(nsims == 1, 1, min(parallel::detectCores(), control$ncores))
  control$ncores <- ncores

  if (is.null(control$resimulate.network)) {
    control$resimulate.network <- FALSE
  }

  s <- NULL
  if (ncores == 1) {
    sout <- lapply(seq_len(control$nsims), function(s) {
      # Run the simulation
      netsim_loop(x, param, init, control, s)
    })
  }

  if (ncores > 1) {
    cluster.type <- if (!is.null(control$cluster.type)) control$cluster.type
                    else "PSOCK"
    cl <- parallel::makeCluster(ncores, type = cluster.type)
    sout <- parallel::parLapply(cl, seq_len(control$nsims), function(s) {
      # Run the simulation
      netsim_loop(x, param, init, control, s)
    })
    parallel::stopCluster(cl)
  }

  # Process the outputs unless `control$raw.output` is `TRUE`
  if (!is.null(control$raw.output) && control$raw.output == TRUE) {
    out <- sout
  } else {
    out <- process_out.net(sout)
  }

  return(out)
}

#' @title Internal function running the network simulation loop
#'
#' @description This function run the initialization and simulation loop for one
#'              simulation. Error, warning and messages are pretty printed using
#'              the `netsim_cond_msg` (utils.R)
#' @inheritParams initialize.net
#' @keywords internal
netsim_loop <- function(x, param, init, control, s) {
  dat <- withCallingHandlers(
    expr = {

      ## Instantiate random parameters
      param <- generate_random_params(param, verbose = FALSE)

      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        current_mod <- "initialize.FUN"
        at <- paste0("`Initialization Step` (", control$start, ")")
        dat <- do.call(control[[current_mod]], list(x, param, init, control, s))
      }

      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {

          ## Module order
          morder <- control$module.order
          if (is.null(morder)) {
            lim.bi.mods <- control$bi.mods[-which(
              control$bi.mods %in% c("initialize.FUN", "verbose.FUN")
            )]
            morder <- c(control$user.mods, lim.bi.mods)
          }

          ## Evaluate modules
          for (i in seq_along(morder)) {
            current_mod <- morder[[i]]
            dat <- do.call(control[[current_mod]], list(dat, at))
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
    message = function(e) message(netsim_cond_msg("MESSAGE",
                                                  current_mod, at, e)),
    warning = function(e) message(netsim_cond_msg("WARNING",
                                                  current_mod, at, e)),
    error = function(e) message(netsim_cond_msg("ERROR", current_mod, at, e))
  )

  return(dat)
}
