
#' @title Stochastic Network Models
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease.
#'
#' @param x Fitted network model object, as an object of class \code{netest}.
#'        Alternatively, if restarting a previous simulation, may be an object of
#'        class \code{netsim}.
#' @param param Model parameters, as an object of class \code{param.net}.
#' @param init Initial conditions, as an object of class \code{init.net}.
#' @param control Control settings, as an object of class
#'        \code{control.net}.
#'
#' @details
#' Stochastic network models explicitly represent phenomena within and across edges
#' (pairs of nodes that remain connected) over time. This enables edges to have duration,
#' allowing for repeated transmission-related acts within the same dyad, specification of
#' edge formation and dissolution rates, control over the temporal sequencing of
#' multiple edges, and specification of network-level features. A detailed
#' description of these models, along with examples, is found in the
#' \href{http://statnet.github.io/tut/BasicNet.html}{Basic Network Models}
#' tutorial.
#'
#' The \code{netsim} function performs modeling of both the base model types
#' and original models. Base model types include one-mode and bipartite models
#' with disease types for Susceptible-Infected (SI), Susceptible-Infected-Recovered
#' (SIR), and Susceptible-Infected-Susceptible (SIS).
#'
#' Original models may be parameterized by writing new process modules that
#' either take the place of existing modules (for example, disease recovery), or
#' supplement the set of existing processes with a new one contained in a new
#' module. This functionality is documented in the
#' \href{http://statnet.github.io/tut/NewNet.html}{Solving New Network Models}
#' tutorial. The list of modules within \code{netsim} available for modification
#' is listed in \code{\link{modules.net}}.
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
#'  \item \strong{network:} a list of \code{networkDynamic} objects (or
#'        \code{network} objects if \code{delete.nodes} was set to \code{TRUE}),
#'        one for each model simulation.
#' }
#'
#' @references
#' Goodreau SM, Carnegie NB, et al. What drives the US and Peruvian HIV epidemics
#' in men who have sex with men (MSM)? PloS One. 2012; 7(11): e50522.
#'
#' @seealso Extract the model results with \code{\link{as.data.frame.netsim}}.
#'          Summarize the time-specific model results with \code{\link{summary.netsim}}.
#'          Plot the model results with \code{\link{plot.netsim}}.
#'
#' @keywords model
#' @export
#'
#' @examples
#' \dontrun{
#' ## Example 1: Independent SI Model
#' # Network model estimation
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Epidemic model
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
#' mod1 <- netsim(est1, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod1
#' plot(mod1)
#' summary(mod1, at = 50)
#'
#' ## Example 2: Dependent SIR Model
#' # Recalculate dissolution coefficient with death rate
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
#'                                d.rate = 0.0021)
#'
#' # Reestimate the model with new coefficient
#' est2 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Reset parameters to include demographic rates
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15,
#'                    rec.rate = 0.02, rec.rate.m2 = 0.02,
#'                    b.rate = 0.002, b.rate.m2 = NA,
#'                    ds.rate = 0.001, ds.rate.m2 = 0.001,
#'                    di.rate = 0.001, di.rate.m2 = 0.001,
#'                    dr.rate = 0.001, dr.rate.m2 = 0.001)
#' init <- init.net(i.num = 10, i.num.m2 = 10,
#'                  r.num = 0, r.num.m2 = 0)
#' control <- control.net(type = "SIR", nsteps = 100, nsims = 5)
#'
#' # Simulate the model with new network fit
#' mod2 <- netsim(est2, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod2
#' plot(mod2)
#' summary(mod2, at = 100)
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

  if (ncores == 1) {
    for (s in 1:control$nsims) {

      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(x, param, init, control, s))
      }


      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {

          ## Module order
          morder <- control$module.order
          if (is.null(morder)) {
            lim.bi.mods <- control$bi.mods[-which(control$bi.mods %in%
                                                    c("initialize.FUN", "verbose.FUN"))]
            morder <- c(control$user.mods, lim.bi.mods)
          }

          ## Evaluate modules
          for (i in seq_along(morder)) {
            dat <- do.call(control[[morder[i]]], list(dat, at))
          }

          ## Verbose module
          if (!is.null(control[["verbose.FUN"]])) {
            do.call(control[["verbose.FUN"]], list(dat, type = "progress", s, at))
          }

        }
      }

      # Set output
      if (s == 1) {
        out <- saveout.net(dat, s)
      } else {
        out <- saveout.net(dat, s, out)
      }
      class(out) <- "netsim"
    }
  }

  if (ncores > 1) {
    doParallel::registerDoParallel(ncores)

    sout <- foreach(s = 1:nsims) %dopar% {

      control$nsims <- 1
      control$currsim <- s

      ## Initialization Module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(x, param, init, control, s))
      }


      ### TIME LOOP
      if (control$nsteps > 1) {
        for (at in max(2, control$start):control$nsteps) {

          ## Module order
          morder <- control$module.order
          if (is.null(morder)) {
            lim.bi.mods <- control$bi.mods[-which(control$bi.mods %in%
                                                    c("initialize.FUN", "verbose.FUN"))]
            morder <- c(control$user.mods, lim.bi.mods)
          }

          ## Evaluate modules
          for (i in seq_along(morder)) {
            dat <- do.call(control[[morder[i]]], list(dat, at))
          }

          ## Verbose module
          if (!is.null(control[["verbose.FUN"]])) {
            do.call(control[["verbose.FUN"]], list(dat, type = "progress", s, at))
          }

        }
      }

      # Set output
      out <- saveout.net(dat, s = 1)
      class(out) <- "netsim"
      return(out)
    }

    merged.out <- sout[[1]]
    for (i in 2:length(sout)) {
      merged.out <- merge(merged.out, sout[[i]], param.error = FALSE)
    }
    out <- merged.out
    class(out) <- "netsim"
  }

  return(out)
}

