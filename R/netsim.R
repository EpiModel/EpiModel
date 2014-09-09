
#' @title Stochastic Network Models
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease.
#'
#' @param x fitted network model object, as an object of class
#'        \code{\link{netest}}.
#' @param param model parameters, as an object of class \code{\link{param.net}}.
#' @param init initial conditions, as an object of class \code{\link{init.net}}.
#' @param control control settings, as an object of class
#'        \code{\link{control.net}}.
#'
#' @details
#' Stochastic network models move beyond stochastic individual contact models by
#' explicitly modeling phenomena within and across edges (pairs of nodes that
#' remain connected) over time. This enables edges to have duration, allowing for
#' repeated transmission-related acts within the same dyad, specification of edge
#' formation and dissolution rates, control over the temporal sequencing of
#' multiple edges, and specification of network-level features. A detailed
#' description of these models, along with examples, is found in Section 4 of
#' the \href{http://statnet.org/EpiModel/vignette/Tutorial.pdf}{EpiModel Tutorial}.
#'
#' The \code{netsim} function performs modeling of both the built-in model types
#' and original models. Built-in model types include one-mode and bipartite models
#' with disease types for Susceptible-Infected (SI), Susceptible-Infected-Recovered
#' (SIR), and Susceptible-Infected-Susceptible (SIS).
#'
#' Original models may be parameterized by writing new process modules that
#' either take the place of existing modules (for example, disease recovery), or
#' supplement the set of existing processes with a new one contained in an new
#' module. This functionality is new to \code{EpiModel} and further documentation
#' will be posted through tutorial vignettes at the
#' \href{http://statnet.org/trac/wiki/EpiModel}{EpiModel website} shortly. This
#' modular approach is the same as building new ICMs, for which there is an
#' existing \href{http://statnet.org/EpiModel/vignette/NewICMs.html}{Solving New
#' ICMs with EpiModel} tutorial. Finally, the list of modules within \code{netsim}
#' available for modification is listed in \code{\link{modules.net}}.
#'
#' @return
#' A list of class \code{netsim} with the following elements:
#' \itemize{
#'  \item \strong{param:} the epidemic parameters passed into the model through
#'        \code{param}, with additional parameters added as necessary.
#'  \item \strong{control:} the control settings passed into the model through,
#'        \code{control}, with additional controls added as necessary.
#'  \item \strong{epi:} a list of of data frames, one for each epidemiological
#'        output from the model. Outputs for built-in models always include the
#'        size of each compartment, as well as flows in, out, and between
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
#' # Initialize network and set network model parameters
#' nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' duration <- 20
#' coef.diss <- dissolution_coefs(dissolution, duration)
#'
#' # Estimate the ERGM models (see help for netest)
#' # Skipping model diagnostics for this, but one should always run these
#' est1 <- netest(nw,
#'                formation,
#'                dissolution,
#'                target.stats,
#'                coef.diss,
#'                verbose = FALSE)
#'
#' # Parameters, initial conditions, and controls for model
#' param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
#' init <- init.net(i.num = 10, i.num.m2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        verbose.int = 0)
#'
#' # Run the model simulation
#' mod1 <- netsim(est1, param, init, control)
#'
#' # Print, plot, and summarize the results
#' mod1
#' plot(mod1)
#' summary(mod1, at = 50)
#'
#' ## Example 2: Dependent SIR Model
#' # Recalculate dissolution coefficient with death rate
#' coef.diss <- dissolution_coefs(dissolution, duration, d.rate = 0.0021)
#'
#' # Reestimate the model with new coefficient
#' est2 <- netest(nw,
#'                formation,
#'                dissolution,
#'                target.stats,
#'                coef.diss,
#'                verbose = FALSE)
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
netsim <- function(x,
                   param,
                   init,
                   control
                   ) {

  crosscheck.net(x, param, init, control)
  do.call(control[["verbose.FUN"]], list(control, type = "startup"))


  ### SIMULATION LOOP
  for (s in 1:control$nsims) {

    ## Initialization Module
    all <- do.call(control[["initialize.FUN"]], list(x, param, init, control))


    ### TIME LOOP
    for (at in 2:control$nsteps) {

      ## User Modules
      um <- control$user.mods
      if (length(um) > 0) {
        for (i in seq_along(um)) {
          all <- do.call(control[[um[i]]], list(all, at))
        }
      }

      ## Demographics Modules
      all <- do.call(control[["deaths.FUN"]], list(all, at))
      all <- do.call(control[["births.FUN"]], list(all, at))


      ## Recovery Module
      all <- do.call(control[["recovery.FUN"]], list(all, at))


      ## Resimulate network
      all <- do.call(control[["edges_correct.FUN"]], list(all, at))
      all <- do.call(control[["resim_nets.FUN"]], list(all, at))


      ## Infection Module
      all <- do.call(control[["infection.FUN"]], list(all, at))


      ## Save Prevalence
      all <- do.call(control[["get_prev.FUN"]], list(all, at))


      ## Progress Console
      do.call(control[["verbose.FUN"]], list(all, type = "progress", s, at))

    }

    # Set output
    if (s == 1) {
      out <- saveout.net(all, s)
    } else {
      out <- saveout.net(all, s, out)
    }

  } # end sim loop


  class(out) <- "netsim"
  invisible(out)

}


#' @title Stochastic Network Models in Parallel
#'
#' @description Simulates stochastic network epidemic models for infectious
#'              disease in parallel.
#'
#' @inheritParams netsim
#'
#' @details
#' This is an experimental implementation of the \code{\link{netsim}} function
#' that runs model simulations in parallel, using the \code{foreach} and
#' \code{doParallel} libraries.
#'
#' To run models in parallel, add an argument to the control settings called
#' \code{ncores} that is equal to the number of parallel cores the simulations
#' should be initiated on. Use \code{\link{detectCores}} to find the maximum on
#' a system.
#'
#' This has been tested on Linux, Mac, and Windows but no guarantees are made
#' that it will work on every platform. It is best-suited to be run in batch
#' mode. Memory management errors have been encounted when running large simulations
#' (large networks, long time steps, saving networkDynamic objects) in
#' interactive environments like Rstudio server.
#'
#' Note that this function may be folded into \code{\link{netsim}} and deprecated
#' in the future.
#'
#' @keywords model
#' @export
#'
#' @examples
#' \dontrun{
#' nw <- network.initialize(n = 1000, directed = FALSE)
#' formation <- ~ edges
#' target.stats <- 500
#' dissolution <- ~ offset(edges)
#' duration <- 50
#' coef.diss <- dissolution_coefs(dissolution, duration)
#'
#' est <- netest(nw,
#'               formation,
#'               dissolution,
#'               target.stats,
#'               coef.diss,
#'               verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.25)
#' init <- init.net(i.num = 50)
#' control <- control.net(type = "SI", nsteps = 250,
#'                        nsims = 4, ncores = 4)
#'
#' sims <- netsim_parallel(est, param, init, control)
#' plot(sims)
#' }
#'
netsim_parallel <- function(x,
                            param,
                            init,
                            control
                            ) {

  nsims <- control$nsims
  ncores <- control$ncores
  partype <- control$partype
  if (is.null(partype)) {
    partype <- "parallel"
  }
  if (partype == "snow") {
    ncores <- sum(cores.per.node)
  }

  if (nsims > 1 && ncores > 1) {
    suppressPackageStartupMessages(require(foreach))

    if (partype == "parallel") {
      suppressPackageStartupMessages(require(doParallel))
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)
    }
    if (partype == "snow") {
      suppressPackageStartupMessages(require(doSNOW))
      nodes <- control$nodes
      cores.per.node <- control$cores.per.node
      cl <- makeSOCKcluster(rep(nodes, cores.per.node))
      registerDoSNOW(cl)
    }

    out <- foreach(i = 1:nsims) %dopar% {

      require(EpiModel)
      control$verbose = FALSE
      control$nsims = 1

      netsim(x, param, init, control)

    }

    all <- out[[1]]
    for (i in 2:length(out)) {
      all <- merge(all, out[[i]])
    }

    if (partype == "snow") {
      stopCluster(cl)
    }

  } else {
    all <- netsim(x, param, init, control)
  }

  return(all)
}
