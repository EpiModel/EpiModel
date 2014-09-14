
#' @title Dynamic Network Model Diagnostics
#'
#' @description Runs dynamic diagnostics on an ERGM/STERGM estimated through
#'              \code{netest}
#'
#' @param x an \code{EpiModel} object of class \code{netest}.
#' @param nsims number of simulations to run.
#' @param nsteps number of time steps per simulation.
#' @param nwstats.formula a right-hand sided ERGM formula with the network
#'        statistics of interest. The default is the formation formula of the
#'        network model contained in \code{x}.
#' @param set.control.ergm control arguments passed to simulate.ergm (see
#'        details).
#' @param set.control.stergm control arguments passed to simulate.stergm (see
#'        details).
#' @param verbose print progress to the console.
#'
#' @details
#' The \code{netdx} function handles dynamic network diagnostics for network
#' models fit with the \code{netest} function. Given the fitted model, \code{netdx}
#' simulates a specified number of dynamic networks for a specified number of
#' time steps per simulation. The network statistics in \code{nwstats.formula}
#' are saved for each time step. Summary statistics for the formation model terms,
#' as well as dissolution model statistics, are then calculated for access when
#' printing or plotting the \code{netdx} object.
#'
#' @section Control Arguments:
#' Models fit with the full STERGM method in \code{netest} (setting \code{edapprox}
#' argument to \code{FALSE}) require only a call to \code{simulate.stergm}.
#' Control parameters for those simulations may be set using
#' \code{set.control.stergm} in \code{netdx}. The parameters should be input
#' through the \code{control.simulate.stergm()} function, with the available
#' parameters listed in the \code{\link{control.simulate.stergm}} help
#' page in the \code{tergm} package.
#'
#' Models fit with the ERGM method with the edges dissolution approximation
#' (setting \code{edapprox} to \code{TRUE}) require a call first to
#' \code{simulate.ergm} for simulating an initial network and second to
#' \code{simulate.network} for simulating that static network forward through
#' time. Control parameters may be set for both processes in \code{netdx}.
#' For the first, the parameters should be input through the
#' \code{control.simulate.ergm()} function, with the available parameters listed
#' in the \code{\link[ergm:control.simulate.ergm]{control.simulate.ergm}} help
#' page in the \code{ergm} package. For the second, parameters should be input
#' through the \code{control.simulate.network()} function, with the available
#' parameters listed in the \code{\link{control.simulate.network}} help page in
#' the \code{tergm} package. An example is shown below.
#'
#' @seealso Plot these model diagnostics with \code{\link{plot.netdx}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Network initialization and model parameterization
#' nw <- network.initialize(100, directed = FALSE)
#' formation <- ~ edges
#' dissolution <- ~ offset(edges)
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution, duration = 10)
#'
#' # Estimate the model
#' est <- netest(
#'   nw,
#'   formation,
#'   dissolution,
#'   target.stats,
#'   coef.diss,
#'   verbose = FALSE)
#'
#' # Run diagnostics
#' dx <- netdx(est,
#'   nsims = 5,
#'   nsteps = 500,
#'   nwstats.formula = ~ edges + meandeg + concurrent,
#'   set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6),
#'   set.control.stergm = control.simulate.network(MCMC.burnin.min = 1e5))
#' dx
#' plot(dx)
#' }
#'
netdx <- function(x,
                  nsims = 1,
                  nsteps,
                  nwstats.formula,
                  set.control.ergm,
                  set.control.stergm,
                  verbose = TRUE) {

  if (class(x$fit) == "network") {
    nw <- x$fit
  } else {
    nw <- x$fit$network
    fit <- x$fit
  }
  formation <- x$formation
  coef.form <- x$coef.form
  dissolution <- x$dissolution
  coef.diss <- x$coef.diss
  constraints <- x$constraints
  if (is.null(constraints)) {
    constraints <- ~ .
  }
  target.stats <- x$target.stats
  edapprox <- x$edapprox

  if (verbose == TRUE) {
    cat("======================")
    cat("\nRunning Diagnostics")
    cat("\n======================\n")
  }

  if (missing(nwstats.formula)) {
    nwstats.formula <- formation
  }

  if (verbose == TRUE) {
    cat("- Simulating", nsims, "networks")
  }

  if (edapprox == FALSE) {

    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.stergm()
    }

    diag.sim <- list()
    if (verbose == TRUE) {
      cat("\n  |")
    }
    for (i in 1:nsims) {
      diag.sim[[i]] <- simulate(fit,
                           time.slices = nsteps,
                           monitor = nwstats.formula,
                           nsim = 1,
                           control = set.control.stergm)
      if (verbose == TRUE) {
        cat("*")
      }
    }
    if (verbose == TRUE) {
      cat("|")
    }
  }

  if (edapprox == TRUE) {

    if (missing(set.control.ergm)) {
      set.control.ergm <- control.simulate.ergm()
    }
    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.network()
    }

    diag.sim <- list()
    if (verbose == TRUE) {
      cat("\n  |")
    }
    for (i in 1:nsims) {
      if (class(x$fit) == "network") {
        fit.sim <- simulate(formation,
                            basis = nw,
                            coef = x$coef.form.crude,
                            constraints = constraints)
      } else {
        fit.sim <- simulate(fit, control = set.control.ergm)
      }
      diag.sim[[i]] <- simulate(fit.sim,
                           formation = formation,
                           dissolution = dissolution,
                           coef.form = coef.form,
                           coef.diss = coef.diss$coef.crude,
                           time.slices = nsteps,
                           constraints = constraints,
                           monitor = nwstats.formula,
                           nsim = 1,
                           control = set.control.stergm)
      if (verbose == TRUE) {
        cat("*")
      }
    }
    if (verbose == TRUE) {
      cat("|")
    }
  }

  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  stats <- list()
  for (i in 1:length(diag.sim)) {
    stats[[i]] <- as.matrix(attributes(diag.sim[[i]])$stats)[1:nsteps, , drop = FALSE]
  }

  ## Merged stats across all simulations
  if (nsims > 1) {
    merged.stats <- matrix(NA, nrow = nrow(stats[[1]])*nsims,
                           ncol = ncol(stats[[1]]))
    for (i in 1:ncol(stats[[1]])) {
      merged.stats[,i] <- as.numeric(sapply(stats, function(x) c(x[,i])))
    }
    colnames(merged.stats) <- colnames(stats[[1]])
  } else {
    merged.stats <- stats[[1]]
  }

  ## Calculate mean/sd from merged stats
  stats.means <- colMeans(merged.stats)
  stats.sd <- apply(merged.stats, 2, sd)
  stats.table <- data.frame(sorder = 1:length(names(stats.means)),
                            names = names(stats.means),
                            stats.means, stats.sd)

  # Which formation terms are offsets?
  is.offset.term <- grep(pattern = "offset",
                         strsplit(as.character(formation), "[+]")[[2]])


  ## Get stats from for target statistics, removing offsets
  ts.attr.names <- names(coef.form)
  if (length(is.offset.term > 0)) {
    ts.attr.names <- ts.attr.names[-is.offset.term]
  }
  ts.out <- data.frame(names = ts.attr.names,
                       targets = target.stats)


  ## Create stats.formation table for output
  stats.table <- merge(ts.out, stats.table, all = TRUE)
  stats.table <- stats.table[do.call("order",
                                     stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
  rownames(stats.table) <- stats.table$names
  stats.table.formation <- stats.table[, -c(1, 3)]


  if (verbose == TRUE) {
    cat("\n- Calculating duration statistics")
  }


  ## Duration calculations
  sim.df <- list()
  for (i in 1:length(diag.sim)) {
    sim.df[[i]] <- as.data.frame(diag.sim[[i]])
  }


  ## Create a merged vector of durations
  ncens <- which(sim.df[[1]]$onset.censored == FALSE &
                 sim.df[[1]]$terminus.censored == FALSE)
  durVec <- sim.df[[1]]$duration[ncens]
  if (nsims > 1) {
    for (i in 1:length(diag.sim)) {
      ncens <- which(sim.df[[i]]$onset.censored == FALSE &
                     sim.df[[i]]$terminus.censored == FALSE)
      durVec <- c(durVec, sim.df[[i]]$duration[ncens])
    }
  }

  # Create duration table for "dissolution = ~ offset(edges)"
  if (dissolution == ~offset(edges)) {
    duration.mean <- mean(durVec)
    duration.sd <- sd(durVec)
    duration.expected <- exp(coef.diss$coef.crude[1]) + 1
    stats.table.duration <- c(target = duration.expected,
                              sim.mean = duration.mean,
                              sim.sd = duration.sd)
  } else {
    stop('Only ~offset(edges) dissolution models currently supported')
  }

  # Calculate mean partnership age from edgelist
  pages <- list()
  if (verbose == TRUE) {
    cat("\n  |")
  }
  for (i in 1:length(diag.sim)) {
    pages[[i]] <- edgelist_meanage(el = sim.df[[i]])
    if (verbose == TRUE) {
      cat("*")
    }
  }
  if (verbose == TRUE) {
    cat("|")
  }


  ## Save output
  out <- list()
  out$nw <- nw
  out$formation <- formation
  out$coef.form <- coef.form
  out$dissolution <- dissolution
  out$coef.diss <- coef.diss
  out$constraints <- constraints
  out$edapprox <- edapprox
  out$target.stats <- ts.out
  out$nsims <- nsims
  out$nsteps <- nsteps

  out$stats <- stats
  out$stats.table.formation <- stats.table.formation
  out$stats.table.duration <- stats.table.duration
  out$edgelist <- sim.df
  out$pages <- pages

  class(out) <- "netdx"
  return(out)
}
