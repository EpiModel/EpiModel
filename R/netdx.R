
#' @title Dynamic Network Model Diagnostics
#'
#' @description Runs dynamic diagnostics on an ERGM/STERGM estimated through
#'              \code{netest}
#'
#' @param x An \code{EpiModel} object of class \code{netest}.
#' @param nsims Number of simulations to run.
#' @param dynamic If \code{TRUE}, runs dynamic diagnostics. If \code{FALSE} and
#'        the \code{netest} object was fit with the Edges Dissolution approximation
#'        method, simulates from the static ERGM fit.
#' @param nsteps Number of time steps per simulation (dynamic simulations only).
#' @param nwstats.formula A right-hand sided ERGM formula with the network
#'        statistics of interest. The default is the formation formula of the
#'        network model contained in \code{x}.
#' @param set.control.ergm Control arguments passed to simulate.ergm (see
#'        details).
#' @param set.control.stergm Control arguments passed to simulate.stergm (see
#'        details).
#' @param keep.tedgelist If \code{TRUE}, keep the timed edgelist generated from
#'        the dynamic simulations, for further analysis on edge durations.
#' @param verbose Print progress to the console.
#' @param ncores Number of processor cores to run multiple simulations
#'        on, using the \code{foreach} and \code{doParallel} implementations.
#'
#' @details
#' The \code{netdx} function handles dynamic network diagnostics for network
#' models fit with the \code{netest} function. Given the fitted model, \code{netdx}
#' simulates a specified number of dynamic networks for a specified number of
#' time steps per simulation. The network statistics in \code{nwstats.formula}
#' are saved for each time step. Summary statistics for the formation model terms,
#' as well as dissolution model and relational duration statistics, are then
#' calculated for access when printing or plotting the \code{netdx} object.
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
#' coef.diss <- dissolution_coefs(dissolution, duration = 25)
#'
#' # Estimate the model
#' est <- netest(nw, formation, dissolution,
#'               target.stats, coef.diss, verbose = FALSE)
#'
#' # Static diagnostics on the ERGM fit
#' dx1 <- netdx(est, nsims = 1e4, dynamic = FALSE,
#'              nwstats.formula = ~ edges + meandeg + concurrent)
#' dx1
#' plot(dx1, method = "b", stats = c("edges", "concurrent"))
#'
#' # Dynamic diagnostics on the STERGM approximation
#' dx2 <- netdx(est, nsims = 5, nsteps = 500,
#'              nwstats.formula = ~ edges + meandeg + concurrent,
#'              set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
#' dx2
#' plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#' plot(dx2, type = "duration")
#' plot(dx2, type = "dissolution", method = "b", col = "bisque")
#' }
#'
netdx <- function(x,
                  nsims = 1,
                  dynamic = TRUE,
                  nsteps,
                  nwstats.formula = "formation",
                  set.control.ergm,
                  set.control.stergm,
                  keep.tedgelist = FALSE,
                  verbose = TRUE,
                  ncores = 1) {

  if (class(x) != "netest") {
    stop("x must be an object of class netest", call. = FALSE)
  }

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

  if (dynamic == TRUE && missing(nsteps)) {
    stop("Specify number of time steps with nsteps", call. = FALSE)
  }

  if (dynamic == FALSE && nwstats.formula == "formation") {
    nwstats.formula <- x$formation
  }

  if (verbose == TRUE) {
    cat("\n======================")
    cat("\nRunning Diagnostics")
    cat("\n======================\n")
  }

  if (verbose == TRUE) {
    cat("- Simulating", nsims, "networks")
  }

  if (edapprox == FALSE) {

    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.stergm()
    }

    if (nsims == 1 || ncores == 1) {
      diag.sim <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      for (i in 1:nsims) {
        diag.sim[[i]] <- simulate(fit,
                                  time.slices = nsteps,
                                  monitor = nwstats.formula,
                                  nsim = 1,
                                  control = set.control.stergm)
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    } else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)

      diag.sim <- foreach(i = 1:nsims) %dopar% {
        simulate(fit,
                 time.slices = nsteps,
                 monitor = nwstats.formula,
                 nsim = 1,
                 control = set.control.stergm)
      }
    }
  }

  if (edapprox == TRUE) {

    if (missing(set.control.ergm)) {
      set.control.ergm <- control.simulate.ergm()
    }
    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.network()
    }

    if (dynamic == TRUE) {
      if (nsims == 1 || ncores == 1) {
        diag.sim <- list()
        if (verbose == TRUE & nsims > 1) {
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
          if (verbose == TRUE & nsims > 1) {
            cat("*")
          }
        }
        if (verbose == TRUE & nsims > 1) {
          cat("|")
        }
      } else {
        cluster.size <- min(nsims, ncores)
        registerDoParallel(cluster.size)

        diag.sim <- foreach(i = 1:nsims) %dopar% {
          if (class(x$fit) == "network") {
            fit.sim <- simulate(formation,
                                basis = nw,
                                coef = x$coef.form.crude,
                                constraints = constraints)
          } else {
            fit.sim <- simulate(fit, control = set.control.ergm)
          }
          simulate(fit.sim,
                   formation = formation,
                   dissolution = dissolution,
                   coef.form = coef.form,
                   coef.diss = coef.diss$coef.crude,
                   time.slices = nsteps,
                   constraints = constraints,
                   monitor = nwstats.formula,
                   nsim = 1,
                   control = set.control.stergm)
        }
      }
    }
    if (dynamic == FALSE) {
      if (class(x$fit) == "network") {
        diag.sim <- simulate(formation,
                             basis = nw,
                             coef = x$coef.form.crude,
                             constraints = constraints,
                             nsim = nsims,
                             statsonly = TRUE,
                             monitor = nwstats.formula)
      } else {
        diag.sim <- simulate(fit, nsim = nsims,
                             statsonly = TRUE,
                             control = set.control.ergm,
                             monitor = nwstats.formula)
      }
    }
  } # end edapprox = TRUE condition

  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  if (dynamic == TRUE) {
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
  } else {
    stats <- list(diag.sim[,!duplicated(colnames(diag.sim)), drop = FALSE])
    merged.stats <- diag.sim[,!duplicated(colnames(diag.sim)), drop = FALSE]
  }


  ## Calculate mean/sd from merged stats
  stats.means <- colMeans(merged.stats)
  stats.sd <- apply(merged.stats, 2, sd)
  stats.table <- data.frame(sorder = 1:length(names(stats.means)),
                            names = names(stats.means),
                            stats.means, stats.sd)


  ## Get stats from for target statistics
  ts.attr.names <- x$target.stats.names
  target.stats <- target.stats[which(target.stats > 0)]
  ts.out <- data.frame(names = ts.attr.names,
                       targets = target.stats)


  ## Create stats.formation table for output
  stats.table <- merge(ts.out, stats.table, all = TRUE)
  stats.table <- stats.table[do.call("order",
                                     stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
  rownames(stats.table) <- stats.table$names

  stats.table$reldiff <- (stats.table$stats.means-stats.table$targets)/stats.table$targets
  stats.table.formation <- stats.table[, c(2, 4, 6, 5)]
  colnames(stats.table.formation) <- c("Target", "Sim Mean", "Pct Diff", "Sim SD")


  if (dynamic == TRUE) {
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
      for (i in 2:length(diag.sim)) {
        ncens <- which(sim.df[[i]]$onset.censored == FALSE &
                         sim.df[[i]]$terminus.censored == FALSE)
        durVec <- c(durVec, sim.df[[i]]$duration[ncens])
      }
    }



    # Calculate mean partnership age from edgelist
    if (nsims == 1 || ncores == 1) {
      pages <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      for (i in 1:length(diag.sim)) {
        pages[[i]] <- edgelist_meanage(el = sim.df[[i]])
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    } else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)

      pages <- foreach(i = 1:nsims) %dopar% {
        edgelist_meanage(el = sim.df[[i]])
      }
    }

    ## Dissolution calculations
    if (verbose == TRUE) {
      cat("\n- Calculating dissolution statistics")
    }

    ## Create a list of dissolution proportions (i.e. dissolutions/edges)
    if (nsims == 1 || ncores == 1) {
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      prop.diss <- list()
      for (i in 1:length(diag.sim)) {
        prop.diss[[i]] <- sapply(1:nsteps, function(x) sum(sim.df[[i]]$terminus==x) /
                                   sum(sim.df[[i]]$onset < x &
                                         sim.df[[i]]$terminus>=x))
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    } else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)

      prop.diss <- foreach(i = 1:nsims) %dopar% {
        sapply(1:nsteps, function(x) sum(sim.df[[i]]$terminus==x) /
                 sum(sim.df[[i]]$onset < x &
                       sim.df[[i]]$terminus>=x))
      }
    }

    if (verbose == TRUE) {
      cat("\n ")
    }


    # Create dissolution tables
    duration.mean <- mean(durVec)
    duration.sd <- sd(durVec)
    duration.expected <- exp(coef.diss$coef.crude[1]) + 1
    duration.pctdiff <- (duration.mean-duration.expected)/duration.expected


    dissolution.mean <- mean(unlist(prop.diss))
    dissolution.sd <- sd(unlist(prop.diss))
    dissolution.expected <- 1/(exp(coef.diss$coef.crude[1]) + 1)
    dissolution.pctdiff <- (dissolution.mean-dissolution.expected)/dissolution.expected

    stats.table.dissolution <- data.frame(Targets = c(duration.expected,
                                                      dissolution.expected),
                                          Sim_Means = c(duration.mean,
                                                        dissolution.mean),
                                          Pct_Diff = c(duration.pctdiff,
                                                       dissolution.pctdiff),
                                          Sim_SD = c(duration.sd,
                                                     dissolution.sd))
    colnames(stats.table.dissolution) <- c("Target", "Sim Mean", "Pct Diff", "Sim SD")
    rownames(stats.table.dissolution) <- c("Edge Duration", "Pct Edges Diss")

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
  out$dynamic <- dynamic

  out$stats <- stats
  out$stats.table.formation <- stats.table.formation
  if (dynamic == TRUE) {
    out$nsteps <- nsteps
    out$stats.table.dissolution <- stats.table.dissolution
    out$edgelist <- sim.df
    out$pages <- pages
    out$prop.diss <- prop.diss
    if (keep.tedgelist == TRUE) {
      out$tedgelist <- sim.df
    }
  }

  class(out) <- "netdx"
  return(out)
}
