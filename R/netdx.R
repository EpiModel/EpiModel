#' @title Dynamic Network Model Diagnostics
#'
#' @description Runs dynamic diagnostics on an ERGM/STERGM estimated through
#'              \code{\link{netest}}.
#'
#' @param x An \code{EpiModel} object of class \code{netest}.
#' @param nsims Number of simulations to run.
#' @param dynamic If \code{TRUE}, runs dynamic diagnostics. If \code{FALSE} and
#'        the \code{netest} object was fit with the Edges Dissolution
#'        approximation method, simulates from the static ERGM fit.
#' @param nsteps Number of time steps per simulation (dynamic simulations only).
#' @param nwstats.formula A right-hand sided ERGM formula with the network
#'        statistics of interest. The default is the formation formula of the
#'        network model contained in \code{x}.
#' @param set.control.ergm Control arguments passed to \code{ergm}'s
#'        \code{simulate_formula.network} (see details).
#' @param set.control.tergm Control arguments passed to \code{tergm}'s
#'        \code{simulate_formula.network} (see details).
#' @param sequential For static diagnostics (\code{dynamic=FALSE}): if
#'        \code{FALSE}, each of the \code{nsims} simulated Markov chains begins
#'        at the initial network; if \code{TRUE}, the end of one simulation is
#'        used as the start of the next.
#' @param keep.tedgelist If \code{TRUE}, keep the timed edgelist generated from
#'        the dynamic simulations. Returned in the form of a list of matrices,
#'        with one entry per simulation. Accessible at \code{$edgelist}.
#' @param keep.tnetwork If \code{TRUE}, keep the full networkDynamic objects
#'        from the dynamic simulations. Returned in the form of a list of nD
#'        objects, with one entry per simulation. Accessible at \code{$network}.
#' @param verbose If \code{TRUE}, print progress to the console.
#' @param ncores Number of processor cores to run multiple simulations
#'        on, using the \code{foreach} and \code{doParallel} implementations.
#' @param skip.dissolution If \code{TRUE}, skip over the calculations of
#'        duration and dissolution stats in \code{netdx}.
#'
#' @details
#' The \code{netdx} function handles dynamic network diagnostics for network
#' models fit with the \code{\link{netest}} function. Given the fitted model,
#' \code{netdx} simulates a specified number of dynamic networks for a specified
#' number of time steps per simulation. The network statistics in
#' \code{nwstats.formula} are saved for each time step. Summary statistics for
#' the formation model terms, as well as dissolution model and relational
#' duration statistics, are then calculated and can be accessed when printing or
#' plotting the \code{netdx} object.  See \code{\link{print.netdx}} and
#' \code{\link{plot.netdx}} for details on printing and plotting.
#'
#' @section Control Arguments:
#' Models fit with the full STERGM method in \code{netest} (setting the
#' \code{edapprox} argument to \code{FALSE}) require only a call to
#' \code{tergm}'s \code{simulate_formula.network}. Control parameters for those
#' simulations may be set using \code{set.control.tergm} in \code{netdx}.
#' The parameters should be input through the
#' \code{control.simulate.formula.tergm} function, with the available
#' parameters listed in the \code{\link{control.simulate.formula.tergm}} help
#' page in the \code{tergm} package.
#'
#' Models fit with the ERGM method with the edges dissolution approximation
#' (setting \code{edapprox} to \code{TRUE}) require a call first to
#' \code{ergm}'s \code{simulate_formula.network} for simulating an initial
#' network, and second to \code{tergm}'s \code{simulate_formula.network} for
#' simulating that static network forward through time. Control parameters may
#' be set for both processes in \code{netdx}. For the first, the parameters
#' should be input through the \code{control.simulate.formula()} function, with
#' the available parameters listed in the
#' \code{\link[ergm:control.simulate.formula]{control.simulate.formula}} help
#' page in the \code{ergm} package. For the second, parameters should be input
#' through the \code{control.simulate.formula.tergm()} function, with the
#' available parameters listed in the
#' \code{\link{control.simulate.formula.tergm}} help page in the \code{tergm}
#' package. An example is shown below.
#'
#' @return A list of class \code{netdx}.
#'
#' @seealso Plot these model diagnostics with \code{\link{plot.netdx}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Network initialization and model parameterization
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 25)
#'
#' # Estimate the model
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Static diagnostics on the ERGM fit
#' dx1 <- netdx(est,
#'   nsims = 1e4, dynamic = FALSE,
#'   nwstats.formula = ~ edges + meandeg + concurrent
#' )
#' dx1
#' plot(dx1, method = "b", stats = c("edges", "concurrent"))
#'
#' # Dynamic diagnostics on the STERGM approximation
#' dx2 <- netdx(est,
#'   nsims = 5, nsteps = 500,
#'   nwstats.formula = ~ edges + meandeg + concurrent,
#'   set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6)
#' )
#' dx2
#' plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#' plot(dx2, type = "duration")
#' plot(dx2, type = "dissolution", qnts.col = "orange2")
#' plot(dx2, type = "dissolution", method = "b", col = "bisque")
#'
#' # Dynamic diagnostics on a more complex model
#' nw <- network_initialize(n = 1000)
#' nw <- set_vertex_attribute(nw, "neighborhood", rep(1:10, 100))
#' formation <- ~edges + nodematch("neighborhood", diff = TRUE)
#' target.stats <- c(800, 45, 81, 24, 16, 32, 19, 42, 21, 24, 31)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges) +
#'                     offset(nodematch("neighborhood", diff = TRUE)),
#'                     duration = c(52, 58, 61, 55, 81, 62, 52, 64, 52, 68, 58))
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#' dx11 <- netdx(est, nsims = 5, nsteps = 100)
#' plot(dx11)
#' plot(dx11, type = "duration", plots.joined = TRUE, qnts = 0.2)
#' plot(dx11, type = "dissolution", mean.smooth = FALSE, mean.col = "red")
#' }
#'
netdx <- function(x, nsims = 1, dynamic = TRUE, nsteps,
                  nwstats.formula = "formation",
                  set.control.ergm = control.simulate.formula(),
                  set.control.tergm = control.simulate.formula.tergm(),
                  sequential = TRUE, keep.tedgelist = FALSE,
                  keep.tnetwork = FALSE, verbose = TRUE, ncores = 1,
                  skip.dissolution = FALSE) {

  if (!inherits(x, "netest")) {
    stop("x must be an object of class netest", call. = FALSE)
  }

  ncores <- ifelse(nsims == 1, 1, min(parallel::detectCores(), ncores))

  formation <- x$formation
  coef.form <- x$coef.form
  dissolution <- x$coef.diss$dissolution
  coef.diss <- x$coef.diss
  constraints <- x$constraints
  target.stats <- x$target.stats
  edapprox <- x$edapprox
  nw <- x$newnetwork

  if (dynamic == TRUE && missing(nsteps)) {
    stop("Specify number of time steps with nsteps", call. = FALSE)
  }

  if (any(x$coef.diss$duration == 1) && dynamic == TRUE) {
    stop("Running dynamic diagnostics on a cross-sectional ERGM (duration = 1)
         is not possible. \nSet netdx parameter 'dynamic' to 'FALSE'",
      call. = FALSE
    )
  }

  if (dynamic == FALSE && nwstats.formula == "formation") {
    nwstats.formula <- x$formation
  }

  if (dynamic == FALSE && edapprox == FALSE) {
    stop("cannot perform static simulation without edapprox")
  }

  if (verbose == TRUE) {
    cat("\nNetwork Diagnostics")
    cat("\n-----------------------\n")
  }

  if (verbose == TRUE) {
    if (nsims == 1) {
      cat("- Simulating 1 network")
    } else {
      cat("- Simulating", nsims, "networks")
    }
  }

  dosim <- function() {
    if (edapprox == TRUE) {
      init <- simulate(x$formula,
                       coef = x$coef.form.crude,
                       basis = x$newnetwork,
                       constraints = x$constraints,
                       control = set.control.ergm,
                       dynamic = FALSE,
                       nsim = if (dynamic == TRUE) 1L else nsims,
                       output = if (dynamic == TRUE) "network" else "stats",
                       sequential = if (dynamic == TRUE) FALSE else sequential,
                       monitor = if (dynamic == TRUE) NULL else nwstats.formula)
    } else {
      init <- x$newnetwork
    }

    if (dynamic == FALSE) {
      stats <- init
      attr(stats, "ess") <- ess(stats)
      return(list(stats = stats))
    }

    if (keep.tedgelist == TRUE || keep.tnetwork == TRUE) {
      output <- "networkDynamic"
    } else {
      output <- "changes"
    }

    diag.sim <- simulate(init ~ Form(x$formation) + Persist(x$coef.diss$dissolution),
                         coef = c(x$coef.form, x$coef.diss$coef.crude),
                         constraints = x$constraints,
                         time.slices = nsteps,
                         monitor = nwstats.formula,
                         time.start = 0L,
                         nsim = 1L,
                         output = output,
                         control = set.control.tergm,
                         dynamic = TRUE)

    stats <- attr(diag.sim, "stats")
    attr(stats, "ess") <- ess(stats)
    out <- list(stats = stats)

    if (output == "networkDynamic") {
      sim.df <- as.data.frame(diag.sim)
      toggles <- tedgelist_to_toggles(sim.df)
    } else {
      changes <- diag.sim
      if (network.edgecount(init) > 0L) {
        changes <- rbind(cbind(0L, as.edgelist(init), 1L),
                         changes)
      }
      toggles <- changes[, -4L, drop = FALSE]
    }

    if (keep.tnetwork == TRUE) {
      out$tnetwork <- diag.sim
    }

    if (keep.tedgelist == TRUE) {
      out$tedgelist <- sim.df
    }

    if (skip.dissolution == FALSE) {
      out <- c(out, toggles_to_diss_stats(toggles, x$coef.diss, nsteps, init))
    }

    out
  }

  if (dynamic == FALSE || nsims == 1) {
    diag.sim <- list(dosim())
  } else if (ncores == 1) {
    diag.sim <- list()
    if (verbose == TRUE) {
      cat("\n |")
    }
    for (i in seq_len(nsims)) {
      diag.sim[[i]] <- dosim()
      if (verbose == TRUE) {
        cat("*")
      }
    }
    if (verbose == TRUE) {
      cat("|")
    }
  } else {
    cluster.size <- min(nsims, ncores)
    registerDoParallel(cluster.size)
    diag.sim <- foreach(i = seq_len(nsims)) %dopar% dosim()
  }

  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  stats <- lapply(diag.sim, function(x) x$stats[, !duplicated(colnames(x$stats)), drop = FALSE])

  ts.attr.names <- x$target.stats.names
  if (length(ts.attr.names) != length(target.stats)) {
    target.stats <- target.stats[which(target.stats > 0)]
  }
  ts.out <- data.frame(
    names = ts.attr.names,
    targets = target.stats
  )
  names(target.stats) <- ts.attr.names

  ## Calculate mean/sd from stats
  stats.table.formation <- make_stats_table(stats, target.stats)

  # Calculate dissolution / duration stats
  if (skip.dissolution == FALSE && dynamic == TRUE) {
    durs <- coef.diss$duration

    dims <- c(nsteps, length(durs), length(diag.sim))
    dissolution.stats <- list("stats.table.duration" = make_stats_table(lapply(diag.sim, `[[`, "meanageimputed"), durs),
                              "stats.table.dissolution" = make_stats_table(lapply(diag.sim, `[[`, "propdiss"),
                                                                           1 / durs),
                              "pages" = array(unlist(lapply(diag.sim, `[[`, "meanage")), dim = dims),
                              "pages_imptd" = array(unlist(lapply(diag.sim, `[[`, "meanageimputed")), dim = dims),
                              "prop.diss" = array(unlist(lapply(diag.sim, `[[`, "propdiss")), dim = dims),
                              "anyNA" = any(unlist(lapply(diag.sim, `[[`, "anyNA"))))
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
    if (skip.dissolution == FALSE) {
      out <- c(out, dissolution.stats)
    }
    if (keep.tedgelist == TRUE) {
      out$tedgelist <- lapply(diag.sim, `[[`, "tedgelist")
    }
    if (keep.tnetwork == TRUE) {
      out$network <- lapply(diag.sim, `[[`, "tnetwork")
    }
  }
  out$anyNA <- NVL(out$anyNA, FALSE)

  class(out) <- "netdx"
  return(out)
}

## internal wrapper around coda::effectiveSize, returning a vector of NAs
##   rather than throwing an error when there are 0 or 1 observations; the
##   argument x should be a matrix with column names
ess <- function(x) {
  if (NROW(x) <= 1L) {
    structure(rep(NA, length.out = NCOL(x)), names = colnames(x))
  } else {
    coda::effectiveSize(x)
  }
}

#' @title Create a Summary Table of Simulation Statistics
#'
#' @param stats A list of simulated statistics matrices, of length equal to the
#'   number of simulations performed.  Each matrix should have one row for each
#'   simulated network if \code{dynamic == FALSE}, one row for each time step
#'   if \code{dynamic == TRUE}, and one column for each statistic.  The columns
#'   should be named for the statistics they correspond to, with all matrices
#'   having the same statistics, in the same order.  Each matrix may have an
#'   \code{attr}-style attribute named \code{"ess"} attached, giving the
#'   effective sample sizes for the columns of the matrix; if this attribute is
#'   \code{NULL}, then the effective sample sizes will be computed within the
#'   call to \code{make_stats_table}.
#' @param targets A vector of target values for the statistics in \code{stats}.
#'   May be named (in which case targets will be matched to statistics based on
#'   column names in matrices in \code{stats}) or unnamed (in which case
#'   targets will be matched to statistics based on position, and the number of
#'   targets must equal the number of columns).
#'
#' @return A \code{data.frame} summarizing the simulated statistics.
#' @keywords internal
make_stats_table <- function(stats, targets) {
  ess_list <- lapply(stats, function(x) NVL(attr(x, "ess"), ess(x)))
  ess_sum <- colSums(do.call(rbind, ess_list), na.rm = TRUE)

  stats.onesim.sd <- apply(do.call(rbind, lapply(stats, colMeans, na.rm = TRUE)), 2, sd, na.rm = TRUE)

  stats <- do.call(rbind, stats)
  stats.means <- colMeans(stats, na.rm = TRUE)
  stats.sd <- apply(stats, 2L, sd, na.rm = TRUE)
  stats.se <- stats.sd / sqrt(ess_sum)

  if (!is.null(names(targets))) {
    stats.targets <- rep(NA, length.out = length(stats.means))
    matches <- match(names(targets), names(stats.means))
    stats.targets[na.omit(matches)] <- targets[!is.na(matches)]
  } else {
    stats.targets <- targets
  }

  stats.table <- data.frame("Target" = stats.targets,
                            "Sim Mean" = stats.means,
                            "Pct Diff" = 100 * (stats.means - stats.targets) / stats.targets,
                            "Sim SE" = stats.se,
                            "Z Score" = (stats.means - stats.targets) / stats.se,
                            "SD(Sim Means)" = stats.onesim.sd,
                            "SD(Statistic)" = stats.sd)
  colnames(stats.table) <- c("Target", "Sim Mean", "Pct Diff", "Sim SE", "Z Score", "SD(Sim Means)", "SD(Statistic)")
  rownames(stats.table) <- names(stats.means)

  return(stats.table)
}

#' @title Convert Timed Edgelist to Matrix of Toggles
#'
#' @param tedgelist A timed edgelist, as produced by
#'   \code{\link{as.data.frame.networkDynamic}}.
#'
#' @return The matrix of toggles corresponding to \code{tedgelist}.
#' @keywords internal
tedgelist_to_toggles <- function(tedgelist) {
  tedgelist <- as.matrix(tedgelist)
  toggles <- rbind(tedgelist[, c(1L, 3L, 4L), drop = FALSE],
                   tedgelist[!tedgelist[, 6L], c(2L, 3L, 4L), drop = FALSE])
  colnames(toggles) <- c("time", "tail", "head")
  toggles
}

#' @title Convert Matrix of Toggles to Dissolution and Duration Statistics
#'
#' @param toggles A matrix of toggles, as produced by
#'   \code{\link{tedgelist_to_toggles}}.
#' @param coef.diss Dissolution coefficients used in the simulation.
#' @param nsteps Number of time steps in the simulation.
#' @param nw Network used in the simulation.
#' @param time.start Starting time for the simulation.
#'
#' @return Named list containing dissolution and duration statistics matrices
#'   and other related information.
#' @keywords internal
toggles_to_diss_stats <- function(toggles, coef.diss,
                                  nsteps, nw, time.start = 0L) {
  nw <- as.network(nw) # drop nwd
  delete.network.attribute(nw, "time")
  delete.network.attribute(nw, "lasttoggle")
  nw[, ] <- FALSE

  diss_formula <- coef.diss$dissolution
  durs <- coef.diss$duration

  changestats <- as.matrix(tergm.godfather(nw ~ Passthrough(diss_formula)
                                                + EdgeAges(diss_formula)
                                                + Persist(diss_formula),
                                           toggles = toggles,
                                           start = time.start - 1L,
                                           end = time.start + nsteps,
                                           stats.start = FALSE))

  edgecounts <- changestats[, seq_along(durs), drop = FALSE]

  # drop offset() from names
  colnames(edgecounts) <- substr(colnames(edgecounts), 8L,
                                 nchar(colnames(edgecounts)) - 1L)

  edgeages <- changestats[, seq_along(durs) + length(durs), drop = FALSE]
  colnames(edgeages) <- colnames(edgecounts)

  edgepers <- changestats[, seq_along(durs) + 2L * length(durs), drop = FALSE]

  if (length(durs) > 1L) {
    edgecounts[, 1L] <- edgecounts[, 1L] - rowSums(edgecounts[, -1L, drop = FALSE])
    edgeages[, 1L] <- edgeages[, 1L] - rowSums(edgeages[, -1L, drop = FALSE])
    edgepers[, 1L] <- edgepers[, 1L] - rowSums(edgepers[, -1L, drop = FALSE])
  }

  edgediss <- edgecounts[-NROW(edgecounts), , drop = FALSE] -
    edgepers[-1L, , drop = FALSE]
  edgeages <- edgeages[-1L, , drop = FALSE]

  edgeagesimputed <- edgeages
  toggles <- toggles[order(toggles[, 2L], toggles[, 3L],
                           toggles[, 1L]), , drop = FALSE]
  w <- which(toggles[, 1] == time.start)
  if (length(w) > 0L) {
    # imputation
    changestats <- as.matrix(tergm.godfather(nw ~ Passthrough(diss_formula),
                                             toggles = cbind(seq_along(w),
                                                             toggles[w, -1L, drop = FALSE]),
                                             stats.start = TRUE))

    for (i in seq_along(w)) {
      dyad_type <- max(which(changestats[i, ] != changestats[i + 1L, ]))
      index <- w[i]
      if (index < NROW(toggles) &&
            toggles[index, 2L] == toggles[index + 1L, 2L] &&
            toggles[index, 3L] == toggles[index + 1L, 3L]) {
        terminus_time <- toggles[index + 1L, 1L]
      } else {
        terminus_time <- time.start + nsteps + 1L
      }
      if (terminus_time > time.start + 1L) {
        edgeagesimputed[seq_len(terminus_time - time.start - 1L), dyad_type] <-
          edgeagesimputed[seq_len(terminus_time - time.start - 1L), dyad_type] +
          rgeom(1L, 1 / durs[dyad_type])
      }
    }
  }

  ## 0/0 is possible, resulting in NaN, which we set to 0 for the time being...
  meanage <- edgeages / edgecounts[-1L, , drop = FALSE]
  meanageimputed <- edgeagesimputed / edgecounts[-1L, , drop = FALSE]
  propdiss <- edgediss / edgecounts[-NROW(edgecounts), , drop = FALSE]

  if (any(is.na(meanage)) || any(is.na(meanageimputed)) || any(is.na(propdiss))) {
    meanage[is.na(meanage)] <- 0
    meanageimputed[is.na(meanageimputed)] <- 0
    propdiss[is.na(propdiss)] <- 0
    anyNA <- TRUE
  } else {
    anyNA <- FALSE
  }

  attr(meanage, "ess") <- ess(meanage)
  attr(meanageimputed, "ess") <- ess(meanageimputed)
  attr(propdiss, "ess") <- ess(propdiss)

  return(list(meanage = meanage,
              meanageimputed = meanageimputed,
              propdiss = propdiss,
              anyNA = anyNA))
}
