#' @title Dynamic Network Model Diagnostics
#'
#' @description Runs diagnostic simulations on an ERGM/STERGM estimated with
#'              [`netest`] to assess whether the fitted model reproduces the
#'              intended network features. Both static (cross-sectional) and
#'              dynamic (temporal) diagnostics are supported. This is the
#'              recommended second step in the network modeling pipeline, after
#'              estimation with [`netest`] and before epidemic simulation with
#'              [`netsim`].
#'
#' @param x An `EpiModel` object of class `netest`.
#' @param nsims Number of simulations to run. For dynamic diagnostics, 5--10
#'        simulations are usually sufficient to assess model fit. For static
#'        diagnostics, use 10,000+ draws to obtain stable estimates.
#' @param dynamic If `TRUE`, runs dynamic diagnostics that simulate the
#'        temporal network forward in time, checking both formation targets
#'        and partnership duration/dissolution. If `FALSE`, draws from the
#'        static ERGM fit to check cross-sectional network structure only
#'        (faster, but does not verify dissolution dynamics). Static
#'        diagnostics are only available when the model was fit with the
#'        edges dissolution approximation (`edapprox = TRUE` in [`netest`]).
#' @param nsteps Number of time steps per simulation (dynamic simulations
#'        only). Should be at least several multiples of the longest target
#'        partnership duration to allow the duration and dissolution statistics
#'        to stabilize. For example, if the target duration is 50, running for
#'        500 time steps is a reasonable starting point.
#' @param nwstats.formula A right-hand sided ERGM formula with the network
#'        statistics of interest. The default is the formation formula of the
#'        network model contained in `x`. You may track additional network
#'        statistics beyond the formation terms by specifying them here, such
#'        as `~ edges + meandeg + concurrent + degree(0:4)`. This is useful
#'        for verifying that the model produces reasonable values for network
#'        features that were not directly targeted in the formation model.
#' @param set.control.ergm Control arguments passed to `ergm`'s
#'        `simulate_formula.network` (see details).
#' @param set.control.tergm Control arguments passed to `tergm`'s
#'        `simulate_formula.network` (see details).
#' @param sequential For static diagnostics (`dynamic=FALSE`): if
#'        `FALSE`, each of the `nsims` simulated Markov chains begins
#'        at the initial network; if `TRUE`, the end of one simulation is
#'        used as the start of the next.
#' @param keep.tedgelist If `TRUE`, keep the timed edgelist generated from
#'        the dynamic simulations. Returned in the form of a list of matrices,
#'        with one entry per simulation. Accessible at `$edgelist`.
#' @param keep.tnetwork If `TRUE`, keep the full networkDynamic objects
#'        from the dynamic simulations. Returned in the form of a list of nD
#'        objects, with one entry per simulation. Accessible at `$network`.
#' @param verbose If `TRUE`, print progress to the console.
#' @param ncores Number of processor cores to run multiple simulations
#'        on, using the `future` framework.
#' @param skip.dissolution If `TRUE`, skip over the calculations of
#'        duration and dissolution stats in `netdx`.
#' @param future.use.plan If `FALSE`, `netdx` will use `multisession` is used with `workers = ncores for its
#'        parallelization. If `TRUE`, `netdx` will use the user defined plan from `globalEnv`. Finally, it can
#'        take the output of a `future::tweak()` call to setup a user defined temporary plan within `netdx`.
#'        Which can be useful for distributed computation (HPC).
#'
#' @details
#' The `netdx` function handles dynamic network diagnostics for network
#' models fit with the [`netest`] function. Given the fitted model,
#' `netdx` simulates a specified number of dynamic networks for a specified
#' number of time steps per simulation. The network statistics in
#' `nwstats.formula` are saved for each time step. Summary statistics for
#' the formation model terms, as well as dissolution model and relational
#' duration statistics, are then calculated and can be accessed when printing or
#' plotting the `netdx` object. See [`print.netdx`] and [`plot.netdx`]
#' for details on printing and plotting.
#'
#' @section Control Arguments:
#' Models fit with the full STERGM method in `netest` (setting the
#' `edapprox` argument to `FALSE`) require only a call to
#' `tergm`'s `simulate_formula.network`. Control parameters for those
#' simulations may be set using `set.control.tergm` in `netdx`.
#' The parameters should be input through the `control.simulate.formula.tergm`
#' function, with the available parameters listed in the
#' [`tergm::control.simulate.formula.tergm`] help page in the `tergm`
#' package.
#'
#' Models fit with the ERGM method with the edges dissolution approximation
#' (setting `edapprox` to `TRUE`) require a call first to
#' `ergm`'s `simulate_formula.network` for simulating an initial
#' network, and second to `tergm`'s `simulate_formula.network` for
#' simulating that static network forward through time. Control parameters may
#' be set for both processes in `netdx`. For the first, the parameters
#' should be input through the `control.simulate.formula()` function, with
#' the available parameters listed in the
#' [`ergm::control.simulate.formula`] help
#' page in the `ergm` package. For the second, parameters should be input
#' through the `control.simulate.formula.tergm()` function, with the
#' available parameters listed in the [`tergm::control.simulate.formula.tergm`]
#' help page in the `tergm` package. An example is shown below.
#'
#' @section Static vs. Dynamic Diagnostics:
#' Static diagnostics (`dynamic = FALSE`) draw many independent networks from
#' the fitted ERGM and compare the resulting statistics to the target values.
#' This is fast and checks whether the cross-sectional structure is correct,
#' but it does not verify partnership durations or dissolution rates. Dynamic
#' diagnostics (`dynamic = TRUE`) simulate the full temporal network forward
#' in time, checking both formation targets and dissolution/duration dynamics.
#' Dynamic diagnostics are slower but more comprehensive, and are required to
#' verify models that will be used with vital dynamics (arrivals/departures).
#'
#' @section Interpreting Diagnostics:
#' After running `netdx`, use `print()` and [`plot.netdx`] to inspect the
#' results. Key indicators of a good model fit include:
#'
#'  * **Formation statistics:** The "Sim Mean" should be close to the
#'    "Target" value. A small "Pct Diff" (< 5\%) and a "Z Score" near 0
#'    indicate good fit.
#'  * **Duration statistics** (dynamic only): The simulated mean edge
#'    durations should match the values passed to [`dissolution_coefs`].
#'  * **Dissolution statistics** (dynamic only): The simulated dissolution
#'    rates should be approximately `1 / duration`.
#'
#' Common problems: If formation statistics are off, the ERGM may need
#' increased burn-in (via `set.control.ergm`), or the target statistics
#' may be incompatible (e.g., specifying more edges than the network can
#' support). If durations are off but formation is correct, verify that
#' `d.rate` was correctly specified in [`dissolution_coefs`] for models
#' with vital dynamics.
#'
#' @return
#' A list of class `netdx`. Use `print()` to view summary tables of
#' formation statistics, duration, and dissolution diagnostics. Use
#' [`plot.netdx`] to visualize these diagnostics over time. Use
#' [as.data.frame.netdx()] to extract timed edgelists (if
#' `keep.tedgelist = TRUE`).
#'
#' @seealso
#' Estimate the network model with [`netest`] before running diagnostics.
#' Plot diagnostics with [`plot.netdx`] and print summary tables with
#' [`print.netdx`]. After diagnostics confirm a good fit, simulate the
#' epidemic with [`netsim`].
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Static diagnostics on a simple model
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#' dx <- netdx(est, nsims = 1e4, dynamic = FALSE, verbose = FALSE)
#' dx
#' plot(dx)
#' }
#'
#' \dontrun{
#' # Static diagnostics with additional network statistics
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
#' est2 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#' dx3 <- netdx(est2, nsims = 5, nsteps = 100)
#' print(dx3)
#' plot(dx3)
#' plot(dx3, type = "duration", plots.joined = TRUE, qnts = 0.2, legend = TRUE)
#' plot(dx3, type = "dissolution", mean.smooth = FALSE, mean.col = "red")
#' }
#'
netdx <- function(x, nsims = 1, dynamic = TRUE, nsteps = NULL,
                  nwstats.formula = "formation",
                  set.control.ergm = control.simulate.formula(),
                  set.control.tergm = control.simulate.formula.tergm(MCMC.maxchanges = .Machine$integer.max),
                  sequential = TRUE, keep.tedgelist = FALSE,
                  keep.tnetwork = FALSE, verbose = TRUE, ncores = 1,
                  skip.dissolution = FALSE, future.use.plan = FALSE) {

  if (!inherits(x, "netest")) {
    stop("x must be an object of class netest")
  }

  ncores <- if (nsims == 1) 1 else min(parallelly::availableCores(), ncores)

  formation <- x$formation
  coef.form <- x$coef.form
  dissolution <- x$coef.diss$dissolution
  coef.diss <- x$coef.diss
  constraints <- x$constraints
  target.stats <- x$target.stats
  edapprox <- x$edapprox
  nw <- x$newnetwork

  if (dynamic && is.null(nsteps)) {
    stop("Specify number of time steps with nsteps")
  }

  if (x$coef.diss$duration[1] == 1 && dynamic == TRUE) {
    stop("Running dynamic diagnostics on a cross-sectional ERGM (duration = 1) is not possible.
         \nSet netdx parameter 'dynamic' to 'FALSE'"
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

    init %n% "time" <- 0L
    init %n% "lasttoggle" <- cbind(as.edgelist(init), 0L)

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
  } else if (ncores == 1 && isFALSE(future.use.plan)) {
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
    if (inherits(future.use.plan, c("tweaked", "future"))) {
      with(future::plan(future.use.plan), local = TRUE)
    } else if (ncores > 1 && isFALSE(future.use.plan)) {
      ncores_eff <- min(nsims, ncores)
      with(future::plan("multisession", workers = ncores_eff), local = TRUE)
    } # else if (isTRUE(future.use.plan)) - use plan defined by user
    diag.sim <- future.apply::future_replicate(
      nsims, dosim(), simplify = FALSE, future.seed = TRUE
    )
  }

  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  stats <- lapply(diag.sim, function(x) x$stats[, !duplicated(colnames(x$stats)), drop = FALSE])
  ts.attr.names <- x$target.stats.names

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
##   if coda::effectiveSize throws an error; argument x should be a matrix
##   with column names
ess <- function(x) {
  tryCatch(coda::effectiveSize(x),
           error = function(e) {
             structure(rep(NA, length.out = NCOL(x)), names = colnames(x))
           })
}

#' @title Create a Summary Table of Simulation Statistics
#'
#' @param stats A list of simulated statistics matrices, of length equal to the
#'   number of simulations performed.  Each matrix should have one row for each
#'   simulated network if `dynamic == FALSE`, one row for each time step
#'   if `dynamic == TRUE`, and one column for each statistic.  The columns
#'   should be named for the statistics they correspond to, with all matrices
#'   having the same statistics, in the same order.  Each matrix may have an
#'   `attr`-style attribute named `"ess"` attached, giving the
#'   effective sample sizes for the columns of the matrix; if this attribute is
#'   `NULL`, then the effective sample sizes will be computed within the
#'   call to `make_stats_table`.
#' @param targets A vector of target values for the statistics in `stats`.
#'   May be named (in which case targets will be matched to statistics based on
#'   column names in matrices in `stats`) or unnamed (in which case
#'   targets will be matched to statistics based on position, and the number of
#'   targets must equal the number of columns).
#'
#' @return A `data.frame` summarizing the simulated statistics.
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
#'   [`networkDynamic::as.data.frame.networkDynamic`].
#'
#' @return The matrix of toggles corresponding to `tedgelist`.
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
#'   [tedgelist_to_toggles()].
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

  # nolint start
  diss_formula <- coef.diss$dissolution
  # nolint end
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
