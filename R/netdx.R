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
#' @param set.control.stergm Deprecated control argument of class
#'        \code{control.simulate.network}; use \code{set.control.tergm}
#'        instead.
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
                  set.control.stergm = control.simulate.network(),
                  set.control.tergm = control.simulate.formula.tergm(),
                  sequential = TRUE, keep.tedgelist = FALSE,
                  keep.tnetwork = FALSE, verbose = TRUE, ncores = 1,
                  skip.dissolution = FALSE) {

  if (!inherits(x, "netest")) {
    stop("x must be an object of class netest", call. = FALSE)
  }

  STERGM <- !missing(set.control.stergm)
  if (STERGM == TRUE) {
    warning("set.control.stergm is deprecated and will be removed in a future
             version; use set.control.tergm instead.")
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

  if (any(x$coef.diss$duration == 1) & dynamic == TRUE) {
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
      return(list(stats = init))
    }

    if (keep.tedgelist == TRUE || keep.tnetwork == TRUE) {
      output <- "networkDynamic"
    } else {
      output <- "changes"
    }

    if (STERGM == TRUE) {
      diag.sim <- simulate(init,
                           formation = x$formation,
                           dissolution = x$coef.diss$dissolution,
                           coef.form = x$coef.form,
                           coef.diss = x$coef.diss$coef.crude,
                           constraints = x$constraints,
                           time.slices = nsteps,
                           monitor = nwstats.formula,
                           time.start = 0L,
                           nsim = 1L,
                           output = output,
                           control = set.control.stergm)
    } else {
      diag.sim <- simulate(init ~ Form(x$formation) +
                             Persist(x$coef.diss$dissolution),
                           coef = c(x$coef.form, x$coef.diss$coef.crude),
                           constraints = x$constraints,
                           time.slices = nsteps,
                           monitor = nwstats.formula,
                           time.start = 0L,
                           nsim = 1L,
                           output = output,
                           control = set.control.tergm,
                           dynamic = TRUE)
    }

    out <- list(stats = attr(diag.sim, "stats"))

    if (output == "networkDynamic") {
      sim.df <- as.data.frame(diag.sim)
      toggles <- tedgelist_to_toggles(sim.df)
    } else {
      changes <- diag.sim
      if (network.edgecount(init) > 0L) {
        changes <- rbind(cbind(0L, as.edgelist(init), 1L),
                         changes)
      }
      toggles <- changes[, -4L, drop=FALSE]
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

  ## Calculate mean/sd from stats
  stats.table.formation <- make_formation_table(stats, ts.out)

  # Calculate dissolution / duration stats
  if (skip.dissolution == FALSE) {
    if (dynamic == TRUE) {
      dissolution.stats <- make_dissolution_stats(
        diag.sim,
        x$coef.diss,
        nsteps,
        verbose
      )
    }
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

  class(out) <- "netdx"
  return(out)
}

#' @title Calculate the Formation Statistics of a Network
#'
#' @param stats A list of formation statistic matrices, one for each 
#'   simulation, with one row for each simulated network if 
#'   \code{dynamic == FALSE}, one row for each time step if 
#'   \code{dynamic == TRUE}, and one column for each statistic.
#' @param targets A \code{data.frame} of the formation targets with two columns:
#'   "names" and "targets".
#'
#' @return A \code{data.frame} of the formation statistics.
#' @keywords internal
make_formation_table <- function(stats, targets) {
  ess <- lapply(stats, function(x) apply(x, 2L, function(y) if (sum(!is.na(y)) <= 1L) NA else effectiveSize(na.omit(y))))
  ess <- colSums(do.call(rbind, ess), na.rm = TRUE)

  stats.onesim.sd <- apply(do.call(rbind, lapply(stats, colMeans, na.rm = TRUE)), 2, sd, na.rm = TRUE)
  
  stats <- do.call(rbind, stats)
  stats.means <- colMeans(stats, na.rm = TRUE)
  stats.sd <- apply(stats, 2L, sd, na.rm = TRUE)
  stats.se <- stats.sd/sqrt(ess)
  
  stats.table <- data.frame(
    sorder = seq_along(names(stats.means)),
    names = names(stats.means),
    stats.means,
    stats.se,
    stats.onesim.sd,
    stats.sd
  )

  ## Create stats.formation table for output
  stats.table <- merge(targets, stats.table, all = TRUE)
  stats.table <- stats.table[order(stats.table[["sorder"]]), , drop = FALSE]
  rownames(stats.table) <- stats.table$names

  stats.table$reldiff <- (stats.table$stats.means - stats.table$targets) /
    stats.table$targets * 100
  stats.table.formation <- stats.table[, c(2, 4, 8, 5)]
  stats.table.formation <- cbind(stats.table.formation, 
                                 zscore = (stats.table.formation[,2] - stats.table.formation[,1]) / 
                                           stats.table.formation[,4],
                                 onesim.sd = stats.table[,6],
                                 stats.sd = stats.table[,7])
  colnames(stats.table.formation) <- c(
    "Target",
    "Sim Mean",
    "Pct Diff",
    "Sim SE",
    "Z Score",
    "SD(1-Sim Mean)",
    "SD(Statistic)"
  )

  return(stats.table.formation)
}

#' @title Calculate the Dissolution Statistics of a Network
#'
#' @param diag.sim A list of dissolution statistics (as created by 
#'   \code{toggles_to_diss_stats}), of length equal to the number of 
#'   simulations.
#' @param coef.diss The \code{coef.diss} element of \code{nwparam}.
#' @param nsteps The number of simulated steps.
#' @param verbose A verbosity toggle (default = TRUE).
#'
#' @return A \code{list} of dissolution statistics, combined across 
#'   simulations.
#' @keywords internal
make_dissolution_stats <- function(diag.sim, coef.diss,
                                   nsteps, verbose = TRUE) {
  if (verbose == TRUE) {
    cat("\n- Calculating duration statistics")
  }

  if (any(unlist(lapply(diag.sim, `[[`, "anyNA")))) {
    warning("duration/dissolution data contains undefined values due to",
            " having zero edges of some dissolution dyad type(s) on some time",
            " step(s); these undefined values will be set to 0 when",
            " processing the data; this behavior, which introduces a bias",
            " towards 0, may be changed in the future")
  }

  ## exclude nodefactor from heterogeneous dissolution calculation
  if (coef.diss$diss.model.type == "nodefactor") {
    durs <- mean(coef.diss$duration)
  } else {
    durs <- coef.diss$duration
  }
  
  meanage_list <- lapply(diag.sim, `[[`, "meanage")
  pages <- array(unlist(meanage_list),
                 dim = c(nsteps, length(durs), length(diag.sim)))

  meanage_imptd_list <- lapply(diag.sim, `[[`, "meanageimputed")
  pages_imptd <- array(unlist(meanage_imptd_list),
                       dim = c(nsteps, length(durs), length(diag.sim)))

  propdiss_list <- lapply(diag.sim, `[[`, "propdiss")
  propdiss <- array(unlist(propdiss_list),
                    dim = c(nsteps, length(durs), length(diag.sim)))

  meanage_imptd_mean <- apply(pages_imptd, 2, mean, na.rm = TRUE)
  meanage_imptd_sd <- apply(pages_imptd, 2, sd, na.rm = TRUE)
  meanage_imptd_ess <- lapply(meanage_imptd_list, function(x) apply(x, 2L, function(y) if (sum(!is.na(y)) <= 1L) NA else effectiveSize(na.omit(y))))
  meanage_imptd_ess <- colSums(do.call(rbind, meanage_imptd_ess), na.rm = TRUE)
  meanage_imptd_se <- meanage_imptd_sd/sqrt(meanage_imptd_ess)

  meanage_imptd_onesim_sd <- apply(do.call(rbind, lapply(diag.sim, `[[`, "meanmeanageimputed")), 2, sd, na.rm = TRUE)

  propdiss_mean <- apply(propdiss, 2, mean, na.rm = TRUE)
  propdiss_sd <- apply(propdiss, 2, sd, na.rm = TRUE)
  propdiss_ess <- lapply(propdiss_list, function(x) apply(x, 2L, function(y) if (sum(!is.na(y)) <= 1L) NA else effectiveSize(na.omit(y))))
  propdiss_ess <- colSums(do.call(rbind, propdiss_ess), na.rm = TRUE)
  propdiss_se <- propdiss_sd/sqrt(propdiss_ess)

  propdiss_onesim_sd <- apply(do.call(rbind, lapply(diag.sim, `[[`, "meanpropdiss")), 2, sd, na.rm = TRUE)
  
  stats.table.duration <- data.frame("Target" = durs,
                                     "Sim Mean" = meanage_imptd_mean,
                                     "Pct Diff" = 100*(meanage_imptd_mean - durs)/durs,
                                     "Sim SE" = meanage_imptd_se,
                                     "Z Score" = (meanage_imptd_mean - durs)/meanage_imptd_se,
                                     "SD(1-Sim Mean)" = meanage_imptd_onesim_sd,
                                     "SD(Statistic)" = meanage_imptd_sd)
  colnames(stats.table.duration) <- c("Target", "Sim Mean", "Pct Diff", "Sim SE", "Z Score", "SD(1-Sim Mean)", "SD(Statistic)")
  
  stats.table.dissolution <- data.frame("Target" = 1/durs,
                                        "Sim Mean" = propdiss_mean,
                                        "Pct Diff" = 100*(propdiss_mean - 1/durs)/(1/durs),
                                        "Sim SE" = propdiss_se,
                                        "Z Score" = (propdiss_mean - 1/durs)/propdiss_se,
                                        "SD(1-Sim Mean)" = propdiss_onesim_sd,
                                        "SD(Statistic)" = propdiss_sd)
  colnames(stats.table.dissolution) <- c("Target", "Sim Mean", "Pct Diff", "Sim SE", "Z Score", "SD(1-Sim Mean)", "SD(Statistic)")

  # Construct return list
  return(
    list(
      "stats.table.duration" = stats.table.duration,
      "stats.table.dissolution" = stats.table.dissolution,
      "pages" = pages,
      "pages_imptd" = pages_imptd,
      "prop.diss" = propdiss
    )
  )
}

tedgelist_to_toggles <- function(tedgelist) {
  tedgelist <- as.matrix(tedgelist)
  toggles <- rbind(tedgelist[, c(1L, 3L, 4L), drop = FALSE],
                   tedgelist[!tedgelist[, 6L], c(2L, 3L, 4L), drop = FALSE])
  colnames(toggles) <- c("time", "tail", "head")
  toggles
}

toggles_to_diss_stats <- function(toggles, coef.diss,
                                  nsteps, nw, time.start = 0L) {
  nw <- as.network(nw) # drop nwd
  delete.network.attribute(nw, "time")
  delete.network.attribute(nw, "lasttoggle")
  nw[, ] <- FALSE

  ## exclude nodefactor from heterogeneous dissolution calculation
  if (coef.diss$diss.model.type == "nodefactor") {
    diss_formula <- ~offset(edges)
    durs <- mean(coef.diss$duration)
  } else {
    diss_formula <- coef.diss$dissolution
    durs <- coef.diss$duration
  }

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

  edgeages <- changestats[, seq_along(durs)+length(durs), drop = FALSE]
  colnames(edgeages) <- colnames(edgecounts)

  edgepers <- changestats[, seq_along(durs) + 2L * length(durs), drop = FALSE]

  if (length(durs) > 1L) {
    edgecounts[,1L] <- edgecounts[, 1L] - rowSums(edgecounts[, -1L, drop = FALSE])
    edgeages[, 1L] <- edgeages[, 1L] - rowSums(edgeages[, -1L, drop = FALSE])
    edgepers[, 1L] <- edgepers[, 1L] - rowSums(edgepers[, -1L, drop = FALSE])
  }

  edgediss <- edgecounts[-NROW(edgecounts), , drop = FALSE] -
    edgepers[-1L, , drop = FALSE]
  edgeages <- edgeages[-1L, , drop = FALSE]

  edgeagesimputed <- edgeages
  toggles <- toggles[order(toggles[, 2L], toggles[, 3L],
                           toggles[, 1L]), , drop = FALSE]
  w <- which(toggles[,1] == time.start)
  if (length(w) > 0L) {
    # imputation
    changestats <- as.matrix(tergm.godfather(nw ~ Passthrough(diss_formula),
                                             toggles = cbind(seq_along(w),
                                                             toggles[w, -1L, drop = FALSE]),
                                             stats.start = TRUE))

    for (i in seq_along(w)) {
      dyad_type <- max(which(changestats[i, ] != changestats[i+1L, ]))
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
          rgeom(1L, 1/durs[dyad_type])
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

  meanmeanageimputed <- colMeans(meanageimputed, na.rm = TRUE)
  meanpropdiss <- colMeans(propdiss, na.rm = TRUE)

  return(list(edgecounts = edgecounts,
              edgeages = edgeages,
              edgeagesimputed = edgeagesimputed,
              edgediss = edgediss,
              meanage = meanage,
              meanageimputed = meanageimputed,
              propdiss = propdiss,
              meanpropdiss = meanpropdiss,
              meanmeanageimputed = meanmeanageimputed,
              anyNA = anyNA))
}
