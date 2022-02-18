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
#' @param set.control.ergm Control arguments passed to \code{simulate.ergm} (see
#'        details).
#' @param set.control.stergm Control arguments passed to \code{simulate.stergm}
#'        (see details).
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
#' plotting the \code{netdx} object.
#'
#' @section Control Arguments:
#' Models fit with the full STERGM method in \code{netest} (setting the
#' \code{edapprox} argument to \code{FALSE}) require only a call to
#' \code{simulate.stergm}. Control parameters for those simulations may be set
#' using \code{set.control.stergm} in \code{netdx}. The parameters should be
#' input through the \code{control.simulate.stergm()} function, with the
#' available parameters listed in the \code{\link{control.simulate.stergm}} help
#' page in the \code{tergm} package.
#'
#' Models fit with the ERGM method with the edges dissolution approximation
#' (setting \code{edapprox} to \code{TRUE}) require a call first to
#' \code{simulate.ergm} for simulating an initial network, and second to
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
#'   set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6)
#' )
#' dx2
#' plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#' plot(dx2, type = "duration")
#' plot(dx2, type = "dissolution", qnts.col = "orange2")
#' plot(dx2, type = "dissolution", method = "b", col = "bisque")
#' }
#'
netdx <- function(x, nsims = 1, dynamic = TRUE, nsteps,
                  nwstats.formula = "formation", set.control.ergm,
                  set.control.stergm, sequential = TRUE, keep.tedgelist = FALSE,
                  keep.tnetwork = FALSE, verbose = TRUE, ncores = 1,
                  skip.dissolution = FALSE) {
  if (class(x) != "netest") {
    stop("x must be an object of class netest", call. = FALSE)
  }

  ncores <- ifelse(nsims == 1, 1, min(parallel::detectCores(), ncores))

  fit <- x$fit
  formation <- x$formation
  coef.form <- x$coef.form
  dissolution <- x$coef.diss$dissolution
  coef.diss <- x$coef.diss
  constraints <- x$constraints
  target.stats <- x$target.stats
  edapprox <- x$edapprox
  if (edapprox == TRUE) {
    nw <- x$fit$newnetwork
  } else {
    nw <- x$fit$network
  }

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

  if (edapprox == FALSE) {
    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.stergm()
    }

    if (nsims == 1 || ncores == 1) {
      diag.sim <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n |")
      }
      for (i in seq_len(nsims)) {
        diag.sim[[i]] <- simulate(fit,
          time.slices = nsteps,
          monitor = nwstats.formula,
          nsim = 1,
          control = set.control.stergm
        )
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

      diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
        simulate(fit,
          time.slices = nsteps,
          monitor = nwstats.formula,
          nsim = 1,
          control = set.control.stergm
        )
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
        for (i in seq_len(nsims)) {
          fit.sim <- simulate(fit,
            basis = fit$newnetwork,
            control = set.control.ergm, dynamic = FALSE
          )
          diag.sim[[i]] <- simulate(fit.sim,
            formation = formation,
            dissolution = dissolution,
            coef.form = coef.form,
            coef.diss = coef.diss$coef.crude,
            time.slices = nsteps,
            constraints = constraints,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.stergm
          )
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

        diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
          fit.sim <- simulate(fit,
            basis = fit$newnetwork,
            control = set.control.ergm, dynamic = FALSE
          )
          simulate(fit.sim,
            formation = formation,
            dissolution = dissolution,
            coef.form = coef.form,
            coef.diss = coef.diss$coef.crude,
            time.slices = nsteps,
            constraints = constraints,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.stergm
          )
        }
      }
    }
    if (dynamic == FALSE) {
      diag.sim <- simulate(fit,
        nsim = nsims,
        output = "stats",
        control = set.control.ergm,
        sequential = sequential,
        monitor = nwstats.formula,
        dynamic = FALSE
      )
    }
  } # end edapprox = TRUE condition

  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  if (dynamic == TRUE) {
    stats <- lapply(diag.sim, function(x) attributes(x)$stats)
    merged.stats <- Reduce(function(a, x) rbind(a, x), stats, init = c())
  } else {
    stats <- list(diag.sim[, !duplicated(colnames(diag.sim)), drop = FALSE])
    merged.stats <- diag.sim[, !duplicated(colnames(diag.sim)), drop = FALSE]
  }

  ts.attr.names <- x$target.stats.names
  if (length(ts.attr.names) != length(target.stats)) {
    target.stats <- target.stats[which(target.stats > 0)]
  }
  ts.out <- data.frame(
    names = ts.attr.names,
    targets = target.stats
  )

  ## Calculate mean/sd from merged stats
  stats.table.formation <- make_formation_table(merged.stats, ts.out)

  if (skip.dissolution == FALSE) {
    if (dynamic == TRUE) {
      sim.df <- lapply(diag.sim, as.data.frame)
      dissolution.stats <- make_dissolution_stats(
        sim.df,
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
      out$stats.table.dissolution <- dissolution.stats$stats.table.dissolution
      out$pages <- dissolution.stats$pages
      out$pages_imptd <- dissolution.stats$pages_imptd
      out$prop.diss <- dissolution.stats$prop.diss
    }
    if (keep.tedgelist == TRUE) {
      out$tedgelist <- sim.df
    }
    if (keep.tnetwork == TRUE) {
      out$network <- diag.sim
    }
  }

  class(out) <- "netdx"
  return(out)
}

#' @title Calculate the Formation Statistics of a Network
#'
#' @param merged.stats A matrix of \code{nsims * nsteps} rows, with a column for
#'   each of the formation targets.
#' @param targets A \code{data.frame} of the formation targets with two columns:
#'   "names" and "targets".
#'
#' @return A \code{data.frame} of the formation statistics.
#' @keywords internal
make_formation_table <- function(merged.stats, targets) {

  ## Calculate mean/sd from merged stats
  stats.means <- colMeans(merged.stats)
  stats.sd <- apply(merged.stats, 2, sd)

  stats.table <- data.frame(
    sorder = seq_along(names(stats.means)),
    names = names(stats.means),
    stats.means,
    stats.sd
  )

  ## Create stats.formation table for output
  stats.table <- merge(targets, stats.table, all = TRUE)
  stats.table <- stats.table[order(stats.table[["sorder"]]), , drop = FALSE]
  rownames(stats.table) <- stats.table$names

  stats.table$reldiff <- (stats.table$stats.means - stats.table$targets) /
    stats.table$targets * 100
  stats.table.formation <- stats.table[, c(2, 4, 6, 5)]
  colnames(stats.table.formation) <- c(
    "Target",
    "Sim Mean",
    "Pct Diff",
    "Sim SD"
  )

  return(stats.table.formation)
}

#' @title Calculate the Dissolution Statistics of a Network
#'
#' @param sim.df A list of network objects (one per simulation).
#' @param coef.diss The \code{coef.diss} element of \code{nwparam}.
#' @param nsteps The number of simulated steps.
#' @param verbose A verbosity toggle (default = TRUE).
#'
#' @return A \code{list} of dissolution statistics.
#' @keywords internal
make_dissolution_stats <- function(sim.df, coef.diss, nsteps, verbose = TRUE) {
  if (verbose == TRUE) {
    cat("\n- Calculating duration statistics")
  }

  nsims <- length(sim.df)

  # Calculate mean partnership age from edgelist and ensure that a value is
  # provided for each timestep
  pages <- lapply(sim.df, function(x) {
    meanage <- edgelist_meanage(el = x)
    l <- nsteps - length(meanage)
    if (l > 0) {
      meanage <- c(meanage, rep(NA, l))
    }
    return(meanage)
  })

  # TODO: imputation currently averaged for heterogeneous models
  if (coef.diss$model.type == "hetero") {
    coef_dur <- mean(coef.diss$duration)
  } else {
    coef_dur <- coef.diss$duration
  }
  pages_imptd <- coef_dur^2 * dgeom(2:(nsteps + 1), 1 / coef_dur)

  ## Dissolution calculations
  if (verbose == TRUE) {
    cat("\n- Calculating dissolution statistics")
  }

  ## Create a list of dissolution proportions (i.e. dissolutions/edges)
  prop.diss <- lapply(sim.df, function(d) {
    vapply(seq_len(nsteps), function(x) {
      sum(d$terminus == x) / sum(d$onset < x & d$terminus >= x)
    }, 0)
  })


  if (verbose == TRUE) {
    cat("\n ")
  }

  # Create dissolution tables
  duration.obs <- matrix(unlist(pages), nrow = nsteps)
  duration.imputed <- duration.obs + pages_imptd
  duration.mean.by.sim <- colMeans(duration.imputed)
  duration.mean <- mean(duration.mean.by.sim, na.rm = TRUE)
  if (nsims > 1) {
    duration.sd <- sd(duration.mean.by.sim, na.rm = TRUE)
  } else {
    duration.sd <- NA
  }

  duration.expected <- exp(coef.diss$coef.crude[1]) + 1
  duration.pctdiff <- (duration.mean - duration.expected) /
    duration.expected * 100

  dissolution.mean <- mean(unlist(prop.diss), na.rm = TRUE)

  if (nsims > 1) {
    dissolution.sd <- sd(sapply(prop.diss, mean, na.rm = TRUE))
  } else {
    dissolution.sd <- NA
  }

  dissolution.expected <- 1 / (exp(coef.diss$coef.crude[1]) + 1)
  dissolution.pctdiff <- (dissolution.mean - dissolution.expected) /
    dissolution.expected * 100

  stats.table.dissolution <- data.frame(
    Targets = c(duration.expected, dissolution.expected),
    Sim_Means = c(duration.mean, dissolution.mean),
    Pct_Diff = c(duration.pctdiff, dissolution.pctdiff),
    Sim_SD = c(duration.sd, dissolution.sd)
  )
  colnames(stats.table.dissolution) <- c(
    "Target", "Sim Mean", "Pct Diff", "Sim SD"
  )
  rownames(stats.table.dissolution) <- c("Edge Duration", "Pct Edges Diss")

  return(
    list(
      "stats.table.dissolution" = stats.table.dissolution,
      "pages" = pages,
      "pages_imptd" = pages_imptd,
      "prop.diss" = prop.diss
    )
  )
}
