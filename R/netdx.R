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
#' plotting the \code{netdx} object.
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
                  nwstats.formula = "formation", set.control.ergm,
                  set.control.stergm, set.control.tergm,
                  sequential = TRUE, keep.tedgelist = FALSE,
                  keep.tnetwork = FALSE, verbose = TRUE, ncores = 1,
                  skip.dissolution = FALSE) {

  if (!inherits(x, "netest")) {
    stop("x must be an object of class netest", call. = FALSE)
  }

  if (!missing(set.control.stergm)) {
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
    if (missing(set.control.stergm) && missing(set.control.tergm)) {
      set.control.tergm <- control.simulate.formula.tergm()
    }

    if (nsims == 1 || ncores == 1) {
      diag.sim <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n |")
      }
      for (i in seq_len(nsims)) {
        if (!missing(set.control.stergm)) {
          diag.sim[[i]] <- simulate(x$newnetwork,
            formation = x$formation,
            dissolution = x$dissolution,
            coef.form = x$coef.form,
            coef.diss = x$coef.diss$coef.crude,
            constraints = x$constraints,
            time.slices = nsteps,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.stergm
          )
        } else {
          diag.sim[[i]] <- simulate(x$formula,
            coef = c(x$coef.form, x$coef.diss$coef.crude),
            constraints = x$constraints,
            basis = x$newnetwork,
            time.slices = nsteps,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.tergm,
            dynamic = TRUE
          )
        }
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

      if (!missing(set.control.stergm)) {
        diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
          simulate(x$newnetwork,
            formation = x$formation,
            dissolution = x$dissolution,
            coef.form = x$coef.form,
            coef.diss = x$coef.diss$coef.crude,
            constraints = x$constraints,
            time.slices = nsteps,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.stergm
          )
        }
      } else {
        diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
          simulate(x$formula,
            coef = c(x$coef.form, x$coef.diss$coef.crude),
            constraints = x$constraints,
            basis = x$newnetwork,
            time.slices = nsteps,
            monitor = nwstats.formula,
            nsim = 1,
            control = set.control.tergm,
            dynamic = TRUE
          )
        }
      }
    }
  }

  if (edapprox == TRUE) {
    if (missing(set.control.ergm)) {
      set.control.ergm <- control.simulate.formula()
    }
    if (missing(set.control.stergm) && missing(set.control.tergm)) {
      set.control.tergm <- control.simulate.formula.tergm()
    }

    if (dynamic == TRUE) {
      if (nsims == 1 || ncores == 1) {
        diag.sim <- list()
        if (verbose == TRUE & nsims > 1) {
          cat("\n  |")
        }
        for (i in seq_len(nsims)) {
          fit.sim <- simulate(x$formula,
            coef = x$coef.form.crude,
            basis = x$newnetwork,
            constraints = constraints,
            control = set.control.ergm, dynamic = FALSE
          )
          if (!missing(set.control.stergm)) {
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
          } else {
            diag.sim[[i]] <- simulate(fit.sim ~ Form(formation) +
                                                Persist(dissolution),
              coef = c(coef.form, coef.diss$coef.crude),
              time.slices = nsteps,
              constraints = constraints,
              monitor = nwstats.formula,
              nsim = 1,
              control = set.control.tergm,
              dynamic = TRUE
            )
          }
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

        if (!missing(set.control.stergm)) {
          diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
            fit.sim <- simulate(x$formula,
              coef = x$coef.form.crude,
              basis = x$newnetwork,
              constraints = x$constraints,
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
        } else {
          diag.sim <- foreach(i = seq_len(nsims)) %dopar% {
            fit.sim <- simulate(x$formula,
              coef = x$coef.form.crude,
              basis = x$newnetwork,
              constraints = x$constraints,
              control = set.control.ergm, dynamic = FALSE
            )

            simulate(fit.sim ~ Form(formation) +
                               Persist(dissolution),
              coef = c(coef.form, coef.diss$coef.crude),
              time.slices = nsteps,
              constraints = constraints,
              monitor = nwstats.formula,
              nsim = 1,
              control = set.control.tergm,
              dynamic = TRUE
            )
          }
        }
      }
    }
    if (dynamic == FALSE) {
      diag.sim <- simulate(x$formula,
        coef = x$coef.form.crude,
        basis = x$newnetwork,
        constraints = x$constraints,
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
      out$stats.table.duration <- dissolution.stats$stats.table.duration
      out$stats.table.dissolution <- dissolution.stats$stats.table.dissolution
      out$pages <- dissolution.stats$pages
      out$pages_imptd <- dissolution.stats$pages_imptd
      out$prop.diss <- dissolution.stats$prop.diss
    }
    if (keep.tedgelist == TRUE) {
      out$tedgelist <- lapply(diag.sim, as.data.frame)
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
#' @param diag.sim A list of network objects (one per simulation).
#' @param coef.diss The \code{coef.diss} element of \code{nwparam}.
#' @param nsteps The number of simulated steps.
#' @param verbose A verbosity toggle (default = TRUE).
#'
#' @return a \code{list} of dissolution statistics
#' @keywords internal
make_dissolution_stats <- function(diag.sim, coef.diss,
                                   nsteps, verbose = TRUE) {
  if (verbose == TRUE) {
    cat("\n- Calculating duration statistics")
  }

  sim.df <- lapply(diag.sim, as.data.frame)
  nsims <- length(sim.df)
  dissolution <- coef.diss$dissolution
  diss_term <- if (coef.diss$diss.model.type == "edgesonly") {
    NULL
  } else {
    coef.diss$diss.model.type
  }

  # Check form of dissolution formula and extract attribute name, if any
  # Code adapted from dissolution_coefs and diss_check
  # TODO: consider moving this into dissolution_coefs, and saving the attribute
  # name there as an additional element in coef.diss.  (It now saves the term
  # as "diss.model.type") That would allow us to need to replicate some of this
  # code here. Alternative plan: consider having diss_check return values for
  # both the dissolution model term and model attribute, which then get saved in
  # the netest object, allowing *both* this code and the similar piece in
  # dissolution_coefs to both be removed.
  diss.terms <- strsplit(as.character(dissolution)[2], "[+]")[[1]]
  diss.terms <- gsub("\\s", "", diss.terms)
  offpos.d <- grep("offset(", diss.terms, fixed = TRUE)
  diss.terms[offpos.d] <- substr(diss.terms[offpos.d], nchar("offset(") + 1,
                                 nchar(diss.terms[offpos.d]) - 1)
  argpos.d <- regexpr("\\(", diss.terms)
  diss.terms <- vapply(regmatches(diss.terms, argpos.d, invert = TRUE),
                       function(x) {
                         if (length(x) < 2) {
                           x <- c(x, "")
                         } else {
                           x[2] <- substr(x[2], 1, nchar(x[2]) - 1)
                         }
                         x
                       },
                       c(term = "", args = ""))
  diss.terms <- gsub("\"", "", diss.terms)
  if (ncol(diss.terms) == 2) {
    if (grepl(",", diss.terms[2, 2]) == TRUE) {
      if (grepl("=", diss.terms[2, 2]) == TRUE) {
        diss_arg <- as.logical(strsplit(diss.terms[2, 2], "=")$args[2])
      } else {
        (stop("Dissolution model does not conform to expected format."))
      }
      diss.terms[2, 2] <- strsplit(diss.terms[2, 2], ",")$args[1]
    }    # Used to remove diff argument if present, regardless of its value
    if (diss_term == "nodematch" & !exists("diss_arg"))
      diss_arg <- FALSE
    diss_attr_name <- diss.terms[2, 2]
  } else {
    diss_attr_name <- NULL
  }

  # Calculate mean partnership age from edgelist
  pages <- sapply(seq_along(sim.df), function(x) {
                      meanage <- edgelist_meanage(
                        el = sim.df[[x]],
                        diss_term = diss_term,
                        diss_attr = if (is.null(diss_term)) NULL else
                            get_vertex_attribute(diag.sim[[x]], diss_attr_name),
                        diss_arg = if (!exists("diss_arg")) NULL else
                          diss_arg)
                      l <- nsteps - nrow(meanage)
                      if (l > 0) {
                        meanage <- rbind(meanage,
                                         matrix(rep(NA, l *
                                                      ncol(meanage)),
                                                nrow = l))
                      }
                      return(meanage)
              }, simplify = "array")

  # when 1 time step and 1 stat (edgesonly or nodefactor)
  if (is.vector(pages)) pages <- array(pages, dim = c(1, 1, nsims))

  # calculate expected time prior to simulation
  # TODO: remove nodefactor in future release
  if (coef.diss$diss.model.type == "nodefactor") {
    coef_dur <- mean(coef.diss$duration)
  } else {
    coef_dur <- coef.diss$duration
  }

  pages_imptd <- sapply(seq_along(coef_dur), function(x)
    coef_dur[x]^2 * dgeom(2:(nsteps + 1), 1 / coef_dur[x]))
  if (nsteps == 1)
    pages_imptd <- matrix(pages_imptd, nrow = 1)

  ## Dissolution calculations
  if (verbose == TRUE) {
    cat("\n- Calculating dissolution statistics")
  }

  if (is.null(diss_term) || diss_term == "nodefactor") {
    # TODO: remove nodefactor in future release
    if (!is.null(diss_term) && diss_term == "nodefactor")
      warning("Support for dissolution models containing a nodefactor term is
              deprecated, and will be removed in a future release.",
              call. = FALSE)
    prop.diss <- sapply(seq_along(sim.df), function(d) {
      matrix(sapply(seq_len(nsteps), function(x) {
        sum(sim.df[[d]]$terminus == x) / sum(sim.df[[d]]$onset < x &
                                               sim.df[[d]]$terminus >= x)
      }), ncol = 1)}, simplify = "array")
    if (nsteps == 1) prop.diss <- array(prop.diss, dim = c(1, 1, nsims))
  } else {
    if (diss_term == "nodematch") {
      # assumes same attribute values across sims -- appropriate for netdx
      attribute <- get_vertex_attribute(diag.sim[[1]], diss_attr_name)

      if (diss_arg == TRUE) {
          attrvalues <- sort(unique(attribute))
          prop.diss <- sapply(seq_along(sim.df), function(d) {
            t(sapply(seq_len(nsteps), function(x) {
              heterogs <- sum(sim.df[[d]]$terminus == x &
                                attribute[sim.df[[d]]$head] !=
                                attribute[sim.df[[d]]$tail]) /
                          sum(sim.df[[d]]$onset < x &
                                sim.df[[d]]$terminus >= x &
                                attribute[sim.df[[d]]$head] !=
                                attribute[sim.df[[d]]$tail])
              homogs <- sapply(seq_along(attrvalues), function(y)
                          sum(sim.df[[d]]$terminus == x &
                                attribute[sim.df[[d]]$head] == attrvalues[y] &
                                attribute[sim.df[[d]]$tail] == attrvalues[y]) /
                          sum(sim.df[[d]]$onset<x & sim.df[[d]]$terminus >= x &
                                attribute[sim.df[[d]]$head] == attrvalues[y] &
                                attribute[sim.df[[d]]$tail] == attrvalues[y])
              )
              return(c(heterogs, homogs))
            }))
          }, simplify = "array")
      } else {
          prop.diss <- sapply(seq_along(sim.df), function(d) {
            t(sapply(seq_len(nsteps), function(x) {
              c(sum(sim.df[[d]]$terminus == x &
                      attribute[sim.df[[d]]$head] !=
                      attribute[sim.df[[d]]$tail]) /
                  sum(sim.df[[d]]$onset<x & sim.df[[d]]$terminus >= x &
                        attribute[sim.df[[d]]$head] !=
                        attribute[sim.df[[d]]$tail]),
                sum(sim.df[[d]]$terminus == x &
                      attribute[sim.df[[d]]$head] ==
                      attribute[sim.df[[d]]$tail]) /
                  sum(sim.df[[d]]$onset < x & sim.df[[d]]$terminus >= x &
                        attribute[sim.df[[d]]$head] ==
                        attribute[sim.df[[d]]$tail]))
            }))
          }, simplify = "array")
      }
    } else {
      if (diss_term == "nodemix") {
        # assumes same attribute values across sims -- appropriate for netdx
        attribute <- get_vertex_attribute(diag.sim[[1]], diss_attr_name)
        attrvalues <- sort(unique(attribute))
        n.attrvalues <- length(attrvalues)
        n.attrcombos <- n.attrvalues*(n.attrvalues+1)/2
        indices2.grid <- expand.grid(row = 1:n.attrvalues, col = 1:n.attrvalues)
        rowleqcol <- indices2.grid$row <= indices2.grid$col #assumes undirected
        indices2.grid <- indices2.grid[rowleqcol, ]
        prop.diss <- sapply(seq_along(sim.df), function(d) {
          t(sapply(seq_len(nsteps), function(x) {
            sapply(seq_len(nrow(indices2.grid)), function(y) {
                ingroup <- (attribute[sim.df[[d]]$head] ==
                              attribute[indices2.grid$row[y]] &
                            attribute[sim.df[[d]]$tail] ==
                              attribute[indices2.grid$col[y]]) |
                           (attribute[sim.df[[d]]$head] ==
                              attribute[indices2.grid$col[y]] &
                            attribute[sim.df[[d]]$tail] ==
                              attribute[indices2.grid$row[y]])
                sum(sim.df[[d]]$terminus==x & ingroup) /
                  sum(sim.df[[d]]$onset<x & sim.df[[d]]$terminus>=x & ingroup)
              })
          }))
        }, simplify = "array")
      } else {
        prop.diss <- NULL
      }
    }
  }
  if (verbose == TRUE) {
    cat("\n ")
  }

  # Create duration table
  duration.imputed <- simplify2array(lapply(1:nsims,
                              function(x)pages[,,x]+pages_imptd))
  if (is.vector(duration.imputed) & nsteps == 1) {
    duration.imputed <- array(duration.imputed, dim = c(1, 1, nsims))
  }
  duration.mean.by.sim <- apply(duration.imputed, 2:3, mean)
  duration.mean <- rowMeans(duration.mean.by.sim, na.rm = TRUE)

  if (nsims > 1) {
    duration.sd <- apply(duration.mean.by.sim, 1, sd, na.rm = TRUE)
  } else {
    duration.sd <- NA
  }
  duration.expected <- coef_dur
  duration.pctdiff <- (duration.mean - duration.expected) /
    duration.expected * 100

  stats.table.duration <- data.frame(
    Targets = c(duration.expected),
    Sim_Means = c(duration.mean),
    Pct_Diff = c(duration.pctdiff),
    Sim_SD = c(duration.sd)
  )

  # Create dissolution table
  dissolution.mean.by.sim <- apply(prop.diss, 2:3, mean)
  dissolution.mean <- rowMeans(dissolution.mean.by.sim, na.rm = TRUE)

  if (nsims > 1) {
    dissolution.sd <- apply(dissolution.mean.by.sim, 1, sd, na.rm = TRUE)
  } else {
    dissolution.sd <- NA
  }

  dissolution.expected <- 1/coef_dur
  dissolution.pctdiff <- (dissolution.mean - dissolution.expected) /
    dissolution.expected * 100

  stats.table.dissolution <- data.frame(
    Targets = c(dissolution.expected),
    Sim_Means = c(dissolution.mean),
    Pct_Diff = c(dissolution.pctdiff),
    Sim_SD = c(dissolution.sd)
  )

  # Set column names for both duration and dissolution tables
  colnames(stats.table.duration) <- colnames(stats.table.dissolution) <- c(
    "Target", "Sim Mean", "Pct Diff", "Sim SD"
  )

  # Set row names for both duration and dissolution tables
  if (is.null(diss_term) || diss_term == "nodefactor") {
  # TODO: remove nodefactor in future release
    rownames(stats.table.duration) <- rownames(stats.table.dissolution) <-
        c("edges")
  } else {
    if (diss_term=="nodematch") {
      if(diss_arg==TRUE) {
        rownames(stats.table.duration) <-
          rownames(stats.table.dissolution) <-
          c(paste("match",diss_attr_name, "FALSE", sep="."),
            sapply(seq_along(attrvalues), function(z) paste("match",
                                                            diss_attr_name,
                                                            "TRUE",
                                                            attrvalues[z],
                                                            sep = '.')))
      } else {
        rownames(stats.table.duration) <-
          rownames(stats.table.dissolution) <-
          c(paste("match", diss_attr_name, "FALSE", sep = "."),
            paste("match", diss_attr_name, "TRUE ", sep = "."))
      }
    } else {
      if (diss_term=="nodemix") {
        # assumes same attribute values across sims -- appropriate for netdx
        attribute <- get_vertex_attribute(diag.sim[[1]], diss_attr_name)
        attrvalues <- sort(unique(attribute))
        n.attrvalues <- length(attrvalues)
        n.attrcombos <- n.attrvalues*(n.attrvalues+1)/2
        indices2.grid <- expand.grid(row = 1:n.attrvalues, col = 1:n.attrvalues)
        rowleqcol <- indices2.grid$row <= indices2.grid$col #assumes undirected
        uun <- as.vector(outer(attrvalues, attrvalues, paste, sep = "."))
        uun <- uun[rowleqcol]
        rownames(stats.table.duration) <-
          rownames(stats.table.dissolution) <-
            paste("mix", diss_attr_name, uun, sep = ".")
      }
    }
  }

  # Construct return list
  return(
    list(
      "stats.table.duration" = stats.table.duration,
      "stats.table.dissolution" = stats.table.dissolution,
      "pages" = pages,
      "pages_imptd" = pages_imptd,
      "prop.diss" = prop.diss
    )
  )
}
