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
#' \code{ergm}'s \code{simulate_formula.network} for simulating an initial network, and second to
#' \code{tergm}'s \code{simulate_formula.network} for simulating that static network forward through
#' time. Control parameters may be set for both processes in \code{netdx}.
#' For the first, the parameters should be input through the
#' \code{control.simulate.formula()} function, with the available parameters listed
#' in the \code{\link[ergm:control.simulate.formula]{control.simulate.formula}} help
#' page in the \code{ergm} package. For the second, parameters should be input
#' through the \code{control.simulate.formula.tergm()} function, with the available
#' parameters listed in the \code{\link{control.simulate.formula.tergm}} help page in
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
#'   set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6)
#' )
#' dx2
#' plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#' plot(dx2, type = "duration")
#' plot(dx2, type = "dissolution", qnts.col = "orange2")
#' plot(dx2, type = "dissolution", method = "b", col = "bisque")
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
  if (STERGM) {
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

  if (!dynamic && !edapprox) {
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

  dosims <- function() {
    if (edapprox == TRUE) {
      init <- simulate(x$formula,
                       coef = x$coef.form.crude,
                       basis = x$newnetwork,
                       constraints = x$constraints,
                       control = set.control.ergm, 
                       dynamic = FALSE,
                       nsim = if(dynamic) 1L else nsims,
                       output = if(dynamic) "network" else "stats",
                       sequential = if(dynamic) FALSE else sequential,
                       monitor = if(dynamic) NULL else nwstats.formula)
    } else {
      init <- x$newnetwork
    }

    if(!dynamic) {
      return(list(stats = init))
    }

    if (STERGM) {
      diag.sim <- simulate(init,
                           formation = x$formation,
                           dissolution = x$coef.diss$dissolution,
                           coef.form = x$coef.form,
                           coef.diss = x$coef.diss$coef.crude,
                           constraints = x$constraints,
                           time.slices = nsteps,
                           monitor = nwstats.formula,
                           time.start = 0,
                           nsim = 1,
                           output = "changes",
                           control = set.control.stergm)
    } else {
      diag.sim <- simulate(init ~ Form(x$formation) + Persist(x$coef.diss$dissolution),
                           coef = c(x$coef.form, x$coef.diss$coef.crude),
                           constraints = x$constraints,
                           time.slices = nsteps,
                           monitor = nwstats.formula,
                           time.start = 0,
                           nsim = 1,
                           output = "changes",
                           control = set.control.tergm,
                           dynamic = TRUE)
    }
    
    stats <- attr(diag.sim, "stats")
    
    changes <- rbind(cbind(0L, as.edgelist(init), 1L),
                     diag.sim)

    if(skip.dissolution) {
      return(list(stats = stats, changes = changes))
    }
    
    tedgelist <- changes_to_tedgelist(changes, nsteps)
    
    diss_stats <- tedgelist_to_diss_stats(tedgelist, x$coef.diss, nsteps, init)
    
    return(c(diss_stats, list(stats = stats, changes = changes)))
  }
  
  if (!dynamic || nsims == 1) {
    diag.sim <- list(dosims())
  } else if (ncores == 1) {
    diag.sim <- list()
    if (verbose == TRUE) {
      cat("\n |")
    }
    for (i in seq_len(nsims)) {    
      diag.sim[[i]] <- dosims()      
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
    diag.sim <- foreach(i = seq_len(nsims)) %dopar% dosims()  
  }
  
  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }

  ## List for stats for each simulation
  stats <- lapply(diag.sim, function(x) { y <- x$stats; y[,!duplicated(colnames(y)),drop=FALSE] })
  merged.stats <- do.call(rbind, stats)

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
      out$tedgelist <- lapply(diag.sim, `[[`, "tedgelist")
    }
    if (keep.tnetwork == TRUE) {
      nw[,] <- FALSE
      out$network <- lapply(diag.sim, function(sim) networkDynamic(base.net = nw, edge.toggles = sim$changes[,1:3,drop=FALSE]))
    }
  }

  class(out) <- "netdx"
  return(out)
}

## preconditions: changes includes on-toggles of initial edges at time.start, and
##                time.start == 0L
changes_to_tedgelist <- function(changes, nsteps) {
  if(NROW(changes) == 0L) {
    ## handle case of zero changes...
    tedgelist <- data.frame(onset = integer(0),
                            terminus = integer(0),
                            tail = integer(0),
                            head = integer(0),
                            onset.censored = logical(0),
                            terminus.censored = logical(0),
                            duration = integer(0))
  } else {    
    changes <- changes[order(changes[,2L], changes[,3L], changes[,1L]),,drop=FALSE]
    
    is_onset <- logical(NROW(changes))
    has_terminus <- logical(NROW(changes))
    
    samenext <- changes[-NROW(changes),2L] == changes[-1L,2L] & changes[-NROW(changes),3L] == changes[-1L,3L]
    
    is_onset[1L] <- TRUE
    i <- 1L
    while(i < NROW(changes)) {
      is_onset[i] <- TRUE
      has_terminus[i] <- samenext[i]
      i <- i + 1L + as.integer(samenext[i])
    }
    is_onset[NROW(changes)] <- i == NROW(changes)
    
    onset.times <- changes[is_onset,1L]
    terminus.times <- rep(nsteps + 1L, length(onset.times))
    terminus.times[has_terminus[is_onset]] <- changes[which(has_terminus) + 1L,1L]
    
    duration <- terminus.times - onset.times
    onset.censored <- onset.times == 0L
    terminus.censored <- !has_terminus[is_onset]
    
    tails <- changes[is_onset,2L]
    heads <- changes[is_onset,3L]
    
    tedgelist <- data.frame(onset = onset.times,
                            terminus = terminus.times,
                            tail = tails,
                            head = heads,
                            onset.censored = onset.censored,
                            terminus.censored = terminus.censored,
                            duration = duration)
  }
  
  tedgelist
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
make_dissolution_stats <- function(diag.sim, coef.diss, nsteps, verbose = TRUE) {
  if (verbose == TRUE) {
    cat("\n- Calculating duration statistics")
  }
  
  ## exclude nodefactor from hetergeneous dissolution calculation
  if(coef.diss$diss.model.type == "nodefactor") {
    durs <- mean(coef.diss$duration)
  } else {
    durs <- coef.diss$duration  
  }
  
  pages <- array(unlist(lapply(diag.sim, `[[`, "meanage")),
                 dim = c(nsteps, length(durs), length(diag.sim)))

  pages_imptd <- array(unlist(lapply(diag.sim, `[[`, "meanageimputed")),
                 dim = c(nsteps, length(durs), length(diag.sim)))

  propdiss <- array(unlist(lapply(diag.sim, `[[`, "propdiss")),
                    dim = c(nsteps, length(durs), length(diag.sim)))

  combinedmeanageimputed <- do.call(rbind, lapply(diag.sim, `[[`, "meanmeanageimputed"))
  meanagesimputed <- colMeans(combinedmeanageimputed)
  meanagesd <- apply(combinedmeanageimputed, 2, sd)

  combinedpropdiss <- do.call(rbind, lapply(diag.sim, `[[`, "propdiss"))
  meanpropdiss <- colMeans(combinedpropdiss)
  propdisssd <- apply(combinedpropdiss, 2, sd)
  
  stats.table.duration <- data.frame("Target" = durs,
                                     "Sim Mean" = meanagesimputed,
                                     "Pct Diff" = 100*(meanagesimputed - durs)/durs,
                                     "Sim SD" = meanagesd)
  colnames(stats.table.duration) <- c("Target", "Sim Mean", "Pct Diff", "Sim SD")
  
  stats.table.dissolution <- data.frame("Target" = 1/durs,
                                        "Sim Mean" = meanpropdiss,
                                        "Pct Diff" = 100*(meanpropdiss - 1/durs)/(1/durs),
                                        "Sim SD" = propdisssd)
  colnames(stats.table.dissolution) <- c("Target", "Sim Mean", "Pct Diff", "Sim SD")

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

tedgelist_to_diss_stats <- function(tedgelist, coef.diss, nsteps, nw) {
  if(coef.diss$diss.model.type == "nodefactor") {
    diss_formula <- ~offset(edges)
    durs <- mean(coef.diss$duration)
  } else {
    diss_formula <- coef.diss$dissolution
    durs <- coef.diss$duration
  }

  diss_formula <- nonsimp_update.formula(diss_formula, nw ~ ., from.new = "nw")
  
  togglemat <- cbind(seq_len(NROW(tedgelist)), tedgelist[,3L], tedgelist[,4L])

  changestats <- tergm.godfather(diss_formula,
                                 toggles = togglemat,
                                 stats.start = TRUE)
                                     
  changemat <- changestats[-1L,,drop=FALSE] != changestats[-NROW(changestats),,drop=FALSE]

  dyad_types <- integer(NROW(tedgelist))
  for(i in seq_len(NROW(tedgelist))) {
    dyad_types[i] <- max(which(changemat[i,]))
  }

  imputed_corrections <- integer(NROW(tedgelist))
  for(i in seq_len(NROW(tedgelist))) {
    if(tedgelist[i,"onset.censored"]) {
      imputed_corrections[i] <- rgeom(1L, 1/durs[dyad_types[i]])
    }
  }
  
  tedgelist <- cbind(tedgelist,
                     imputed_corrections = imputed_corrections,
                     dyad_types = dyad_types)
    
  edgecounts <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(edgecounts) <- colnames(changestats)
  
  edgeages <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(edgeages) <- colnames(changestats)
  
  edgeagesimputed <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(edgeagesimputed) <- colnames(changestats)
  
  edgediss <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(edgediss) <- colnames(changestats)
  
  nform <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(nform) <- colnames(changestats)
  ndiss <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(ndiss) <- colnames(changestats)
  nedges <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(nedges) <- colnames(changestats)
  ageform <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(ageform) <- colnames(changestats)
  agediss <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(agediss) <- colnames(changestats)
  ageform_imputed <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(ageform_imputed) <- colnames(changestats)
  agediss_imputed <- matrix(0L, nrow = nsteps, ncol = NCOL(changestats), byrow = TRUE)
  colnames(agediss_imputed) <- colnames(changestats)

  init.edges <- integer(NCOL(changestats))
  init.age <- integer(NCOL(changestats))
  init.age_imputed <- integer(NCOL(changestats))
  
  onset.censored <- tedgelist$onset.censored
  terminus.censored <- tedgelist$terminus.censored
  onset <- tedgelist$onset
  terminus <- tedgelist$terminus
  dyad_types <- dyad_types
  imputed_corrections <- imputed_corrections
  
  for(i in seq_len(NROW(tedgelist))) {
    type <- dyad_types[i]
    ot <- onset[i]
    tt <- terminus[i]
    ic <- imputed_corrections[i]
    
    if(!onset.censored[i]) {
      nform[ot,type] <- nform[ot,type] + 1L
      ageform[ot,type] <- ageform[ot,type] + 1L
      ageform_imputed[ot,type] <- ageform_imputed[ot,type] + 1L
    } else {
      init.edges[type] <- init.edges[type] + 1L
      init.age[type] <- init.age[type] + 1L
      init.age_imputed[type] <- init.age_imputed[type] + 1L + ic
    }
    
    if(!terminus.censored[i]) {
      ndiss[tt,type] <- ndiss[tt,type] + 1L
      agediss[tt,type] <- agediss[tt,type] + tt - ot + 1L
      agediss_imputed[tt,type] <- agediss_imputed[tt,type] + tt - ot + 1L + ic
    }
  }
  
  edgediss <- ndiss
  
  for(j in seq_len(NCOL(edgecounts))) {
    edgecounts[,j] <- init.edges[j] + cumsum(nform[,j] - ndiss[,j])
    edgeages[,j] <- init.age[j] + cumsum(ageform[,j] - agediss[,j] + c(init.edges[j], edgecounts[-NROW(edgecounts),j]))
    edgeagesimputed[,j] <- init.age_imputed[j] + cumsum(ageform_imputed[,j] - agediss_imputed[,j] + c(init.edges[j], edgecounts[-NROW(edgecounts),j]))
  }
  
  init_stats <- changestats[1L,]
  if(length(init_stats) > 1L) {
    init_stats[1] <- init_stats[1] - sum(init_stats[-1])
  }
  edgecounts <- rbind(init_stats, edgecounts)
  
  meanage <- edgeages/edgecounts[-1L,]
  meanage[is.nan(meanage)] <- 0
  meanageimputed <- edgeagesimputed/edgecounts[-1L,]
  meanageimputed[is.nan(meanageimputed)] <- 0
  
  propdiss <- edgediss/edgecounts[-NROW(edgecounts),]
  propdiss[is.nan(propdiss)] <- 0
  
  meanmeanageimputed <- colMeans(meanageimputed)
  
  return(list(tedgelist = tedgelist, 
              dyad_types = dyad_types,
              edgecounts = edgecounts,
              edgeages = edgeages,
              edgeagesimputed = edgeagesimputed,
              edgediss = edgediss,
              meanage = meanage,
              meanageimputed = meanageimputed,
              propdiss = propdiss,
              meanmeanageimputed = meanmeanageimputed))
}
