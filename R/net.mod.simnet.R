
#' @title Simulate Dynamic Network at Time 1
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps, to be used in \code{\link{netsim}} models.
#'
#' @param x an \code{EpiModel} object of class \code{\link{netest}}.
#' @param nw a \code{networkDynamic} object.
#' @param nsteps number of time steps to simulate the network over.
#' @param nsims number of independent network simulations.
#' @param control an \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets <- function(x, nw, nsteps, control) {

  if (x$edapprox == FALSE) {
    suppressWarnings(
      sim <- simulate(nw,
                      time.slices = nsteps,
                      monitor = control$nwstats.formula,
                      nsim = 1,
                      time.start = 1,
                      time.offset = 0,
                      control = control$set.control.stergm))
  } else {
    suppressWarnings(
      sim <- simulate(nw,
                      formation = x$formation,
                      dissolution = x$dissolution,
                      coef.form = x$coef.form,
                      coef.diss = x$coef.diss$coef.crude,
                      time.slices = nsteps,
                      time.start = 1,
                      time.offset = 0,
                      constraints = x$constraints,
                      monitor = control$nwstats.formula,
                      nsim = 1,
                      control = control$set.control.stergm))
  }

  return(sim)
}


#' @title Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @param x a master object passed through \code{\link{netsim}}.
#' @param at current time step.
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at) {

  idsActive <- which(dat$attr$active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$modes == 2) {
    nActiveM1 <- length(intersect(modeids(dat$nw, mode = 1), idsActive))
    nActiveM2 <- length(intersect(modeids(dat$nw, mode = 2), idsActive))
    anyActive <- ifelse(nActiveM1 > 0 & nActiveM2 > 0, TRUE, FALSE)
  }

  # Pull network model parameters
  nwparam <- get_nwparam(dat)

  # Serosorting model check
  statOnNw <- ("status" %in% get_formula_terms(nwparam$formation))
  status <- dat$attr$status
  if (statOnNw == TRUE && length(unique(status)) == 1) {
    stop("Stopping simulation because status in formation formula and no longer any discordant nodes",
         call. = TRUE)
  }

  # Set up nwstats df
  if (dat$control$save.nwstats == TRUE) {
    if (at == 2) {
      nwstats <- attributes(dat$nw)$stats
      dat$stats$nwstats <- as.data.frame(nwstats)
    }
  }

  # Network simulation
  if (anyActive > 0 & dat$control$depend == TRUE) {
    suppressWarnings(
      dat$nw <- simulate(dat$nw,
                         formation = nwparam$formation,
                         dissolution = nwparam$dissolution,
                         coef.form = nwparam$coef.form,
                         coef.diss = nwparam$coef.diss$coef.adj,
                         constraints = nwparam$constraints,
                         time.start = at,
                         time.slices = 1,
                         time.offset = 0,
                         monitor = dat$control$nwstats.formula,
                         control = dat$control$set.control.stergm))

    # Set up nwstats df
    if (dat$control$save.nwstats == TRUE) {
      dat$stats$nwstats[at, ] <- tail(attributes(dat$nw)$stats, 1)
    }

    if (dat$control$delete.nodes == TRUE) {
      dat$nw <- network.extract(dat$nw, at = at)
      inactive <- which(dat$attr$active == 0)
      dat$attr <- deleteAttr(dat$attr, inactive)
    }

  }

  return(dat)
}
