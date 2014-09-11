
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
resim_nets <- function(all, at) {

  idsActive <- which(all$attr$active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (all$param$modes == 2) {
    nActiveM1 <- length(intersect(modeids(all$nw, mode = 1), idsActive))
    nActiveM2 <- length(intersect(modeids(all$nw, mode = 2), idsActive))
    anyActive <- ifelse(nActiveM1 > 0 & nActiveM2 > 0, TRUE, FALSE)
  }

  # Pull network model parameters
  nwparam <- get_nwparam(all)

  # Serosorting model check
  statOnNw <- ("status" %in% get_formula_terms(nwparam$formation))
  status <- all$attr$status
  if (statOnNw == TRUE && length(unique(status)) == 1) {
    stop("Stopping simulation because status in formation formula and no longer any discordant nodes",
         call. = TRUE)
  }

  # Set up nwstats df
  if (all$control$save.nwstats == TRUE) {
    if (at == 2) {
      nwstats <- attributes(all$nw)$stats
      all$stats$nwstats <- as.data.frame(nwstats)
    }
  }

  # Network simulation
  if (anyActive > 0 & all$control$depend == TRUE) {
    suppressWarnings(
      all$nw <- simulate(all$nw,
                         formation = nwparam$formation,
                         dissolution = nwparam$dissolution,
                         coef.form = nwparam$coef.form,
                         coef.diss = nwparam$coef.diss$coef.adj,
                         constraints = nwparam$constraints,
                         time.start = at,
                         time.slices = 1,
                         time.offset = 0,
                         monitor = all$control$nwstats.formula,
                         control = all$control$set.control.stergm))

    # Set up nwstats df
    if (all$control$save.nwstats == TRUE) {
      all$stats$nwstats[at, ] <- tail(attributes(all$nw)$stats, 1)
    }

    if (all$control$delete.nodes == TRUE) {
      all$nw <- network.extract(all$nw, at = at)
      inactive <- which(all$attr$active == 0)
      all$attr <- deleteAttr(all$attr, inactive)
    }

  }



  return(all)

}
