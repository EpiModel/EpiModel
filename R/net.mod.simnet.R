
#' @title Simulate Dynamic Network at Time 1
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps, to be used in \code{\link{netsim}} models.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param nw A \code{networkDynamic} object.
#' @param nsteps Number of time steps to simulate the network over.
#' @param nsims Number of independent network simulations.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets <- function(x, nw, nsteps, control) {

  suppressWarnings({
    sim <- simulate(nw,
                    formation = x$formation,
                    dissolution = x$coef.diss$dissolution,
                    coef.form = x$coef.form,
                    coef.diss = x$coef.diss$coef.crude,
                    time.slices = nsteps,
                    time.start = 1,
                    time.offset = 0,
                    constraints = x$constraints,
                    monitor = control$nwstats.formula,
                    nsim = 1,
                    control = control$set.control.stergm)
  })

  return(sim)
}


#' @title Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @param x A master object passed through \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at) {

  # Variables
  tergmLite <- get_control(dat, "tergmLite")
  save.nwstats <- get_control(dat, "save.nwstats")
  resimulate.network <- get_control(dat, "resimulate.network")

  # Edges Correction
  dat <- edges_correct(dat, at)

  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1)
    groupids.2 <- which(group == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  # Pull network model parameters
  nwparam <- get_nwparam(dat)

  # Full tergm/network Method
  if (tergmLite == FALSE) {

    # Set up nwstats df
    if (save.nwstats == TRUE) {
      if (at == 2) {
        nwstats <- attributes(dat$nw[[1]])$stats
        dat$stats$nwstats <- as.data.frame(nwstats)
      }
    }

    # Network simulation
    if (anyActive > 0 & resimulate.network == TRUE) {
      suppressWarnings(
        dat$nw[[1]] <- simulate(dat$nw[[1]],
                           formation = nwparam$formation,
                           dissolution = nwparam$coef.diss$dissolution,
                           coef.form = nwparam$coef.form,
                           coef.diss = nwparam$coef.diss$coef.adj,
                           constraints = nwparam$constraints,
                           time.start = at,
                           time.slices = 1,
                           time.offset = 0,
                           monitor = dat$control$nwstats.formula,
                           control = dat$control$set.control.stergm))

      # Update nwstats df
      if (save.nwstats == TRUE) {
        dat$stats$nwstats <- rbind(dat$stats$nwstats,
                                   tail(attributes(dat$nw[[1]])$stats, 1)[,])
      }
    }
  }

  # networkLite/tergmLite Method
  if (tergmLite == TRUE & resimulate.network == TRUE) {
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE)
    dat <- tergmLite::updateModelTermInputs(dat)
    if (isTERGM == TRUE) {
      dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                                 el = dat$el[[1]],
                                                 coef.form = nwparam$coef.form,
                                                 coef.diss = nwparam$coef.diss$coef.adj,
                                                 save.changes = TRUE)
    } else {
      dat$el[[1]] <- tergmLite::simulate_ergm(p = dat$p[[1]],
                                              el = dat$el[[1]],
                                              coef = nwparam$coef.form)
    }
  }

  return(dat)
}


#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficient in a dynamic network model
#'              simulated in \code{\link{netsim}} to preserve the mean
#'              degree of nodes in the network.
#'
#' @param dat Master object in \code{netsim} simulations.
#' @param at Current time step.
#'
#' @keywords internal
#' @export
#'
edges_correct <- function(dat, at) {

  resimulate.network <- get_control(dat, "resimulate.network")
  tergmLite <- get_control(dat, "tergmLite")
  groups <- get_param(dat, "groups")
  active <- get_attr(dat, "active")

  if (resimulate.network == TRUE) {

    if (groups == 1) {
      index <- at - 1
      old.num <- get_epi(dat, "num", index)
      new.num <- sum(active == 1)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(old.num) -
        log(new.num)
    }
    if (groups == 2) {
      index <- at - 1
      group <- get_attr(dat, "group")
      old.num.g1 <- get_epi(dat, "num", index)
      old.num.g2 <- get_epi(dat, "num.g2", index)
      new.num.g1 <- sum(active == 1 & group == 1)
      new.num.g2 <- sum(active == 1 & group == 2)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
        log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
    }
  }
  return(dat)
}
