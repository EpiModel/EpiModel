
#' @title Simulate Initial Network at Time 1 for Model Initialization
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps for TERGMs or one or multiple cross-sectional network
#'              panels for ERGMs, for use in \code{\link{netsim}} modeling.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @inheritParams recovery.net
#' @param nsteps For TERGMs, the number of time steps to simulate the network
#'        over; for ERGMs, the number of independent network panels to simulate.
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets_t1 <- function(x, dat, nsteps) {

  # control settings
  nwparam <- get_nwparam(dat)
  isTERGM <- all(nwparam$coef.diss$duration > 1)
  dat <- set_control(dat, "isTERGM", isTERGM)

  # Reset default formula for nwstats.formula
  if (get_control(dat, "nwstats.formula") == "formation") {
    dat <- set_control(dat, "nwstats.formula", x$formation)
  }
  nwstats.formula <- get_control(dat, "nwstats.formula")

  save.nwstats <- get_control(dat, "save.nwstats")

  # if TERGM, then use tergm package simulation for dynamic network
  if (isTERGM == TRUE) {

    set.control.ergm <- get_control(dat, "set.control.ergm")

    # Simulate t0 basis network
    if (x$edapprox == TRUE) {
      nw <- simulate(x$formation,
                     coef = coef(x$fit),
                     basis = x$fit$newnetwork,
                     constraints = x$constraints,
                     control = set.control.ergm,
                     dynamic = FALSE)
    } else {
      nw <- x$fit$network
    }

    set.control.stergm <- get_control(dat, "set.control.stergm")
    set.control.tergm <- get_control(dat, "set.control.tergm")

    if (!is.null(set.control.stergm)) {
      # Simulate dynamic network
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
                        monitor = nwstats.formula,
                        nsim = 1,
                        control = set.control.stergm)
        sim <- networkDynamic::activate.vertices(sim, onset = 1, terminus = Inf)
      })
    } else {
      # Simulate dynamic network
      suppressWarnings({
        sim <- simulate(nw ~ Form(x$formation) + 
                             Persist(x$coef.diss$dissolution),
                        coef = c(x$coef.form, x$coef.diss$coef.crude),
                        time.slices = nsteps,
                        time.start = 1,
                        time.offset = 0,
                        constraints = x$constraints,
                        monitor = nwstats.formula,
                        nsim = 1,
                        control = set.control.tergm,
                        dynamic = TRUE)
        sim <- networkDynamic::activate.vertices(sim, onset = 1, terminus = Inf)
      })
    }
    dat$nw[[1]] <- sim

  # If ERGM, then use ergm package simulation for x-sectional network panels
  } else {

    set.control.ergm <- get_control(dat, "set.control.ergm")

    sim <- simulate(x$formation,
                    coef = coef(x$fit),
                    basis = x$fit$newnetwork,
                    constraints = x$constraints,
                    control = set.control.ergm,
                    dynamic = FALSE,
                    monitor = nwstats.formula,
                    nsim = nsteps)

    # save working network and temp network list
    if (nsteps == 1) {
      dat$nw[[1]] <- sim
      dat$temp$nw_list <- list(sim)
    } else {
      dat$nw[[1]] <- sim[[1]]
      dat$temp$nw_list <- sim
    }
  }

  # Set up nwstats df
  if (save.nwstats == TRUE) {
    if (isTERGM == FALSE & nsteps > 1) {
      nwstats <- as.data.frame(
        simulate(x$formation,
                 coef = coef(x$fit),
                 basis = x$fit$newnetwork,
                 constraints = x$constraints,
                 control = set.control.ergm,
                 dynamic = FALSE,
                 monitor = nwstats.formula,
                 nsim = nsteps, output = "stats"))
    } else {
      nwstats <- attributes(dat$nw[[1]])$stats
    }
    keep.cols <- which(!duplicated(colnames(nwstats)))
    nwstats <- nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[1]] <- nwstats
  }

  return(dat)
}


#' @title Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at) {

  # Control settings for resimulation
  tergmLite <- get_control(dat, "tergmLite")
  isTERGM <- get_control(dat, "isTERGM")
  save.nwstats <- get_control(dat, "save.nwstats")
  resimulate.network <- get_control(dat, "resimulate.network")
  nwstats.formula <- get_control(dat, "nwstats.formula")
  set.control.stergm <- get_control(dat, "set.control.stergm")
  set.control.tergm <- get_control(dat, "set.control.tergm")
  set.control.ergm <- get_control(dat, "set.control.ergm")
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")

  # Edges Correction
  dat <- edges_correct(dat, at)

  # active attribute (all models)
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)

  # group attribute (built-in models)
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

  # Network resimulation
  if (anyActive == TRUE & resimulate.network == TRUE) {

    # Full tergm/network Method
    if (tergmLite == FALSE) {

      # TERGM simulation
      if (isTERGM == TRUE) {
        if (!is.null(set.control.stergm)) {
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
                                    monitor = nwstats.formula,
                                    control = set.control.stergm))
        } else {
          suppressWarnings(
            dat$nw[[1]] <- simulate(dat$nw[[1]] ~ Form(nwparam$formation) + 
                                                  Persist(nwparam$coef.diss$dissolution),
                                    coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
                                    constraints = nwparam$constraints,
                                    time.start = at,
                                    time.slices = 1,
                                    time.offset = 0,
                                    monitor = nwstats.formula,
                                    control = set.control.tergm,
                                    dynamic = TRUE))
        }
      } else {
        dat$nw[[1]] <- simulate(object = nwparam$formation,
                                basis = dat$nw[[1]],
                                coef = nwparam$coef.form,
                                constraints = nwparam$constraints,
                                control = set.control.ergm,
                                dynamic = FALSE,
                                monitor = nwstats.formula,
                                nsim = 1)
      }

      # Update nwstats data frame
      if (save.nwstats == TRUE) {
        new.nwstats <- tail(attributes(dat$nw[[1]])$stats, 1)
        keep.cols <- which(!duplicated(colnames(new.nwstats)))
        new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], new.nwstats)
      }
    }

    # networkLite/tergmLite Method
    if (tergmLite == TRUE) {
      nwL <- networkLite(dat$el[[1]], dat$attr)
      if (tergmLite.track.duration == TRUE) {
        nwL %n% "time" <- dat$nw[[1]] %n% "time"
        nwL %n% "lasttoggle" <- dat$nw[[1]] %n% "lasttoggle"
      }

      if (isTERGM == TRUE) {
        if (!is.null(set.control.stergm)) {
          dat$nw[[1]] <- simulate(nwL,
                                  formation = nwparam$formation,
                                  dissolution = nwparam$coef.diss$dissolution,
                                  coef.form = nwparam$coef.form,
                                  coef.diss = nwparam$coef.diss$coef.adj,
                                  constraints = nwparam$constraints,
                                  time.start = at - 1,
                                  time.slices = 1,
                                  time.offset = 1, # default value
                                  control = set.control.stergm,
                                  output = "final")
        } else {
          dat$nw[[1]] <- simulate(nwL ~ Form(nwparam$formation) + 
                                        Persist(nwparam$coef.diss$dissolution),
                                  coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
                                  constraints = nwparam$constraints,
                                  time.start = at - 1,
                                  time.slices = 1,
                                  time.offset = 1, # default value
                                  control = set.control.tergm,
                                  output = "final",
                                  dynamic = TRUE)
        }
      } else {
        dat$nw[[1]] <- simulate(object = nwparam$formation,
                                basis = nwL,
                                coef = nwparam$coef.form,
                                constraints = nwparam$constraints,
                                control = set.control.ergm,
                                dynamic = FALSE,
                                nsim = 1,
                                output = "network")
      }

      dat$el[[1]] <- as.edgelist(dat$nw[[1]])

      if (save.nwstats == TRUE) {
        term.options <- if (isTERGM == TRUE) {
          if (!is.null(set.control.stergm)) {
            set.control.stergm$term.options
          } else {
            set.control.tergm$term.options          
          }
        } else {
          set.control.ergm$term.options
        }
        nwstats <- summary(dat$control$nwstats.formulas[[1]],
                           basis = dat$nw[[1]],
                           term.options = term.options,
                           dynamic = isTERGM)
        keep.cols <- which(!duplicated(names(nwstats)))
        dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]],
                                        nwstats[keep.cols])
      }
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
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @keywords internal
#' @export
#'
edges_correct <- function(dat, at) {

  resimulate.network <- get_control(dat, "resimulate.network")
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
