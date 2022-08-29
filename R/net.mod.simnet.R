
#' @title Initialize Network Used in netsim
#'
#' @description This function initializes the network used in
#'              \code{\link{netsim}}, simulating a dynamic network over one or
#'              multiple time steps for TERGMs or one or multiple
#'              cross-sectional network panels for ERGMs.
#'
#' @inheritParams recovery.net
#' @param network index of the network to initialize
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets_t1 <- function(dat, network = 1L) {

  nwparam <- get_nwparam(dat, network = network)

  # Simulate t0 basis network
  if (nwparam$edapprox == TRUE) {
    dat$nw[[network]] <- simulate(nwparam$formula,
                                  coef = nwparam$coef.form.crude,
                                  basis = dat$nw[[network]],
                                  constraints = nwparam$constraints,
                                  control = get_network_control(dat, "set.control.ergm", network = network),
                                  dynamic = FALSE)
  }

  if (get_control(dat, "tergmLite") == TRUE) {
    ## set up el
    dat$el[[network]] <- as.edgelist(dat$nw[[network]])
    if (get_control(dat, "tergmLite.track.duration") == TRUE) {
      ## set up time, lasttoggle
      dat$nw[[network]] %n% "time" <- 0L
      dat$nw[[network]] %n% "lasttoggle" <- cbind(dat$el[[network]], 0L)
    }
    ## copy over network attributes
    for (netattrname in setdiff(list.network.attributes(dat$nw[[network]]),
                                names(attributes(dat$el[[network]])))) {
      attr(dat$el[[network]], netattrname) <-
        get.network.attribute(dat$nw[[network]], netattrname)
    }
  }

  if (get_control(dat, "resimulate.network") == TRUE) {
    nsteps <- 1L
  } else {
    nsteps <- get_control(dat, "nsteps")
  }

  dat <- simulate_dat(dat, at = 1L, network = network, nsteps = nsteps)

  if (get_control(dat, "tergmLite") == FALSE) {
    dat$nw[[network]] <- networkDynamic::activate.vertices(dat$nw[[network]],
                                                           onset = 0L,
                                                           terminus = Inf)
  }

  return(dat)
}

#' @title Construct Network from dat Object
#'
#' @description This function returns the network object representing the
#'              current state of the simulation.
#'
#' @inheritParams recovery.net
#' @param network index of the network to construct
#'
#' @return The network.
#'
#' @export
#' @keywords netUtils internal
#'
make_sim_network <- function(dat, network = 1L) {
  if (get_control(dat, "tergmLite") == FALSE) {
    ## networkDynamic
    nw <- dat$nw[[network]]
  } else {
    ## networkLite
    nw <- networkLite(dat$el[[network]], dat$attr)
    if (get_control(dat, "tergmLite.track.duration") == TRUE) {
      nw %n% "time" <- dat$nw[[network]] %n% "time"
      nw %n% "lasttoggle" <- dat$nw[[network]] %n% "lasttoggle"
    }
  }
  return(nw)
}

#' @title Set Network on dat Object
#'
#' @description This function updates the dat object given the network
#'              representing the current state of the simulation.
#'
#' @inheritParams recovery.net
#' @param nw the network
#' @param network index of the network to set
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
set_sim_network <- function(dat, nw, network = 1L) {
  dat$nw[[network]] <- nw
  if (get_control(dat, "tergmLite") == TRUE) {
    dat$el[[network]] <- as.edgelist(nw)
  }
  return(dat)
}

#' @title Simulate a Network for a Specified Number of Time Steps
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps for TERGMs or one or multiple cross-sectional network
#'              panels for ERGMs, for use in \code{\link{netsim}} modeling.
#'              Network statistics are also extracted and saved if
#'              \code{save.nwstats == TRUE}.
#'
#' @inheritParams recovery.net
#' @param network index of the network to simulate
#' @param nsteps number of time steps to simulate
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
simulate_dat <- function(dat, at, network = 1L, nsteps = 1L) {
  ## get/construct network for (re)simulation
  nw <- make_sim_network(dat, network = network)

  ## determine output type
  if (get_control(dat, "tergmLite") == FALSE) {
    output <- "networkDynamic"
  } else {
    output <- "final"
  }

  nwparam <- get_nwparam(dat, network = network)

  if (all(nwparam$coef.diss$duration > 1)) {
    formula <- ~Form(nwparam$formation) +
                Persist(nwparam$coef.diss$dissolution)
    coef <- c(nwparam$coef.form, nwparam$coef.diss$coef.adj)
  } else {
    formula <- nwparam$formation
    coef <- nwparam$coef.form
  }

  if (get_control(dat, "save.nwstats") == TRUE &&
      get_control(dat, "resimulate.network") == FALSE) {
    monitor <- get_network_control(dat, "nwstats.formula", network = network)
  } else {
    monitor <- NULL # will be handled by summary_nets, if needed
  }

  # TERGM simulation
  nw <- simulate(formula,
                 coef = coef,
                 basis = nw,
                 constraints = nwparam$constraints,
                 time.start = at - 1,
                 time.slices = nsteps,
                 output = output,
                 control = get_network_control(dat, "set.control.tergm", network = network),
                 monitor = monitor,
                 dynamic = TRUE)

  dat <- set_sim_network(dat, nw, network = network)

  if (!is.null(monitor)) {
    new.nwstats <- attributes(nw)$stats
    keep.cols <- which(!duplicated(colnames(new.nwstats)))
    new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[network]] <- rbind(dat$stats$nwstats[[network]], new.nwstats)
    rownames(dat$stats$nwstats[[network]]) <- NULL
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

  # Network resimulation
  if (anyActive == TRUE && get_control(dat, "resimulate.network") == TRUE) {
    dat.updates <- get_control(dat, "dat.updates")
    if (!is.null(dat.updates)) {
      dat <- dat.updates(dat = dat, at = at, network = 0L)
    }
    for (network in seq_along(dat$nwparam)) {
      dat <- simulate_dat(dat, at, network = network)
      if (!is.null(dat.updates)) {
        dat <- dat.updates(dat = dat, at = at, network = network)
      }
    }
  }

  cumulative.edgelist <- get_control(
    dat, "cumulative.edgelist", override.null.error = TRUE)
  # Cummulative edgelist
  if (!is.null(cumulative.edgelist) && cumulative.edgelist == TRUE) {

    truncate.el.cuml <- get_control(
      dat, "truncate.el.cuml", override.null.error = TRUE)
    truncate.el.cuml <- if (is.null(truncate.el.cuml)) 1 else truncate.el.cuml

    for (network in seq_along(dat[["nwparam"]])) {
      dat <- update_cumulative_edgelist(dat, network, truncate.el.cuml)
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
      adjustment <- log(old.num) - log(new.num)
    }
    if (groups == 2) {
      index <- at - 1
      group <- get_attr(dat, "group")
      old.num.g1 <- get_epi(dat, "num", index)
      old.num.g2 <- get_epi(dat, "num.g2", index)
      new.num.g1 <- sum(active == 1 & group == 1)
      new.num.g2 <- sum(active == 1 & group == 2)
      adjustment <-
        log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
        log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
    }

    for (network in seq_along(dat$nwparam)) {
      dat$nwparam[[network]]$coef.form[1] <-
        dat$nwparam[[network]]$coef.form[1] + adjustment
    }
  }
  return(dat)
}
