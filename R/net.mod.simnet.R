
#' @title Simulate Initial Network at Time 1 for Model Initialization
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps for TERGMs or one or multiple cross-sectional network
#'              panels for ERGMs, for use in \code{\link{netsim}} modeling.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @inheritParams recovery.net
#' @param network Index of network to simulate.
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets_t1 <- function(x, dat, network = 1L) {

  # Simulate t0 basis network
  if (x$edapprox == TRUE) {
    dat$nw[[network]] <- simulate(x$formula,
                                  coef = x$coef.form.crude,
                                  basis = x$newnetwork,
                                  constraints = x$constraints,
                                  control = get_control(dat, "set.control.ergm"),
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
    for (netattrname in setdiff(list.network.attributes(dat$nw[[network]]), names(attributes(dat$el[[network]])))) {
      attr(dat$el[[network]], netattrname) <- get.network.attribute(dat$nw[[network]], netattrname)
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
                                                           onset = 0,
                                                           terminus = Inf)
  }

  return(dat)
}

simulate_dat <- function(dat, at, network = 1L, nsteps = 1L) {
  ## get network for (re)simulation
  if (get_control(dat, "tergmLite") == FALSE) {
    # Full tergm/network Method
    nw <- dat$nw[[network]]
    output <- "networkDynamic"
  } else {
    # tergmLite/networkLite Method
    nw <- networkLite(dat$el[[network]], dat$attr)
    output <- "final"
    if (get_control(dat, "tergmLite.track.duration") == TRUE) {
      nw %n% "time" <- dat$nw[[network]] %n% "time"
      nw %n% "lasttoggle" <- dat$nw[[network]] %n% "lasttoggle"
    }
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
  
  # TERGM simulation
  dat$nw[[network]] <- simulate(formula,
                                coef = coef,
                                basis = nw,
                                constraints = nwparam$constraints,
                                time.start = at - 1,
                                time.slices = nsteps,
                                output = output,
                                monitor = get_control(dat, "nwstats.formula"),
                                control = get_control(dat, "set.control.tergm"),
                                dynamic = TRUE)

  # Update nwstats data frame
  if (get_control(dat, "save.nwstats") == TRUE) {
    new.nwstats <- attributes(dat$nw[[network]])$stats
    keep.cols <- which(!duplicated(colnames(new.nwstats)))
    new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
    dat$stats$nwstats[[network]] <- rbind(dat$stats$nwstats[[network]], new.nwstats)
  }
  
  if (get_control(dat, "tergmLite") == TRUE) {
    dat$el[[network]] <- as.edgelist(dat$nw[[network]])
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
    dat <- simulate_dat(dat, at)
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
