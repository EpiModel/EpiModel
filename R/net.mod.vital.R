#' @title Departures: netsim Module
#'
#' @description This function simulates departure for use in \link{netsim} simulations.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
departures.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  type <- dat$control$type

  #Book-keeping

  idsDpt <- list()

  # Susceptible departures ------------------------------------------------------

  # Initialize counts and pull rates
  nDepartures.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)

  if (nElig.sus > 0) {

    # Departure rates by group
    rates.sus <- dat$param$ds.rate

    # Stochastic exits
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, rates.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsDpt$sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- length(idsDpt.sus)
      dat$attr$active[idsDpt.sus] <- 0
      dat$attr$exitTime[idsDpt.sus] <- at
    }
  }


  # Infected departures ---------------------------------------------------------

  # Initialize counts and query rates
  nDepartures.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)

  if (nElig.inf > 0) {

    # Departure rates by group
    rates.inf <- dat$param$di.rate

    # Stochastic exits
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, rates.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsDpt$inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- length(idsDpt.inf)
      dat$attr$active[idsDpt.inf] <- 0
      dat$attr$exitTime[idsDpt.inf] <- at
    }
  }


  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {

    # Initialize counts and query rates
    nDepartures.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)

    if (nElig.rec > 0) {

      # Departure rates by group
      rates.rec <- dat$param$dr.rate

      # Stochastic exits
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, rates.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsDpt$rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- length(idsDpt.rec)
        dat$attr$active[idsDpt.rec] <- 0
        dat$attr$exitTime[idsDpt.rec] <- at
      }
    }
  }


  # Output ------------------------------------------------------------------

  dat$nw.touch$dpt$idsDpt <- idsDpt

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDepartures.sus)
    dat$epi$di.flow <- c(0, nDepartures.inf)
    if (type == "SIR") {
      dat$epi$dr.flow <- c(0, nDepartures.rec)
    }
  } else {
    dat$epi$ds.flow[at] <- nDepartures.sus
    dat$epi$di.flow[at] <- nDepartures.inf
    if (type == "SIR") {
      dat$epi$dr.flow[at] <- nDepartures.rec
    }
  }

  return(dat)
}


#' @title Arrivals: netsim Module
#'
#' @description This function simulates new arrivals into the network
#'   for use in \code{\link{netsim}} simulations.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
arrivals.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  tea.status <- dat$control$tea.status
  nOld <- dat$epi$num[at - 1]
  nCurr <- network.size(dat$nw)

  nArrivals <- 0
  newNodes <- NULL


  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {
    nArrivals <- sum(rbinom(nOld, 1, a.rate))
    if (nArrivals > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nArrivals)
      newNodes <- (nCurr + 1):(nCurr + nArrivals)
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
    }
  }

  # Update Nodal Attributes -------------------------------------------------
  if (length(newNodes) > 0) {

    # Set attributes on nw
    fterms <- dat$temp$fterms
    curr.tab <- get_attr_prop(dat$nw, fterms)
    if (length(curr.tab) > 0) {
      dat$nw <- update_nwattr(dat$nw, newNodes, dat$control$attr.rules,
                              curr.tab, dat$temp$t1.tab)
    }

    # Save any val on attr
    dat <- copy_toall_attr(dat, at, fterms)

    if (tea.status == TRUE) {
      if ("status" %in% fterms) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = dat$attr$status[newNodes],
                                            onset = at, terminus = Inf,
                                            v = newNodes)
      } else {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = "s", onset = at, terminus = Inf,
                                            v = newNodes)
      }
    }
    if (!("status" %in% fterms)) {
      dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
    }
    dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))

    ## Handles infTime when incoming nodes are infected
    newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
    dat$attr$infTime[newNodesInf] <- at

    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.")
    }

  }


  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$a.flow <- c(0, nArrivals)
  } else {
    dat$epi$a.flow[at] <- nArrivals
  }

  return(dat)
}


#' @title Departures: netsim Module
#'
#' @description This function simulates departure for use in \link{netsim} simulations.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'

departures.net.grp <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  group <- get.vertex.attribute(dat$nw, "group")
  idsDpt <-

    # Bookkeeping
    type <- dat$control$type

  # Susceptible departures ------------------------------------------------------

  # Initialize counts and pull rates
  nDepartures.sus <- nDeparturesG2.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    # Departure rates by group
    gElig.sus <- group[idsElig.sus]
    rates.sus <- c(dat$param$ds.rate, dat$param$ds.rate.g2)
    ratesElig.sus <- rates.sus[gElig.sus]
    # Stochastic exits
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, ratesElig.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsDpt$sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- sum(group[idsDpt.sus] == 1)
      nDeparturesG2.sus <- sum(group[idsDpt.sus] == 2)
      dat$attr$active[idsDpt.sus] <- 0
      dat$attr$exitTime[idsDpt.sus] <- at
    }
  }

  # Infected departures ---------------------------------------------------------
  # Initialize counts and query rates
  nDepartures.inf <- nDeparturesG2.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)

  if (nElig.inf > 0) {
    # Departure rates by mode
    gElig.inf <- group[idsElig.inf]
    rates.inf <- c(dat$param$di.rate, dat$param$di.rate.g2)
    ratesElig.inf <- rates.inf[gElig.inf]
    # Stochastic exits
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, ratesElig.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsDpt$inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- sum(group[idsDpt.inf] == 1)
      nDeparturesG2.inf <- sum(group[idsDpt.inf] == 2)
      dat$attr$active[idsDpt.inf] <- 0
      dat$attr$exitTime[idsDpt.inf] <- at
    }
  }

  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {
    # Initialize counts and query rates
    nDepartures.rec <- nDeparturesG2.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)
    if (nElig.rec > 0) {
      # Departure rates by group
      gElig.rec <- group[idsElig.rec]
      rates.rec <- c(dat$param$dr.rate, dat$param$dr.rate.g2)
      ratesElig.rec <- rates.rec[gElig.rec]
      # Stochastic exits
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, ratesElig.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsDpt$rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- sum(group[idsDpt.rec] == 1)
        nDeparturesG2.rec <- sum(group[idsDpt.rec] == 2)
        dat$attr$active[idsDpt.rec] <- 0
        dat$attr$exitTime[idsDpt.rec] <- at
      }
    }
  }

  # Output ------------------------------------------------------------------
  dat$nw.update$idsDpt <- idsDpt

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDepartures.sus)
    dat$epi$di.flow <- c(0, nDepartures.inf)
    if (type == "SIR") {
      dat$epi$dr.flow <- c(0, nDepartures.rec)
    }
    dat$epi$ds.flow.g2 <- c(0, nDeparturesG2.sus)
    dat$epi$di.flow.g2 <- c(0, nDeparturesG2.inf)
    if (type == "SIR") {
      dat$epi$dr.flow.g2 <- c(0, nDeparturesG2.rec)
    }
  } else {
    dat$epi$ds.flow[at] <- nDepartures.sus
    dat$epi$di.flow[at] <- nDepartures.inf
    if (type == "SIR") {
      dat$epi$dr.flow[at] <- nDepartures.rec
    }
    dat$epi$ds.flow.g2[at] <- nDeparturesG2.sus
    dat$epi$di.flow.g2[at] <- nDeparturesG2.inf
    if (type == "SIR") {
      dat$epi$dr.flow.g2[at] <- nDeparturesG2.rec
    }
  }
  return(dat)
}

#' @title Arrivals: netsim Module
#'
#' @description This function simulates new arrivals into the network
#'   for use in \code{\link{netsim}} simulations.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'

arrivals.net.grp <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  a.rate.g2 <- dat$param$a.rate.g2
  tea.status <- dat$control$tea.status
  nOld <- dat$epi$num[at - 1]
  nOldG2 <- dat$epi$num.g2[at - 1]
  a.rand <- dat$control$a.rand

  nArrivals <- nArrivalsG2 <- 0
  newNodes <- newNodesG2 <- NULL
  nCurr <- network.size(dat$nw)


  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {

    if (a.rand == TRUE) {
      if (is.na(a.rate.g2)) {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivalsG2 <- sum(rbinom(nOld, 1, a.rate))
      } else {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivalsG2 <- sum(rbinom(nOldG2, 1, a.rate.g2))
      }
    }
    nArrivals.tot <- nArrivals + nArrivalsG2

    if (nArrivals.tot > 0) {
      newNodes <- (nCurr + 1):(nCurr + nArrivals.tot)
      dat$nw <- add.vertices(dat$nw, nv = nArrivals.tot)
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
      dat$nw <- set.vertex.attribute(dat$nw, "group",
                                     rep(1:2, c(nArrivals, nArrivalsG2)),
                                     newNodes)
    }
  }


  # Update Nodal Attributes -------------------------------------------------
  if (length(newNodes) > 0) {

    # Set attributes on nw
    fterms <- dat$temp$fterms
    curr.tab <- get_attr_prop(dat$nw, fterms)
    if (length(curr.tab) > 0) {
      dat$nw <- update_nwattr(dat$nw, newNodes, dat$control$attr.rules,
                              curr.tab, dat$temp$t1.tab)
    }

    # Save any val on attr
    dat <- copy_toall_attr(dat, at, fterms)

    if (tea.status == TRUE) {
      if ("status" %in% fterms) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = dat$attr$status[newNodes],
                                            onset = at, terminus = Inf,
                                            v = newNodes)
      } else {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = "s", onset = at, terminus = Inf,
                                            v = newNodes)
      }
    }


    if (!("status" %in% fterms)) {
      dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
    }
    dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))

    ## Handles infTime when incoming nodes are infected
    newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
    dat$attr$infTime[newNodesInf] <- at

    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.")
    }

  }


  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$a.flow <- c(0, nArrivals)
    dat$epi$a.flow.g2 <- c(0, nArrivalsG2)
  } else {
    dat$epi$a.flow[at] <- nArrivals
    dat$epi$a.flow.g2[at] <- nArrivalsG2
  }

  return(dat)
}
