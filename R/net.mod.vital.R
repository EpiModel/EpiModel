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
  type <- dat$control$type

  # Susceptible departures ------------------------------------------------------
  nDepartures.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    rates.sus <- dat$param$ds.rate
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, rates.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- length(idsDpt.sus)
      dat$attr$active[idsDpt.sus] <- 0
      dat$attr$exitTime[idsDpt.sus] <- at
    }
  }

  # Infected departures ---------------------------------------------------------
  nDepartures.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)
  if (nElig.inf > 0) {
    rates.inf <- dat$param$di.rate
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, rates.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- length(idsDpt.inf)
      dat$attr$active[idsDpt.inf] <- 0
      dat$attr$exitTime[idsDpt.inf] <- at
    }
  }

  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {
    nDepartures.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)
    if (nElig.rec > 0) {
      rates.rec <- dat$param$dr.rate
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, rates.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- length(idsDpt.rec)
        dat$attr$active[idsDpt.rec] <- 0
        dat$attr$exitTime[idsDpt.rec] <- at
      }
    }
  }

  # Output ------------------------------------------------------------------

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
  nOld <- dat$epi$num[at - 1]
  tergmLite <- dat$control$tergmLite
  nCurr <- length(dat$attr$active)

  nArrivals <- 0

  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {
    nArrivals <- sum(rbinom(nOld, 1, a.rate))
    if (nArrivals > 0) {
      if (tergmLite == FALSE) {
        nwterms <- dat$temp$nwterms
        if (!("status" %in% nwterms)) {
          dat$attr$status <- c(dat$attr$status, rep("s", nArrivals))
        }
        dat$attr$active <- c(dat$attr$active, rep(1, nArrivals))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrivals))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nArrivals))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nArrivals))

        ## Handles infTime when incoming nodes are infected
        newNodes <- c((nCurr+1):(nCurr+nArrivals))
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1 & is.null(nwterms)) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }

      if (tergmLite == TRUE) {
        dat$attr$status <- c(dat$attr$status, rep("s", nArrivals))
        dat$attr$active <- c(dat$attr$active, rep(1, nArrivals))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrivals))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nArrivals))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nArrivals))

        ## Handles infTime when incoming nodes are infected
        newNodes <- c((nCurr+1):(nCurr+nArrivals))
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }
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
departures.2g.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  if (dat$control$tergmLite == FALSE){
    group <- get.vertex.attribute(dat$nw, "group")
  } else {
    group <- dat$attr$group
  }

  type <- dat$control$type

  # Susceptible departures ------------------------------------------------------
  nDepartures.sus <- nDeparturesG2.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    gElig.sus <- group[idsElig.sus]
    rates.sus <- c(dat$param$ds.rate, dat$param$ds.rate.g2)
    ratesElig.sus <- rates.sus[gElig.sus]
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, ratesElig.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- sum(group[idsDpt.sus] == 1)
      nDeparturesG2.sus <- sum(group[idsDpt.sus] == 2)
      dat$attr$active[idsDpt.sus] <- 0
      dat$attr$exitTime[idsDpt.sus] <- at
    }
  }

  # Infected departures ---------------------------------------------------------
  nDepartures.inf <- nDeparturesG2.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)
  if (nElig.inf > 0) {
    gElig.inf <- group[idsElig.inf]
    rates.inf <- c(dat$param$di.rate, dat$param$di.rate.g2)
    ratesElig.inf <- rates.inf[gElig.inf]
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, ratesElig.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- sum(group[idsDpt.inf] == 1)
      nDeparturesG2.inf <- sum(group[idsDpt.inf] == 2)
      dat$attr$active[idsDpt.inf] <- 0
      dat$attr$exitTime[idsDpt.inf] <- at
    }
  }

  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {
    nDepartures.rec <- nDeparturesG2.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)
    if (nElig.rec > 0) {
      gElig.rec <- group[idsElig.rec]
      rates.rec <- c(dat$param$dr.rate, dat$param$dr.rate.g2)
      ratesElig.rec <- rates.rec[gElig.rec]
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, ratesElig.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- sum(group[idsDpt.rec] == 1)
        nDeparturesG2.rec <- sum(group[idsDpt.rec] == 2)
        dat$attr$active[idsDpt.rec] <- 0
        dat$attr$exitTime[idsDpt.rec] <- at
      }
    }
  }

  # Output ------------------------------------------------------------------
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
arrivals.2g.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  a.rate.g2 <- dat$param$a.rate.g2
  nOld <- dat$epi$num[at - 1]
  nOldG2 <- dat$epi$num.g2[at - 1]
  a.rand <- dat$control$a.rand
  tergmLite <- dat$control$tergmLite

  totArr <- nArrivals <- nArrivalsG2 <- 0
  newNodes <- newNodesG2 <- NULL

  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {

    if (is.na(a.rate.g2)) {
      nArrivals <- sum(rbinom(nOld, 1, a.rate))
      nArrivalsG2 <- sum(rbinom(nOld, 1, a.rate))
      totArr <- nArrivals + nArrivalsG2
    } else {
      nArrivals <- sum(rbinom(nOld, 1, a.rate))
      nArrivalsG2 <- sum(rbinom(nOldG2, 1, a.rate.g2))
      totArr <- nArrivals + nArrivalsG2
    }

    if (totArr > 0) {
      if (tergmLite == FALSE) {
        nCurr <- length(dat$attr$active)
        newNodes <- (nCurr + 1):(nCurr + totArr)
        nwterms <- dat$temp$nwterms
        if (!("status" %in% nwterms)) {
          dat$attr$status <- c(dat$attr$status, rep("s", totArr))
        }
        dat$attr$group <- c(dat$attr$group, c(rep(1, nArrivals),
                                              rep(2, nArrivalsG2)))
        dat$attr$active <- c(dat$attr$active, rep(1, totArr))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArr))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, totArr))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, totArr))

        ## Handles infTime when incoming nodes are infected
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }

      if (tergmLite == TRUE) {
        nCurr <- length(which(dat$attr$active == 1))
        newNodes <- (nCurr + 1):(nCurr + totArr)

        dat$attr$group <- c(dat$attr$group, c(rep(1, nArrivals),
                                              rep(2, nArrivalsG2)))
        dat$attr$status <- c(dat$attr$status, rep("s", sum(totArr)))
        dat$attr$active <- c(dat$attr$active, rep(1, sum(totArr)))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, sum(totArr)))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, sum(totArr)))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, sum(totArr)))

        ## Handles infTime when incoming nodes are infected
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }
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
