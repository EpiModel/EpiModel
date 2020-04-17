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
  if (get_param(dat, "vital") == FALSE) {
    return(dat)
  }
  type <- get_control(dat, "type")

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  exitTime <- get_attr(dat, "exitTime")
  rates.sus <- get_param(dat, "ds.rate")
  rates.inf <- get_param(dat, "di.rate")
  if (type == "SIR") {
    rates.rec <- get_param(dat, "dr.rate")
  }

  # Susceptible departures ------------------------------------------------------
  nDepartures.sus <- 0
  idsElig.sus <- which(active == 1 & status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, rates.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- length(idsDpt.sus)
      active[idsDpt.sus] <- 0
      exitTime[idsDpt.sus] <- at
    }
  }

  # Infected departures ---------------------------------------------------------
  nDepartures.inf <- 0
  idsElig.inf <- which(active == 1 & status == "i")
  nElig.inf <- length(idsElig.inf)
  if (nElig.inf > 0) {
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, rates.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- length(idsDpt.inf)
      active[idsDpt.inf] <- 0
      exitTime[idsDpt.inf] <- at
    }
  }

  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {
    nDepartures.rec <- 0
    idsElig.rec <- which(active == 1 & status == "r")
    nElig.rec <- length(idsElig.rec)
    if (nElig.rec > 0) {
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, rates.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- length(idsDpt.rec)
        active[idsDpt.rec] <- 0
        exitTime[idsDpt.rec] <- at
      }
    }
  }

  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  if (at == 2) {
    dat <- set_epi(dat, "di.flow", 1, 0)
    dat <- set_epi(dat, "di.flow", 2, nDepartures.inf)
    dat <- set_epi(dat, "dr.flow", 1, 0)
    dat <- set_epi(dat, "ds.flow", 2, nDepartures.sus)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow", 1, 0)
      dat <- set_epi(dat, "dr.flow", 2, nDepartures.rec)
    }
  } else {
    dat <- set_epi(dat, "ds.flow", at, nDepartures.sus)
    dat <- set_epi(dat, "di.flow", at, nDepartures.inf)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow", at, nDepartures.rec)
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
  if (get_param(dat, "vital") == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- get_param(dat, "a.rate")
  index <- at-1
  nOld <- get_epi(dat, "num", index)
  nArrivals <- 0

  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {
    nArrivals <- sum(rbinom(nOld, 1, a.rate))
    if (nArrivals > 0) {
      dat <- append_attr(dat, "status", "s", nArrivals)
      dat <- append_attr(dat, "active", 1, nArrivals)
      dat <- append_attr(dat, "infTime", NA, nArrivals)
      dat <- append_attr(dat, "entrTime", at, nArrivals)
      dat <- append_attr(dat, "exitTime", NA, nArrivals)
    }
  }

  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat <- set_epi(dat, "a.flow", 1, 0)
    dat <- set_epi(dat, "a.flow", 2, nArrivals)
  } else {
    dat <- set_epi(dat, "a.flow", at, nArrivals)
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
  if (get_param(dat, "vital") == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------

  type <- get_control(dat, "type")

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  exitTime <- get_attr(dat, "exitTime")
  group <- get_attr(dat, "group")
  rates.sus <- c(get_param(dat, "ds.rate"), get_param(dat, "ds.rate.g2"))
  rates.inf <- c(get_param(dat, "di.rate"), get_param(dat, "di.rate.g2"))
  if (type == "SIR")
  rates.rec <- c(get_param(dat, "dr.rate"), get_param(dat, "dr.rate.g2"))

  # Susceptible departures ------------------------------------------------------
  nDepartures.sus <- nDeparturesG2.sus <- 0
  idsElig.sus <- which(active == 1 & status == "s")
  nElig.sus <- length(idsElig.sus)
  if (nElig.sus > 0) {
    gElig.sus <- group[idsElig.sus]
    ratesElig.sus <- rates.sus[gElig.sus]
    vecDepartures.sus <- which(rbinom(nElig.sus, 1, ratesElig.sus) == 1)
    if (length(vecDepartures.sus) > 0) {
      idsDpt.sus <- idsElig.sus[vecDepartures.sus]
      nDepartures.sus <- sum(group[idsDpt.sus] == 1)
      nDeparturesG2.sus <- sum(group[idsDpt.sus] == 2)
      active[idsDpt.sus] <- 0
      exitTime[idsDpt.sus] <- at
    }
  }

  # Infected departures ---------------------------------------------------------
  nDepartures.inf <- nDeparturesG2.inf <- 0
  idsElig.inf <- which(active == 1 & status == "i")
  nElig.inf <- length(idsElig.inf)
  if (nElig.inf > 0) {
    gElig.inf <- group[idsElig.inf]
    ratesElig.inf <- rates.inf[gElig.inf]
    vecDepartures.inf <- which(rbinom(nElig.inf, 1, ratesElig.inf) == 1)
    if (length(vecDepartures.inf) > 0) {
      idsDpt.inf <- idsElig.inf[vecDepartures.inf]
      nDepartures.inf <- sum(group[idsDpt.inf] == 1)
      nDeparturesG2.inf <- sum(group[idsDpt.inf] == 2)
      active[idsDpt.inf] <- 0
      exitTime[idsDpt.inf] <- at
    }
  }

  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {
    nDepartures.rec <- nDeparturesG2.rec <- 0
    idsElig.rec <- which(active == 1 & status == "r")
    nElig.rec <- length(idsElig.rec)
    if (nElig.rec > 0) {
      gElig.rec <- group[idsElig.rec]
      ratesElig.rec <- rates.rec[gElig.rec]
      vecDepartures.rec <- which(rbinom(nElig.rec, 1, ratesElig.rec) == 1)
      if (length(vecDepartures.rec) > 0) {
        idsDpt.rec <- idsElig.rec[vecDepartures.rec]
        nDepartures.rec <- sum(group[idsDpt.rec] == 1)
        nDeparturesG2.rec <- sum(group[idsDpt.rec] == 2)
        active[idsDpt.rec] <- 0
        exitTime[idsDpt.rec] <- at
      }
    }
  }

  # Output ------------------------------------------------------------------
  dat <- set_attr(dat, "active", active)
  if (at == 2) {
    dat <- set_epi(dat, "di.flow", 1, 0)
    dat <- set_epi(dat, "di.flow", 2, nDepartures.inf)
    dat <- set_epi(dat, "ds.flow", 1, 0)
    dat <- set_epi(dat, "ds.flow", 2, nDepartures.sus)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow", 1, 0)
      dat <- set_epi(dat, "dr.flow", 2, nDepartures.rec)
    }
    dat <- set_epi(dat, "di.flow.g2", 1, 0)
    dat <- set_epi(dat, "di.flow.g2", 2, nDeparturesG2.inf)
    dat <- set_epi(dat, "ds.flow.g2", 1, 0)
    dat <- set_epi(dat, "ds.flow.g2", 2, nDeparturesG2.sus)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow.g2", 1, 0)
      dat <- set_epi(dat, "dr.flow.g2", 2, nDeparturesG2.rec)
    }
  } else {
    dat <- set_epi(dat, "ds.flow", at, nDepartures.sus)
    dat <- set_epi(dat, "di.flow", at, nDepartures.inf)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow", at, nDepartures.rec)
    }
    dat <- set_epi(dat, "ds.flow.g2", at, nDeparturesG2.sus)
    dat <- set_epi(dat, "di.flow.g2", at, nDeparturesG2.inf)
    if (type == "SIR") {
      dat <- set_epi(dat, "dr.flow.g2", at, nDeparturesG2.rec)
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
  if (get_param(dat, "vital") == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- get_param(dat, "a.rate")
  a.rate.g2 <- get_param(dat, "a.rate.g2")
  index <- at-1
  nCurr <- length(get_attr(dat, "active"))
  nOld <- get_epi(dat, "num", index)
  nOldG2 <- get_epi(dat, "num.g2", index)
  tergmLite <- get_control(dat, "tergmLite")
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
      newNodes <- (nCurr + 1):(nCurr + totArr)
      dat <- append_attr(dat, "group", 1, nArrivals)
      dat <- append_attr(dat, "group", 2, nArrivalsG2)
      dat <- append_attr(dat, "status", "s", totArr)
      dat <- append_attr(dat, "active", 1, totArr)
      dat <- append_attr(dat, "infTime", NA, totArr)
      dat <- append_attr(dat, "entrTime", at, totArr)
      dat <- append_attr(dat, "exitTime", NA, totArr)
    }
  }

  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat <- set_epi(dat, "a.flow", 1, 0)
    dat <- set_epi(dat, "a.flow", 2, nArrivals)
    dat <- set_epi(dat, "a.flow.g2", 1, 0)
    dat <- set_epi(dat, "a.flow.g2", 2, nArrivalsG2)
  } else {
    dat <- set_epi(dat, "a.flow", at, nArrivals)
    dat <- set_epi(dat, "a.flow.g2",at, nArrivalsG2)
  }

  return(dat)
}
