
#' @title Departure: icm Module
#'
#' @description This function simulates departure for use in `\link{icm}`
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso `\link{icm}`
#'
#' @export
#' @keywords internal
#'
departures.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Susceptible departures ------------------------------------------------------
  nDepartures <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig <- length(idsElig)
  if (nElig > 0) {

    rates <- dat$param$ds.rate

    vecDepartures <- which(rbinom(nElig, 1, rates) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- length(idsDpt)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDepartures)
  } else {
    dat$epi$ds.flow[at] <- nDepartures
  }


  # Infected Departures ---------------------------------------------------------
  nDepartures <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig <- length(idsElig)
  if (nElig > 0) {

    rates <- dat$param$di.rate

    vecDepartures <- which(rbinom(nElig, 1, rates) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- length(idsDpt)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$di.flow <- c(0, nDepartures)
  } else {
    dat$epi$di.flow[at] <- nDepartures
  }


  # Recovered Departures --------------------------------------------------------
  nDepartures <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "r")
  nElig <- length(idsElig)
  if (nElig > 0) {

    rates <- dat$param$dr.rate

    vecDepartures <- which(rbinom(nElig, 1, rates) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- length(idsDpt)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$dr.flow <- c(0, nDepartures)
  } else {
    dat$epi$dr.flow[at] <- nDepartures
  }

  return(dat)
}


#' @title Arrivals: icm Module
#'
#' @description This function simulates arrival for use in `\link{icm}`
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso `\link{icm}`
#'
#' @export
#' @keywords internal
#'
arrivals.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  nOld <- dat$epi$num[at - 1]


  # Process -----------------------------------------------------------------
  nArrivals <- sum(rbinom(nOld, 1, a.rate))


  ## Set attributes
  totArrivals <- nArrivals
  dat$attr$active <- c(dat$attr$active, rep(1, totArrivals))
  dat$attr$group <- c(dat$attr$group, c(rep(1, nArrivals)))
  dat$attr$status <- c(dat$attr$status, rep("s", totArrivals))
  dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArrivals))


  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$a.flow <- c(0, nArrivals)
  } else {
    dat$epi$a.flow[at] <- nArrivals
  }

  return(dat)
}


#' @title Departure: bipartite icm Module
#'
#' @description This function simulates departure for use in `\link{icm}`
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso `\link{icm}`
#'
#' @export
#' @keywords internal
#'
departures.icm.bip <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  group <- dat$attr$group


  # Susceptible departures ------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$ds.rate, dat$param$ds.rate.g2)
    ratesElig <- rates[gElig]

    vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- sum(group[idsDpt] == 1)
      nDeparturesG2 <- sum(group[idsDpt] == 2)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDepartures)
    dat$epi$ds.flow.g2 <- c(0, nDeparturesG2)
  } else {
    dat$epi$ds.flow[at] <- nDepartures
    dat$epi$ds.flow.g2[at] <- nDeparturesG2
  }


  # Infected Departures ---------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$di.rate, dat$param$di.rate.g2)
    ratesElig <- rates[gElig]

    vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- sum(group[idsDpt] == 1)
      nDeparturesG2 <- sum(group[idsDpt] == 2)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$di.flow <- c(0, nDepartures)
    dat$epi$di.flow.g2 <- c(0, nDeparturesG2)
  } else {
    dat$epi$di.flow[at] <- nDepartures
    dat$epi$di.flow.g2[at] <- nDeparturesG2
  }


  # Recovered Departures --------------------------------------------------------
  nDepartures <- nDeparturesG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "r")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$dr.rate, dat$param$dr.rate.g2)
    ratesElig <- rates[gElig]

    vecDepartures <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecDepartures) > 0) {
      idsDpt <- idsElig[vecDepartures]
      nDepartures <- sum(group[idsDpt] == 1)
      nDeparturesG2 <- sum(group[idsDpt] == 2)
      dat$attr$active[idsDpt] <- 0
    }
  }

  if (at == 2) {
    dat$epi$dr.flow <- c(0, nDepartures)
    dat$epi$dr.flow.g2 <- c(0, nDeparturesG2)
  } else {
    dat$epi$dr.flow[at] <- nDepartures
    dat$epi$dr.flow.g2[at] <- nDeparturesG2
  }

  return(dat)
}

#' @title Arrivals: bipartite icm Module
#'
#' @description This function simulates arrival for use in `\link{icm}`
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso `\link{icm}`
#'
#' @export
#' @keywords internal
#'
arrivals.icm.bip <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }
  # Variables ---------------------------------------------------------------
  a.rate <- dat$param$a.rate
  a.rate.g2 <- dat$param$a.rate.g2
  nOld <- dat$epi$num[at - 1]

  # Process -----------------------------------------------------------------
  nArrivals <- nArrivalsG2 <- 0
  nOldG2 <- dat$epi$num.g2[at - 1]

  if (is.na(a.rate.g2)) {
    nArrivals <- sum(rbinom(nOld, 1, a.rate))
    nArrivalsG2 <- sum(rbinom(nOld, 1, a.rate))
  } else {
    nArrivals <- sum(rbinom(nOld, 1, a.rate))
    nArrivalsG2 <- sum(rbinom(nOldG2, 1, a.rate.g2))
  }

  ## Set attributes
  totArrivals <- nArrivals + nArrivalsG2
  dat$attr$active <- c(dat$attr$active, rep(1, totArrivals))
  dat$attr$group <- c(dat$attr$group, c(rep(1, nArrivals), rep(2, nArrivalsG2)))
  dat$attr$status <- c(dat$attr$status, rep("s", totArrivals))
  dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totArrivals))

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
