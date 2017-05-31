
#' @title Deaths: icm Module
#'
#' @description This function simulates death for use in \code{\link{icm}}
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso \code{\link{icm}}
#'
#' @export
#' @keywords internal
#'
deaths.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  groups <- dat$param$groups
  group <- dat$attr$group


  # Susceptible deaths ------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$ds.rate, dat$param$ds.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        dat$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeathsG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDeaths)
    if (groups == 2) {
      dat$epi$ds.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    dat$epi$ds.flow[at] <- nDeaths
    if (groups == 2) {
      dat$epi$ds.flow.g2[at] <- nDeathsG2
    }
  }


  # Infected Deaths ---------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$di.rate, dat$param$di.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        dat$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeathsG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$di.flow <- c(0, nDeaths)
    if (groups == 2) {
      dat$epi$di.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    dat$epi$di.flow[at] <- nDeaths
    if (groups == 2) {
      dat$epi$di.flow.g2[at] <- nDeathsG2
    }
  }


  # Recovered Deaths --------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "r")
  nElig <- length(idsElig)
  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(dat$param$dr.rate, dat$param$dr.rate.g2)
    ratesElig <- rates[gElig]

    if (dat$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        dat$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      dat$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        dat$attr$active[ssample(idsElig[gElig == 2], nDeathsG2)] <- 0
      }
    }
  }

  if (at == 2) {
    dat$epi$dr.flow <- c(0, nDeaths)
    if (groups == 2) {
      dat$epi$dr.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    dat$epi$dr.flow[at] <- nDeaths
    if (groups == 2) {
      dat$epi$dr.flow.g2[at] <- nDeathsG2
    }
  }

  return(dat)
}



#' @title Births: icm Module
#'
#' @description This function simulates birth for use in \code{\link{icm}}
#'              simulations.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @seealso \code{\link{icm}}
#'
#' @export
#' @keywords internal
#'
births.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  b.rate <- dat$param$b.rate
  b.rate.g2 <- dat$param$b.rate.g2
  b.rand <- dat$control$b.rand
  groups <- dat$param$groups
  nOld <- dat$epi$num[at - 1]


  # Process -----------------------------------------------------------------
  nBirths <- nBirthsG2 <- 0

  if (groups == 1) {
    if (b.rand == TRUE) {
      nBirths <- sum(rbinom(nOld, 1, b.rate))
    }
    if (b.rand == FALSE) {
      nBirths <- round(nOld * b.rate)
    }
  }
  if (groups == 2) {
    nOldG2 <- dat$epi$num.g2[at - 1]
    if (b.rand == TRUE) {
      if (is.na(b.rate.g2)) {
        nBirths <- sum(rbinom(nOld, 1, b.rate))
        nBirthsG2 <- sum(rbinom(nOld, 1, b.rate))
      } else {
        nBirths <- sum(rbinom(nOld, 1, b.rate))
        nBirthsG2 <- sum(rbinom(nOldG2, 1, b.rate.g2))
      }
    }
    if (b.rand == FALSE) {
      if (is.na(b.rate.g2)) {
        nBirths <- round(nOld * b.rate)
        nBirthsG2 <- round(nOld * b.rate)
      } else {
        nBirths <- round(nOld * b.rate)
        nBirthsG2 <- round(nOldG2 * b.rate.g2)
      }
    }
  }

  ## Set attributes
  totBirths <- nBirths + nBirthsG2
  dat$attr$active <- c(dat$attr$active, rep(1, totBirths))
  dat$attr$group <- c(dat$attr$group, c(rep(1, nBirths), rep(2, nBirthsG2)))
  dat$attr$status <- c(dat$attr$status, rep("s", totBirths))
  dat$attr$infTime <- c(dat$attr$infTime, rep(NA, totBirths))


  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$b.flow <- c(0, nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }
  if (dat$param$groups == 2) {
    if (at == 2) {
      dat$epi$b.flow.g2 <- c(0, nBirthsG2)
    } else {
      dat$epi$b.flow.g2[at] <- nBirthsG2
    }
  }

  return(dat)
}
