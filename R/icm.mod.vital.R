
#' @title Deaths: icm Module
#'
#' @description This function simulates death for use in \code{\link{icm}} 
#'              simulations.
#' 
#' @param all master data list object.
#' @param at current time step.
#' 
#' @seealso \code{\link{icm}}
#' 
#' @export
#' @keywords internal
#'
deaths.icm <- function(all, at) {
  
  # Conditions --------------------------------------------------------------
  if (all$param$vital == FALSE) {
    return(all)
  }
  
  
  # Variables ---------------------------------------------------------------
  groups <- all$param$groups
  group <- all$attr$group
  
  
  
  # Susceptible deaths ------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(all$attr$active == 1 & all$attr$status == 0)
  nElig <- length(idsElig)
  if (nElig > 0) {
    
    gElig <- group[idsElig]
    rates <- c(all$param$ds.rate, all$param$ds.rate.g2)
    ratesElig <- rates[gElig]
    
    if (all$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        all$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      all$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        all$attr$active[ssample(idsElig[gElig == 2], nDeaths)] <- 0
      }
    }
  }
    
  if (at == 2) {
    all$out$ds.flow <- c(0, nDeaths)
    if (groups == 2) {
      all$out$ds.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    all$out$ds.flow[at] <- nDeaths
    if (groups == 2) {
      all$out$ds.flow.g2[at] <- nDeathsG2
    }
  }
  

  # Infected Deaths ---------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(all$attr$active == 1 & all$attr$status == 1)
  nElig <- length(idsElig)
  if (nElig > 0) {
    
    gElig <- group[idsElig]
    rates <- c(all$param$di.rate, all$param$di.rate.g2)
    ratesElig <- rates[gElig]
    
    if (all$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        all$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      all$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        all$attr$active[ssample(idsElig[gElig == 2], nDeaths)] <- 0
      }
    }
  }
  
  if (at == 2) {
    all$out$di.flow <- c(0, nDeaths)
    if (groups == 2) {
      all$out$di.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    all$out$di.flow[at] <- nDeaths
    if (groups == 2) {
      all$out$di.flow.g2[at] <- nDeathsG2
    }
  }

  
  # Recovered Deaths --------------------------------------------------------
  nDeaths <- nDeathsG2 <- 0
  idsElig <- which(all$attr$active == 1 & all$attr$status == 2)
  nElig <- length(idsElig)
  if (nElig > 0) {
    
    gElig <- group[idsElig]
    rates <- c(all$param$dr.rate, all$param$dr.rate.g2)
    ratesElig <- rates[gElig]
    
    if (all$control$d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(group[idsDth] == 1)
        nDeathsG2 <- sum(group[idsDth] == 2)
        all$attr$active[idsDth] <- 0
      }
    } else {
      nDeaths <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      all$attr$active[ssample(idsElig[gElig == 1], nDeaths)] <- 0
      if (groups == 2) {
        nDeathsG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        all$attr$active[ssample(idsElig[gElig == 2], nDeaths)] <- 0
      }
    }
  }
  
  if (at == 2) {
    all$out$dr.flow <- c(0, nDeaths)
    if (groups == 2) {
      all$out$dr.flow.g2 <- c(0, nDeathsG2)
    }
  } else {
    all$out$dr.flow[at] <- nDeaths
    if (groups == 2) {
      all$out$dr.flow.g2[at] <- nDeathsG2
    }
  }
  
  return(all)
}



#' @title Births: icm Module
#'
#' @description This function simulates birth for use in \code{\link{icm}} 
#'              simulations.
#' 
#' @param all master data list object.
#' @param at current time step.
#' 
#' @seealso \code{\link{icm}}
#' 
#' @export
#' @keywords internal
#'
births.icm <- function(all, at) {
  
  # Conditions --------------------------------------------------------------
  if (all$param$vital == FALSE) {
    return(all)
  }
  
  # Variables ---------------------------------------------------------------
  b.rate <- all$param$b.rate
  b.rate.g2 <- all$param$b.rate.g2
  b.rand <- all$control$b.rand
  groups <- all$param$groups
  nOld <- all$out$num[at-1]
  
  
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
    nOldG2 <- all$out$num.g2[at-1]
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
  all$attr$active <- c(all$attr$active, rep(1, totBirths))
  all$attr$group <- c(all$attr$group, c(rep(1, nBirths), rep(2, nBirthsG2)))
  all$attr$status <- c(all$attr$status, rep(0, totBirths))
  all$attr$infTime <- c(all$attr$infTime, rep(NA, totBirths))
  
  
  # Output ------------------------------------------------------------------
  if (at == 2) {
    all$out$b.flow <- c(0, nBirths)
  } else {
    all$out$b.flow[at] <- nBirths
  }
  if (all$param$groups == 2) {
    if (at == 2) {
      all$out$b.flow.g2 <- c(0, nBirthsG2)
    } else {
      all$out$b.flow.g2[at] <- nBirthsG2
    }
  }
  
  return(all)
}
