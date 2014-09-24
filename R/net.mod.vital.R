
#' @title Deaths: netsim Module
#'
#' @description This function simulates death for use in \link{netsim} simulations.
#'
#' @param dat a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
deaths.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  modes <- dat$param$modes
  mode <- idmode(dat$nw)

  type <- dat$control$type
  d.rand <- dat$control$d.rand


  # Susceptible deaths ------------------------------------------------------

  # Initialize counts and query rates
  nDeaths <- nDeathsM2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig <- length(idsElig)

  if (nElig > 0) {

    # Pull rates by mode
    mElig <- mode[idsElig]
    rates <- c(dat$param$ds.rate, dat$param$ds.rate.m2)
    ratesElig <- rates[mElig]

    # Stochastic deaths
    if (d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(mode[idsDth] == 1)
        nDeathsM2 <- sum(mode[idsDth] == 2)
        dat$attr$active[idsDth] <- 0
        dat$nw <- deactivate.vertices(dat$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = idsDth,
                                      deactivate.edges = TRUE)
      }
    }

    # Deterministic deaths
    if (d.rand == FALSE) {
      idsDth <- idsDthM2 <- NULL
      nDeaths <- min(round(sum(ratesElig[mElig == 1])), sum(mElig == 1))
      dat$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
      if (modes == 2) {
        nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
        dat$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
      }
      totDth <- nDeaths + nDeathsM2
      if (totDth > 0) {
        allids <- c(idsDth, idsDthM2)
        dat$nw <- deactivate.vertices(dat$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = allids,
                                      deactivate.edges = TRUE)
      }
    }
  }

  # Output
  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDeaths)
    if (modes == 2) {
      dat$epi$ds.flow.m2 <- c(0, nDeathsM2)
    }
  } else {
    dat$epi$ds.flow[at] <- nDeaths
    if (modes == 2) {
      dat$epi$ds.flow.m2[at] <- nDeathsM2
    }
  }


  # Infected deaths ---------------------------------------------------------

  # Initialize counts and query rates
  nDeaths <- nDeathsM2 <- 0
  idsElig <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig <- length(idsElig)

  if (nElig > 0) {

    # Pull rates by mode
    mElig <- mode[idsElig]
    rates <- c(dat$param$di.rate, dat$param$di.rate.m2)
    ratesElig <- rates[mElig]

    # Stochastic deaths
    if (d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(mode[idsDth] == 1)
        nDeathsM2 <- sum(mode[idsDth] == 2)
        dat$attr$active[idsDth] <- 0
        dat$nw <- deactivate.vertices(dat$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = idsDth,
                                      deactivate.edges = TRUE)
      }
    }

    # Deterministic deaths
    if (d.rand == FALSE) {
      idsDth <- idsDthM2 <- NULL
      nDeaths <- min(round(sum(ratesElig[mElig == 1])), sum(mElig == 1))
      dat$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
      if (modes == 2) {
        nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
        dat$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
      }
      totDth <- nDeaths + nDeathsM2
      if (totDth > 0) {
        allids <- c(idsDth, idsDthM2)
        dat$nw <- deactivate.vertices(dat$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = allids,
                                      deactivate.edges = TRUE)
      }
    }
  }

  # Output
  if (at == 2) {
    dat$epi$di.flow <- c(0, nDeaths)
    if (modes == 2) {
      dat$epi$di.flow.m2 <- c(0, nDeathsM2)
    }
  } else {
    dat$epi$di.flow[at] <- nDeaths
    if (modes == 2) {
      dat$epi$di.flow.m2[at] <- nDeathsM2
    }
  }


  # Recovered deaths --------------------------------------------------------
  if (type == "SIR") {

    # Initialize counts and query rates
    nDeaths <- nDeathsM2 <- 0
    idsElig <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig <- length(idsElig)

    if (nElig > 0) {

      # Pull rates by mode
      mElig <- mode[idsElig]
      rates <- c(dat$param$dr.rate, dat$param$dr.rate.m2)
      ratesElig <- rates[mElig]

      # Stochastic deaths
      if (d.rand == TRUE) {
        vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecDeaths) > 0) {
          idsDth <- idsElig[vecDeaths]
          nDeaths <- sum(mode[idsDth] == 1)
          nDeathsM2 <- sum(mode[idsDth] == 2)
          dat$attr$active[idsDth] <- 0
          dat$nw <- deactivate.vertices(dat$nw,
                                        onset = at,
                                        terminus = Inf,
                                        v = idsDth,
                                        deactivate.edges = TRUE)
        }
      }

      # Deterministic deaths
      if (d.rand == FALSE) {
        idsDth <- idsDthM2 <- NULL
        nDeaths <- min(round(sum(ratesElig[mElig == 1])), sum(mElig == 1))
        dat$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
        if (modes == 2) {
          nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
          dat$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
        }
        totDth <- nDeaths + nDeathsM2
        if (totDth > 0) {
          allids <- c(idsDth, idsDthM2)
          dat$nw <- deactivate.vertices(dat$nw,
                                        onset = at,
                                        terminus = Inf,
                                        v = allids,
                                        deactivate.edges = TRUE)
        }
      }
    }

    # Output
    if (at == 2) {
      dat$epi$dr.flow <- c(0, nDeaths)
      if (modes == 2) {
        dat$epi$dr.flow.m2 <- c(0, nDeathsM2)
      }
    } else {
      dat$epi$dr.flow[at] <- nDeaths
      if (modes == 2) {
        dat$epi$dr.flow.m2[at] <- nDeathsM2
      }
    }
  }


  return(dat)
}



#' @title Births: netsim Module
#'
#' @description This function simulates new births into the network
#'   for use in \code{\link{netsim}} simulations.
#'
#' @param dat a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
births.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (dat$param$vital == FALSE) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  b.rate <- dat$param$b.rate
  b.rate.m2 <- dat$param$b.rate.m2
  formation <- get_nwparam(dat)$formation
  modes <- dat$param$modes
  tea.status <- dat$control$tea.status
  nOld <- dat$epi$num[at - 1]
  nCurr <- network.size(dat$nw)
  b.rand <- dat$control$b.rand
  delete.nodes <- dat$control$delete.nodes

  nBirths <- nBirthsM2 <- 0
  newNodes <- newNodesM2 <- NULL


  # Add Nodes ---------------------------------------------------------------
  if (modes == 1 && nOld > 0) {
    if (b.rand == TRUE) {
      nBirths <- sum(rbinom(nOld, 1, b.rate))
    }
    if (b.rand == FALSE) {
      nBirths <- round(nOld * b.rate)
    }
    if (nBirths > 0) {
      dat$nw <- add.vertices(dat$nw,
                             nv = nBirths)
      newNodes <- (nCurr + 1):(nCurr + nBirths)
      dat$nw <- activate.vertices(dat$nw,
                                  onset = at,
                                  terminus = Inf,
                                  v = newNodes)
    }
  }
  if (modes == 2 && nOld > 0) {
    nOldM2 <- dat$epi$num.m2[at - 1]
    if (b.rand == TRUE) {
      if (is.na(b.rate.m2)) {
        nBirths <- sum(rbinom(nOld, 1, b.rate))
        nBirthsM2 <- sum(rbinom(nOld, 1, b.rate))
      } else {
        nBirths <- sum(rbinom(nOld, 1, b.rate))
        nBirthsM2 <- sum(rbinom(nOldM2, 1, b.rate.m2))
      }
    }
    if (b.rand == FALSE) {
      if (is.na(b.rate.m2)) {
        nBirths <- round(nOld * b.rate)
        nBirthsM2 <- round(nOld * b.rate)
      } else {
        nBirths <- round(nOld * b.rate)
        nBirthsM2 <- round(nOldM2 * b.rate.m2)
      }
    }

    nCurrM1 <- length(modeids(dat$nw, 1))
    nCurrM2 <- length(modeids(dat$nw, 2))
    prefixes <- unique(substr(dat$nw %v% "vertex.names", 1, 1))

    if (nBirths > 0) {
      newNodeIDs <- (nCurrM1 + 1):(nCurrM1 + nBirths)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[1], newNodeIDs)
        dat$nw <- add.vertices(dat$nw,
                               nv = nBirths,
                               last.mode = FALSE,
                               vertex.pid = newPids)
      } else {
        dat$nw <- add.vertices(dat$nw,
                               nv = nBirths,
                               last.mode = FALSE)
      }
      newNodes <- newNodeIDs
    }
    if (nBirthsM2 > 0) {
      newNodeIDs <- (nCurrM2 + 1):(nCurrM2 + nBirthsM2)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[2], newNodeIDs)
        dat$nw <- add.vertices(dat$nw,
                               nv = nBirthsM2,
                               last.mode = TRUE,
                               vertex.pid = newPids)
      } else {
        dat$nw <- add.vertices(dat$nw,
                               nv = nBirthsM2,
                               last.mode = TRUE)
      }
      newSize <- network.size(dat$nw)
      newNodesM2 <- (newSize - nBirthsM2 + 1):newSize
    }
    newNodes <- c(newNodes, newNodesM2)
    if (!is.null(newNodes)) {
      dat$nw <- activate.vertices(dat$nw,
                                  onset = at,
                                  terminus = Inf,
                                  v = newNodes)
    }
  }


  # Update Nodal Attributes -------------------------------------------------
  if (length(newNodes) > 0) {

    # Set attributes on nw
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_terms(form)
    curr.tab <- get_attr_prop(dat$nw, fterms)
    if (length(curr.tab) > 0) {
      dat$nw <- update_nwattr(dat$nw,
                              newNodes,
                              dat$control$attr.rules,
                              curr.tab,
                              dat$temp$t1.tab)
    }

    # Save any val on attr
    dat <- copy_toall_attr(dat, at, fterms)

    if (tea.status == TRUE) {
      if ("status" %in% fterms) {
        dat$nw <- activate.vertex.attribute(dat$nw,
                                            prefix = "testatus",
                                            value = dat$attr$status[newNodes],
                                            onset = at,
                                            terminus = Inf,
                                            v = newNodes)
      } else {
        dat$nw <- activate.vertex.attribute(dat$nw,
                                            prefix = "testatus",
                                            value = "s",
                                            onset = at,
                                            terminus = Inf,
                                            v = newNodes)
      }
    }
    if (modes == 1) {
      if (!("status" %in% fterms)) {
        dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
      }
      dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
    }
    if (modes == 2) {
      if (!("status" %in% fterms)) {
        dat <- split_bip(dat, "status", "s",
                         nCurrM1, nCurrM2, nBirths, nBirthsM2)
      }
      dat <- split_bip(dat, "active", 1,
                       nCurrM1, nCurrM2, nBirths, nBirthsM2)
      dat <- split_bip(dat, "infTime", NA,
                       nCurrM1, nCurrM2, nBirths, nBirthsM2)
    }

    ## Handles infTime when incoming nodes are infected
    newNodesInf <- intersect(newNodes, which(dat$attr$status == "s"))
    dat$attr$infTime[newNodesInf] <- at
  }


  # Output ------------------------------------------------------------------
  if (at == 2) {
    dat$epi$b.flow <- c(0, nBirths)
    if (modes == 2) {
      dat$epi$b.flow.m2 <- c(0, nBirthsM2)
    }
  } else {
    dat$epi$b.flow[at] <- nBirths
    if (modes == 2) {
      dat$epi$b.flow.m2[at] <- nBirthsM2
    }
  }

  return(dat)
}

