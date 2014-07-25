
#' @title Deaths: netsim Module
#'
#' @description This function simulates death for use in \link{netsim} simulations.
#'
#' @param all a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
deaths.net <- function(all, at) {

  # Conditions --------------------------------------------------------------
  if (all$param$vital == FALSE) {
    return(all)
  }


  # Variables ---------------------------------------------------------------
  modes <- all$param$modes
  mode <- idmode(all$nw)

  type <- all$control$type
  d.rand <- all$control$d.rand


  # Susceptible deaths ------------------------------------------------------

  # Initialize counts and query rates
  nDeaths <- nDeathsM2 <- 0
  idsElig <- which(all$attr$active == 1 & all$attr$status == 0)
  nElig <- length(idsElig)

  if (nElig > 0) {

    # Pull rates by mode
    mElig <- mode[idsElig]
    rates <- c(all$param$ds.rate, all$param$ds.rate.m2)
    ratesElig <- rates[mElig]

    # Stochastic deaths
    if (d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(mode[idsDth] == 1)
        nDeathsM2 <- sum(mode[idsDth] == 2)
        all$attr$active[idsDth] <- 0
        all$nw <- deactivate.vertices(all$nw,
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
      all$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
      if (modes == 2) {
        nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
        all$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
      }
      totDth <- nDeaths + nDeathsM2
      if (totDth > 0) {
        allids <- c(idsDth, idsDthM2)
        all$nw <- deactivate.vertices(all$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = allids,
                                      deactivate.edges = TRUE)
      }
    }
  }

  # Output
  if (at == 2) {
    all$out$ds.flow <- c(0, nDeaths)
    if (modes == 2) {
      all$out$ds.flow.m2 <- c(0, nDeathsM2)
    }
  } else {
    all$out$ds.flow[at] <- nDeaths
    if (modes == 2) {
      all$out$ds.flow.m2[at] <- nDeathsM2
    }
  }


  # Infected deaths ---------------------------------------------------------

  # Initialize counts and query rates
  nDeaths <- nDeathsM2 <- 0
  idsElig <- which(all$attr$active == 1 & all$attr$status == 1)
  nElig <- length(idsElig)

  if (nElig > 0) {

    # Pull rates by mode
    mElig <- mode[idsElig]
    rates <- c(all$param$di.rate, all$param$di.rate.m2)
    ratesElig <- rates[mElig]

    # Stochastic deaths
    if (d.rand == TRUE) {
      vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecDeaths) > 0) {
        idsDth <- idsElig[vecDeaths]
        nDeaths <- sum(mode[idsDth] == 1)
        nDeathsM2 <- sum(mode[idsDth] == 2)
        all$attr$active[idsDth] <- 0
        all$nw <- deactivate.vertices(all$nw,
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
      all$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
      if (modes == 2) {
        nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
        all$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
      }
      totDth <- nDeaths + nDeathsM2
      if (totDth > 0) {
        allids <- c(idsDth, idsDthM2)
        all$nw <- deactivate.vertices(all$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = allids,
                                      deactivate.edges = TRUE)
      }
    }
  }

  # Output
  if (at == 2) {
    all$out$di.flow <- c(0, nDeaths)
    if (modes == 2) {
      all$out$di.flow.m2 <- c(0, nDeathsM2)
    }
  } else {
    all$out$di.flow[at] <- nDeaths
    if (modes == 2) {
      all$out$di.flow.m2[at] <- nDeathsM2
    }
  }


  # Recovered deaths --------------------------------------------------------
  if (type == "SIR") {

    # Initialize counts and query rates
    nDeaths <- nDeathsM2 <- 0
    idsElig <- which(all$attr$active == 1 & all$attr$status == 2)
    nElig <- length(idsElig)

    if (nElig > 0) {

      # Pull rates by mode
      mElig <- mode[idsElig]
      rates <- c(all$param$dr.rate, all$param$dr.rate.m2)
      ratesElig <- rates[mElig]

      # Stochastic deaths
      if (d.rand == TRUE) {
        vecDeaths <- which(rbinom(nElig, 1, ratesElig) == 1)
        if (length(vecDeaths) > 0) {
          idsDth <- idsElig[vecDeaths]
          nDeaths <- sum(mode[idsDth] == 1)
          nDeathsM2 <- sum(mode[idsDth] == 2)
          all$attr$active[idsDth] <- 0
          all$nw <- deactivate.vertices(all$nw,
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
        all$attr$active[ssample(idsElig[mElig == 1], nDeaths)] <- 0
        if (modes == 2) {
          nDeathsM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
          all$attr$active[ssample(idsElig[mElig == 2], nDeaths)] <- 0
        }
        totDth <- nDeaths + nDeathsM2
        if (totDth > 0) {
          allids <- c(idsDth, idsDthM2)
          all$nw <- deactivate.vertices(all$nw,
                                        onset = at,
                                        terminus = Inf,
                                        v = allids,
                                        deactivate.edges = TRUE)
        }
      }
    }

    # Output
    if (at == 2) {
      all$out$dr.flow <- c(0, nDeaths)
      if (modes == 2) {
        all$out$dr.flow.m2 <- c(0, nDeathsM2)
      }
    } else {
      all$out$dr.flow[at] <- nDeaths
      if (modes == 2) {
        all$out$dr.flow.m2[at] <- nDeathsM2
      }
    }
  }


  return(all)
}



#' @title Births: netsim Module
#'
#' @description This function simulates new births into the network
#'   for use in \code{\link{netsim}} simulations.
#'
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
births.net <- function(all, at) {

  # Conditions --------------------------------------------------------------
  if (all$param$vital == FALSE) {
    return(all)
  }


  # Variables ---------------------------------------------------------------
  b.rate <- all$param$b.rate
  b.rate.m2 <- all$param$b.rate.m2
  formation <- all$nwparam$formation
  modes <- all$param$modes
  tea.status <- all$control$tea.status
  nOld <- all$out$num[at - 1]
  nCurr <- network.size(all$nw)
  b.rand <- all$control$b.rand
  delete.nodes <- all$control$delete.nodes

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
      all$nw <- add.vertices(all$nw,
                             nv = nBirths)
      newNodes <- (nCurr + 1):(nCurr + nBirths)
      all$nw <- activate.vertices(all$nw,
                                  onset = at,
                                  terminus = Inf,
                                  v = newNodes)
    }
  }
  if (modes == 2 && nOld > 0) {
    nOldM2 <- all$out$num.m2[at - 1]
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

    nCurrM1 <- length(modeids(all$nw, 1))
    nCurrM2 <- length(modeids(all$nw, 2))
    prefixes <- unique(substr(all$nw %v% "vertex.names", 1, 1))

    if (nBirths > 0) {
      newNodeIDs <- (nCurrM1 + 1):(nCurrM1 + nBirths)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[1], newNodeIDs)
        all$nw <- add.vertices(all$nw,
                               nv = nBirths,
                               last.mode = FALSE,
                               vertex.pid = newPids)
      } else {
        all$nw <- add.vertices(all$nw,
                               nv = nBirths,
                               last.mode = FALSE)
      }
      newNodes <- newNodeIDs
    }
    if (nBirthsM2 > 0) {
      newNodeIDs <- (nCurrM2 + 1):(nCurrM2 + nBirthsM2)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[2], newNodeIDs)
        all$nw <- add.vertices(all$nw,
                               nv = nBirthsM2,
                               last.mode = TRUE,
                               vertex.pid = newPids)
      } else {
        all$nw <- add.vertices(all$nw,
                               nv = nBirthsM2,
                               last.mode = TRUE)
      }
      newSize <- network.size(all$nw)
      newNodesM2 <- (newSize - nBirthsM2 + 1):newSize
    }
    newNodes <- c(newNodes, newNodesM2)
    if (!is.null(newNodes)) {
      all$nw <- activate.vertices(all$nw,
                                  onset = at,
                                  terminus = Inf,
                                  v = newNodes)
    }
  }


  # Update Nodal Attributes -------------------------------------------------
  if (length(newNodes) > 0) {

    # Set attributes on nw
    t <- get_formula_terms(all$nwparam$formation)
    curr.tab <- get_attr_prop(all$nw, t)
    if (length(curr.tab) > 0) {
      all$nw <- update_nwattr(all$nw,
                              newNodes,
                              all$control$attr.rules,
                              curr.tab,
                              all$temp$t1.tab)
    }

    # Save any val on attr
    all <- copy_toall_attr(all, at)

    if (tea.status == TRUE) {
      if ("status" %in% t) {
        all$nw <- activate.vertex.attribute(all$nw,
                                            prefix = "testatus",
                                            value = all$attr$status[newNodes],
                                            onset = at,
                                            terminus = Inf,
                                            v = newNodes)
      } else {
        all$nw <- activate.vertex.attribute(all$nw,
                                            prefix = "testatus",
                                            value = 0,
                                            onset = at,
                                            terminus = Inf,
                                            v = newNodes)
      }
    }
    if (modes == 1) {
      if (!("status" %in% t)) {
        all$attr$status <- c(all$attr$status, rep(0, length(newNodes)))
      }
      all$attr$active <- c(all$attr$active, rep(1, length(newNodes)))
      all$attr$infTime <- c(all$attr$infTime, rep(NA, length(newNodes)))
    }
    if (modes == 2) {
      if (!("status" %in% t)) {
        all <- split_bip(all, "status", 0,
                         nCurrM1, nCurrM2, nBirths, nBirthsM2)
      }
      all <- split_bip(all, "active", 1,
                       nCurrM1, nCurrM2, nBirths, nBirthsM2)
      all <- split_bip(all, "infTime", NA,
                       nCurrM1, nCurrM2, nBirths, nBirthsM2)
    }

    ## Handles infTime when incoming nodes are infected
    newNodesInf <- intersect(newNodes, which(all$attr$status == 1))
    all$attr$infTime[newNodesInf] <- at
  }


  # Output ------------------------------------------------------------------
  if (at == 2) {
    all$out$b.flow <- c(0, nBirths)
    if (modes == 2) {
      all$out$b.flow.m2 <- c(0, nBirthsM2)
    }
  } else {
    all$out$b.flow[at] <- nBirths
    if (modes == 2) {
      all$out$b.flow.m2[at] <- nBirthsM2
    }
  }

  return(all)
}

