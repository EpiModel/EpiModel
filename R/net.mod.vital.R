
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
  modes <- dat$param$modes
  mode <- idmode(dat$nw)

  type <- dat$control$type
  d.rand <- dat$control$d.rand


  # Susceptible departures ------------------------------------------------------

  # Initialize counts and pull rates
  nDepartures.sus <- nDeparturesM2.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)

  if (nElig.sus > 0) {

    # Departure rates by mode
    mElig.sus <- mode[idsElig.sus]
    rates.sus <- c(dat$param$ds.rate, dat$param$ds.rate.m2)
    ratesElig.sus <- rates.sus[mElig.sus]

    # Stochastic exits
    if (d.rand == TRUE) {
      vecDepartures.sus <- which(rbinom(nElig.sus, 1, ratesElig.sus) == 1)
      if (length(vecDepartures.sus) > 0) {
        idsDpt.sus <- idsElig.sus[vecDepartures.sus]
        nDepartures.sus <- sum(mode[idsDpt.sus] == 1)
        nDeparturesM2.sus <- sum(mode[idsDpt.sus] == 2)
        dat$attr$active[idsDpt.sus] <- 0
        dat$attr$exitTime[idsDpt.sus] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDpt.sus, deactivate.edges = TRUE)
      }
    }

    # Deterministic exits
    if (d.rand == FALSE) {
      idsDpt.sus <- idsDptM2.sus <- NULL
      nDepartures.sus <- min(round(sum(ratesElig.sus[mElig.sus == 1])), sum(mElig.sus == 1))
      idsDpt.sus <- ssample(idsElig.sus[mElig.sus == 1], nDepartures.sus)
      if (modes == 2) {
        nDeparturesM2.sus <- min(round(sum(ratesElig.sus[mElig.sus == 2])), sum(mElig.sus == 2))
        idsDptM2.sus <- ssample(idsElig.sus[mElig.sus == 2], nDeparturesM2.sus)
      }
      totDpt.sus <- nDepartures.sus + nDeparturesM2.sus
      if (totDpt.sus > 0) {
        idsDptAll.sus <- c(idsDpt.sus, idsDptM2.sus)
        dat$attr$active[idsDptAll.sus] <- 0
        dat$attr$exitTime[idsDptAll.sus] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDptAll.sus, deactivate.edges = TRUE)
      }
    }
  }


  # Infected departures ---------------------------------------------------------

  # Initialize counts and query rates
  nDepartures.inf <- nDeparturesM2.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)

  if (nElig.inf > 0) {

    # Departure rates by mode
    mElig.inf <- mode[idsElig.inf]
    rates.inf <- c(dat$param$di.rate, dat$param$di.rate.m2)
    ratesElig.inf <- rates.inf[mElig.inf]

    # Stochastic exits
    if (d.rand == TRUE) {
      vecDepartures.inf <- which(rbinom(nElig.inf, 1, ratesElig.inf) == 1)
      if (length(vecDepartures.inf) > 0) {
        idsDpt.inf <- idsElig.inf[vecDepartures.inf]
        nDepartures.inf <- sum(mode[idsDpt.inf] == 1)
        nDeparturesM2.inf <- sum(mode[idsDpt.inf] == 2)
        dat$attr$active[idsDpt.inf] <- 0
        dat$attr$exitTime[idsDpt.inf] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDpt.inf, deactivate.edges = TRUE)
      }
    }

    # Deterministic exits
    if (d.rand == FALSE) {
      idsDpt.inf <- idsDptM2.inf <- NULL
      nDepartures.inf <- min(round(sum(ratesElig.inf[mElig.inf == 1])), sum(mElig.inf == 1))
      idsDpt.inf <- ssample(idsElig.inf[mElig.inf == 1], nDepartures.inf)
      if (modes == 2) {
        nDeparturesM2.inf <- min(round(sum(ratesElig.inf[mElig.inf == 2])), sum(mElig.inf == 2))
        idsDptM2.inf <- ssample(idsElig.inf[mElig.inf == 2], nDeparturesM2.inf)
      }
      totDpt.inf <- nDepartures.inf + nDeparturesM2.inf
      if (totDpt.inf > 0) {
        idsDptAll.inf <- c(idsDpt.inf, idsDptM2.inf)
        dat$attr$active[idsDptAll.inf] <- 0
        dat$attr$exitTime[idsDptAll.inf] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDptAll.inf, deactivate.edges = TRUE)
      }
    }
  }


  # Recovered departures --------------------------------------------------------
  if (type == "SIR") {

    # Initialize counts and query rates
    nDepartures.rec <- nDeparturesM2.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)

    if (nElig.rec > 0) {

      # Departure rates by mode
      mElig.rec <- mode[idsElig.rec]
      rates.rec <- c(dat$param$dr.rate, dat$param$dr.rate.m2)
      ratesElig.rec <- rates.rec[mElig.rec]

      # Stochastic exits
      if (d.rand == TRUE) {
        vecDepartures.rec <- which(rbinom(nElig.rec, 1, ratesElig.rec) == 1)
        if (length(vecDepartures.rec) > 0) {
          idsDpt.rec <- idsElig.rec[vecDepartures.rec]
          nDepartures.rec <- sum(mode[idsDpt.rec] == 1)
          nDeparturesM2.rec <- sum(mode[idsDpt.rec] == 2)
          dat$attr$active[idsDpt.rec] <- 0
          dat$attr$exitTime[idsDpt.rec] <- at
          dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                        v = idsDpt.rec, deactivate.edges = TRUE)
        }
      }

      # Deterministic exits
      if (d.rand == FALSE) {
        idsDpt.rec <- idsDptM2.rec <- NULL
        nDepartures.rec <- min(round(sum(ratesElig.rec[mElig.rec == 1])), sum(mElig.rec == 1))
        idsDpt.rec <- ssample(idsElig.rec[mElig.rec == 1], nDepartures.rec)
        if (modes == 2) {
          nDeparturesM2.rec <- min(round(sum(ratesElig.rec[mElig.rec == 2])), sum(mElig.rec == 2))
          idsDptM2.rec <- ssample(idsElig.rec[mElig.rec == 2], nDeparturesM2.rec)
        }
        totDpt.rec <- nDepartures.rec + nDeparturesM2.rec
        if (totDpt.rec > 0) {
          idsDptAll.rec <- c(idsDpt.rec, idsDptM2.rec)
          dat$attr$active[idsDptAll.rec] <- 0
          dat$attr$exitTime[idsDptAll.rec] <- at
          dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                        v = idsDptAll.rec, deactivate.edges = TRUE)
        }
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
    if (modes == 2) {
      dat$epi$ds.flow.m2 <- c(0, nDeparturesM2.sus)
      dat$epi$di.flow.m2 <- c(0, nDeparturesM2.inf)
      if (type == "SIR") {
        dat$epi$dr.flow.m2 <- c(0, nDeparturesM2.rec)
      }
    }
  } else {
    dat$epi$ds.flow[at] <- nDepartures.sus
    dat$epi$di.flow[at] <- nDepartures.inf
    if (type == "SIR") {
      dat$epi$dr.flow[at] <- nDepartures.rec
    }
    if (modes == 2) {
      dat$epi$ds.flow.m2[at] <- nDeparturesM2.sus
      dat$epi$di.flow.m2[at] <- nDeparturesM2.inf
      if (type == "SIR") {
        dat$epi$dr.flow.m2[at] <- nDeparturesM2.rec
      }
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
  a.rate.m2 <- dat$param$a.rate.m2
  modes <- dat$param$modes
  tea.status <- dat$control$tea.status
  nOld <- dat$epi$num[at - 1]
  nCurr <- network.size(dat$nw)
  a.rand <- dat$control$a.rand
  delete.nodes <- dat$control$delete.nodes

  nArrivals <- nArrivalsM2 <- 0
  newNodes <- newNodesM2 <- NULL


  # Add Nodes ---------------------------------------------------------------
  if (modes == 1 && nOld > 0) {
    if (a.rand == TRUE) {
      nArrivals <- sum(rbinom(nOld, 1, a.rate))
    }
    if (a.rand == FALSE) {
      nArrivals <- round(nOld * a.rate)
    }
    if (nArrivals > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nArrivals)
      newNodes <- (nCurr + 1):(nCurr + nArrivals)
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
    }
  }
  if (modes == 2 && nOld > 0) {
    nOldM2 <- dat$epi$num.m2[at - 1]
    if (a.rand == TRUE) {
      if (is.na(a.rate.m2)) {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivalsM2 <- sum(rbinom(nOld, 1, a.rate))
      } else {
        nArrivals <- sum(rbinom(nOld, 1, a.rate))
        nArrivalsM2 <- sum(rbinom(nOldM2, 1, a.rate.m2))
      }
    }
    if (a.rand == FALSE) {
      if (is.na(a.rate.m2)) {
        nArrivals <- round(nOld * a.rate)
        nArrivalsM2 <- round(nOld * a.rate)
      } else {
        nArrivals <- round(nOld * a.rate)
        nArrivalsM2 <- round(nOldM2 * a.rate.m2)
      }
    }

    nCurrM1 <- length(modeids(dat$nw, 1))
    nCurrM2 <- length(modeids(dat$nw, 2))
    prefixes <- unique(substr(dat$nw %v% "vertex.names", 1, 1))

    if (nArrivals > 0) {
      newNodeIDs <- (nCurrM1 + 1):(nCurrM1 + nArrivals)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[1], newNodeIDs)
        dat$nw <- add.vertices(dat$nw,
                               nv = nArrivals,
                               last.mode = FALSE,
                               vertex.pid = newPids)
      } else {
        dat$nw <- add.vertices(dat$nw,
                               nv = nArrivals,
                               last.mode = FALSE)
      }
      newNodes <- newNodeIDs
    }
    if (nArrivalsM2 > 0) {
      newNodeIDs <- (nCurrM2 + 1):(nCurrM2 + nArrivalsM2)
      if (delete.nodes == FALSE) {
        newPids <- paste0(prefixes[2], newNodeIDs)
        dat$nw <- add.vertices(dat$nw,
                               nv = nArrivalsM2,
                               last.mode = TRUE,
                               vertex.pid = newPids)
      } else {
        dat$nw <- add.vertices(dat$nw,
                               nv = nArrivalsM2,
                               last.mode = TRUE)
      }
      newSize <- network.size(dat$nw)
      newNodesM2 <- (newSize - nArrivalsM2 + 1):newSize
    }
    newNodes <- c(newNodes, newNodesM2)
    if (!is.null(newNodes)) {
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
    }
  }


  # Update Nodal Attributes -------------------------------------------------
  if (length(newNodes) > 0) {

    # Set attributes on nw
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_terms(form)
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
    if (modes == 1) {
      if (!("status" %in% fterms)) {
        dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
      }
      dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
      dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
      dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))
    }
    if (modes == 2) {
      if (!("status" %in% fterms)) {
        dat <- split_bip(dat, "status", "s", nCurrM1, nCurrM2, nArrivals, nArrivalsM2)
      }
      dat <- split_bip(dat, "active", 1, nCurrM1, nCurrM2, nArrivals, nArrivalsM2)
      dat <- split_bip(dat, "infTime", NA, nCurrM1, nCurrM2, nArrivals, nArrivalsM2)
      dat <- split_bip(dat, "entrTime", at, nCurrM1, nCurrM2, nArrivals, nArrivalsM2)
      dat <- split_bip(dat, "exitTime", NA, nCurrM1, nCurrM2, nArrivals, nArrivalsM2)
    }

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
    if (modes == 2) {
      dat$epi$a.flow.m2 <- c(0, nArrivalsM2)
    }
  } else {
    dat$epi$a.flow[at] <- nArrivals
    if (modes == 2) {
      dat$epi$a.flow.m2[at] <- nArrivalsM2
    }
  }

  return(dat)
}

