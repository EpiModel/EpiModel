
#' @title Deaths: netsim Module
#'
#' @description This function simulates death for use in \link{netsim} simulations.
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

  # Initialize counts and pull rates
  nDeaths.sus <- nDeathsM2.sus <- 0
  idsElig.sus <- which(dat$attr$active == 1 & dat$attr$status == "s")
  nElig.sus <- length(idsElig.sus)

  if (nElig.sus > 0) {

    # Mortality rates by mode
    mElig.sus <- mode[idsElig.sus]
    rates.sus <- c(dat$param$ds.rate, dat$param$ds.rate.m2)
    ratesElig.sus <- rates.sus[mElig.sus]

    # Stochastic exits
    if (d.rand == TRUE) {
      vecDeaths.sus <- which(rbinom(nElig.sus, 1, ratesElig.sus) == 1)
      if (length(vecDeaths.sus) > 0) {
        idsDth.sus <- idsElig.sus[vecDeaths.sus]
        nDeaths.sus <- sum(mode[idsDth.sus] == 1)
        nDeathsM2.sus <- sum(mode[idsDth.sus] == 2)
        dat$attr$active[idsDth.sus] <- 0
        dat$attr$exitTime[idsDth.sus] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDth.sus, deactivate.edges = TRUE)
      }
    }

    # Deterministic exits
    if (d.rand == FALSE) {
      idsDth.sus <- idsDthM2.sus <- NULL
      nDeaths.sus <- min(round(sum(ratesElig.sus[mElig.sus == 1])), sum(mElig.sus == 1))
      idsDth.sus <- ssample(idsElig.sus[mElig.sus == 1], nDeaths.sus)
      if (modes == 2) {
        nDeathsM2.sus <- min(round(sum(ratesElig.sus[mElig.sus == 2])), sum(mElig.sus == 2))
        idsDthM2.sus <- ssample(idsElig.sus[mElig.sus == 2], nDeathsM2.sus)
      }
      totDth.sus <- nDeaths.sus + nDeathsM2.sus
      if (totDth.sus > 0) {
        idsDthAll.sus <- c(idsDth.sus, idsDthM2.sus)
        dat$attr$active[idsDthAll.sus] <- 0
        dat$attr$exitTime[idsDthAll.sus] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDthAll.sus, deactivate.edges = TRUE)
      }
    }
  }


  # Infected deaths ---------------------------------------------------------

  # Initialize counts and query rates
  nDeaths.inf <- nDeathsM2.inf <- 0
  idsElig.inf <- which(dat$attr$active == 1 & dat$attr$status == "i")
  nElig.inf <- length(idsElig.inf)

  if (nElig.inf > 0) {

    # Mortality rates by mode
    mElig.inf <- mode[idsElig.inf]
    rates.inf <- c(dat$param$di.rate, dat$param$di.rate.m2)
    ratesElig.inf <- rates.inf[mElig.inf]

    # Stochastic exits
    if (d.rand == TRUE) {
      vecDeaths.inf <- which(rbinom(nElig.inf, 1, ratesElig.inf) == 1)
      if (length(vecDeaths.inf) > 0) {
        idsDth.inf <- idsElig.inf[vecDeaths.inf]
        nDeaths.inf <- sum(mode[idsDth.inf] == 1)
        nDeathsM2.inf <- sum(mode[idsDth.inf] == 2)
        dat$attr$active[idsDth.inf] <- 0
        dat$attr$exitTime[idsDth.inf] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDth.inf, deactivate.edges = TRUE)
      }
    }

    # Deterministic exits
    if (d.rand == FALSE) {
      idsDth.inf <- idsDthM2.inf <- NULL
      nDeaths.inf <- min(round(sum(ratesElig.inf[mElig.inf == 1])), sum(mElig.inf == 1))
      idsDth.inf <- ssample(idsElig.inf[mElig.inf == 1], nDeaths.inf)
      if (modes == 2) {
        nDeathsM2.inf <- min(round(sum(ratesElig.inf[mElig.inf == 2])), sum(mElig.inf == 2))
        idsDthM2.inf <- ssample(idsElig.inf[mElig.inf == 2], nDeathsM2.inf)
      }
      totDth.inf <- nDeaths.inf + nDeathsM2.inf
      if (totDth.inf > 0) {
        idsDthAll.inf <- c(idsDth.inf, idsDthM2.inf)
        dat$attr$active[idsDthAll.inf] <- 0
        dat$attr$exitTime[idsDthAll.inf] <- at
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDthAll.inf, deactivate.edges = TRUE)
      }
    }
  }


  # Recovered deaths --------------------------------------------------------
  if (type == "SIR") {

    # Initialize counts and query rates
    nDeaths.rec <- nDeathsM2.rec <- 0
    idsElig.rec <- which(dat$attr$active == 1 & dat$attr$status == "r")
    nElig.rec <- length(idsElig.rec)

    if (nElig.rec > 0) {

      # Mortality rates by mode
      mElig.rec <- mode[idsElig.rec]
      rates.rec <- c(dat$param$dr.rate, dat$param$dr.rate.m2)
      ratesElig.rec <- rates.rec[mElig.rec]

      # Stochastic exits
      if (d.rand == TRUE) {
        vecDeaths.rec <- which(rbinom(nElig.rec, 1, ratesElig.rec) == 1)
        if (length(vecDeaths.rec) > 0) {
          idsDth.rec <- idsElig.rec[vecDeaths.rec]
          nDeaths.rec <- sum(mode[idsDth.rec] == 1)
          nDeathsM2.rec <- sum(mode[idsDth.rec] == 2)
          dat$attr$active[idsDth.rec] <- 0
          dat$attr$exitTime[idsDth.rec] <- at
          dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                        v = idsDth.rec, deactivate.edges = TRUE)
        }
      }

      # Deterministic exits
      if (d.rand == FALSE) {
        idsDth.rec <- idsDthM2.rec <- NULL
        nDeaths.rec <- min(round(sum(ratesElig.rec[mElig.rec == 1])), sum(mElig.rec == 1))
        idsDth.rec <- ssample(idsElig.rec[mElig.rec == 1], nDeaths.rec)
        if (modes == 2) {
          nDeathsM2.rec <- min(round(sum(ratesElig.rec[mElig.rec == 2])), sum(mElig.rec == 2))
          idsDthM2.rec <- ssample(idsElig.rec[mElig.rec == 2], nDeathsM2.rec)
        }
        totDth.rec <- nDeaths.rec + nDeathsM2.rec
        if (totDth.rec > 0) {
          idsDthAll.rec <- c(idsDth.rec, idsDthM2.rec)
          dat$attr$active[idsDthAll.rec] <- 0
          dat$attr$exitTime[idsDthAll.rec] <- at
          dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                        v = idsDthAll.rec, deactivate.edges = TRUE)
        }
      }
    }
  }


  # Output ------------------------------------------------------------------

  if (at == 2) {
    dat$epi$ds.flow <- c(0, nDeaths.sus)
    dat$epi$di.flow <- c(0, nDeaths.inf)
    if (type == "SIR") {
      dat$epi$dr.flow <- c(0, nDeaths.rec)
    }
    if (modes == 2) {
      dat$epi$ds.flow.m2 <- c(0, nDeathsM2.sus)
      dat$epi$di.flow.m2 <- c(0, nDeathsM2.inf)
      if (type == "SIR") {
        dat$epi$dr.flow.m2 <- c(0, nDeathsM2.rec)
      }
    }
  } else {
    dat$epi$ds.flow[at] <- nDeaths.sus
    dat$epi$di.flow[at] <- nDeaths.inf
    if (type == "SIR") {
      dat$epi$dr.flow[at] <- nDeaths.rec
    }
    if (modes == 2) {
      dat$epi$ds.flow.m2[at] <- nDeathsM2.sus
      dat$epi$di.flow.m2[at] <- nDeathsM2.inf
      if (type == "SIR") {
        dat$epi$dr.flow.m2[at] <- nDeathsM2.rec
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
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
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
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
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
        dat <- split_bip(dat, "status", "s", nCurrM1, nCurrM2, nBirths, nBirthsM2)
      }
      dat <- split_bip(dat, "active", 1, nCurrM1, nCurrM2, nBirths, nBirthsM2)
      dat <- split_bip(dat, "infTime", NA, nCurrM1, nCurrM2, nBirths, nBirthsM2)
      dat <- split_bip(dat, "entrTime", at, nCurrM1, nCurrM2, nBirths, nBirthsM2)
      dat <- split_bip(dat, "exitTime", NA, nCurrM1, nCurrM2, nBirths, nBirthsM2)
    }

    ## Handles infTime when incoming nodes are infected
    newNodesInf <- intersect(newNodes, which(dat$attr$status == "s"))
    dat$attr$infTime[newNodesInf] <- at

    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check births.net module.")
    }

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

