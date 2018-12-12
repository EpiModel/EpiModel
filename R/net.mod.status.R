
#' @title Primary Infection Module for netsim
#'
#' @description This function simulates the main infection process given the
#'              current state of the partnerships and disease in the system.
#'
#' @param dat A list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @details
#' The main steps in this infection module are as follows:
#' \enumerate{
#'  \item Get IDs for current infected and susceptibles given the current disease
#'        status.
#'  \item Call \code{\link{discord_edgelist}} to get the current discordant edgelist
#'        given step 1.
#'  \item Determine the transmission rates (e.g., as a function of mode).
#'  \item Pull the number of acts per partnership in a time step from the
#'        \code{act.rate} parameter.
#'  \item Calculate the final transmission probabilities given the transmission
#'        rates and act rates.
#'  \item Randomly transmit on the discordant edgelist.
#'  \item Conduct bookkeeping for new infections to update status on the nodes
#'        and calculate disease incidence.
#' }
#'
#' @return
#' The main \code{dat} object is returned with updated disease status and summary
#' incidence measures.
#'
#' @export
#' @keywords netMod internal
#'
#' @seealso \code{\link{discord_edgelist}} is used within \code{infection.net}
#' to obtain a discordant edgelist.
#'
infection.net <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  modes <- dat$param$modes
  mode <- idmode(dat$nw)

  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate

  nw <- dat$nw
  tea.status <- dat$control$tea.status

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- totInf <- 0


  # Process -----------------------------------------------------------------
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist
    del <- discord_edgelist(dat, at)

    # If some discordant edges, then proceed
    if (!(is.null(del))) {

      # Infection duration to at
      del$infDur <- at - dat$attr$infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
        del$transProb <- ifelse(del$infDur <= linf.prob,
                                inf.prob[del$infDur],
                                inf.prob[linf.prob])

      # Interventions
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$transProb <- del$transProb * (1 - dat$param$inter.eff)
      }

      # Calculate infection-stage act/contact rates
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate,
                            act.rate[del$infDur],
                            act.rate[lact.rate])

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate

      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      nInf <- sum(mode[idsNewInf] == 1)
      totInf <- nInf

      # Update nw attributes
      if (totInf > 0) {
        if (tea.status == TRUE) {
          nw <- activate.vertex.attribute(nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
        }
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at

        if ("status" %in% dat$temp$fterms) {
          nw <- set.vertex.attribute(nw, "status", dat$attr$status)
        }
      }

      # Substitute PIDs for vital bipartite sims
      if (any(names(nw$gal) %in% "vertex.pid")) {
        del$sus <- get.vertex.pid(nw, del$sus)
        del$inf <- get.vertex.pid(nw, del$inf)
      }

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    del <- del[!duplicated(del$sus), ]
    if (at == 2) {
      dat$stats$transmat <- del
    } else {
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
    }
  }

  ## Save incidence vector
  if (at == 2) {
    dat$epi$si.flow <- c(0, nInf)
  } else {
    dat$epi$si.flow[at] <- nInf
  }

  dat$nw <- nw
  return(dat)
}

#' @title Primary Infection Module for netsim.bip
#'
#' @description This function simulates the main infection process given the
#'              current state of the partnerships and disease in the system.
#'
#' @param dat A list object containing a \code{networkDynamic.bip} object and other
#'        initialization information passed from \code{\link{netsim.bip}}.
#' @param at Current time step.
#'
#' @details
#' The main steps in this infection module are as follows:
#' \enumerate{
#'  \item Get IDs for current infected and susceptibles given the current disease
#'        status.
#'  \item Call \code{\link{discord_edgelist.bip}} to get the current discordant edgelist
#'        given step 1.
#'  \item Determine the transmission rates (e.g., as a function of mode).
#'  \item Pull the number of acts per partnership in a time step from the
#'        \code{act.rate} parameter.
#'  \item Calculate the final transmission probabilities given the transmission
#'        rates and act rates.
#'  \item Randomly transmit on the discordant edgelist.
#'  \item Conduct bookkeeping for new infections to update status on the nodes
#'        and calculate disease incidence.
#' }
#'
#' @return
#' The main \code{dat} object is returned with updated disease status and summary
#' incidence measures.
#'
#' @export
#' @keywords netMod internal
#'
#' @seealso \code{\link{discord_edgelist.bip}} is used within \code{infection.net.bip}
#' to obtain a discordant edgelist.
#'
infection.net.bip <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  mode <- idmode(dat$nw)

  inf.prob <- dat$param$inf.prob
  inf.prob.m2 <- dat$param$inf.prob.m2
  act.rate <- dat$param$act.rate

  nw <- dat$nw
  tea.status <- dat$control$tea.status

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- nInfM2 <- totInf <- 0


  # Process -----------------------------------------------------------------
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist
    del <- discord_edgelist(dat, at)

    # If some discordant edges, then proceed
    if (!(is.null(del))) {

      # Infection duration to at
      del$infDur <- at - dat$attr$infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
      if (is.null(inf.prob.m2)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
                                inf.prob[del$infDur],
                                inf.prob[linf.prob])
      } else {
        del$transProb <- ifelse(del$sus <= nw %n% "bipartite",
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob[del$infDur],
                                       inf.prob[linf.prob]),
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob.m2[del$infDur],
                                       inf.prob.m2[linf.prob]))
      }

      # Interventions
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$transProb <- del$transProb * (1 - dat$param$inter.eff)
      }

      # Calculate infection-stage act/contact rates
      lact.rate <- length(act.rate)
      del$actRate <- ifelse(del$infDur <= lact.rate,
                            act.rate[del$infDur],
                            act.rate[lact.rate])

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - del$transProb) ^ del$actRate

      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections vector
      idsNewInf <- unique(del$sus)
      nInf <- sum(mode[idsNewInf] == 1)
      nInfM2 <- sum(mode[idsNewInf] == 2)
      totInf <- nInf + nInfM2

      # Update nw attributes
      if (totInf > 0) {
        if (tea.status == TRUE) {
          nw <- activate.vertex.attribute(nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
        }
        dat$attr$status[idsNewInf] <- "i"
        dat$attr$infTime[idsNewInf] <- at

        if ("status" %in% dat$temp$fterms) {
          nw <- set.vertex.attribute(nw, "status", dat$attr$status)
        }
      }

      # Substitute PIDs for vital bipartite sims
      if (any(names(nw$gal) %in% "vertex.pid")) {
        del$sus <- get.vertex.pid(nw, del$sus)
        del$inf <- get.vertex.pid(nw, del$inf)
      }

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    del <- del[!duplicated(del$sus), ]
    if (at == 2) {
      dat$stats$transmat <- del
    } else {
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
    }
  }

  ## Save incidence vector
  if (at == 2) {
    dat$epi$si.flow <- c(0, nInf)
    dat$epi$si.flow.m2 <- c(0, nInfM2)

  } else {
    dat$epi$si.flow[at] <- nInf
    dat$epi$si.flow.m2[at] <- nInfM2
  }

  dat$nw <- nw
  return(dat)
}

#' @title Discordant Edgelist from NetworkDynamic Object
#'
#' @description This function returns a \code{data.frame} with a discordant
#'              edgelist, defined as the set of edges in which the status of the
#'              two partners is one susceptible and one infected.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @details
#' This internal function works within the parent \code{\link{infection.net}} function
#' to pull the current edgelist from the dynamic network object, look up the disease
#' status of the head and tails on the edge, and subset the list to those edges
#' with one susceptible and one infected node.
#'
#' @return
#' This function returns a \code{data.frame} with the following columns:
#' \itemize{
#'  \item \strong{time:} time step queried
#'  \item \strong{sus:} ID number for the susceptible partner
#'  \item \strong{inf:} ID number for the infected partner
#' }
#' The output from this function is added to the transmission \code{data.frame}
#' object that is requested as output in \code{netsim} simulations with
#' the \code{save.trans=TRUE} argument.
#'
#' @seealso \code{\link{netsim}}, \code{\link{infection.net}}
#'
#' @export
#' @keywords netMod internal
#'
discord_edgelist <- function(dat, at) {

  status <- dat$attr$status
  el <- get.dyads.active(dat$nw, at = at)

  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% "i", ncol = 2)
    isSus <- matrix(stat %in% "s", ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(at, sus, inf)
    }
  }

  return(del)
}


#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to an distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              \code{\link{netsim}}.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
recovery.net <- function(dat, at) {

  ## Only run with SIR/SIS
  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  tea.status <- dat$control$tea.status

  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate

  nRecov <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)


  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1

  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    ratesElig <- rec.rate[idsElig]
  } else {
    #mElig <- mode[idsElig]
    rateseLIG <- ifelse(infDur <= lrec.rate, rec.rate[infDur], rec.rate[lrec.rate])
  }



  # Process -----------------------------------------------------------------
  if (nElig > 0) {

    if (rec.rand == TRUE) {
      vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecRecov) > 0) {
        idsRecov <- idsElig[vecRecov]
        nRecov <- sum(mode[idsRecov])
        status[idsRecov] <- recovState
        if (tea.status == TRUE) {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = recovState, onset = at,
                                              terminus = Inf, v = idsRecov)
        }
      }
    } else {
      idsRecov <- NULL
      nRecov <- min(round(sum(ratesElig)), length(idsElig))
      if (nRecov > 0) {
        idsRecov <- ssample(idsElig, nRecov)
        status[idsRecov] <- recovState
      }

      if (tea.status == TRUE & nRecov > 0) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = recovState, onset = at,
                                            terminus = Inf, v = idsRecov)
      }
    }
  }

  dat$attr$status <- status
  if ("status" %in% dat$temp$fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }

  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")

  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  } else {
    dat$epi[[outName[1]]][at] <- nRecov
  }

  return(dat)
}

#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to an distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              \code{\link{netsim}}.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
recovery.net.bip <- function(dat, at) {

  ## Only run with SIR/SIS
  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  tea.status <- dat$control$tea.status

  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.m2 <- dat$param$rec.rate.m2

  nRecov <- nRecovM2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)


  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1

  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    mElig <- mode[idsElig]
    rates <- c(rec.rate, rec.rate.m2)
    ratesElig <- rates[mElig]
  } else {
    mElig <- mode[idsElig]
    if (is.null(rec.rate.m2)) {
      #Is this equiv. to a check of whether modes = 2?
      rates <- ifelse(infDur <= lrec.rate, rec.rate[infDur], rec.rate[lrec.rate])
    } else {
      rates <- ifelse(mElig == 1, ifelse(infDur <= lrec.rate,
                                         rec.rate[infDur],
                                         rec.rate[lrec.rate]),
                      ifelse(infDur <= lrec.rate,
                             rec.rate.m2[infDur],
                             rec.rate.m2[lrec.rate]))
    }
    ratesElig <- rates
  }


  # Process -----------------------------------------------------------------
  if (nElig > 0) {

    if (rec.rand == TRUE) {
      vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecRecov) > 0) {
        idsRecov <- idsElig[vecRecov]
        nRecov <- sum(mode[idsRecov] == 1)
        nRecovM2 <- sum(mode[idsRecov] == 2)
        status[idsRecov] <- recovState
        if (tea.status == TRUE) {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = recovState, onset = at,
                                              terminus = Inf, v = idsRecov)
        }
      }
    } else {
      idsRecov <- idsRecovM2 <- NULL
      nRecov <- min(round(sum(ratesElig[mElig == 1])), sum(mElig == 1))
      if (nRecov > 0) {
        idsRecov <- ssample(idsElig[mElig == 1], nRecov)
        status[idsRecov] <- recovState
      }

      nRecovM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
      if (nRecovM2 > 0) {
        idsRecovM2 <- ssample(idsElig[mElig == 2], nRecovM2)
        status[idsRecovM2] <- recovState
      }

      totRecov <- nRecov + nRecovM2
      if (tea.status == TRUE & totRecov > 0) {
        allids <- c(idsRecov, idsRecovM2)
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = recovState, onset = at,
                                            terminus = Inf, v = allids)
      }
    }
  }

  dat$attr$status <- status
  if ("status" %in% dat$temp$fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }

  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".m2")

  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  } else {
    dat$epi[[outName[1]]][at] <- nRecov
  }
  if (at == 2) {
    dat$epi[[outName[2]]] <- c(0, nRecovM2)
  } else {
    dat$epi[[outName[2]]][at] <- nRecovM2
  }


  return(dat)
}
