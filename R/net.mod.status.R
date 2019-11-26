
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
#'  \item Determine the transmission rates (e.g., as a function of group).
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

  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate

  nw <- dat$nw
  tea.status <- dat$control$tea.status

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- 0



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
      nInf <- length(idsNewInf)

      #Output to nw.update
      dat$nw.update$inf$idsNewInf <- idsNewInf

      # Substitute PIDs for vital two-group sims
      if (any(names(nw$gal) %in% "vertex.pid")) {
        del$sus <- get.vertex.pid(nw, del$sus)
        del$inf <- get.vertex.pid(nw, del$inf)
      }

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix

  if (nInf > 0) {
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
#'  \item Determine the transmission rates (e.g., as a function of group).
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

infection.net.grp <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  nw <- dat$nw
  group <- get.vertex.attribute(nw, "group")

  inf.prob <- dat$param$inf.prob
  inf.prob.g2 <- dat$param$inf.prob.g2
  act.rate <- dat$param$act.rate

  tea.status <- dat$control$tea.status

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- nInfG2 <- totInf <- 0


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
      if (is.null(inf.prob.g2)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
                                inf.prob[del$infDur],
                                inf.prob[linf.prob])
      } else {
        #FLAG
        del$transProb <- ifelse(del$sus <= nw %n% "bipartite",
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob[del$infDur],
                                       inf.prob[linf.prob]),
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob.g2[del$infDur],
                                       inf.prob.g2[linf.prob]))
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
      nInf <- sum(group[idsNewInf] == 1)
      nInfG2 <- sum(group[idsNewInf] == 2)

      #Out to network upate
      dat$nw.update$inf$nInf <- nInf + nInfG2
      dat$nw.update$inf$idsNewInf

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
    dat$epi$si.flow.g2 <- c(0, nInfG2)

  } else {
    dat$epi$si.flow[at] <- nInf
    dat$epi$si.flow.g2[at] <- nInfG2
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

  rec.rate <- dat$param$rec.rate

  nRecov <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)


  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1
  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    ratesElig <- rec.rate
  } else {
    ratesElig <- ifelse(infDur <= lrec.rate, rec.rate[infDur], rec.rate[lrec.rate])
  }


  # Process -----------------------------------------------------------------
  if (nElig > 0) {
    vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecRecov) > 0) {
      idsRecov <- idsElig[vecRecov]
      nRecov <- length(idsRecov)
      status[idsRecov] <- recovState
      dat$nw.update$rec$idsRecov <- idsRecov
      dat$nw.update$rec$recovState <- recovState
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
recovery.net.grp <- function(dat, at) {

  ## Only run with SIR/SIS
  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime
  tea.status <- dat$control$tea.status

  group <- get.vertex.attribute(dat$nw, "group")

  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rate <- dat$param$rec.rate
  rec.rate.g2 <- dat$param$rec.rate.g2

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1
  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    gElig <- group[idsElig]
    rates <- c(rec.rate, rec.rate.g2)
    ratesElig <- rates[gElig]
  } else {
    gElig <- group[idsElig]
    if (is.null(rec.rate.g2)) {
      rates <- ifelse(infDur <= lrec.rate, rec.rate[infDur], rec.rate[lrec.rate])
    } else {
      rates <- rep(NA, length(infDur))
      rates[gElig == 1] <- ifelse(infDur[gElig == 1] <= lrec.rate,
                                  rec.rate[infDur[gElig == 1]], rec.rate[lrec.rate])
      rates[gElig == 2] <- ifelse(infDur[gElig == 2] <= length(rec.rate.g2),
                                  rec.rate.g2[infDur[gElig == 2]], rec.rate.g2[length(rec.rate.g2)])
    }
    ratesElig <- rates
  }


  # Process -----------------------------------------------------------------
  if (nElig > 0) {
    vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecRecov) > 0) {
      dat$nw.update$rec$idsRecov <- idsRecov <- idsElig[vecRecov]
      nRecov <- sum(group[idsRecov] == 1)
      nRecovG2 <- sum(group[idsRecov] == 2)
      dat$nw.update$rec$recovState <- recovState <- status[idsRecov]
    }
  }

  dat$attr$status <- status
  if ("status" %in% dat$temp$fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }

  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".g2")

  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  } else {
    dat$epi[[outName[1]]][at] <- nRecov
  }
  if (at == 2) {
    dat$epi[[outName[2]]] <- c(0, nRecovG2)
  } else {
    dat$epi[[outName[2]]][at] <- nRecovG2
  }


  return(dat)
}
