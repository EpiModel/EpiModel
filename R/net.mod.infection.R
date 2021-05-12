
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
#'  \item Get IDs for current infected and susceptibles given the current
#'        disease status.
#'  \item Call \code{\link{discord_edgelist}} to get the current discordant
#'        edgelist given step 1.
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
#' The main \code{dat} object is returned with updated disease status and
#' summary incidence measures.
#'
#' @export
#' @keywords netMod internal
#'
#' @seealso \code{\link{discord_edgelist}} is used within \code{infection.net}
#' to obtain a discordant edgelist.
#'
infection.net <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)

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
    del <- discord_edgelist(dat, at, network = 1)

    # If some discordant edges, then proceed
    if (!(is.null(del))) {

      # Infection duration to at
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
      del$transProb <- ifelse(del$infDur <= linf.prob,
                              inf.prob[del$infDur],
                              inf.prob[linf.prob])

      # Interventions
      if (!is.null(inter.eff) && at >= inter.start) {
        del$transProb <- del$transProb * (1 - inter.eff)
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
      status <- get_attr(dat, "status")
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- length(idsNewInf)

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix

  if (nInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow", at, nInf)

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
#'  \item Get IDs for current infected and susceptibles given the current
#'        disease status.
#'  \item Call \code{\link{discord_edgelist}} to get the current discordant
#'        edgelist given step 1.
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
#' The main \code{dat} object is returned with updated disease status and
#' summary incidence measures.
#'
#' @export
#' @keywords netMod internal
#'
#' @seealso \code{\link{discord_edgelist}} is used within \code{infection.net}
#' to obtain a discordant edgelist.
#'
infection.2g.net <- function(dat, at) {

  # Variables ---------------------------------------------------------------

  active <- get_attr(dat, "active")
  infTime <- get_attr(dat, "infTime")
  status <- get_attr(dat, "status")
  group <- get_attr(dat, "group")

  inf.prob <- get_param(dat, "inf.prob")
  inf.prob.g2 <- get_param(dat, "inf.prob.g2")
  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)


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
    del <- discord_edgelist(dat, at, network = 1)

    # If some discordant edges, then proceed
    if (!(is.null(del))) {

      # Infection duration to at
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf.prob)
      if (is.null(inf.prob.g2)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
                                inf.prob[del$infDur],
                                inf.prob[linf.prob])
      } else {
        #FLAG
        del$transProb <- ifelse(group[del$sus] == 1,
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob[del$infDur],
                                       inf.prob[linf.prob]),
                                ifelse(del$infDur <= linf.prob,
                                       inf.prob.g2[del$infDur],
                                       inf.prob.g2[linf.prob]))
      }

      # Interventions
      if (!is.null(inter.eff) && at >= inter.start) {
        del$transProb <- del$transProb * (1 - inter.eff)
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
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)
      nInf <- sum(group[idsNewInf] == 1)
      nInfG2 <- sum(group[idsNewInf] == 2)
      totInf <- nInf + nInfG2

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow", at, nInf)
  dat <- set_epi(dat, "si.flow.g2", at, nInfG2)

  return(dat)
}

#' @title Discordant Edgelist
#'
#' @description This function returns a \code{data.frame} with a discordant
#'              edgelist, defined as the set of edges in which the status of the
#'              two partners is one susceptible and one infected.
#'
#' @param dat Master list object containing a \code{networkDynamic} object or
#'        edgelist (if tergmLite is used) and other initialization information
#'        passed from \code{\link{netsim}}.
#' @param at Current time step.
#' @param network In case of models with multiple networks, the network to pull
#'        the current edgelist from. Default of \code{network = 1}.
#' @param infstat Character vector of disease status values that are considered
#'        infectious, defining the SI pairs.
#'
#' @details
#' This internal function works within the parent \code{\link{infection.net}}
#' function to pull the current edgelist from the dynamic network object, look
#' up the disease status of the head and tails on the edge, and subset the list
#' to those edges with one susceptible and one infected node.
#'
#' EpiModel v2.0.3 extended the function by allowing flexibility in the
#' definition what disease status counts as infectious, with the \code{infstat}
#' parameter. For extension models with multiple infectious states, this can be
#' a vector of length greater than 1: \code{infstat = c("i", "a")}.
#'
#' @return
#' This function returns a \code{data.frame} with the following columns:
#' \itemize{
#'  \item \strong{time:} time step queried
#'  \item \strong{sus:} ID number for the susceptible partner
#'  \item \strong{inf:} ID number for the infectious partner
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
discord_edgelist <- function(dat, at, network = 1, infstat = "i") {

  status <- get_attr(dat, "status")
  active <- get_attr(dat, "active")
  tergmLite <- get_control(dat, "tergmLite")

  if (tergmLite == TRUE) {
    el <- dat$el[[network]]
  } else {
    el <- get.dyads.active(dat$nw[[network]], at = at)
  }

  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% infstat, ncol = 2)
    isSus <- matrix(stat %in% "s", ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(at, sus, inf)

      # Check for active status
      keep <- rowSums(matrix(c(active[del$sus], active[del$inf]),
                             ncol = 2)) == 2
      del <- del[keep, ]
      if (nrow(del) < 1) {
        del <- NULL
      }
    }
  }

  return(del)
}

#' @title Save Transmission Matrix
#'
#' @description This function appends the transmission matrix created during
#'              \code{infection.net} and \code{infection.2g.net}.
#'
#' @param dat Master list object containing a \code{networkDynamic} object or
#'        edgelist (if tergmLite is used) and other initialization information
#'        passed from \code{\link{netsim}}.
#' @param at Current time step.
#' @param del Discordant edgelist created within \code{infection.net} and
#'        \code{infection.2g.net}.
#'
#' @details
#' This internal function works within the parent \code{\link{infection.net}}
#' functions to save the transmission matrix created at time step \code{at} to
#' the master list object \code{dat}.
#'
#' @export
#'
set_transmat <- function(dat, del, at) {
  del <- del[!duplicated(del$sus), ]
  # Convert the Discordant Edglist indexes to the corresponding unique_ids
  del[["sus"]] <- get_unique_ids(dat, del[["sus"]])
  del[["inf"]] <- get_unique_ids(dat, del[["inf"]])

  if (at != 2) {
    del <- rbind(dat$stats$transmat, del)
  }

  dat$stats$transmat <- del

  return(dat)
}
