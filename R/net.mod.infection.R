
#' @title Primary Infection Module for netsim
#'
#' @description This function simulates the main infection process given the
#'              current state of the partnerships and disease in the system.
#'
#' @inheritParams recovery.net
#'
#' @details
#' The main steps in this infection module are as follows:
#' \enumerate{
#'  \item Get IDs for current infected and susceptible nodes given the current
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
#' @inherit recovery.net return
#'
#' @export
#' @keywords netMod internal
#'
#' @seealso \code{\link{discord_edgelist}} is used within \code{infection.net}
#' to obtain a discordant edgelist.
#'
infection.net <- function(dat, at) {
  infection_with_ngroups(dat, at, 1)
}

#' @rdname infection.net
#' @export
#' @keywords netMod internal
#'
infection.2g.net <- function(dat, at) {
  infection_with_ngroups(dat, at, 2)
}

infection_with_ngroups <- function(dat, at, ngroups) {

  # Variables ---------------------------------------------------------------

  active <- get_attr(dat, "active")
  infTime <- get_attr(dat, "infTime")
  status <- get_attr(dat, "status")

  if (ngroups > 1) {
    group <- get_attr(dat, "group")
  } else {
    group <- rep(1, length(active))
  }
  suffixes <- c("", if (ngroups > 1) paste0(".g", seq_len(ngroups)[-1]))

  inf.probs <- lapply(seq_along(suffixes), function(i) get_param(dat, paste0("inf.prob", suffixes[i])))

  act.rate <- get_param(dat, "act.rate")
  inter.eff <- get_param(dat, "inter.eff", override.null.error = TRUE)
  inter.start <- get_param(dat, "inter.start", override.null.error = TRUE)

  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  nInf <- integer(ngroups)
  totInf <- 0

  # Process -----------------------------------------------------------------
  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {

    # Get discordant edgelist
    del_list <- lapply(seq_len(dat$num.nw), discord_edgelist, dat = dat, at = at, include.network = TRUE)
    del <- dplyr::bind_rows(del_list)

    # If some discordant edges, then proceed
    if (NROW(del) > 0) {

      # Infection duration to at
      del$infDur <- pmax(1, at - infTime[del$inf])

      # Calculate infection-stage transmission rates
      for (g in seq_len(ngroups)) {
        ## allow NULL inf.prob.g2, for backwards compatiblity
        inf.prob <- NVL(inf.probs[[g]], inf.probs[[1]])
        del$transProb[group[del$sus] == g] <- inf.prob[pmin(del$infDur[group[del$sus] == g], length(inf.prob))]
      }

      # Interventions
      if (!is.null(inter.eff) && at >= inter.start) {
        del$transProb <- del$transProb * (1 - inter.eff)
      }

      # Calculate infection-stage act/contact rates
      del$actRate <- act.rate[pmin(del$infDur, length(act.rate))]

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
      nInf <- tabulate(group[idsNewInf], nbins = ngroups)
      totInf <- sum(nInf)

    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  for (g in seq_len(ngroups)) {
    dat <- set_epi(dat, paste0("si.flow", suffixes[g]), at, nInf[g])
  }

  return(dat)
}

#' @title Discordant Edgelist
#'
#' @description This function returns a \code{data.frame} with a discordant
#'              edgelist, defined as the set of edges in which the status of the
#'              two partners is one susceptible and one infected.
#'
#' @inheritParams recovery.net
#' @param network In case of models with multiple networks, the network to pull
#'        the current edgelist from. Default of \code{network = 1}.
#' @param infstat Character vector of disease status values that are considered
#'        infectious, defining the SI pairs.
#' @param include.network Should the \code{network} value be included as the
#'        final column of the discordant edgelist?
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
#'  \item \strong{time:} time step queried.
#'  \item \strong{sus:} ID number for the susceptible partner.
#'  \item \strong{inf:} ID number for the infectious partner.
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
discord_edgelist <- function(dat, at, network = 1, infstat = "i", include.network = FALSE) {

  status <- get_attr(dat, "status")
  active <- get_attr(dat, "active")

  el <- get_edgelist(dat, network)

  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(seq_len(nrow(el))), , drop = FALSE]
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

  if (NROW(del) > 0 && include.network == TRUE) {
    del <- dplyr::bind_cols(del, network = network)
  }

  return(del)
}

#' @title Save Transmission Matrix
#'
#' @description This function appends the transmission matrix created during
#'              \code{infection.net} and \code{infection.2g.net}.
#'
#' @inheritParams recovery.net
#' @param del Discordant edgelist created within \code{\link{infection.net}} and
#'        \code{\link{infection.2g.net}}.
#'
#' @details
#' This internal function works within the parent \code{\link{infection.net}}
#' functions to save the transmission matrix created at time step \code{at} to
#' the main list object \code{dat}.
#'
#' @inherit recovery.net return
#'
#' @export
#'
set_transmat <- function(dat, del, at) {
  if (is.null(dat$stats$transmat)) {
    dat$stats$transmat <- list()
  }

  if (length(dat$stats$transmat) < at) {
    nsteps <- get_control(dat, "nsteps")
    dat$stats$transmat <- padded_vector(dat$stats$transmat, nsteps)
  }

  del <- del[!duplicated(del$sus), ]
  # Convert the Discordant Edglist indexes to the corresponding unique_ids
  del[["sus"]] <- get_unique_ids(dat, del[["sus"]])
  del[["inf"]] <- get_unique_ids(dat, del[["inf"]])

  if (length(dat$stats$transmat) < at) {
    dat$stats$transmat <- padded_vector(
      dat$stats$transmat,
      get_control(dat, "nsteps")
    )
  }

  dat$stats$transmat[[at]] <- del

  return(dat)
}
