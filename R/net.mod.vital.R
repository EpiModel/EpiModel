#' @title Departures: netsim Module
#'
#' @description This function simulates departure for use in \link{netsim}
#'        simulations.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
departures.net <- function(dat, at) {
  departures_with_ngroups(dat, at, 1)
}

#' @rdname departures.net
#' @export
#' @keywords netMod internal
#'
departures.2g.net <- function(dat, at) {
  departures_with_ngroups(dat, at, 2)
}

departures_with_ngroups <- function(dat, at, ngroups) {
  # Conditions --------------------------------------------------------------
  vital <- get_param(dat, "vital")
  if (vital == FALSE) {
    return(dat)
  }

  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  exitTime <- get_attr(dat, "exitTime")

  statuses <- c("s", "i", if (type == "SIR") "r")
  suffixes <- c("", if (ngroups > 1) paste0(".g", seq_len(ngroups)[-1L]))
  if (ngroups > 1) {
    group <- get_attr(dat, "group")
  } else {
    group <- rep(1, length.out = length(active))
  }

  for (status_value in statuses) {
    gDepartures <- integer(ngroups)
    idsElig <- which(active == 1 & status == status_value)
    nElig <- length(idsElig)
    if (nElig > 0) {
      rates <- unlist(get_param_list(dat, paste0("d", status_value, ".rate", suffixes)))
      gElig <- group[idsElig]
      gRates <- rates[gElig]
      vecDepartures <- which(rbinom(nElig, 1, gRates) == 1)
      if (length(vecDepartures) > 0) {
        idsDpt <- idsElig[vecDepartures]
        gDepartures <- tabulate(group[idsDpt], nbins = ngroups)
        active[idsDpt] <- 0
        exitTime[idsDpt] <- at
      }
    }
    for (i in seq_along(suffixes)) {
      dat <- set_epi(dat, paste0("d", status_value, ".flow", suffixes[i]), at, gDepartures[i])
    }
  }

  # Output ------------------------------------------------------------------

  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  return(dat)
}


#' @title Arrivals: netsim Module
#'
#' @description This function simulates new arrivals into the network
#'   for use in \code{\link{netsim}} simulations.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
arrivals.net <- function(dat, at) {
  arrivals_with_ngroups(dat, at, 1)
}

#' @rdname arrivals.net
#' @export
#' @keywords netMod internal
#'
arrivals.2g.net <- function(dat, at) {
  arrivals_with_ngroups(dat, at, 2)
}

arrivals_with_ngroups <- function(dat, at, ngroups) {
  # Conditions --------------------------------------------------------------
  vital <- get_param(dat, "vital")
  if (vital == FALSE) {
    return(dat)
  }

  suffixes <- c("", if (ngroups > 1) paste0(".g", seq_len(ngroups)[-1L]))

  # Variables ---------------------------------------------------------------
  a.rate <- unlist(get_param_list(dat, paste0("a.rate", suffixes)))
  a.rate[is.na(a.rate)] <- a.rate[1] # allow a.rate.g2 to be NA, for backwards compatibility

  index <- at - 1
  nOld <- integer(ngroups)
  for (i in seq_len(ngroups)) {
    nOld[i] <- get_epi(dat, paste0("num", suffixes[i]), index)
  }

  # Add Nodes ---------------------------------------------------------------
  nArrivals <- integer(ngroups)
  if (nOld[1] > 0) { # should this be any(nOld > 0)?
    for (i in seq_len(ngroups)) {
      nArrivals[i] <- rbinom(1, nOld[i], a.rate[i])
    }
    totArr <- sum(nArrivals)

    if (totArr > 0) {
      dat <- append_core_attr(dat, at, totArr)
      if (ngroups > 1) {
        for (i in seq_len(ngroups)) {
          dat <- append_attr(dat, "group", i, nArrivals[i])
        }
      }
      dat <- append_attr(dat, "status", "s", totArr)
      dat <- append_attr(dat, "infTime", NA, totArr)
    }
  }

  # Output ------------------------------------------------------------------

  for (i in seq_len(ngroups)) {
    dat <- set_epi(dat, paste0("a.flow", suffixes[i]), at, nArrivals[i])
  }

  return(dat)
}
