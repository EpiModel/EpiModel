
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
  rec.rate <- get_param(dat, "rec.rate", override.null.error = TRUE)
  if (is.null(rec.rate)) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  recov.flag <- get_init(dat, "r.num", override.null.error = TRUE)

  recovState <- ifelse(is.null(recov.flag), "s", "r")

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
    }
  }
  dat <- set_attr(dat, "status", status)

  # Output ------------------------------------------------------------------
  outName <- ifelse(is.null(recov.flag), "is.flow", "ir.flow")
  dat <- set_epi(dat, outName, at, nRecov)

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
recovery.2g.net <- function(dat, at) {

  ## Only run with SIR/SIS
  rec.rate <- get_param(dat, "rec.rate", override.null.error = TRUE)
  if (is.null(rec.rate)) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  group <- get_attr(dat, "group")

  recov.flag <- get_init(dat, "r.num", override.null.error = TRUE)
  rec.rate.g2 <- get_param(dat, "rec.rate.g2")

  recovState <- ifelse(is.null(recov.flag), "s", "r")



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
      idsRecov <- idsElig[vecRecov]
      nRecov <- sum(group[idsRecov] == 1)
      nRecovG2 <- sum(group[idsRecov] == 2)
      status[idsRecov] <- recovState
    }
  }
  dat <- set_attr(dat, "status", status)

  # Output ------------------------------------------------------------------
  outName <- ifelse(is.null(recov.flag), "is.flow", "ir.flow")
  outName[2] <- paste0(outName, ".g2")

  dat <- set_epi(dat, outName[1], at, nRecov)
  dat <- set_epi(dat, outName[2], at, nRecovG2)


  return(dat)
}
