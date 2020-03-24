
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
      dat$temp$idsRecov <- idsRecov
    }
  }
  dat$attr$status <- status

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
recovery.2g.net <- function(dat, at) {

  ## Only run with SIR/SIS
  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status
  infTime <- dat$attr$infTime

  if (dat$control$tergmLite == FALSE) {
    group <- get.vertex.attribute(dat$nw, "group")
  } else {
    group <- dat$attr$group
  }

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
      dat$temp$idsRecov <- idsRecov <- idsElig[vecRecov]
      nRecov <- sum(group[idsRecov] == 1)
      nRecovG2 <- sum(group[idsRecov] == 2)
      status[idsRecov] <- recovState
    }
  }
  dat$attr$status <- status

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
