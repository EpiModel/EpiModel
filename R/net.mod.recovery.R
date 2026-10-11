
#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to a distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              [netsim()].
#'
#' @param dat Main `netsim_dat` object containing a `networkDynamic`
#'        object and other initialization information passed from
#'        [netsim()].
#' @param at Current time step.
#'
#' @return The updated `netsim_dat` main list object.
#'
#' @export
#' @keywords internal
#'
recovery.net <- function(dat, at) {

  ## Only run with SIR/SIS
  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  if (!(type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rate <- get_param(dat, "rec.rate")

  nRecov <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)


  # Recovery Rate by Duration of Infection ----------------------------------
  infDur <- at - infTime[idsElig]
  infDur[infDur == 0] <- 1
  ratesElig <- infection_duration_value(rec.rate, infDur, "rec.rate")


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
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  dat <- set_epi(dat, outName, at, nRecov)

  return(dat)
}

#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to a distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              [netsim()].
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords internal
#'
recovery.2g.net <- function(dat, at) {

  ## Only run with SIR/SIS
  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  if (!(type %in% c("SIR", "SIS"))) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  group <- get_attr(dat, "group")
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rate <- get_param(dat, "rec.rate")
  rec.rate.g2 <- get_param(dat, "rec.rate.g2")

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  # Recovery Rate by Duration of Infection ----------------------------------
  infDur <- at - infTime[idsElig]
  infDur[infDur == 0] <- 1
  gElig <- group[idsElig]
  ratesElig <- infection_duration_value(rec.rate, infDur, "rec.rate")
  if (!is.null(rec.rate.g2)) {
    g2 <- which(gElig == 2)
    ratesElig[g2] <- infection_duration_value(rec.rate.g2, infDur[g2],
                                              "rec.rate.g2")
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
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".g2")

  dat <- set_epi(dat, outName[1], at, nRecov)
  dat <- set_epi(dat, outName[2], at, nRecovG2)


  return(dat)
}
