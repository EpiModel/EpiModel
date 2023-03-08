
#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to a distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              \code{\link{netsim}}.
#'
#' @param dat Main list object containing a \code{networkDynamic} object and
#'        other initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @return The updated \code{dat} main list object.
#'
#' @export
#' @keywords internal
#'
recovery.net <- function(dat, at) {
  recovery_with_ngroups(dat, at, 1)
}

#' @rdname recovery.net
#' @export
#' @keywords internal
#'
recovery.2g.net <- function(dat, at) {
  recovery_with_ngroups(dat, at, 2)
}

recovery_with_ngroups <- function(dat, at, ngroups) {
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
  if (ngroups > 1) {
    group <- get_attr(dat, "group")
  } else {
    group <- rep(1, length.out = length(active))
  }
  recovState <- ifelse(type == "SIR", "r", "s")

  suffixes <- c("", if (ngroups > 1) paste0(".g", seq_len(ngroups)[-1L]))

  rec.rates <- lapply(seq_len(ngroups), function(i) get_param(dat, paste0("rec.rate", suffixes[i])))

  nRecov <- integer(ngroups)
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- pmax(1, at - infTime[active == 1 & status == "i"])

  gElig <- group[idsElig]
  ratesElig <- numeric(length(gElig))
  for(g in seq_len(ngroups)) {
    rec.rate <- NVL(rec.rates[[g]], rec.rates[[1]]) # allow rec.rates.g2 to be NULL, for backwards compatibility
    ratesElig[gElig == g] <- rec.rate[pmin(infDur[gElig == g], length(rec.rate))]
  }

  # Process -----------------------------------------------------------------
  if (nElig > 0) {
    vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecRecov) > 0) {
      idsRecov <- idsElig[vecRecov]
      nRecov <- tabulate(group[idsRecov], nbins = ngroups)
      status[idsRecov] <- recovState
    }
  }
  dat <- set_attr(dat, "status", status)

  # Output ------------------------------------------------------------------
  if (type == "SIR") {
    outName <- "ir.flow"
  } else {
    outName <- "is.flow"
  }

  for (g in seq_len(ngroups)) {
    dat <- set_epi(dat, paste0(outName, suffixes[g]), at, nRecov[g])
  }

  return(dat)
}
