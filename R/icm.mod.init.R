
#' @title Initialization: icm Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, and simulates disease status and other attributes.
#'
#' @param param An \code{EpiModel} object of class \code{\link{param.icm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.icm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.icm}}.
#'
#' @export
#' @keywords internal
#'
initialize.icm <- function(param, init, control) {

  ## Master List for Data ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control


  # Set attributes
  dat$attr <- list()
  numeric.init <- init[which(sapply(init, class) == "numeric")]
  n <- do.call("sum", numeric.init)
  dat$attr$active <- rep(1, n)
  if (dat$param$groups == 1) {
    dat$attr$group <- rep(1, n)
  } else {
    g2inits <- grep(".g2", names(numeric.init))
    g1inits <- setdiff(1:length(numeric.init), g2inits)
    nG1 <- sum(sapply(g1inits, function(x) init[[x]]))
    nG2 <- sum(sapply(g2inits, function(x) init[[x]]))
    dat$attr$group <- c(rep(1, nG1), rep(2, max(0, nG2)))
  }

  # Initialize status and infection time
  dat <- init_status.icm(dat)


  # Summary out list
  dat <- get_prev.icm(dat, at = 1)

  return(dat)
}


#' @title Disease Status Initialization Module for icm
#'
#' @description This function sets the initial disease status on the
#'              network given the specified initial conditions.
#'
#' @param dat Master data list object.
#'
#' @seealso This is an initialization module for \code{\link{icm}}.
#'
#' @export
#' @keywords internal
#'
init_status.icm <- function(dat) {

  # Variables ---------------------------------------------------------------
  type <- dat$control$type
  active <- dat$attr$active
  group <- dat$attr$group
  nGroups <- dat$param$groups

  nG1 <- sum(group == 1)
  nG2 <- sum(group == 2)

  i.num <- dat$init$i.num
  r.num <- dat$init$r.num
  i.num.g2 <- dat$init$i.num.g2
  r.num.g2 <- dat$init$r.num.g2

  status.rand <- dat$init$status.rand


  # Status ------------------------------------------------------------------
  if (status.rand == FALSE) {
    status <- rep("s", nG1 + nG2)
    status[sample(which(group == 1), size = i.num)] <- "i"
    if (nGroups == 2) {
      status[sample(which(group == 2), size = i.num.g2)] <- "i"
    }
    if (type == "SIR") {
      status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
      if (nGroups == 2) {
        status[sample(which(group == 2 & status == "s"), size = r.num.g2)] <- "r"
      }
    }
  } else {
    status <- rep(NA, nG1 + nG2)
    if (type == "SIR") {
      status[which(group == 1)] <- sample(
        x = c("s", "i", "r"),
        size = nG1,
        replace = TRUE,
        prob = c(1-(i.num/nG1)-(r.num/nG1), i.num/nG1, r.num/nG1))
      if (sum(status == "i" & group == 1) == 0 & i.num > 0) {
        status[sample(which(group == 1), size = i.num)] <- "i"
      }
      if (sum(status == "r" & group == 1) == 0 & r.num > 0) {
        status[sample(which(group == 1), size = r.num)] <- "r"
      }
      if (nGroups == 2) {
        status[which(group == 2)] <- sample(
          x = c("s", "i", "r"),
          size = nG2,
          replace = TRUE,
          prob = c(1-(i.num.g2/nG2)-(r.num.g2/nG2), i.num.g2/nG2, r.num.g2/nG2))
        if (sum(status == "i" & group == 2) == 0 & i.num.g2 > 0) {
          status[sample(which(group == 2), size = i.num.g2)] <- "i"
        }
        if (sum(status == "r" & group == 2) == 0 & r.num.g2 > 0) {
          status[sample(which(group == 2), size = r.num.g2)] <- "r"
        }
      }
    } else {
      status[which(group == 1)] <- sample(
        x = c("s", "i"),
        size = nG1,
        replace = TRUE,
        prob = c(1-(i.num/nG1), i.num/nG1))
      if (sum(status == "i" & group == 1) == 0 & i.num > 0) {
        status[sample(which(group == 1), size = i.num)] <- "i"
      }
      if (nGroups == 2) {
        status[which(group == 2)] <- sample(
          x = c("s", "i"),
          size = nG2,
          replace = TRUE,
          prob = c(1-(i.num.g2/nG2), i.num.g2/nG2))
        if (sum(status == "i" & group == 2) == 0 & i.num.g2 > 0) {
          status[sample(which(group == 2), size = i.num.g2)] <- "i"
        }
      }
    }
  }
  dat$attr$status <- status


  # Infection Time ----------------------------------------------------------
  idsInf <- which(status == "i")
  infTime <- rep(NA, length(status))

  # If vital=TRUE, infTime is a uniform draw over the duration of infection
  if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
    infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate)+2
  } else {
    if (dat$control$type == "SI" || dat$param$rec.rate == 0) {
      # infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                 length(idsInf), replace = TRUE)
    } else {
      if (nGroups == 1) {
        infTime[idsInf] <- ssample(1:(-round(1/dat$param$rec.rate) + 2),
                                   length(idsInf), replace = TRUE)
      }
      if (nGroups == 2) {
        infG1 <- which(status == "i" & group == 1)
        infTime[infG1] <- ssample(1:(-round(1/dat$param$rec.rate) + 2),
                                  length(infG1), replace = TRUE)
        infG2 <- which(status == "i" & group == 2)
        infTime[infG2] <- ssample(1:(-round(1/dat$param$rec.rate.g2) + 2),
                                  length(infG2), replace = TRUE)
      }
    }
  }
  dat$attr$infTime <- infTime

  return(dat)
}
