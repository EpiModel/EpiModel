
#' @title Primary Infection Module for icm
#'
#' @description This function simulates the main infection process given the
#'              current state of the actors in the system.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
infection.icm <- function(dat, at) {

  ## Expected acts
  if (dat$param$groups == 1) {
    acts <- round(dat$param$act.rate * dat$epi$num[at-1] / 2)
  }
  if (dat$param$groups == 2) {
    if (dat$param$balance == "g1") {
      acts <- round(dat$param$act.rate * (dat$epi$num[at-1] + dat$epi$num.g2[at-1]) / 2)
    }
    if (dat$param$balance == "g2") {
      acts <- round(dat$param$act.rate.g2 * (dat$epi$num[at-1] + dat$epi$num.g2[at-1]) / 2)
    }
  }


  ## Edgelist
  if (dat$param$groups == 1) {
    p1 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1), acts, replace = TRUE)
  } else {
    p1 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 1), acts, replace = TRUE)
    p2 <- ssample(which(dat$attr$active == 1 & dat$attr$group == 2), acts, replace = TRUE)
  }

  del <- NULL
  if (length(p1) > 0 & length(p2) > 0) {
    del <- data.frame(p1, p2)
    if (dat$param$groups == 1) {
      while (any(del$p1 == del$p2)) {
        del$p2 <- ifelse(del$p1 == del$p2, ssample(which(dat$attr$active == 1), 1), del$p2)
      }
    }


    ## Discordant edgelist
    del$p1.stat <- dat$attr$status[del$p1]
    del$p2.stat <- dat$attr$status[del$p2]
    serodis <- (del$p1.stat == "s" & del$p2.stat == "i") |
               (del$p1.stat == "i" & del$p2.stat == "s")
    del <- del[serodis == TRUE, ]


    ## Transmission on edgelist
    if (nrow(del) > 0) {
      if (dat$param$groups == 1) {
        del$tprob <- dat$param$inf.prob
      } else {
        del$tprob <- ifelse(del$p1.stat == "s", dat$param$inf.prob,
                                                dat$param$inf.prob.g2)
      }
      if (!is.null(dat$param$inter.eff) && at >= dat$param$inter.start) {
        del$tprob <- del$tprob * (1 - dat$param$inter.eff)
      }
      del$trans <- rbinom(nrow(del), 1, del$tprob)
      del <- del[del$trans == TRUE, ]
      if (nrow(del) > 0) {
        if (dat$param$groups == 1) {
          newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
          nInf <- length(newIds)
        }
        if (dat$param$groups == 2) {
          newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
          newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
          nInf <- length(newIdsg1)
          nInfg2 <- length(newIdsg2)
          newIds <- c(newIdsg1, newIdsg2)
        }
        dat$attr$status[newIds] <- "i"
        dat$attr$infTime[newIds] <- at
      } else {
        nInf <- nInfg2 <- 0
      }
    } else {
      nInf <- nInfg2 <- 0
    }
  } else {
    nInf <- nInfg2 <- 0
  }


  ## Output
  if (at == 2) {
    dat$epi$si.flow <- c(0, nInf)
  } else {
    dat$epi$si.flow[at] <- nInf
  }
  if (dat$param$groups == 2) {
    if (at == 2) {
      dat$epi$si.flow.g2 <- c(0, nInfg2)
    } else {
      dat$epi$si.flow.g2[at] <- nInfg2
    }
  }

  return(dat)

}


#' @title Recovery: icm Module
#'
#' @description This function simulates recovery from the infected state
#'              either to an distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              \code{\link{icm}}.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
recovery.icm <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  if (!(dat$control$type %in% c("SIR", "SIS"))) {
    return(dat)
  }


  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  status <- dat$attr$status

  groups <- dat$param$groups
  group <- dat$attr$group

  type <- dat$control$type
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rand <- dat$control$rec.rand
  rec.rate <- dat$param$rec.rate
  rec.rate.g2 <- dat$param$rec.rate.g2

  nRecov <- nRecovG2 <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)

  if (nElig > 0) {

    gElig <- group[idsElig]
    rates <- c(rec.rate, rec.rate.g2)
    ratesElig <- rates[gElig]

    if (rec.rand == TRUE) {
      vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecRecov) > 0) {
        idsRecov <- idsElig[vecRecov]
        nRecov <- sum(group[idsRecov] == 1)
        nRecovG2 <- sum(group[idsRecov] == 2)
        status[idsRecov] <- recovState
      }
    } else {
      nRecov <- min(round(sum(ratesElig[gElig == 1])), sum(gElig == 1))
      status[ssample(idsElig[gElig == 1], nRecov)] <- recovState
      if (groups == 2) {
        nRecovG2 <- min(round(sum(ratesElig[gElig == 2])), sum(gElig == 2))
        status[ssample(idsElig[gElig == 2], nRecov)] <- recovState
      }
    }
  }
  dat$attr$status <- status


  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".g2")

  ## Summary statistics
  if (at == 2) {
    dat$epi[[outName[1]]] <- c(0, nRecov)
  } else {
    dat$epi[[outName[1]]][at] <- nRecov
  }
  if (groups == 2) {
    if (at == 2) {
      dat$epi[[outName[2]]] <- c(0, nRecovG2)
    } else {
      dat$epi[[outName[2]]][at] <- nRecovG2
    }
  }

  return(dat)
}
