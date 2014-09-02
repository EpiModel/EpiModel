
#' @title Primary Infection Module for icm
#'
#' @description This function simulates the main infection process given the
#'              current state of the actors in the system.
#'
#' @param all master data list object.
#' @param at current time step.
#'
#' @export
#' @keywords internal
#'
infection.icm <- function(all, at) {

  ## Expected acts
  if (all$param$groups == 1) {
    acts <- round(all$param$act.rate * all$out$num[at-1] / 2)
  }
  if (all$param$groups == 2) {
    if (all$param$balance == "g1") {
      acts <- round(all$param$act.rate * (all$out$num[at-1] + all$out$num.g2[at-1]) / 2)
    }
    if (all$param$balance == "g2") {
      acts <- round(all$param$act.rate.g2 * (all$out$num[at-1] + all$out$num.g2[at-1]) / 2)
    }
  }


  ## Edgelist
  if (all$param$groups == 1) {
    p1 <- ssample(which(all$attr$active == 1), acts, replace = TRUE)
    p2 <- ssample(which(all$attr$active == 1), acts, replace = TRUE)
  } else {
    p1 <- ssample(which(all$attr$active == 1 & all$attr$group == 1), acts, replace = TRUE)
    p2 <- ssample(which(all$attr$active == 1 & all$attr$group == 2), acts, replace = TRUE)
  }

  if (length(p1) > 0 & length(p2) > 0) {
    del <- data.frame(p1, p2)
    if (all$param$groups == 1) {
      while (any(del$p1 == del$p2)) {
        del$p2 <- ifelse(del$p1 == del$p2, ssample(which(all$attr$active == 1), 1), del$p2)
      }
    }


    ## Discordant edgelist
    del$p1.stat <- all$attr$status[del$p1]
    del$p2.stat <- all$attr$status[del$p2]
    serodis <- (del$p1.stat == "s" & del$p2.stat == "i") |
               (del$p1.stat == "i" & del$p2.stat == "s")
    del <- del[serodis == TRUE, ]


    ## Transmission on edgelist
    if (nrow(del) > 0) {
      if (all$param$groups == 1) {
        del$tprob <- all$param$inf.prob
      } else {
        del$tprob <- ifelse(del$p1.stat == 0, all$param$inf.prob, all$param$inf.prob.g2)
      }
      del$trans <- rbinom(nrow(del), 1, del$tprob)
      del <- del[del$trans == TRUE, ]
      if (nrow(del) > 0) {
        if (all$param$groups == 1) {
          newIds <- unique(ifelse(del$p1.stat == "s", del$p1, del$p2))
          nInf <- length(newIds)
        }
        if (all$param$groups == 2) {
          newIdsg1 <- unique(del$p1[del$p1.stat == "s"])
          newIdsg2 <- unique(del$p2[del$p2.stat == "s"])
          nInf <- length(newIdsg1)
          nInfg2 <- length(newIdsg2)
          newIds <- c(newIdsg1, newIdsg2)
        }
        all$attr$status[newIds] <- "i"
        all$attr$infTime[newIds] <- at
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
    all$out$si.flow <- c(0, nInf)
  } else {
    all$out$si.flow[at] <- nInf
  }
  if (all$param$groups == 2) {
    if (at == 2) {
      all$out$si.flow.g2 <- c(0, nInfg2)
    } else {
      all$out$si.flow.g2[at] <- nInfg2
    }
  }

  return(all)

}


#' @title Recovery: icm Module
#'
#' @description This function simulates recovery from the infected state
#'              either to an distinct recovered state (SIR model type) or back
#'              to a susceptible state (SIS model type), for use in
#'              \code{\link{icm}}.
#'
#' @param all master data list object.
#' @param at current time step.
#'
#' @export
#' @keywords internal
#'
recovery.icm <- function(all, at) {

  # Conditions --------------------------------------------------------------
  if (!(all$control$type %in% c("SIR", "SIS"))) {
    return(all)
  }


  # Variables ---------------------------------------------------------------
  active <- all$attr$active
  status <- all$attr$status

  groups <- all$param$groups
  group <- all$attr$group

  type <- all$control$type
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rand <- all$control$rec.rand
  rec.rate <- all$param$rec.rate
  rec.rate.g2 <- all$param$rec.rate.g2

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
  all$attr$status <- status


  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".g2")

  ## Summary statistics
  if (at == 2) {
    all$out[[outName[1]]] <- c(0, nRecov)
  } else {
    all$out[[outName[1]]][at] <- nRecov
  }
  if (groups == 2) {
    if (at == 2) {
      all$out[[outName[2]]] <- c(0, nRecovG2)
    } else {
      all$out[[outName[2]]][at] <- nRecovG2
    }
  }

  return(all)
}
