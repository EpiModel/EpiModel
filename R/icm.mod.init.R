
#' @title Initialization: icm Module
#'
#' @description This function initializes the master \code{all} object on which 
#'              data are stored, and simulates disease status and other attributes.
#' 
#' @param param an \code{EpiModel} object of class \code{\link{param.icm}}.
#' @param init an \code{EpiModel} object of class \code{\link{init.icm}}.
#' @param control an \code{EpiModel} object of class \code{\link{control.icm}}.
#' 
#' @export
#' @keywords internal
#'
initialize.icm <- function(param, init, control) {
  
  ## Master List for Data ##
  all <- list()
  all$param <- param
  all$init <- init
  all$control <- control
  
  
  # Set attributes
  all$attr <- list()
  numeric.init <- init[which(sapply(init, class) == "numeric")]
  n <- do.call("sum", numeric.init)
  all$attr$active <- rep(1, n)
  if (all$param$groups == 1) {
    all$attr$group <- rep(1, n)
  } else {
    g2inits <- grep(".g2", names(numeric.init))
    g1inits <- setdiff(1:length(numeric.init), g2inits)
    nG1 <- sum(sapply(g1inits, function(x) init[[x]]))
    nG2 <- sum(sapply(g2inits, function(x) init[[x]]))
    all$attr$group <- c(rep(1, nG1), rep(2, max(0, nG2)))
  }
   
  # Initialize status and infection time
  all <- init_status.icm(all)
  
  
  # Summary out list
  all <- get_prev.icm(all, at = 1)
  
  return(all)
}


#' @title Disease Status Initialization Module for icm
#'
#' @description This function sets the initial disease status on the 
#'              network given the specified initial conditions.
#' 
#' @param all master data list object.
#' 
#' @seealso This is an initialization module for \code{\link{icm}}.
#' 
#' @export
#' @keywords internal
#'
init_status.icm <- function(all) {

  # Variables ---------------------------------------------------------------
  type <- all$control$type
  active <- all$attr$active
  group <- all$attr$group
  nGroups <- all$param$groups
  
  nG1 <- sum(group == 1)
  nG2 <- sum(group == 2)
  
  i.num <- all$init$i.num
  r.num <- all$init$r.num
  i.num.g2 <- all$init$i.num.g2
  r.num.g2 <- all$init$r.num.g2
  
  status.rand <- all$init$status.rand
  
  
  # Status ------------------------------------------------------------------
  if (status.rand == FALSE) {
    status <- rep(0, nG1 + nG2)
    status[sample(which(group == 1), size = i.num)] <- 1
    if (nGroups == 2) {
      status[sample(which(group == 2), size = i.num.g2)] <- 1
    }
    if (type == "SIR") {
      status[sample(which(group == 1 & status == 0), size = r.num)] <- 2
      if (nGroups == 2) {
        status[sample(which(group == 2 & status == 0), size = r.num.g2)] <- 2
      }
    }
  } else {
    status <- rep(NA, nG1 + nG2)
    if (type == "SIR") {
      status[which(group == 1)] <- sample(
        x = c(0:2), 
        size = nG1, 
        replace = TRUE, 
        prob = c(1-(i.num/nG1)-(r.num/nG1), i.num/nG1, r.num/nG1))
      if (sum(status == 1 & group == 1) == 0 & i.num > 0) {
        status[sample(which(group == 1), size = i.num)] <- 1
      }
      if (sum(status == 2 & group == 1) == 0 & r.num > 0) {
        status[sample(which(group == 1), size = r.num)] <- 2
      }
      if (nGroups == 2) {
        status[which(group == 2)] <- sample(
          x = c(0:2), 
          size = nG2, 
          replace = TRUE, 
          prob = c(1-(i.num.g2/nG2)-(r.num.g2/nG2), i.num.g2/nG2, r.num.g2/nG2))
        if (sum(status == 1 & group == 2) == 0 & i.num.g2 > 0) {
          status[sample(which(group == 2), size = i.num.g2)] <- 1
        }
        if (sum(status == 2 & group == 2) == 0 & r.num.g2 > 0) {
          status[sample(which(group == 2), size = r.num.g2)] <- 2
        }
      }
    } else {
      status[which(group == 1)] <- sample(
        x = c(0:1), 
        size = nG1, 
        replace = TRUE, 
        prob = c(1-(i.num/nG1), i.num/nG1))
      if (sum(status == 1 & group == 1) == 0 & i.num > 0) {
        status[sample(which(group == 1), size = i.num)] <- 1
      }
      if (nGroups == 2) {
        status[which(group == 2)] <- sample(
          x = c(0:1), 
          size = nG2, 
          replace = TRUE, 
          prob = c(1-(i.num.g2/nG2), i.num.g2/nG2))
        if (sum(status == 1 & group == 2) == 0 & i.num.g2 > 0) {
          status[sample(which(group == 2), size = i.num.g2)] <- 1
        }
      }
    } 
  }
  all$attr$status <- status
  
  
  # Infection Time ----------------------------------------------------------  
  ## Set up inf.time vector
  idsInf <- which(status == 1)
  infTime <- rep(NA, length(status))
  
  # If vital=TRUE, infTime is a uniform draw over the duration of infection
  if (all$param$vital == TRUE && all$param$di.rate > 0) {
    infTime[idsInf] <- -rgeom(n = length(idsInf), 
                              prob = all$param$di.rate)+2
  } else {
    if (all$control$type == "SI" || all$param$rec.rate == 0) {
      # infTime a uniform draw over the number of sim time steps
      infTime[idsInf] <- ssample(1:(-all$control$nsteps + 2), 
                                 length(idsInf), 
                                 replace = TRUE)
    } else {
      infTime[idsInf] <- ssample(1:(-round(1/all$param$rec.rate) + 2), 
                                 length(idsInf), 
                                 replace = TRUE)
      #TODO: divide this by group if rec.rate != rec.rate.g2
    }
  }
  all$attr$infTime <- infTime
  
  return(all)
}


