
##
## Modules that change the status of nodes
##

#' @title Primary Infection Module for netsim
#'
#' @description This function simulates the main infection process given the 
#'              current state of the partnerships and disease in the system.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#' 
#' @details
#' The main steps in this infection module are as follows: 
#' \enumerate{
#'  \item Get IDs for current infected and susceptibles given the current disease
#'        status.
#'  \item Call \code{\link{discord_edgelist}} to get the current discordant edgelist 
#'        given step 1.
#'  \item Determine the transmission rates (e.g., as a function of mode).
#'  \item Pull the number of acts per partnership in a time step from the 
#'        \code{act.rate} parameter.
#'  \item Calculate the final transmission probabilities given the transmission 
#'        rates and act rates.
#'  \item Randomly transmit on the discordant edgelist.
#'  \item Conduct bookkeeping for new infections to update status on the nodes
#'        and calculate disease incidence.
#' }
#' 
#' @return
#' The main \code{all} object is returned with updated disease status and summary
#' incidence measures.
#' 
#' @export
#' @keywords netMod internal
#' 
#' @seealso \code{\link{discord_edgelist}} is used within \code{infection.net} 
#' to obtain a discordant edgelist.
#'
infection.net <- function(all, at) {
  
    # Variables ---------------------------------------------------------------
    active <- all$attr$active
    status <- all$attr$status
    modes <- all$param$modes
    mode <- idmode(all$nw)
    
    trans.rate <- all$param$trans.rate
    trans.rate.m2 <- all$param$trans.rate.m2
    act.rate <- all$param$act.rate
    
    nw <- all$nw
    tea.status <- all$control$tea.status
        
    # Vector of infected and susceptible IDs
    idsSus <- which(active == 1 & status == 0)
    idsInf <- which(active == 1 & status == 1)
    nActive <- sum(active == 1)
    nElig <- length(idsInf)
    
    # Initialize vectors
    nInf <- nInfM2 <- totInf <- 0
    
    
    # Process -----------------------------------------------------------------
    # If some infected AND some susceptible, then proceed
    if (nElig > 0 && nElig < nActive) {
      
      # Get discordant edgelist
      del <- discord_edgelist(all, idsInf, idsSus, at)
      
      # If some discordant edges, then proceed
      if (!(is.null(del))) {
        
        # Infection duration to at
        del$infDur <- at - all$attr$infTime[del$inf]
        del$infDur[del$infDur == 0] <- 1
        
        # Calculate infection-stage transmission rates
        ltrans.rate <- length(trans.rate)
        if (is.null(trans.rate.m2)) {
          del$transProb <- ifelse(del$infDur <= ltrans.rate, 
                                  trans.rate[del$infDur], 
                                  trans.rate[ltrans.rate])
        } else {
          del$transProb <- ifelse(del$sus <= nw %n% "bipartite", 
                                  ifelse(del$infDur <= ltrans.rate, 
                                         trans.rate[del$infDur], 
                                         trans.rate[ltrans.rate]),
                                  ifelse(del$infDur <= ltrans.rate, 
                                         trans.rate.m2[del$infDur], 
                                         trans.rate.m2[ltrans.rate]))
        }
        
        # Calculate infection-stage act/contact rates
        lact.rate <- length(act.rate)
        del$actRate <- ifelse(del$infDur <= lact.rate, 
                              act.rate[del$infDur], 
                              act.rate[lact.rate])
        
        # Calculate final transmission probability per timestep
        del$finalProb <- 1-(1-del$transProb)^del$actRate
        
        # Randomize transmissions and subset df
        transmit <- rbinom(nrow(del), 1, del$finalProb)
        del <- del[which(transmit == 1), ]
        
        # Set new infections vector
        idsNewInf <- unique(del$sus)
        nInf <- sum(mode[idsNewInf] == 1)
        nInfM2 <- sum(mode[idsNewInf] == 2)
        totInf <- nInf + nInfM2
        
        # Update nw attributes
        if (totInf > 0) {
          if (tea.status == TRUE) {
            nw <- activate.vertex.attribute(nw, 
                                            prefix = "testatus", 
                                            value = 1, 
                                            onset = at, 
                                            terminus = Inf, 
                                            v = idsNewInf)
          }
          all$attr$status[idsNewInf] <- 1
          all$attr$infTime[idsNewInf] <- at
          
          t <- get_formula_terms(all$nwparam$formation)
          if ("status" %in% t) {
            nw <- set.vertex.attribute(nw, "status", all$attr$status)
          }
        }
        
        # Substitute PIDs for vital bipartite sims
        if (any(names(nw$gal) %in% "vertex.pid")) {
          del$sus <- get.vertex.pid(nw, del$sus)
          del$inf <- get.vertex.pid(nw, del$inf)
        }
        
      } # end some discordant edges condition
    } # end some active discordant nodes condition
  
    
   # Output ------------------------------------------------------------------
  
   # Save transmission matrix
   if (totInf > 0) {
     if (at == 2) {
       all$stats$transmat <- del
     } else {
       all$stats$transmat <- rbind(all$stats$transmat, del)
     }
   }
  
   ## Save incidence vector
   if (at == 2) {
     all$out$si.flow <- c(0, nInf)
     if (modes == 2) {
       all$out$si.flow.m2 <- c(0, nInfM2)
     }
   } else {
     all$out$si.flow[at] <- nInf 
     if (modes == 2) {
       all$out$si.flow.m2[at] <- nInfM2
     }
   }
  
   all$nw <- nw
   return(all)
}




#' @title Discordant Edgelist from NetworkDynamic Object
#'
#' @description This function returns a \code{data.frame} with a discordant 
#'              edgelist, defined as the set of edges in which the status of the 
#'              two partners is one susceptible and one infected.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param idsInf vector of IDs for currently infecteds.
#' @param idsSus vector of IDs for currently susceptible.
#' @param at current time step.
#' 
#' @details
#' This internal function works within the parent \code{\link{infection.net}} function
#' to pull the current edgelist from the dynamic network object, look up the disease
#' status of the head and tails on the edge, and subset the list to those edges
#' with one susceptible and one infected node.
#' 
#' @return
#' This function returns a \code{data.frame} with the following columns:
#' \itemize{
#'  \item \strong{time:} time step queried
#'  \item \strong{sus:} ID number for the susceptible partner
#'  \item \strong{inf:} ID number for the infected partner
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
discord_edgelist <- function(all, idsInf, idsSus, at) {
  
  status <- all$attr$status
  el <- get.dyads.active(all$nw, at = at)
  
  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% 1, ncol = 2)
    isSus <- matrix(stat %in% 0, ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(at, sus, inf)
    } 
  }
  
  return(del)
}


#' @title Recovery: netsim Module
#'
#' @description This function simulates recovery from the infected state
#'              either to an distinct recovered state (SIR model type) or back 
#'              to a susceptible state (SIS model type), for use in 
#'              \code{\link{netsim}}.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
#' 
#' @export
#' @keywords internal
#'
recovery.net <- function(all, at) {
  
  ## Only run with SIR/SIS
  if (!(all$control$type %in% c("SIR", "SIS"))) {
    return(all)
  }
  
  # Variables ---------------------------------------------------------------
  active <- all$attr$active
  status <- all$attr$status
  tea.status <- all$control$tea.status

  modes <- all$param$modes
  mode <- idmode(all$nw)
  
  type <- all$control$type
  recovState <- ifelse(type == "SIR", 2, 0)
  
  rec.rand <- all$control$rec.rand
  rec.rate <- all$param$rec.rate
  rec.rate.m2 <- all$param$rec.rate.m2
  
  nRecov <- nRecovM2 <- 0
  idsElig <- which(active == 1 & status == 1)
  nElig <- length(idsElig)
  
  
  # Process -----------------------------------------------------------------
  if (nElig > 0) {
    
    mElig <- mode[idsElig]
    rates <- c(rec.rate, rec.rate.m2)
    ratesElig <- rates[mElig]
    
    if (rec.rand == TRUE) {
      vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
      if (length(vecRecov) > 0) {
        idsRecov <- idsElig[vecRecov]
        nRecov <- sum(mode[idsRecov] == 1)
        nRecovM2 <- sum(mode[idsRecov] == 2)
        status[idsRecov] <- recovState
        if (tea.status == TRUE) {
          all$nw <- activate.vertex.attribute(all$nw, 
                                              prefix = "testatus", 
                                              value = recovState, 
                                              onset = at, 
                                              terminus = Inf, 
                                              v = idsRecov)
        }
      }
    } else {
      idsRecov <- idsRecovM2 <- NULL
      nRecov <- min(round(sum(ratesElig[mElig == 1])), sum(mElig == 1))
      status[ssample(idsElig[mElig == 1], nRecov)] <- recovState
      if (modes == 2) {
        nRecovM2 <- min(round(sum(ratesElig[mElig == 2])), sum(mElig == 2))
        status[ssample(idsElig[mElig == 2], nRecov)] <- recovState
      }
      totRecov <- nRecov + nRecovM2
      if (tea.status == TRUE & totRecov > 0) {
        allids <- c(idsRecov, idsRecovM2)
        all$nw <- activate.vertex.attribute(all$nw, 
                                            prefix = "testatus", 
                                            value = recovState, 
                                            onset = at, 
                                            terminus = Inf, 
                                            v = allids)
      }
    }
  }
  
  all$attr$status <- status
  t <- get_formula_terms(all$nwparam$formation)
  if ("status" %in% t) {
    all$nw <- set.vertex.attribute(all$nw, "status", all$attr$status)
  }
  
  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  outName[2] <- paste0(outName, ".m2")
  
  if (at == 2) {
    all$out[[outName[1]]] <- c(0, nRecov)
  } else {
    all$out[[outName[1]]][at] <- nRecov
  }
  if (modes == 2) {
    if (at == 2) {
      all$out[[outName[2]]] <- c(0, nRecovM2)
    } else {
      all$out[[outName[2]]][at] <- nRecovM2
    }
  }
  
  return(all)
}
