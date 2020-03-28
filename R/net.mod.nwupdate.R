
#' @title EpiModel Network Writes
#'
#' @description This function handles all call to network object contained in
#' \code{dat$nw} during simulation, whether a direct or indirect manipulation.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#'
nwupdate.net <- function(dat, at) {

  groups <- dat$param$groups
  tL <- dat$control$tergmLite

  ## Arrivals
  if (groups == 1) {
    nArrivals <- dat$epi$a.flow[at]
  } else {
    nArrivals <- c(dat$epi$a.flow[at], dat$epi$a.flow.g2[at])
  }
  if (sum(nArrivals) > 0) {
    nCurr <- length(dat$attr$active)
    newNodes <- (nCurr + 1):(nCurr + sum(nArrivals))
    nwterms <- dat$temp$nwterms
    if (!is.null(nwterms)) {
      curr.tab <- get_attr_prop(dat, nwterms)
      if (length(curr.tab) > 0) {
        dat <- auto_update_attr(dat, newNodes, dat$control$attr.rules,
                                curr.tab, dat$temp$t1.tab)
      }
    }
    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.")
    }
    if (tL == FALSE) {
      dat$nw <- add.vertices(dat$nw, nv = sum(nArrivals))
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = dat$attr$status[newNodes],
                                          onset = at, terminus = Inf,
                                          v = newNodes)
    }
    if (tL == TRUE) {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = sum(nArrivals))
    }
  }


  ## Departures
  inactive <- which(dat$attr$active == 0 & dat$attr$exitTime == at)
  if (length(inactive) > 0) {
    if (tL == FALSE) {
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = inactive, deactivate.edges = TRUE)
    }
    if (tL == TRUE) {
      dat$attr <- deleteAttr(dat$attr, inactive)
      dat$el[[1]] <- delete_vertices(dat$el[[1]], inactive)
    }
  }

  ## Recovery
  if (tL == FALSE) {
    status <- dat$attr$status
    recovState <- ifelse(dat$control$type == "SIR", "r", "s")
    attr.status <- which(status == recovState)
    nw.status <- get.vertex.attribute(dat$nw, "status")
    idsRecov <- setdiff(attr.status, nw.status)

    if (length(idsRecov) > 0) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = recovState, onset = at,
                                          terminus = Inf, v = idsRecov)
    }
  }


  ## Infection
  if (tL == FALSE) {
    idsNewInf <- which(dat$attr$status == "i" & dat$attr$infTime == at)
    if (length(idsNewInf) > 0) {
      dat$nw <- activate.vertex.attribute(dat$nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
      dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
    }
  }


  ## Output
  return(dat)
}
