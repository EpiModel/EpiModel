
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

  if (dat$control$tergmLite == FALSE) {

    if (dat$param$vital != FALSE) {

      ## Departures
      inactive <- which(dat$attr$active == 0 & dat$attr$exitTime == at)
      if (length(inactive) > 0) {
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = inactive, deactivate.edges = TRUE)
      }

      ## Arrivals
      if (groups == 1) {
        nArrivals <- dat$epi$a.flow[at]
      } else {
        nArrivals <- c(dat$epi$a.flow[at], dat$epi$a.flow.g2[at])
      }
      if (sum(nArrivals) > 0) {
        nCurr <- network.size(dat$nw)
        dat$nw <- add.vertices(dat$nw, nv = sum(nArrivals))
        newNodes <- (nCurr + 1):(nCurr + sum(nArrivals))
        dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
        if (length(nArrivals) > 1) {
          dat$nw <- set.vertex.attribute(dat$nw, "group",
                                         rep(1:2, c(nArrivals[1], nArrivals[2])),
                                         newNodes)
        }

        # Set attributes on nw
        fterms <- dat$temp$fterms
        curr.tab <- get_attr_prop(dat$nw, fterms)
        if (length(curr.tab) > 0) {
          dat$nw <- update_nwattr(dat$nw, newNodes, dat$control$attr.rules,
                                  curr.tab, dat$temp$t1.tab)
        }

        # Save any val on attr
        if (!is.null(fterms)) {
          dat <- copy_toall_attr(dat, at, fterms)
          if (length(unique(sapply(dat$attr, length))) != 1) {
            stop("Attribute list of unequal length. Check arrivals.net module.")
          }
        }
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = dat$attr$status[newNodes],
                                            onset = at, terminus = Inf,
                                            v = newNodes)

      }
    }

    ## Recovery
    idsRecov <- dat$temp$idsRecov
    recovState <- ifelse(dat$control$type == "SIR", "r", "s")
    status <- dat$attr$status

    if (length(idsRecov) > 0) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = recovState, onset = at,
                                          terminus = Inf, v = idsRecov)
    }

    ## Infection
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

  if (dat$control$tergmLite == TRUE) {

    if (dat$param$vital != FALSE) {

      ## Departures
      inactive <- which(dat$attr$active == 0 & dat$attr$exitTime == at)
      if (length(inactive) > 0) {
        dat$attr <- deleteAttr(dat$attr, inactive)
        dat$el[[1]] <- delete_vertices(dat$el[[1]], inactive)
      }

      ## Arrivals
      if (groups == 1) {
        nArrivals <- dat$epi$a.flow[at]
      } else {
        nArrivals <- c(dat$epi$a.flow[at], dat$epi$a.flow.g2[at])
      }
      if (sum(nArrivals) > 0) {
        nCurr <- sum(dat$attr$active == 1)

        dat$el[[1]] <- add_vertices(dat$el[[1]], nv = sum(nArrivals))

      }
    }

  }

  ## Output
  return(dat)
}
