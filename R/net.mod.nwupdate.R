
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
nw.update.net <- function(dat, at) {

  if (dat$control$tgl == FALSE) {

    tea.status <- dat$control$tea.status

    #Resimulate Network----

    #Deactive inactive nodes
    inactive <- dat$nw.update$resim$inactive
    if (length(inactive) > 0) {
      dat$attr <- deleteAttr(dat$attr, inactive)
    }

    if (dat$param$vital != FALSE) {

      #Departures----

      idsDpt <- unlist(dat$nw.update$idsDpt)
      idsDpt <- as.vector(idsDpt)
      #Deactive all departures on the network -

      if (length(idsDpt) > 0) {
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDpt, deactivate.edges = TRUE)
      }

      #Arrivals----
      nArrivals <- dat$nw.update$arr$nArrivals
      if (sum(nArrivals) > 0) {
        nCurr <- network.size(dat$nw)
        #New Arrivals
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
        dat <- copy_toall_attr(dat, at, fterms)

        if (tea.status == TRUE) {
          if ("status" %in% fterms) {
            dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                                value = dat$attr$status[newNodes],
                                                onset = at, terminus = Inf,
                                                v = newNodes)
          } else {
            dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                                value = "s", onset = at, terminus = Inf,
                                                v = newNodes)
          }
        }
        if (!("status" %in% fterms)) {
          dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
        }
        dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))

        ## Handles infTime when incoming nodes are infected
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }
    }

    #Recovery----

    idsRecov <- dat$nw.update$rec$idsRecov
    recovState <- dat$nw.update$rec$recovState
    status <- dat$attr$status

    if (length(idsRecov) > 0) {
      if (tea.status == TRUE) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = recovState, onset = at,
                                            terminus = Inf, v = idsRecov)
      }
    }

    if ("status" %in% dat$temp$fterms) {
      dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
    }

    #Infection----

    #Active and set vertex attribute of infected
    idsNewInf <- dat$nw.update$inf$idsNewInf
    tea.status <- dat$control$tea.status
    if (length(idsNewInf) > 0) {
      if (tea.status == TRUE) {
        dat$nw <- activate.vertex.attribute(dat$nw,
                                            prefix = "testatus",
                                            value = "i",
                                            onset = at,
                                            terminus = Inf,
                                            v = idsNewInf)
      }

      if ("status" %in% dat$temp$fterms) {
        dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
      }
    }

  }

  if (dat$control$tgl == TRUE) {

    tea.status <- dat$control$tea.status

    ## Resimulate Network----

    # Delete inactive nodes
    inactive <- which(dat$attr$active == 0)

    if (length(inactive) > 0) {
      el.temp <- dat$nw$el[[1]]
      el.temp <- delete_vertices(el.temp, inactive)
      dat$nw$el[[1]] <- el.temp
    }


    if (dat$param$vital != FALSE) {

      ## Departures----

      idsDpt <- unlist(dat$nw.update$idsDpt)
      idsDpt <- as.vector(idsDpt)
      # Deactive all departures on the network -

      if (length(idsDpt) > 0) {
        el.temp <- dat$el[[1]]
        el.temp <- delete_vertices(el.temp, idsDpt)
        dat$nw$el[[1]] <- el.temp
      }

      ## Arrivals----
      nArrivals <- dat$nw.update$arr$nArrivals
      if (sum(nArrivals) > 0) {
        el.temp <- dat$el[[1]]
        nCurr <- length(dat$attr$group)
        #New Arrivals
        el.temp <- add_vertices(el.temp, nv = sum(nArrivals))

        if (length(nArrivals) > 1) {
          dat$attr$group <- c(dat$attr$group, c(rep(1, nArrivals[1]),
                                                rep(2, nArrivals[2])))
        } else {
          dat$attr$group <- c(dat$attr$group, rep(1, nArrivals))
        }

        dat$attr$status <- c(dat$attr$status, rep("s", sum(nArrivals)))
        dat$attr$active <- c(dat$attr$active, rep(1, sum(nArrivals)))
        dat$attr$infTime <- c(dat$attr$infTime, rep(NA, sum(nArrivals)))
        dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, sum(nArrivals)))
        dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, sum(nArrivals)))

        ## Handles infTime when incoming nodes are infected
        #When would this happen? - Infection occurs before arrivals
        newNodes <- c((nCurr+1):(nCurr+sum(nArrivals)))
        newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
        dat$attr$infTime[newNodesInf] <- at

        if (length(unique(sapply(dat$attr, length))) != 1) {
          stop("Attribute list of unequal length. Check arrivals.net module.")
        }
      }
    }

    ## Recovery----

    #No network updates

    ## Infection----

    #No network updates needed

  }

  # Output
  return(dat)
}
