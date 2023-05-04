
#' @title Dynamic Network Updates
#'
#' @description This function handles all calls to the network object contained
#'              on the main \code{netsim_dat} object handled in \code{netsim}.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#'
nwupdate.net <- function(dat, at) {

  ## Attributes
  status <- get_attr(dat, "status")
  active <- get_attr(dat, "active")
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")

  ## Controls
  tergmLite <- get_control(dat, "tergmLite")
  resimulate.network <- get_control(dat, "resimulate.network")

  ## Vital Dynamics
  arrivals <- which(active == 1 & entrTime == at)
  departures <- which(active == 0 & exitTime == at)

  nArrivals <- length(arrivals)
  if (nArrivals > 0) {

    ## Arrivals
    nwterms <- dat$temp$nwterms
    if (!is.null(nwterms)) {
      curr.tab <- get_attr_prop(dat, nwterms)
      dat <- auto_update_attr(dat, arrivals, curr.tab)
    }
    if (length(unique(sapply(dat$attr, length))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.\n",
           print(cbind(sapply(get_attr_list(dat), length))))
    }
    if (tergmLite == FALSE) {
      for (network in seq_len(dat$num.nw)) {
        dat$nw[[network]] <- add.vertices(dat$nw[[network]], nv = nArrivals)

        dat$nw[[network]] <- activate.vertices(dat$nw[[network]],
                                               onset = at,
                                               terminus = Inf,
                                               v = arrivals)

        dat$nw[[network]] <- activate.vertex.attribute(dat$nw[[network]],
                                                       prefix = "testatus",
                                                       value = status[arrivals],
                                                       onset = at,
                                                       terminus = Inf,
                                                       v = arrivals)
      }
    }
    if (tergmLite == TRUE) {
      for (network in seq_len(dat$num.nw)) {
        dat$el[[network]] <- add_vertices(dat$el[[network]], nv = nArrivals)
        dat$net_attr[[network]][["n"]] <- dat$net_attr[[network]][["n"]] + nArrivals
      }
    }
  }


  ## Departures
  if (length(departures) > 0) {
    if (tergmLite == FALSE) {
      for (network in seq_len(dat$num.nw)) {
        dat$nw[[network]] <- deactivate.vertices(dat$nw[[network]],
                                                 onset = at,
                                                 terminus = Inf,
                                                 v = departures,
                                                 deactivate.edges = TRUE)
      }
    }
    if (tergmLite == TRUE) {
      dat <- delete_attr(dat, departures)
      for (network in seq_len(dat$num.nw)) {
        dat$el[[network]] <- delete_vertices(dat$el[[network]], departures)
        dat$net_attr[[network]][["n"]] <- dat$net_attr[[network]][["n"]] - length(departures)

        if (get_network_control(dat, network, "tergmLite.track.duration") == TRUE) {
          dat$net_attr[[network]][["lasttoggle"]] <-
            delete_vertices(dat$net_attr[[network]][["lasttoggle"]], departures)
        }
      }
    }
  }

  ## Copy static attributes to network object
  if (tergmLite == FALSE && resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }

  ## Update temporally extended disease status
  if (tergmLite == FALSE) {
    for (network in seq_len(dat$num.nw)) {
      dat$nw[[network]] <- activate.vertex.attribute(dat$nw[[network]],
                                                     prefix = "testatus",
                                                     value = status,
                                                     onset = at,
                                                     terminus = Inf)
    }
  }

  ## Output
  return(dat)
}
