
#' @title Dynamic Network Updates
#'
#' @description This function handles all calls to the network object contained
#'              on the main \code{dat} object handled in \code{netsim}.
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
  cumulative.edgelist <- get_control(dat, "cumulative.edgelist")

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
      dat$nw[[1]] <- add.vertices(dat$nw[[1]], nv = nArrivals)

      dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at,
                                       terminus = Inf, v = arrivals)

      dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                               prefix = "testatus",
                                               value = status[arrivals],
                                               onset = at, terminus = Inf,
                                               v = arrivals)
    }
    if (tergmLite == TRUE) {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = nArrivals)
    }
  }


  ## Departures
  if (length(departures) > 0) {
    if (tergmLite == FALSE) {
      dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at,
                                         terminus = Inf, v = departures,
                                         deactivate.edges = TRUE)
    }
    if (tergmLite == TRUE) {
      dat <- delete_attr(dat, departures)
      dat$el[[1]] <- delete_vertices(dat$el[[1]], departures)

      if (dat$control$tergmLite.track.duration) {
        dat$nw[[1]] %n% "lasttoggle" <-
          delete_vertices(dat$nw[[1]] %n% "lasttoggle", departures)
      }
    }
  }

  ## Copy static attributes to network object
  if (tergmLite == FALSE && resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }

  ## Update temporally extended disease status
  if (tergmLite == FALSE) {
    dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                             prefix = "testatus",
                                             value = status,
                                             onset = at,
                                             terminus = Inf)
  }

  # Cummulative edgelist
  if (cumulative.edgelist) {
    truncate.el.cuml <- get_control(dat, "truncate.el.cuml")

    for (network in seq_along(dat[["nwparam"]])) {
      dat <- update_cumulative_edgelist(dat, network, truncate.el.cuml)
    }
  }

  ## Output
  return(dat)
}
