
#' @title Dynamic Network Updates
#'
#' @description This function handles all calls to the network object contained
#'              on the master dat object handled in \code{netsim}..
#'
#' @param dat Master list object containing a full \code{networkDynamic} object
#'        or networkLite edgelist (if using tergmLite), and other initialization
#'        information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#'
nwupdate.net <- function(dat, at) {

  ## Attributes
  type <- get_control(dat, "type", override.null.error = TRUE)
  tergmLite <- get_control(dat, "tergmLite")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  active <- get_attr(dat, "active")
  entrTime <- get_attr(dat, "entrTime")
  exitTime <- get_attr(dat, "exitTime")

  # statOnNw <- "status" %in% dat$temp$nwterms
  resimulate.network <- get_control(dat, "resimulate.network")
  isTERGM <- get_control(dat, "isTERGM")

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
      if (isTERGM == TRUE) {
        dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at,
                                         terminus = Inf, v = arrivals)
        dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                                 prefix = "testatus",
                                                 value = status[arrivals],
                                                 onset = at, terminus = Inf,
                                                 v = arrivals)
      }
    }
    if (tergmLite == TRUE) {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = nArrivals)
    }
  }


  ## Departures
  if (length(departures) > 0) {
    if (tergmLite == FALSE) {
      if (isTERGM == TRUE) {
        dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at,
                                           terminus = Inf, v = departures,
                                           deactivate.edges = TRUE)
      } else {
        dat$nw[[1]] <- delete.vertices(dat$nw[[1]], vid = departures)
        dat <- delete_attr(dat, departures)
      }
    }
    if (tergmLite == TRUE) {
      dat <- delete_attr(dat, departures)
      dat$el[[1]] <- delete_vertices(dat$el[[1]], departures)

      if (dat$control$tergmLite.track.duration) {
        dat$p[[1]]$state$nw0 %n% "lasttoggle" <-
          delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", departures)
      }
    }
  }

  ## Copy static attributes to network object
  if (tergmLite == FALSE & resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }

  ## Update temporally extended disease status
  if (tergmLite == FALSE & isTERGM == TRUE) {
    dat$nw[[1]] <- activate.vertex.attribute(dat$nw[[1]],
                                             prefix = "testatus",
                                             value = status,
                                             onset = at,
                                             terminus = Inf)
  }

  # Record network in nw_list for x-sect ERGM simulations
  if (tergmLite == FALSE & isTERGM == FALSE & resimulate.network == TRUE) {
    dat$temp$nw_list[[at]] <- dat$nw[[1]]
  }
  if (tergmLite == FALSE & isTERGM == FALSE & resimulate.network == FALSE) {
    dat$temp$nw_list[[at]] <- set_vertex_attribute(dat$temp$nw_list[[at]],
                                                   "status", status)
  }

  ## Output
  return(dat)
}
