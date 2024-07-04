
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
    nwterms <- dat$run$nwterms
    if (!is.null(nwterms)) {
      curr.tab <- get_attr_prop(dat, nwterms)
      dat <- auto_update_attr(dat, arrivals, curr.tab)
    }
    if (length(unique(vapply(get_attr_list(dat), length, 1))) != 1) {
      stop("Attribute list of unequal length. Check arrivals.net module.\n",
           print(cbind(vapply(get_attr_list(dat), length, 1))))
    }
    dat <- arrive_nodes(dat, nArrivals)
  }


  ## Departures
  dat <- depart_nodes(dat, departures)

  ## Copy static attributes to network object
  if (tergmLite == FALSE && resimulate.network == TRUE) {
    dat <- copy_datattr_to_nwattr(dat)
  }

  ## Update temporally extended disease status
  if (tergmLite == FALSE) {
    for (network in seq_len(dat$num.nw)) {
      dat$run$nw[[network]] <- activate.vertex.attribute(dat$run$nw[[network]],
                                                     prefix = "testatus",
                                                     value = status,
                                                     onset = at,
                                                     terminus = Inf)
    }
  }

  ## Output
  return(dat)
}
