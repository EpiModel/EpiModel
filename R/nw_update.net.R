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
nw_update.net <- function(dat, at) {

  #Resimulate Network----

  #Deactive inactive nodes
  inactive <- dat$nw.update$resim$inactive
  if (length(inactive) > 0) {
    dat$attr <- deleteAttr(dat$attr, inactive)
  }

  #Departures----

  if (dat$param$vital != FALSE) {
    idsDpt <- unlist(dat$nw.update$dep)
    #Deactive all departures on the network -
    dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                  v = idsDpt, deactivate.edges = TRUE)
  }




  #Output-----
  return(dat)
}
