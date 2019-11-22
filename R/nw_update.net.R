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
    dat$attr$status[idsNewInf] <- "i"
    dat$attr$infTime[idsNewInf] <- at

    if ("status" %in% dat$temp$fterms) {
      dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
    }
  }

  if (dat$param$vital != FALSE) {
    #Departures----

    idsDpt <- unlist(dat$nw.update$dpt$idsDpt)
    idsDpt <- as.vector(idsDpt)
    #Deactive all departures on the network -

    if (length(idsDpt) > 0) {
    dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                  v = idsDpt, deactivate.edges = TRUE)
    }

    #Arrivals----
  }




  #Output-----
  return(dat)
}
