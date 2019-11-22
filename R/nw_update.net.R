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
      nw <- activate.vertex.attribute(nw,
                                      prefix = "testatus",
                                      value = "i",
                                      onset = at,
                                      terminus = Inf,
                                      v = idsNewInf)
    }
    dat$attr$status[idsNewInf] <- "i"
    dat$attr$infTime[idsNewInf] <- at

    if ("status" %in% dat$temp$fterms) {
      nw <- set.vertex.attribute(nw, "status", dat$attr$status)
    }
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
