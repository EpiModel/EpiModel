
#' @title Extract Summary Statistics of Networks Used in netsim
#'
#' @description This function calls \code{summary} on each network being
#'              simulated in \code{netsim}, provided \code{save.nwstats} and
#'              \code{resimulate.network} are both \code{TRUE}. It records the
#'              statistics represented by \code{nwstats.formula} in
#'              \code{dat$stats$nwstats}, where \code{dat} is the main
#'              \code{netsim_dat} class object.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
summary_nets <- function(dat, at) {
  if (get_control(dat, "save.nwstats") == TRUE &&
      get_control(dat, "resimulate.network") == TRUE) {
    for (network in seq_len(dat$num.nw)) {
      nwstats <- summary(get_network_control(dat, network, "nwstats.formula"),
                         basis = get_network(dat, network = network),
                         at = at, # needed for networkDynamic case
                         dynamic = TRUE,
                         term.options = get_network_control(dat, network, "set.control.tergm")$term.options)
      if (is(nwstats, "matrix")) {
        nwstats <- nwstats[, !duplicated(colnames(nwstats)), drop = TRUE]
      } else {
        nwstats <- nwstats[!duplicated(names(nwstats))]
      }
      start_time <- get_control(dat, "start")
      if (start_time > 1L) {
        loc <- at - start_time + 2L
      } else {
        loc <- at
      }
      dat$stats$nwstats[[network]][[loc]] <- nwstats
    }
  }
  return(dat)
}
