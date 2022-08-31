
#' @title Extract Summary Statistics of Networks Used in netsim
#'
#' @description This function calls \code{summary} on each network being
#'              simulated in \code{netsim}, provided \code{save.nwstats} and
#'              \code{resimulate.network} are both \code{TRUE}. It records the
#'              statistics represented by \code{nwstats.formula} in
#'              \code{dat$stats$nwstats}.
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
    for (network in seq_along(dat$nwparam)) {
      nwstats <- summary(get_network_control(dat, network, "nwstats.formula"),
                         basis = make_sim_network(dat, network = network),
                         at = at, # needed for networkDynamic case
                         term.options = get_network_control(dat, network, "set.control.tergm")$term.options)
      if (is(nwstats, "matrix")) {
        nwstats <- nwstats[, !duplicated(colnames(nwstats)), drop = FALSE]
      } else {
        nwstats <- nwstats[!duplicated(names(nwstats))]
      }
      dat$stats$nwstats[[network]] <- rbind(dat$stats$nwstats[[network]], nwstats)
      rownames(dat$stats$nwstats[[network]]) <- NULL
    }
  }
  return(dat)
}
