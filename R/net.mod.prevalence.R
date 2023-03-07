
#' @title Get Epidemic Output from netsim Model
#'
#' @description Provides all active model state sizes from the network at the
#'              specified time step, output to a list of vectors.
#'
#' @inheritParams recovery.net
#'
#' @details
#' This network utility is used during the \code{\link{netsim}} simulation
#' process to efficiently query the current size of each state or compartment
#' in the model at any given timestep. For a two-group network, the current
#' state size for each group and overall is provided.
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
prevalence.net <- function(dat, at) {
  active <- get_attr(dat, "active")

  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) "None" else type

  groups <- get_param(dat, "groups")
  if (groups > 1) {
    group <- get_attr(dat, "group")[get_attr(dat, "active") == 1]
  } else {
    group <- rep(1, length.out = sum(active == 1))
  }

  suffixes <- c("", if (groups > 1) paste0(".g", seq_len(groups)[-1]))
  statuses <- c("s", "i", if (type == "SIR") "r")

  status <- get_attr(dat, "status")[get_attr(dat, "active") == 1]

  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- get_control(dat, "epi.by")
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    if (ebn %in% c("active", "infTime")) {
      ebvec <- NULL
    } else {
      ebvec <- get_attr(dat, ebn)[get_attr(dat, "active") == 1]
    }
  }

  for (g in seq_len(groups)) {
    for (status_value in statuses) {
      dat <- set_epi(dat, paste0(status_value, ".num", suffixes[g]), at,
                     sum(status == status_value & group == g))
      if (eb == TRUE) {
        for (i in seq_along(ebun)) {
          ebn.temp <- paste0(status_value, ".num", suffixes[g], ebun[i])
          dat <- set_epi(dat, ebn.temp, at, sum(status == status_value &
            group == g & ebvec == ebv[i]))
        }
      }
    }

    dat <- set_epi(dat, paste0("num", suffixes[g]), at, sum(group == g))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("num", suffixes[g], ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(ebvec == ebv[i] & group == g))
      }
    }
  }

  return(dat)
}
