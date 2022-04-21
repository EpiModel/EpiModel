
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

  # Subset attr to active == 1
  l <- lapply(seq_along(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL

  status <- l$status

  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- get_control(dat, "epi.by")
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }

  if (at == 1) {
    dat$epi <- list()
  }

  if (groups == 1) {
    dat <- set_epi(dat, "s.num", at, sum(status == "s"))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("s.num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "s" &
          get(ebn) == ebv[i]))
      }
    }

    dat <- set_epi(dat, "i.num", at, sum(status == "i"))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("i.num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "i" &
          get(ebn) == ebv[i]))
      }
    }

    if (type == "SIR") {
      dat <- set_epi(dat, "r.num", at, sum(status == "r"))
      if (eb == TRUE) {
        for (i in seq_along(ebun)) {
          ebn.temp <- paste0("r.num", ebun[i])
          dat <- set_epi(dat, ebn.temp, at, sum(status == "r" &
            get(ebn) == ebv[i]))
        }
      }
    }
    dat <- set_epi(dat, "num", at, length(status))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(get(ebn) == ebv[i]))
      }
    }
  }

  if (groups == 2) {
    group <- get_attr(dat, "group")
    group <- group[active == 1]

    dat <- set_epi(dat, "s.num", at, sum(status == "s" & group == 1))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("s.num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "s" &
          group == 1 &
          get(ebn) == ebv[i]))
      }
    }
    dat <- set_epi(dat, "i.num", at, sum(status == "i" & group == 1))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("i.num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "i" &
          group == 1 &
          get(ebn) == ebv[i]))
      }
    }
    if (type == "SIR") {
      dat <- set_epi(dat, "r.num", at, sum(status == "r" & group == 1))
      if (eb == TRUE) {
        for (i in seq_along(ebun)) {
          ebn.temp <- paste0("r.num", ebun[i])
          dat <- set_epi(dat, ebn.temp, at, sum(status == "r" &
            group == 1 &
            get(ebn) == ebv[i]))
        }
      }
    }
    dat <- set_epi(dat, "num", at, sum(group == 1))
    dat$epi$num[at] <- sum(group == 1)
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("num", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(group == 1 &
          get(ebn) == ebv[i]))
      }
    }
    dat <- set_epi(dat, "s.num.g2", at, sum(status == "s" & group == 2))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("s.num.g2", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "s" &
          group == 2 &
          get(ebn) == ebv[i]))
      }
    }
    dat <- set_epi(dat, "i.num.g2", at, sum(status == "i" & group == 2))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("i.num.g2", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(status == "i" &
          group == 2 &
          get(ebn) == ebv[i]))
      }
    }
    if (type == "SIR") {
      dat <- set_epi(dat, "r.num.g2", at, sum(status == "r" & group == 2))
      if (eb == TRUE) {
        for (i in seq_along(ebun)) {
          ebn.temp <- paste0("r.num.g2", ebun[i])
          dat <- set_epi(dat, ebn.temp, at, sum(status == "r" &
            group == 2 &
            get(ebn) == ebv[i]))
        }
      }
    }
    dat <- set_epi(dat, "num.g2", at, sum(group == 2))
    if (eb == TRUE) {
      for (i in seq_along(ebun)) {
        ebn.temp <- paste0("num.g2", ebun[i])
        dat <- set_epi(dat, ebn.temp, at, sum(group == 2 &
          get(ebn) == ebv[i]))
      }
    }
  }
  return(dat)
}
