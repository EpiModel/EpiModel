#' @title Create a Minimal netsim_dat Main List Object for a Network Model
#'
#' @description This helper function populates a \code{netsim_dat} main list
#'              object with the minimal required elements. All parameters are
#'              optional. When none are given the resulting object is only a
#'              shell list of class \code{netsim_dat} with the different named
#'              elements defined as empty lists.
#'
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @return A \code{netsim_dat} main list object.
#' @export
create_dat_object <- function(param = list(), init = list(), control = list()) {
  dat <- list(
    "param"     = param,
    "init"      = init,
    "control"   = control,
    "attr"      = list(),
    "epi"       = list(),
    "stats"     = list(),
    "temp"      = list(),
    "run"       = list(),
    "_timestep" = 1
  )

  class(dat) <- c("netsim_dat", class(dat))

  return(dat)
}

#' @title Return the Current Timestep
#'
#' @inheritParams recovery.net
#'
#' @return The current timestep.
#' @export
get_current_timestep <- function(dat) {
  return(dat[["_timestep"]])
}

#' @title Set the Current Timestep
#'
#' @description Changes the current timestep in the \code{netsim_dat} object.
#'              Use with caution. This function exists to work around unforeseen
#'              corner cases. In most situation, \code{increment_timestep} is
#'              preferred.
#'
#' @inheritParams recovery.net
#' @param timestep The new value for the timestep.
#'
#' @inherit recovery.net return
#'
#' @section Mutability:
#' This DOES NOT modify the \code{netsim_dat} object in place. The result must
#' be assigned back to \code{dat} in order to be registered:
#' \code{dat <- increment_timestep(dat)}.
#'
#' @export
set_current_timestep <- function(dat, timestep) {
  dat[["_timestep"]] <- timestep
  return(dat)
}

#' @title Increment the Current Timestep
#'
#' @description This function adds 1 to the timestep counter stored in the
#'              \code{netsim_dat} main list object.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @section Mutability:
#' This DOES NOT modify the \code{netsim_dat} object in place. The result must
#' be assigned back to \code{dat} in order to be registered:
#' \code{dat <- increment_timestep(dat)}.
#'
#' @export
increment_timestep <- function(dat) {
  dat[["_timestep"]] <- dat[["_timestep"]] + 1
  return(dat)
}
