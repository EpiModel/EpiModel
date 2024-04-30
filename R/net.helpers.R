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
#' @param run A \code{list} that will contains the objects created by
#' \code{\link{netsim}} that are required for between step communication. This
#' list must be preserved for restarting models.
#'
#' @return A \code{netsim_dat} main list object.
#' @export
create_dat_object <- function(param = list(), init = list(), control = list(),
                              run = list()) {
  dat <- list(
    "run"       = run,
    "param"     = param,
    "control"   = control,
    "epi"       = list(),
    "temp"      = list(),
    "stats"     = list(),
    "init"      = init
  )

  dat <- validate_run(dat)

  class(dat) <- c("netsim_dat", class(dat))

  return(dat)
}

#' Ensures that the `run` sublist contains all the mandatory elements
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
validate_run <- function(dat) {
  defaults <- list(
    current_timestep = 1L,
    last_unique_id = 0L
  )

  for (i in seq_along(defaults)) {
    if (is.null(dat$run[[ names(defaults)[i] ]])) {
      dat$run[[ names(defaults)[i] ]] <- defaults[[i]]
    }
  }

  return(dat)
}

#' @title Return the Current Timestep
#'
#' @inheritParams recovery.net
#'
#' @return The current timestep.
#' @export
get_current_timestep <- function(dat) {
  return(dat$run$current_timestep)
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
  dat$run$current_timestep <- timestep
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
  dat$run$current_timestep <- dat$run$current_timestep + 1
  return(dat)
}
