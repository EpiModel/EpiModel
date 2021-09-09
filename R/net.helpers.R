#' @title Create a Minimal Master List Object of Network Model
#'
#' @description This helper function populate a \code{dat} Master List object
#'              with the minimal required elements. All parameters are optional.
#'              When none is given the resulting object is only a shell list
#'              with the different named elements defined as empty lists.
#'
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @return A \code{dat} Master list object
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
    "_timestep" = 1
  )

  return(dat)
}

#' @title Returns the Current Timestep
#'
#' @param dat a Master list object of network models
#'
#' @return the current timestep
#' @export
get_current_timestep <- function(dat) {
  return(dat[["_timestep"]])
}

#' @title Sets the Current Timestep
#'
#' @description Changes the current timestep in the \code{dat} object. Use with
#'              caution. This function exists to workaround unforseen corner
#'              cases. In most situation, \code{increment_timestep} should be
#'              prefered
#'
#' @param dat a Master list object of network models
#' @par timestep the new value for the timestep
#'
#' @return A \code{dat} Master list object
#'
#' @section Mutability:
#' This DOES NOT modify the dat object in place. The result must be assigned
#' back to \code{dat} in order to be registered
#' \code{dat <- increment_timestep(dat)}
#'
#' @export
set_current_timestep <- function(dat, timestep) {
  dat[["_timestep"]] <- timestep
  return(dat)
}

#' @title Increment the Current Timestep
#'
#' @description This function adds 1 to the timestep counter stored in the
#'              \code{dat} Master list object.
#'
#' @param dat a Master list object of network models
#'
#' @return A \code{dat} Master list object
#'
#' @section Mutability:
#' This DOES NOT modify the dat object in place. The result must be assigned
#' back to \code{dat} in order to be registered
#' \code{dat <- increment_timestep(dat)}
#'
#' @export
increment_timestep <- function(dat) {
  dat[["_timestep"]] <- dat[["_timestep"]] + 1
  return(dat)
}
