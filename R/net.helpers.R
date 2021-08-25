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
    "param"    = param,
    "init"     = init,
    "control"  = control,

    "attr"     = list(),
    "epi"      = list(),
    "stats"    = list(),
    "temp"     = list()
  )

  return(dat)
}
