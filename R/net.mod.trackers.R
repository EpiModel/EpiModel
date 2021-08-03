#' @title Trackers: netsim Module
#'
#' @description This function apply the user provided epi trackers
#'
#' @param dat Master list object containing a \code{networkDynamic} object and
#'        other initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @return the updated dat Master list object.
#'
#' @section Optional Module:
#' This module is not included by default
#'
#' @section The \code{tracker.list} list:
#' \code{tracker.list} is a list of NAMED functions stored in the \code{param}
#' list of the \code{dat} master list object.
#'
#' @section Tracker Functions:
#' This module will apply the tracker functions present in the parameter list
#' \code{tracker.list}. Each tracker must be a function with EXACTLY two
#' arguments: the \code{dat} Master list object and \code{at} the current time
#' step. They must return a VALUE of length one (numeric, logical or character).
#'
#' @examples
#' \dontrun{
#'
#' # Create some trackers
#' epi_prop_infected <- function(dat, at) {
#'   needed_attributes <- c("status", "active")
#'   output <- with(get_attr_list(dat, needed_attributes), {
#'     pop <- active == 1
#'     cond <- status == "i"
#'
#'     out <- sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
#'
#'     out
#'   })
#'
#'   return(output)
#' }
#'
#' epi_s_num <- function(dat, at) {
#'   needed_attributes <- c("status")
#'   output <- with(get_attr_list(dat, needed_attributes), {
#'     out <- sum(status == "s", na.rm = TRUE)
#'
#'     out
#'   })
#'
#'   return(output)
#' }
#'
#' # Create the `tracker.list` list
#' tracker.list <- list(
#'   prop_infected = epi_prop_infected,
#'   s_num         = epi_s_num
#' )
#'
#'  # Do not forget to add it to `param`
#'  param <- param.net(
#'    inf.prob = 0.3,
#'    act.rate = 0.5,
#'    tracker.list = tracker.list
#'  )
#'
#' # Enable the module in `control`
#'  control <- control.net(
#'    type = NULL, # must be NULL as we use a custom module
#'    nsims = 2,
#'    nsteps = 5,
#'    verbose = FALSE,
#'    trackers.FUN = trackers.net
#'  )
#'
#' nw <- network_initialize(n = 50)
#' nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
#' est <- netest(
#'   nw,
#'   formation = ~edges,
#'   target.stats = 25,
#'   coef.diss = dissolution_coefs(~offset(edges), 10, 0),
#'   verbose = FALSE
#' )
#'
#' init <- init.net(i.num = 10)
#' mod <- netsim(est, param, init, control)
#'
#' df <- as.data.frame(mod)
#'
#' df
#'
#' }
#'
#' @seealso \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
trackers.net <- function(dat, at) {
  tracker.list <- get_param(dat, "tracker.list", override.null.error = TRUE)

  if (is.null(tracker.list)) {
    return(dat)
  }

  tryCatch(
    expr = {
      for (nm in names(tracker.list)) {
        dat <- set_epi(dat, nm, at, tracker.list[[nm]](dat, at))
      }
    },
    message = function(e) message("\nIn tracker '", nm, "':\n", e),
    warning = function(e) warning("\nIn tracker '", nm, "':\n", e),
    error   = function(e) stop("\nIn tracker '", nm, "':\n", e)
  )

  return(dat)
}
