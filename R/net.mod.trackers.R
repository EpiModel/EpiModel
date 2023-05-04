#' @title Function to run the user-provided epi trackers
#'
#' @description
#' see the "Working with Custom Attributes and Summary Statistics in EpiModel"
#' vignette.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @section The \code{tracker.list} list:
#' \code{.tracker.list} is a list of NAMED functions stored in the
#' \code{control} list of the main \code{netsim_dat} class object.
#'
#' @section Tracker Functions:
#' This function will apply the tracker functions present in the control list
#' \code{.tracker.list}. Each tracker must be a function with EXACTLY one
#' argument: the \code{netsim_dat} main list object. They must return a VALUE of
#' length one (numeric, logical or character).
#'
#' @examples
#' \dontrun{
#'
#' # Create some trackers
#' epi_prop_infected <- function(dat) {
#'   # we need two attributes for our calculation: `status` and `active`
#'   needed_attributes <- c("status", "active")
#'   # we use `with` to simplify code
#'   output <- with(EpiModel::get_attr_list(dat, needed_attributes), {
#'     pop <- active == 1             # we only look at active nodes
#'     cond <- status == "i"   # which are infected
#'     # how many are `infected` among the `active`
#'     sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
#'   })
#'   return(output)
#' }
#'
#' epi_s_num <- function(dat) {
#'   needed_attributes <- c("status")
#'   output <- with(get_attr_list(dat, needed_attributes), {
#'     sum(status == "s", na.rm = TRUE)
#'   })
#'   return(output)
#' }
#'
#' # Store the trackers in a named list. The names will be used as column names
#' # for in the `epi` list
#' some.trackers <- list(
#'   prop_infected = epi_prop_infected,
#'   s_num         = epi_s_num
#' )
#'
#' # Make a simple SI model with custom trackers
#' control <- EpiModel::control.net(
#'   type = "SI",
#'   nsims = 1,
#'   nsteps = 50,
#'   verbose = FALSE,
#'   .tracker.list = some.trackers
#' )
#'
#' param <- EpiModel::param.net(
#'   inf.prob = 0.3,
#'   act.rate = 0.1
#' )
#'
#' nw <- network_initialize(n = 50)
#' nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
#' est <- EpiModel::netest(
#'   nw,
#'   formation = ~edges,
#'   target.stats = 25,
#'   coef.diss = dissolution_coefs(~offset(edges), 10, 0),
#'   verbose = FALSE
#' )
#'
#' init <- EpiModel::init.net(i.num = 10)
#' sim <- EpiModel::netsim(est, param, init, control)
#'
#' d <- as.data.frame(sim)
#' d
#' }
#'
#' @seealso \code{\link{netsim}}
#'
#' @keywords internal
epi_trackers <- function(dat) {
  tracker.list <- get_control(dat, ".tracker.list", override.null.error = TRUE)

  if (is.null(tracker.list)) {
    return(dat)
  }

  at <- get_current_timestep(dat)

  tryCatch(
    expr = {
      for (nm in names(tracker.list)) {
        dat <- set_epi(dat, nm, at, tracker.list[[nm]](dat))
      }
    },
    message = function(e) message("\nIn tracker '", nm, "':\n", e),
    warning = function(e) warning("\nIn tracker '", nm, "':\n", e),
    error   = function(e) stop("\nIn tracker '", nm, "':\n", e)
  )

  return(dat)
}
