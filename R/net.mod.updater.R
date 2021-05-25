#' Update list code{x} using the elements of list code{new_x}
#'
#' @param x a list
#' @param new_x a list
#'
#' @return the full code{x} list with the modifications added by code{new_x}
#'
#' @details
#' This function updates list code{x} by name. If code{x} and code{new_x} elements are not
#' named, the function will not work properly.
#' If a function is provided to replace an element that was originaly not a
#' function, this function will be applied to the original value.
#'
#' @keywords internal
update_list <- function(x, new_x) {
  for (nm in names(new_x)) {
    if (is.list(new_x[[nm]])) {
      x[[nm]] <- update_list(x[[nm]], new_x[[nm]])
    } else if (is.function(new_x[[nm]]) && ! is.function(x[[nm]])) {
      x[[nm]] <- new_x[[nm]](x[[nm]])
    } else {
      x[[nm]] <- new_x[[nm]]
    }
  }

  return(x)
}

#' Module to modify the parameters during the simulation
#'
#' @param dat Master list object containing a \code{networkDynamic} object and
#'        other initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @return the updated dat Master list object.
#'
#' @details
#' if a list code{updaters} is present in the parameters, this module will
#' update the code{param} list with new values at given timesteps.
#' An updater is a list containing an code{at} element governing when the
#' changes will happen, an optional code{verbose} boolean controlling whether to
#' output a message when a change is made (default = TRUE) and a code{param}
#' named list with the names being the same as the paramter names and the new
#' value to update with.
#' If the new value is a function but the old one is not, the
#' function will be applied to the current element (see example) .
#'
#' @examples
#' \dontrun{
#'
#' updaters = list(
#'   list(
#'     at = 10,
#'     param = list(
#'       hiv.test.rate = rep(0.0128, 3),
#'       trans.scale = c(1.61, 0.836, 0.622)
#'     )
#'   ),
#'   list(
#'     at = 12,
#'     verbose = TRUE,
#'     param = list(
#'       hiv.test.rate = function(x) x * 3,
#'       trans.scale = function(x) x^2 / 3
#'     )
#'   )
#' )
#'
#'  # Do not forget to add it to `param`
#'  param <- param.net(
#'    inf.prob = 0.3,
#'    act.rate = 0.5,
#'    hiv.test.rate = rep(0.256, 3),
#'    trans.scale = c(1, 2, 3),
#'    updaters = updaters
#'  )
#'
#' # Enable the module in `control`
#'  control <- control.net(
#'    type = NULL, # must be NULL as we use a custom module
#'    nsims = 1,
#'    nsteps = 20,
#'    verbose = FALSE,
#'    updater.FUN = updater.net
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
#' }
#'
updater.net <- function(dat, at) {
  updaters <- get_param(dat, "updaters", override.null.error = TRUE)
  if (is.null(updaters))
    return(dat)

  for (i in seq_along(updaters)) {
    if (updaters[[i]][["at"]] == at) {
      verbose <- updaters[[i]][["verbose"]]
      verbose <- if (is.null(verbose)) TRUE else verbose

      new_params <- updaters[[i]][["param"]]

      if (verbose) {
        message(
          "\n\nAt time step = ", at, " the following parameters where modified:",
          "\n'", paste0(names(new_params), collapse = "', '"), "'"
        )
      }

      old_params <- get_param_list(dat, names(new_params))
      updated_params <- update_list(old_params, new_params)

      for(nm in names(updated_params)) {
        dat <- set_param(dat, nm, updated_params[[nm]])
      }
    }
  }

  return(dat)
}
