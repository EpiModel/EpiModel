#' @title Update List \code{x} Using the Elements of List \code{new.x}.
#'
#' @param x A list.
#' @param new.x A list.
#'
#' @return The full \code{x} list with the modifications added by \code{new.x}.
#'
#' @details
#' This function updates list \code{x} by name. If \code{x} and \code{new.x}
#' elements are not named, the function will not work properly. If a function is
#' provided to replace an element that was originally not a function, this
#' function will be applied to the original value.
#'
#' @keywords internal
update_list <- function(x, new.x) {
  for (nm in names(new.x)) {
    if (is.list(new.x[[nm]])) {
      x[[nm]] <- update_list(x[[nm]], new.x[[nm]])
    } else if (is.function(new.x[[nm]]) && !is.function(x[[nm]])) {
      x[[nm]] <- new.x[[nm]](x[[nm]])
    } else {
      x[[nm]] <- new.x[[nm]]
    }
  }

  return(x)
}

#' @title Function to Modify the Controls or Parameters During the Simulation
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @details
#' If a list \code{.param.updater.list} is present in the parameters, this
#' function will update the \code{param} list with new values at given time
#' steps. Similarly, if a list \code{.control.updater.list} is present in the
#' controls, this function will update the \code{param} list with new values at
#' given time steps. An updater is a list containing an \code{at} element
#' governing when the changes will happen, an optional \code{verbose} Boolean
#' controlling whether to output a message when a change is made (default =
#' TRUE), and a \code{param} or \code{control} named list containing the new
#' values(s) with which to update, with the names being the same as the
#' parameter / control names. If the new value is a function but the old one is
#' not, the function will be applied to the current element (see example).
#'
#' @examples
#' \dontrun{
#'
#' # Create the .param.updater.list
#' .param.updater.list <- list(
#'   # this is one updater
#'   list(
#'     at = 10,
#'     param = list(
#'       hiv.test.rate = rep(0.0128, 3),
#'       trans.scale = c(1.61, 0.836, 0.622)
#'     )
#'   ),
#'   # this is another updater
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
#'  # Add it to params
#'  param <- param.net(
#'    inf.prob = 0.3,
#'    act.rate = 0.5,
#'    hiv.test.rate = rep(0.256, 3),
#'    trans.scale = c(1, 2, 3),
#'    .param.updater.list = param.updater.list
#'  )
#'
#' # Create the .control.updater.list
#' # these updaters will toggle on and off the verbosity of the model
#' control.updater.list <- list(
#'   list(
#'     at = 50,
#'     verbose = TRUE,
#'     control = list(
#'       verbose = TRUE
#'     )
#'   ),
#'   # this is another updater
#'   list(
#'     at = 75,
#'     verbose = TRUE,
#'     control = list(
#'       verbose = FALSE
#'     )
#'   )
#' )
#'
#' Add the control updaters to the `control` object
#'  control <- control.net(
#'    nsims = 1,
#'    nsteps = 20,
#'    verbose = FALSE,
#'    .control.updater.list = control.updater.list
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
#' @keywords internal
#' @noRd
input_updater <- function(dat) {
  for (type in c("param", "control")) {
    dat <- common_updater(dat, type)
  }

  return(dat)
}

#' @title Update Either the "param" or "control" List
#'
#' @inheritParams recovery.net
#' @param type Either \code{"param"} or \code{"control"}.
#'
#' @inherit recovery.net return
#'
#' @keywords internal
common_updater <- function(dat, type) {
  # Set the variables and functions for either `param` or `control`
  if (type == "param") {
    type.label <- "parameters"
    type.set <- set_param
    type.get_list <- get_param_list
    updaters.label <- ".param.updater.list"
    updaters <- get_param(dat, updaters.label, override.null.error = TRUE)
  } else if (type == "control") {
    type.label <- "controls"
    type.set <- set_control
    type.get_list <- get_control_list
    updaters.label <- ".control.updater.list"
    updaters <- get_control(dat, updaters.label, override.null.error = TRUE)
  } else {
    stop("`type` must be either 'param' or 'control'")
  }

  # Common update mechanism
  if (is.null(updaters) || length(updaters) == 0) {
    return(dat)
  }

  used.updaters <- numeric(0)
  at <- get_current_timestep(dat)

  for (i in seq_along(updaters)) {
    if (updaters[[i]][["at"]] == at) {
      verbose <- updaters[[i]][["verbose"]]
      verbose <- if (is.null(verbose)) TRUE else verbose

      new.list <- updaters[[i]][[type]]

      if (verbose) {
        message(
          "\n\nAt timestep = ", at, " the following ", type.label,
          " were modified:",
          "\n'", paste0(names(new.list), collapse = "', '"), "'"
        )
      }

      old.list <- type.get_list(dat, names(new.list))
      updated.list <- update_list(old.list, new.list)

      for (nm in names(updated.list)) {
        dat <- type.set(dat, nm, updated.list[[nm]])
      }

      used.updaters <- c(used.updaters, i)
    }
  }

  # Remove the used updaters from the list
  if (length(used.updaters) > 0) {
    dat <- type.set(
      dat, updaters.label,
      updaters[-used.updaters]
    )
  }

  return(dat)
}
