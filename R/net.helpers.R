#' @title Create a Minimal netsim_dat Main List Object for a Network Model
#'
#' @description This helper function populates a `netsim_dat` main list
#'              object with the minimal required elements. All parameters are
#'              optional. When none are given the resulting object is only a
#'              shell list of class `netsim_dat` with the different named
#'              elements defined as empty lists.
#'
#' @param param An `EpiModel` object of class [param.net()].
#' @param init An `EpiModel` object of class [init.net()].
#' @param control An `EpiModel` object of class [control.net()].
#' @param run A `list` that will contains the objects created by
#' [netsim()] that are required for between step communication. This
#' list must be preserved for restarting models.
#'
#' @return A `netsim_dat` main list object.
#' @export
create_dat_object <- function(param = list(), init = list(), control = list(),
                              run = list()) {
  dat <- list(
    "run"       = validate_run(run),
    "param"     = param,
    "control"   = control,
    "epi"       = list(),
    "temp"      = list(),
    "stats"     = list(),
    "init"      = init
  )

  class(dat) <- c("netsim_dat", class(dat))

  return(dat)
}

#' Ensures that the `run` sublist contains all the mandatory elements
#'
#' @param run A `run` sublist to validate
#' @return A valid `run` sublist
#' @noRd
validate_run <- function(run) {
  defaults <- list(
    attr = list(),
    current_timestep = 1L,
    last_unique_id = 0L,
    tracking_attrs = FALSE
  )

  for (elt_name in names(defaults)) {
    if (is.null(run[[elt_name]])) {
      run[[elt_name]] <- defaults[[elt_name]]
    }
  }

  return(run)
}

#' @title Return the Current Timestep
#'
#' @inheritParams recovery.net
#'
#' @return The current timestep.
#' @export
get_current_timestep <- function(dat) {
  return(as.integer(dat$run$current_timestep))
}

#' @title Set the Current Timestep
#'
#' @description Changes the current timestep in the `netsim_dat` object.
#'              Use with caution. This function exists to work around unforeseen
#'              corner cases. In most situation, `increment_timestep` is
#'              preferred.
#'
#' @inheritParams recovery.net
#' @param timestep The new value for the timestep.
#'
#' @inherit recovery.net return
#'
#' @section Mutability:
#' This DOES NOT modify the `netsim_dat` object in place. The result must
#' be assigned back to `dat` in order to be registered:
#' `dat <- increment_timestep(dat)`.
#'
#' @export
set_current_timestep <- function(dat, timestep) {
  dat$run$current_timestep <- as.integer(timestep)
  return(dat)
}

#' @title Increment the Current Timestep
#'
#' @description This function adds 1 to the timestep counter stored in the
#'              `netsim_dat` main list object.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @section Mutability:
#' This DOES NOT modify the `netsim_dat` object in place. The result must
#' be assigned back to `dat` in order to be registered:
#' `dat <- increment_timestep(dat)`.
#'
#' @export
increment_timestep <- function(dat) {
  dat$run$current_timestep <- dat$run$current_timestep + 1L
  return(dat)
}

# Guards for the recording functions whose signatures changed in v2.7.0
#
# `set_epi()`, `record_attr_history()` and `record_raw_object()` moved `at`
# last and made it optional. A call written for the old argument order still
# matches the new positional slots, so instead of failing it would silently
# store the time step as the value (or the value at the wrong time step). The
# three functions take `...` before `at` so that `at` must be named, and any
# argument captured by `...` is an error pointing to the new signature.
#
# @param fun The name of the calling function.
# @param new_usage,old_usage The v2.7.0 and v2.6.x signatures.
# @param fix The recommended call for the most common case.
# @param dots_names The names of the arguments captured by `...`.
# @noRd
stop_old_signature <- function(fun, new_usage, old_usage, fix,
                               dots_names = NULL) {
  renamed <- c(attribute = "item", values = "value")
  hits <- intersect(dots_names, names(renamed))
  msg <- paste0(
    "`", fun, "()` got unexpected arguments. Its signature changed in ",
    "EpiModel 2.7.0 from `", old_usage, "` to `", new_usage, "`: `at` now ",
    "comes last, defaults to the current time step, and must be passed by ",
    "name. Use `", fix, "` to record at the current time step."
  )
  if (length(hits) > 0) {
    msg <- paste0(
      msg, " Renamed arguments: ",
      paste0("`", hits, "` is now `", renamed[hits], "`", collapse = ", "),
      "."
    )
  }
  stop(simpleError(msg, call = sys.call(-1)))
}

# Check that `x` is a single string, as the name or label of a record
# @noRd
check_record_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(simpleError(
      paste0("`", arg, "` must be a single character string"),
      call = sys.call(-1)
    ))
  }
  invisible(x)
}

# Default `at` to the current time step, or check a user-provided one
# @noRd
resolve_record_at <- function(dat, at) {
  if (is.null(at)) {
    at <- get_current_timestep(dat)
    if (length(at) != 1L || is.na(at)) {
      stop(simpleError(
        paste0("`at` is NULL and `dat` holds no current time step ",
               "(`dat$run$current_timestep`); pass `at` explicitly"),
        call = sys.call(-1)
      ))
    }
    return(at)
  }
  if (length(at) != 1L || !is.numeric(at) || is.na(at)) {
    stop(simpleError("`at` must be numeric and of length one",
                     call = sys.call(-1)))
  }
  at
}
