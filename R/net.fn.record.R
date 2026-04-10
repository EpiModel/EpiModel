#' @title Record Attribute History
#'
#' @description
#' Record values for a set of nodes at a specific time step. Records are stored
#' in `dat[["attr.history"]]` as a `collections::queue` and can be accessed
#' from the `netsim` object with `get_attr_history`.
#'
#' This function is called automatically by the tracked attributes system when
#' `tracked.attributes` is set in `control.net`, but can also be called
#' manually inside custom modules for ad-hoc recording.
#'
#' @inheritParams recovery.net
#' @param item The name of the attribute to record.
#' @param value The values to be recorded. Must be of length 1 (recycled for
#'   all nodes) or the same length as the number of nodes identified by
#'   `posit_ids` or `unique_ids`.
#' @param at The time step at which the recording happens. Defaults to the
#'   current time step if `NULL`.
#' @param posit_ids A numeric vector of positional IDs to which the values
#'   apply. Converted internally to unique IDs. Either `posit_ids` or
#'   `unique_ids` must be provided.
#' @param unique_ids A numeric vector of unique IDs to which the values apply.
#'   When provided, skips the positional-to-unique ID conversion.
#'
#' @inherit recovery.net return
#'
#' @details
#' Exactly one of `posit_ids` or `unique_ids` must be non-`NULL`. When both
#' are `NULL`, the function errors. When `unique_ids` is provided directly,
#' the conversion from positional IDs is skipped, which is more efficient
#' when unique IDs are already known (e.g., from the automatic tracking
#' system).
#'
#' See the "Time-Varying Parameters" section of the "Working With Model
#' Parameters" vignette.
#'
#' @examples
#' \dontrun{
#' # Manual recording inside a custom module
#' dat <- record_attr_history(dat, "status", get_attr(dat, "status"),
#'                            posit_ids = get_posit_ids(dat))
#'
#' # Record for a subset of nodes using unique IDs directly
#' some_uids <- c(10, 25, 42)
#' dat <- record_attr_history(dat, "risk", c(1, 0, 1),
#'                            unique_ids = some_uids)
#' }
#'
#' @export
record_attr_history <- function(dat, item, value, at = NULL, posit_ids = NULL,
                                unique_ids = NULL) {
  if (is.null(dat[["attr.history"]]))
    dat[["attr.history"]] <- collections::queue()

  if (is.null(unique_ids)) {
    if (is.null(posit_ids)) {
      stop("Either `unique_ids` or `posit_ids` must be non NULL")
    } else {
      unique_ids <- get_unique_ids(dat, posit_ids)
    }
  }

  if (is.null(at))
    at <- get_current_timestep(dat)

  if (length(value) != 1 && length(value) != length(unique_ids)) {
    stop(
      "When trying to record a value for `",
      item,
      "` at time ",
      at,
      "The size of the `values` vector is not equal to the number of nodes ",
      "selected by the `posit_ids` or `unique_ids` vector nor of length 1. \n",
      "Expected: ",
      length(unique_ids),
      " or 1 \n",
      "Given: ",
      length(value)
    )
  }

  element <- list(at, item, unique_ids, value)
  names(element) <- c("time", "attribute", "uids", "values")

  dat[["attr.history"]]$push(element)
  return(dat)
}

#' @title Record an Arbitrary Object During a Simulation
#'
#' @description
#' This function records any object during a simulation to allow its
#' inspection afterward. The records are stored in `dat[["raw.records"]]`
#' during the simulation, where `dat` is the main `netsim_dat` class
#' object, and in the `netsim` object under the `raw.records`
#' `collections::queue` object.
#'
#' @inheritParams recovery.net
#' @param at The time where the recording happens.
#' @param label The name to give to the recorded object.
#' @param object The object to be recorded.
#'
#' @inherit recovery.net return
#'
#' @details
#' See the "Time-Varying Parameters" section of the "Working With Model
#' Parameters" vignette.
#'
#' @examples
#' \dontrun{
#'
#' dat <- record_raw_object(dat, at, "a.df", data.frame(x = 2:200))
#' dat <- record_raw_object(dat, at, "a.message", "I recorded something")
#'
#' }
#'
#' @export
record_raw_object <- function(dat, at, label, object) {
  if (is.null(dat[["raw.records"]])) {
    dat[["raw.records"]] <- collections::queue()
  }

  element <- list(at, label, object)
  names(element) <- c("time", "label", "object")

  dat[["raw.records"]]$push(element)

  return(dat)
}
