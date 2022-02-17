#' @title Record Attribute History
#'
#' @description
#' This function records values specific to a time-step and a group of nodes.
#' The nodes are identified by their \code{unique_ids}, which allows the
#' recording of data for nodes that are no longer in the network by the end of
#' the run. The records are stored in \code{dat[["attr.history"]]} and can be
#' accessed from the \code{netsim} object with \code{get_attr_history}.
#'
#' @param dat Main data object passed through \code{netsim} simulations.
#' @param at The time where the recording happens.
#' @param attribute The name of the value to record.
#' @param posit_ids A numeric vector of posit_ids to which the measure applies.
#'   (see \code{get_posit_ids}).
#' @param values The values to be recorded.
#'
#' @return The main list object.
#'
#'
#' @details
#' See the "Time Varying Attributes" vignette.
#'
#' @examples
#' \dontrun{
#'
#' dat <- record_attr_history(dat, at, "attr_1", get_posit_ids(dat), 5)
#' some_nodes <- get_posit_ids(dat)
#' some_nodes <- some_nodes[runif(length(some_nodes)) < 0.2]
#' dat <- record_attr_history(
#'   dat, at,
#'   "attr_2",
#'   some_nodes,
#'   rnorm(length(some_nodes))
#' )
#'
#' }
#'
#' @export
record_attr_history <- function(dat, at, attribute, posit_ids, values) {
  if (is.null(dat[["attr.history"]])) {
    dat[["attr.history"]] <- list()
  }

  if (length(values) != 1 && length(values) != length(posit_ids)) {
    stop(
      "When trying to record a value for `", attribute, "` at time ", at,
      "The size of the `values` vector is not equal to the number of nodes ",
      "selected by the `posit_ids` vector nor of length 1. \n",
      "Expected: ", length(posit_ids), " or 1 \n",
      "Given: ", length(posit_ids)
    )
  }

  element <- list(at, attribute, get_unique_ids(dat, posit_ids), values)
  names(element) <- c("time", "attribute", "uids", "values")

  dat[["attr.history"]] <- append(
    dat[["attr.history"]],
    list(element)
  )

  return(dat)
}

#' @title Record an Arbitrary Object During a Simulation
#'
#' @description
#' This function records any object during a simulation to allow its
#' inspection afterward. The records are stored in \code{dat[["raw.records"]]}
#' during the simulation and in the \code{netsim} object under the
#' \code{raw.records} sublists.
#'
#' @param dat Main data object passed through \code{netsim} simulations.
#' @param at The time where the recording happens.
#' @param label The name to give to the recorded object.
#' @param object The object to be recorded.
#'
#' @return The main list object.
#'
#' @details
#' See the "Time Varying Attributes" vignette.
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
    dat[["raw.records"]] <- list()
  }

  element <- list(at, label, object)
  names(element) <- c("time", "label", "object")

  dat[["raw.records"]] <- append(
    dat[["raw.records"]],
    list(element)
  )

  return(dat)
}
