#' @title Cumulative Edgelist: netsim Module
#'
#' @description This function generate the cumulative edgelists for each network
#'              in the model.
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
#' @section Truncation:
#' To avoid storing a cumulative edgelist to long, the "truncate.el_cuml"
#' control value defines a number of steps after which an edge that is no longer
#' active is truncated out of the cummulative edgelist. (see \code{control.net})
#'
#' @section Accessing the Cumulative Edgelists:
#' This module only computes and store the cumulative edgelists. They can then
#' be accessed using \code{\link{get_cumulative_edgelist}}
#' \code{\link{get_cumulative_edgelists_df}}
#
#' @seealso \code{\link{get_cumulative_edgelist}}
#' \code{\link{get_cumulative_edgelists_df}} \code{\link{netsim}}
#'
#' @export
#' @keywords netMod internal
#'
cumulative_edgelist.net <- function(dat, at) {
  for (network in seq_along(dat[["nwparam"]])) {
    dat <- update_cumulative_edgelist(dat, at, network)
  }

  return(dat)
}
