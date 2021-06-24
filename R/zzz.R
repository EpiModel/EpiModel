.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(c(strwrap(
      paste("EpiModel NOTE: EpiModel 2.0 implements significant changes to the
            EpiModel application programming interface (API), which may require
            updates to model code written under EpiModel 1.x. See https://epimodel.org/ for migration guidance.", sep = "")), ""),
      collapse = "\n"))
}

.onLoad <- function(libname, pkgname) {
  eval(COLLATE_ALL_MY_CONTROLS_EXPR)
}

# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists,
#' reexported from `statnet.common`.
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL

#' @export
snctrl <- statnet.common::snctrl
## BEGIN boilerplate: should be kept in sync with statnet.common.


eval(UPDATE_MY_SCTRL_EXPR)
