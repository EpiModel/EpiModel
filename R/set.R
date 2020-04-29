#' @title Set Network in Master Object List
#'
#' @description Helper function that sets an object of class \code{network}
#'              onto a master data list.
#' @param dat Master data object.
#' @param nw  An object of class \code{network}.
#' @param index Where in the network sublist to save.
#'
#' @return A master data object with specified network saved on to sublist \code{nw}.
#'
#' @keywords network
#' @export
#'
set_network <- function(dat, nw, index) {

  ## Warnings and checks

  if (is.list(dat) == FALSE) {
    stop("Supplied dat object is not of class list",
         call. = FALSE)
  }

  if (is.numeric(index) == FALSE) {
    stop("index not a valid position",
         call. = FALSE)
  }

  ## Full Network or Edgelist
  tergmLite <- get_control(dat, "tergmLite")

  if (tergmLite == FALSE) {

    if (is.network(nw) == FALSE) {
      stop("nw is not of class network",
           call. = FALSE)
    }

    dat[["nw"]][[index]] <- nw
  } else {
    if (is.matrix(nw) == FALSE) {
      stop("nw is not of class matrix",
           call. = FALSE)
    }
    dat[["el"]][[index]] <- nw
  }

  ## Output -----------------------

  return(dat)

}
