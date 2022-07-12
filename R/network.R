
#' @title Initialize Network Object
#'
#' @description Initialize an undirected network object for use in EpiModel
#'              workflows.
#'
#' @param n Network size.
#'
#' @details
#' This function is used in \code{EpiModel} workflows to initialize an empty
#' network object.  The network attributes \code{directed}, \code{bipartite},
#' \code{hyper}, \code{loops}, and \code{multiple} are set to \code{FALSE}.
#'
#' @return
#' Returns an object of class \code{network}.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' nw
#'
network_initialize <- function(n) {
  nw <- network.initialize(n, directed = FALSE, hyper = FALSE, loops = FALSE,
                           multiple = FALSE, bipartite = FALSE)
  return(nw)
}

#' @title Set Vertex Attribute on Network Object
#'
#' @description Sets a vertex attribute on an object of class \code{network}.
#'              This function simplifies the related function in the
#'              \code{network} package.
#'
#' @param x An object of class network.
#' @param attrname The name of the attribute to set.
#' @param value A vector of values of the attribute to be set.
#' @param v IDs for the vertices whose attributes are to be altered.
#'
#' @details
#' This function is used in \code{EpiModel} workflows to set vertex attributes
#' on an initialized empty network object (see \code{\link{network_initialize}}.
#'
#' @return
#' Returns an object of class \code{network}.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
#' nw
#'
set_vertex_attribute <- function(x, attrname, value, v) {
  if (missing(v)) {
    v <- seq_len(network.size(x))
  }
  g <- set.vertex.attribute(x, attrname, value, v)
  return(g)
}

#' @title Get Vertex Attribute on Network Object
#'
#' @description Gets a vertex attribute from an object of class \code{network}.
#'              This functions simplifies the related function in the
#'              \code{network} package.
#'
#' @param x An object of class network.
#' @param attrname The name of the attribute to get.
#'
#' @details
#' This function is used in \code{EpiModel} workflows to query vertex
#' attributes on an initialized empty network object (see
#' \code{\link{network_initialize}}).
#'
#' @return
#' Returns an object of class \code{network}.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
#' get_vertex_attribute(nw, "age")
#'
get_vertex_attribute <- function(x, attrname) {
  attr <- get.vertex.attribute(x, attrname, na.omit = FALSE,
                               null.na = TRUE, unlist = TRUE)
  return(attr)
}
