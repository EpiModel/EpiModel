

#' @title Fast Version of network::add.vertices for Edgelist-formatted Network
#'
#' @description This function performs a simple operation of updating the
#'              edgelist attribute \code{n} that tracks the total network
#'              size implicit in an edgelist representation of the network.
#'
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size.
#' @param nv A integer equal to the number of nodes to add to the network
#'           size at the given time step.
#'
#' @details
#' This function is used in \code{EpiModel} modules to add vertices (nodes) to
#' the edgelist object to account for entries into the population (e.g., births
#' and in-migration).
#'
#' @return
#' Returns the matrix of current edges, \code{el}, with the population size
#' attribute updated based on the number of new vertices specified in \code{nv}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("EpiModel")
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        tergmLite = TRUE)
#'
#' # networkLite representation after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#'
#' # Check current network size
#' attributes(dat$el[[1]])$n
#'
#' # Add 10 vertices
#' dat$el[[1]] <- add_vertices(dat$el[[1]], 10)
#'
#' # Check new network size
#' attributes(dat$el[[1]])$n
#' }
#'
add_vertices <- function(el, nv) {
  attributes(el)$n <- attributes(el)$n + nv
  return(el)
}


#' @title Fast Version of network::delete.vertices for Edgelist-formatted
#'        Network
#'
#' @description Given a current two-column matrix of edges and a vector of IDs
#'              to delete from the matrix, this function first removes any rows
#'              of the edgelist in which the IDs are present and then permutes
#'              downward the index of IDs on the edgelist that were numerically
#'              larger than the IDs deleted.
#'
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size.
#' @param vid A vector of IDs to delete from the edgelist.
#'
#' @details
#' This function is used in \code{EpiModel} modules to remove vertices (nodes)
#' from the edgelist object to account for exits from the population (e.g.,
#' deaths and out-migration).
#'
#' @return
#' Returns an updated edgelist object, \code{el}, with the edges of deleted
#' vertices removed from the edgelist and the ID numbers of the remaining edges
#' permuted downward.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library("EpiModel")
#' set.seed(12345)
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        tergmLite = TRUE)
#'
#' # Set seed for reproducibility
#' set.seed(123456)
#'
#' # networkLite representation structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#'
#' # Current edges
#' head(dat$el[[1]], 20)
#'
#' # Remove nodes 1 and 2
#' nodes.to.delete <- 1:2
#' dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes.to.delete)
#'
#' # Newly permuted edges
#' head(dat$el[[1]], 20)
#' }
#'
delete_vertices <- function(el, vid) {

  vid <- sort(vid)

  new.el <- delete_edges(el, vid)
  if (length(vid) > 0 && NROW(new.el) > 0) {
    o1 <- order(new.el[, 1])
    new.el[, 1] <- shiftVec(new.el[o1, 1], vid)[order(o1)]
    o2 <- order(new.el[, 2])
    new.el[, 2] <- shiftVec(new.el[o2, 2], vid)[order(o2)]
  }
  if (!is.null(attr(el, "n"))) attr(new.el, "n") <- attr(el, "n") - length(vid)

  return(new.el)
}


#' @title Remove Edges That Include Specified Vertices
#'
#' @description Given a current two-column matrix of edges and a vector of
#'              vertex IDs, this function removes any rows of the edgelist in
#'              which the IDs are present.
#'
#' @param el A two-column matrix of current edges (edgelist).
#' @param vid A vector of vertex IDs whose edges are to be deleted from the
#'            edgelist.
#'
#' @return
#' Returns an updated edgelist object, with any edges including the specified
#' vertices removed.
#'
#' @export
#'
delete_edges <- function(el, vid) {
  new.el <- el
  if (length(vid) > 0) {
    el.rows.to.del <- which(el[, 1] %in% vid | el[, 2] %in% vid)
    if (length(el.rows.to.del) > 0) {
      new.el <- el[-el.rows.to.del, , drop = FALSE]
      if (!is.null(attr(el, "n"))) attr(new.el, "n") <- attr(el, "n")
    }
  }
  new.el
}


#' Depart Nodes from the netsim_dat Object
#'
#' @param dat the \code{netsim_dat} object
#' @param departures the vertex ids of nodes to depart
#'
#' @details If \code{tergmLite} is \code{FALSE}, the vertex ids
#' \code{departures} are deactivated (from the current timestep onward) in each
#' \code{networkDynamic} stored in \code{dat$nw}. If \code{tergmLite} is
#' \code{TRUE}, the vertex ids \code{departures} are deleted from \code{dat$el},
#' \code{dat$attr}, and \code{dat$net_attr}.
#'
#' @return the updated \code{netsim_dat} object with the nodes in
#' \code{departures} departed
#'
#' @export
#'
depart_nodes <- function(dat, departures) {
  if (length(departures) > 0) {
    if (get_control(dat, "tergmLite") == FALSE) {
      for (network in seq_len(dat$num.nw)) {
        dat$nw[[network]] <- deactivate.vertices(dat$nw[[network]],
                                                 onset = get_current_timestep(dat),
                                                 terminus = Inf,
                                                 v = departures,
                                                 deactivate.edges = TRUE)
      }
    } else {
      dat <- delete_attr(dat, departures)
      for (network in seq_len(dat$num.nw)) {
        dat$run$el[[network]] <- delete_vertices(dat$run$el[[network]], departures)
        dat$run$net_attr[[network]][["n"]] <- dat$run$net_attr[[network]][["n"]] - length(departures)

        if (get_network_control(dat, network, "tergmLite.track.duration") == TRUE) {
          dat$run$net_attr[[network]][["lasttoggle"]] <-
            delete_vertices(dat$run$net_attr[[network]][["lasttoggle"]], departures)
        }
      }
    }
  }
  return(dat)
}


#' Arrive New Nodes to the netsim_dat Object
#'
#' @param dat the \code{netsim_dat} object
#' @param nArrivals number of new nodes to arrive
#'
#' @details \code{nArrivals} new nodes are added to the network data stored on
#' the \code{netsim_dat} object. If \code{tergmLite} is \code{FALSE}, these
#' nodes are activated from the current timestep onward. Attributes for the new
#' nodes must be set separately.
#'
#' Note that this function only supports arriving new nodes; returning to an
#' active state nodes that were previously active in the network is not
#' supported.
#'
#' @return the updated \code{netsim_dat} object with \code{nArrivals} new nodes
#' added
#'
#'
#' @export
#'
arrive_nodes <- function(dat, nArrivals) {
  if (nArrivals > 0) {
    if (get_control(dat, "tergmLite") == FALSE) {
      for (network in seq_len(dat$num.nw)) {
        dat$nw[[network]] <- add.vertices.active(dat$nw[[network]],
                                                 nv = nArrivals,
                                                 onset = get_current_timestep(dat),
                                                 terminus = Inf)
      }
    } else {
      for (network in seq_len(dat$num.nw)) {
        dat$run$el[[network]] <- add_vertices(dat$run$el[[network]], nv = nArrivals)
        dat$run$net_attr[[network]][["n"]] <- dat$run$net_attr[[network]][["n"]] + nArrivals
      }
    }
  }
  return(dat)
}
