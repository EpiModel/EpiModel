#' @title Get an Edgelist From the Specified Network
#'
#' @description This function outputs an edgelist from the specified network.
#'              it chooses the right method depending on the type of network and
#'              the output is standardly formatted.
#'
#' @param dat a Master list object of network models
#' @param at the current timestep
#' @param network the index of the network to extract the edgelist from
#'
#' @return an edgelist in matrix form with the two columns. Each column
#'          contains the posit_ids (see \code{get_posit_ids})of the nodes in
#'          each edge
#'
#' @export
get_edgelist <- function(dat, at, network) {
  if (get_control(dat, "tergmLite")) {
    el <- dat[["el"]][[network]]
  } else {
    el <- networkDynamic::get.dyads.active(dat[["nw"]][[network]], at = at)
  }

  return(el)
}

#' @title Get a Cumulative Edgelist From the Specified Network
#'
#' @param dat a Master list object of network models
#' @param network the index of the network to extract the cumulative edgelist
#'                from
#'
#' @return a cumulative edgelist in \code{data.frame} form with 4 columns:
#'          head and tail, the unique_ids (see \code{get_unique_ids}) of the
#'          nodes on the edge. Start and stop, the timestep where the edges
#'          started and stopped.
#'
#' @export
get_cumulative_edgelist <- function(dat, network) {
  el_cuml <- dat[["el_cuml"]][[network]]

  if(is.null(el_cuml)) {
    el_cuml <- data.frame(
      head  = numeric(0),
      tail  = numeric(0),
      start = numeric(0),
      stop  = numeric(0)
    )
  }

  return(el_cuml)
}

#' @title Update a Cumulative Edgelist of the Specified Network
#'
#' @param dat a Master list object of network models
#' @param at the current timestep
#' @param network the index of the network whose cumulative edgelist will be
#'                updated
#'
#' @section Truncation:
#' To avoid storing a cumulative edgelist to long, the "truncate.el_cuml"
#' control value defines a number of steps after which an edge that is no longer
#' active is truncated out of the cummulative edgelist.
#'
#' @return an updated Master list object of network models
#'
#' @export
update_cumulative_edgelist <- function(dat, at, network) {
  truncate.el_cuml <- get_control(
    dat, "truncate.el_cuml", override.null.error = TRUE)
  truncate.el_cuml <- if (is.null(truncate.el_cuml)) Inf else truncate.el_cuml

  el <- get_edgelist(dat, at, network)
  el_cuml <- get_cumulative_edgelist(dat, network)

  el <- data.frame(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = TRUE
  )

  el_cuml <- merge(el_cuml, el, by = c("head", "tail"), all = TRUE)

  new_edges <- is.na(el_cuml[["start"]])
  if (any(new_edges)) {
    el_cuml[new_edges, ][["start"]] <- at
  }

  terminated_edges <- is.na(el_cuml[["current"]]) & is.na(el_cuml[["stop"]])
  if (any(terminated_edges)) {
    el_cuml[terminated_edges, ][["stop"]] <- at - 1
  }

  if (truncate.el_cuml != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate.el_cuml, ]
  }

  dat[["el_cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]

  return(dat)
}

#' @title Get the Cumulative Edgelists of a Model in \code{data.frame} format
#'
#' @param dat a Master list object of network models
#' @param networks the indexes of the networks to extract the cumulative
#'                 edgelists from. If NULL (default) extract all the cumulative
#'                 edgelists.
#'
#' @return a \code{data.frame} with 5 columns:
#'          head and tail, the unique_ids (see \code{get_unique_ids}) of the
#'          nodes on the edge. Start and stop, the timestep where the edges
#'          started and stopped. network, the network on which the edge is.
#'
#' @export
get_cumulative_edgelists_df <- function(dat, networks = NULL) {
  networks <- if (is.null(networks)) seq_along(dat[["nwparam"]]) else networks

  el_cuml_df <- Reduce(
    function(a, x) {
      el_cuml_tmp <- get_cumulative_edgelist(dat, network = x)
      el_cuml_tmp[["network"]] <- x
      rbind(a, el_cuml_tmp)
    },
    networks,
    init = data.frame()
  )

  return(el_cuml_df)
}
