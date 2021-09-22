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
    if (!is.null(dat[["temp"]][["nw_list"]])) {
      if (!get_control(dat, "resimulate.network")) {
        el <- network::as.edgelist(dat[["temp"]][["nw_list"]][[at]])
      } else {
        el <- network::as.edgelist(dat[["nw"]][[network]])
      }
    } else {
      el <- networkDynamic::get.dyads.active(dat[["nw"]][[network]], at = at)
    }
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
  if (length(dat[["el_cuml"]]) < network) {
    el_cuml <- NULL
  } else {
    el_cuml <- dat[["el_cuml"]][[network]]
  }

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
#' @param truncate after how many steps an edge that is no longer active should
#'                 be removed from the cumulative edgelist (default = Inf)
#'
#' @section Truncation:
#' To avoid storing a cumulative edgelist to long, the `truncate`
#' parameter defines a number of steps after which an edge that is no longer
#' active is truncated out of the cummulative edgelist. When `truncate == Inf`
#' (default), no edge is ever removed. When `truncate == 0`, only the active
#' edges are kept. You may want this behavior to keep track of the active edges
#' start step.
#'
#' @return an updated Master list object of network models
#'
#' @export
update_cumulative_edgelist <- function(dat, at, network, truncate = Inf) {
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

  if (truncate != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate, ]
  }

  dat[["el_cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]

  return(dat)
}

#' @title Get the Cumulative Edgelists of a Model in \code{data.frame} Format
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

#' @title Get the Previous Partners of a Set of Indexes
#'
#' @param dat a Master list object of network models
#' @param at the current timestep
#' @param index_posit_ids the positional IDs of the indexes of interest
#' @param networks the indexes of the networks to extract the partnerships from.
#'                 If NULL (default) extract from all networks.
#' @param max.age after how many steps a partnershipl that is no longer active
#'                should be removed from the results (default = Inf)
#' @param only.active should the inactive partners be removed from the results
#'                    (default = FALSE)
#'
#' @return a \code{data.frame} with 5 columns:
#'          index and partner, the unique_ids (see \code{get_unique_ids}) of the
#'          indexes and partners respectively. Start and stop, the timestep
#'          where the edges started and stopped. network, the network on which
#'          the partnership is.
#'
#' @export
get_partners <- function(dat, at, index_posit_ids, networks = NULL,
                         max.age = Inf, only.active = FALSE) {
  el_cuml_df <- get_cumulative_edgelists_df(dat, networks)
  index_unique_ids <- get_unique_ids(dat, index_posit_ids)

  partner_head_df <- el_cuml_df[el_cuml_df[["head"]] %in% index_unique_ids, ]
  partner_tail_df <- el_cuml_df[
    el_cuml_df[["tail"]] %in% index_unique_ids,
    c(2, 1, 3:5) # switch the head and tail columns
  ]

  colnames(partner_head_df) <- c("index", "partner", "start", "stop", "network")
  colnames(partner_tail_df) <- colnames(partner_head_df)

  partner_df <- rbind(partner_head_df, partner_tail_df)

  if (only.active) {
    active_partners <- is_active_unique_ids(dat, partner_df[["partner"]])
    partner_df <- partner_df[active_partners, ]
  }

  if (max.age != Inf) {
    rel.age <- at - partner_df[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    partner_df <- partner_df[rel.age <= max.age, ]
  }

  return(partner_df)
}
