#' Convert an object to a `cumulative_edgelist`
#'
#' @param x An object to be converted to a cumulative edgelist
#'
#' @return A `cumulative_edgelist` object, a `data.frame` with at least the
#' following columns: `head`, `tail`, `start`, `stop`.
#'
#' @details
#' The edges are active from time `start` to time `stop` included. If stop is
#' `NA`, the edge was not disolved in the simulation that generated the list.
#'
#' @export
as_cumulative_edgelist <- function(x) {
  UseMethod("as_cumulative_edgelist")
}

#' @export
as_cumulative_edgelist.networkDynamic <- function(x) {
  d <- as.data.frame(x)
  d <- d[c("head", "tail", "onset", "terminus")]
  names(d) <- c("head", "tail", "start", "stop")
  d$stop <- d$stop - 1
  class(d) <- c("cumulative_edgelist", class(d))
  return(d)
}

#' Deduplicate a cumulative edgelist by combining overlapping edges
#'
#' @param el A cumulative edgelist with potentially overlapping edges
#'
#' @return A cumulative edgelist with no overlapping edges
dedup_cumulative_edgelist <- function(el) {
  el_n <- el |>
    dplyr::group_by(head, tail) |>
    dplyr::mutate(n = dplyr::n()) |>
    dplyr::ungroup()

  e_unique <- el_n |>
    dplyr::filter(n == 1) |>
    dplyr::select(-n)

  e_dup <- el_n |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n) |>
    dplyr::arrange(head, tail, start, stop)

  e_dedup <- e_dup |>
    dplyr::group_by(head, tail) |>
    dplyr::mutate(
      lstart = dplyr::lag(start),
      lstop = dplyr::lag(stop),
      overlap = !is.na(lstop) & !is.na(lstart) & start <= lstop,
      stop = ifelse(overlap, max(stop, lstop, na.rm = TRUE), stop),
      start = ifelse(overlap, min(start, lstart, na.rm = TRUE), start)
    ) |>
    dplyr::select(-c(lstart, lstop, overlap)) |>
    dplyr::ungroup() |>
    unique()

  dplyr::bind_rows(e_unique, e_dedup)
}


#' @title Get the Forward or Backward Reachable Nodes for a Set of Nodes
#'
#' @description
#' These functions return the Forward or Backward Reachable Nodes of a set of
#' nodes in a network over a time. It is much faster than iterating
#' \code{tsna::tPath} over all nodes but does not give the distance between each
#' node. Only the reachable status over a time period.
#'
#' @param el_cuml a cumulative edgelist object. That is a data.frame with at
#'   least columns: head, tail, start and stop. Start and stop are inclusive.
#' @param from_step the beginning of the time period.
#' @param to_step the end of the time period.
#' @param nodes the subset of nodes to calculate the FRP for. (default = NULL,
#'        all nodes)
#'
#' @return
#' A named list of reachable nodes for each of the `nodes`.
#'
#' @section Time and Memory Use:
#' These functions may be used to efficiently calculate multiple sets of
#' reachable nodes. As cumulative edgelists are way smaller than full
#' `networkDynamic` objects, theses functions are suited for large and dense
#' networks.
#' Also, as long as the size of the `nodes` set is greater than 5, theses
#' functions are faster than iterating over `tsna::tPath`.
#'
#' @section Displaying Progress:
#' These functions are using the
#' \href{https://progressr.futureverse.org/articles/progressr-intro.html}{progressr package}
#' to display its progression. Use
#' \code{progressr::with_progress({ fwd_reach <- get_forward_reachable(el, from = 1, to = 260) })}
#' to display the progress bar. Or see the
#' \href{https://progressr.futureverse.org/articles/progressr-intro.html}{progressr package}
#' for more information and customization.
#'
#' @section Number of Nodes:
#' As cumulative edgelists do not contain the total number of nodes in the
#' network, if `nodes = NULL`, the total number of node on the network is
#' assumed to be the highest ID recorded in an edge.
#' We can therefore arrive to a situation where there is less elements in the
#' output than node in the network if the last N nodes (by ID) are never
#' connected. And therefore are not recorded in the cumulative edgelist.
#' However, it means that these nodes have no reachable nodes appart themselves.
#'
#' @examples
#' \dontrun{
#'
#' # load a network dynamic object
#' nd <- readRDS("nd_obj.Rds")
#' # convert it to a cumulative edgelist
#' el_cuml <- as_cumulative_edgelist(nd)
#'
#' # sample 100 node indexes
#' nnodes <- max(el_cuml$head, el_cuml$tail)
#' nodes <- sample(nnodes, 100)
#'
#' # `get_forward_reachable` uses steps [from_step, to_step] inclusive
#' el_fwd <- get_forward_reachable(el_cuml, 1, 52, nodes)
#'
#' # check if the results are consistent with `tsna::tPath`
#' for (i in seq_along(el_tp)) {
#'   t_fwd <- tsna::tPath(
#'     nd, v = nodes[[i]],
#'     start = 1, end = 52 + 1, # tPath works from [start, end) right exclusive
#'     direction = "fwd"
#'   )
#'
#'   t_fwd_set <- which(t_fwd$tdist < Inf)
#'   if(!setequal(el_fwd[[i]], t_fwd_set))
#'     stop("Missmatch on node: ", names(nodes)[[i]])
#' }
#'
#' # Backward:
#' el_bkwd <- get_backward_reachable(el_cuml, 1, 52, nodes = 1)
#' t_bkwd <- tsna::tPath(
#'   nd, v = nodes[[i]],
#'   start = 1, end = 52 + 1,
#'   direction = "bkwd", type = "latest.depart"
#' )
#' t_bkwd_set <- which(t_bkwd$tdist < Inf)
#' setequal(el_bkwd[[1]], t_bkwd_set)
#'
#' }
#'
#' @name reachable-nodes
NULL

#' @rdname reachable-nodes
#' @export
get_forward_reachable <- function(el_cuml, from_step, to_step, nodes = NULL) {
  n_nodes <- max(c(el_cuml$head, el_cuml$tail))
  if (is.null(nodes))
    nodes <- seq_len(n_nodes)

  # the initial FRP contains only the vertex itself
  reach_cur <- as.list(nodes)
  names(reach_cur) <- paste0("node_", nodes)

  # Prepare the cumulative edgelist:
  # - set the `stop` time to `Inf` instead of NA
  # - remove the edges that don't exist in the analysis period
  # - set the oldest edges `start` to `to_step`
  #     (QoL to calcuate the `change_times`)

  # nolint start
  el_cuml <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop)) |> # current edges never ends
    dplyr::filter(start <= to_step, stop >= from_step) |> # remove edges before and after analysis period
    dplyr::mutate(start = ifelse(start < from_step, from_step, start)) |> # set older edges to start at beginning of analysis
    dplyr::select(start, stop, head, tail)
  # nolint end

  # Only consider steps where new edges occur
  change_times <- sort(unique(el_cuml$start))
  p <- progressr::progressor(length(change_times))

  for (cur_step in change_times) {
    p()

    # nolint start
    #
    # IN EL_CUML: duration is [start, stop] (inclusive)
    # all active edges a time `cur_step` are needed (not only start).
    # This is because getting connected to a node X means also indirectly
    # connecting to its connections, even the ones that started earlier.
    el_cur <- dplyr::filter(el_cuml, start <= cur_step, stop >= cur_step)
    # nolint end

    # get subnet works with adjacency list
    adj_list <- get_adj_list(el_cur, n_nodes)

    # at time T, the REACH(T) of a node is the subnet connected to the REACH(T-1)
    reach_cur <- lapply(reach_cur, get_connected_subnet, adj_list = adj_list)
  }
  return(reach_cur)
}

#' Returns all the node connected directly or indirectly to a set of nodes
#'
#' @param adj_list The network represented as an adjacency list
#' @param nodes A set of nodes
#'
#' @return A vector of nodes indexes that are connected together with the ones
#'         provided in the `nodes` argument
get_connected_subnet <- function(adj_list, nodes) {
  new <- nodes
  subnet <- nodes
  n_nodes <- length(adj_list)
  while (length(new) > 0 && length(subnet) < n_nodes) {
    new <- unlist(adj_list[new])
    new <- setdiff(new, subnet)
    subnet <- c(subnet, new)
  }
  subnet
}

#' Returns an adjacency list from an edge list
#'
#' @param el An edge list as a data.frame with columns `head` and `tail`
#' @param n_nodes The size number of node in the network
#'
#' @return An adjacency list for the network
#'
#' @details
#' The adjacency list is a `list` of length `n_nodes`. The entry for each node
#' is a integer vector containing the index of all the nodes connected to it.
#' This layout makes it directly subsetable in O(1) at the expanse of memory
#' usage.
#' To get all connections to the nodes 10 and 15 : `unlist(adj_list[c(10, 15)]`
get_adj_list <- function(el, n_nodes) {
  head <- el$head
  tail <- el$tail

  adj_list <- vector(mode = "list", length = n_nodes)
  for (i in seq_len(nrow(el))) {
    e_head <- head[i]
    e_tail <- tail[i]
    adj_list[[e_head]] <- c(adj_list[[e_head]], e_tail)
    adj_list[[e_tail]] <- c(adj_list[[e_tail]], e_head)
  }
  adj_list
}

#' @rdname reachable-nodes
#' @export
get_backward_reachable <- function(el_cuml, from_step, to_step, nodes = NULL) {
  # simply invert the time before calling get_forward_reachable`
  el_cuml$stop <- ifelse(is.na(el_cuml$stop), Inf, el_cuml$stop)
  tmp <- el_cuml$start
  el_cuml$start <- - el_cuml$stop
  el_cuml$stop <- - tmp

  get_forward_reachable(el_cuml, -to_step, -from_step, nodes)
}
