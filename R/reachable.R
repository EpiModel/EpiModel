#' @title Get the Forward or Backward Reachable Nodes for a Set of Nodes
#'
#' @description
#' These functions return the Forward or Backward Reachable Nodes of a set of
#' nodes in a network over a time. Warning, these functions ignore nodes without
#' edges in the period of interest. See the `Number of Nodes` section for
#' details It is much faster than iterating
#' \code{tsna::tPath}. The distance between to each node can be back calculated
#' using the length of the reachable set at each time step and the fact that the
#' reachable sets are ordered by the time to arrival.
#'
#' @param el_cuml a cumulative edgelist object. That is a data.frame with at
#'   least columns: head, tail, start and stop. Start and stop are inclusive.
#' @param from_step the beginning of the time period.
#' @param to_step the end of the time period.
#' @param nodes the subset of nodes to calculate the FRP for. (default = NULL,
#'        all nodes)
#' @param dense_optim pre-process the adjacency list to speed up the
#'        computations on dense networks. "auto" (default), enable the
#'        optimisation when `n_edges > n_nodes`. "yes" always enables and "no"
#'        always disables. The overhead of the optimization is not worth it on
#'        sparse networks.
#'
#' @return
#' A named list containing:
#'  `reached`: the set of reachable nodes for each of the `nodes`.
#'  `lengths`: A matrix of `length(nodes)` rows and one column per timestep + 1
#'      with the length of the reachable set at each step from `from_step - 1`
#'      to `to_step`. The first column is always one as the set of reachables
#'      at the beginning is just the node itself.
#'
#' @section Number of Nodes:
#' To speed up the calculations and lower the memory usage, these functions
#' only take into account nodes with edges in the cumulative edgelist over the
#' period of interest. The nodes are identified in the `reached` and `lengths`
#' sublists by names (e.g. `node_1093`). Nodes without any edges are therefore
#' not calculated as the only node they reach is themselve (length of 1).
#' Take this fact into account when exploring the distribution of Forward
#' Reachable Paths for example. As the nodes with FRP == 1 are not in the
#' output.
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
#' el_fwd <- get_forward_reachable(el_cuml, 1, 52, nodes)[["reached"]]
#'
#' # check if the results are consistent with `tsna::tPath`
#' nodes <- strsplit(names(el_fwd), "_")
#' for (i in seq_along(el_fwd)) {
#'   node <- as.integer(nodes[[i]][2])
#'   t_fwd <- tsna::tPath(
#'     nd, v = node,
#'     start = 1, end = 52 + 1, # tPath works from [start, end) right exclusive
#'     direction = "fwd"
#'   )
#'
#'   t_fwd_set <- which(t_fwd$tdist < Inf)
#'   if(!setequal(el_fwd[[i]], t_fwd_set))
#'     stop("Missmatch on node: ", node)
#' }
#'
#' # Backward:
#' el_bkwd <- get_backward_reachable(el_cuml, 1, 52, nodes = 1)[["reached"]]
#' nodes <- strsplit(names(el_bkwd), "_")
#' t_bkwd <- tsna::tPath(
#'   nd, v = nodes[i][2],
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
get_forward_reachable <- function(el_cuml, from_step, to_step, nodes = NULL,
                                  dense_optim = "auto") {
  # Prepare the cumulative edgelist:
  # - set the `stop` time to `Inf` instead of NA
  # - remove the edges that don't exist in the analysis period
  # - set the oldest edges `start` to `to_step`
  #     (QoL to calcuate the `change_times`)
  # nolint start
  el_cuml <- el_cuml |>
    dplyr::mutate(stop = ifelse(is.na(.data$stop), Inf, .data$stop)) |> # current edges never ends
    dplyr::filter(.data$start <= to_step, .data$stop >= from_step) |> # remove edges before and after analysis period
    dplyr::mutate(start = ifelse(.data$start < from_step, from_step, .data$start)) |> # set older edges to start at beginning of analysis
    dplyr::select("start", "stop", "head", "tail")
  # nolint end

  # Make the node indexes continuous
  # (original indexes are reset at the end)
  orig_indexes <- sort(unique(c(el_cuml$head, el_cuml$tail)))
  el_cuml$head <- match(el_cuml$head, orig_indexes)
  el_cuml$tail <- match(el_cuml$tail, orig_indexes)

  if (is.null(nodes)) {
    nodes <- seq_along(orig_indexes)
  } else {
    nodes <- match(nodes, orig_indexes)
    nodes <- nodes[!is.na(nodes)]
  }

  if (length(nodes) == 0) {
    warning("No nodes have edges in the network over this period")
    return(list(reached = list(), lengths = matrix(0)))
  }

  # Only consider steps where new edges occur
  change_times <- sort(unique(el_cuml$start))
  n_steps <- to_step - from_step + 1

  # the initial FRP contains only the vertex itself
  reach_cur <- as.list(nodes)
  reach_lengths <- matrix(0, ncol = n_steps + 1, nrow = length(nodes))
  reach_lengths[, 1] <- 1

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

    adj_list <- get_adj_list(el_cur, length(orig_indexes))

    # precompute the network components. Useful for dense networks
    if (dense_optim == "yes") {
      adj_list <- get_subnet_adj_list(adj_list)
    } else if (dense_optim == "auto" && nrow(el_cur) > length(orig_indexes)) {
      adj_list <- get_subnet_adj_list(adj_list)
    }

    # at time T, the REACH(T) of a node is the subnet connected to the REACH(T-1)
    new_reached <- lapply(reach_cur, get_connected_nodes, adj_list = adj_list)
    reach_lengths[, cur_step - from_step + 2] <- vapply(new_reached, length, 0)
    reach_cur <- Map(c, reach_cur, new_reached)
  }

  reach_lengths <- t(apply(reach_lengths, 1, cumsum))
  reach_cur <- lapply(reach_cur, \(x) orig_indexes[x])

  names(reach_cur) <- paste0("node_", orig_indexes[nodes])
  rownames(reach_lengths) <- names(reach_cur)
  colnames(reach_lengths) <- paste0("step_", (from_step - 1):to_step)

  return(
    list(
      reached = reach_cur,
      lengths = reach_lengths
    )
  )
}

#' @rdname reachable-nodes
#' @export
get_backward_reachable <- function(el_cuml, from_step, to_step, nodes = NULL,
                                   dense_optim = "auto") {
  # simply invert the time before calling get_forward_reachable`
  el_cuml$stop <- ifelse(is.na(el_cuml$stop), Inf, el_cuml$stop)
  tmp <- el_cuml$start
  el_cuml$start <- - el_cuml$stop
  el_cuml$stop <- - tmp

  get_forward_reachable(el_cuml, -to_step, -from_step, nodes, dense_optim)
}

#' Returns all the node connected directly or indirectly to a set of nodes
#'
#' @param adj_list The network represented as an adjacency list
#' @param nodes A set of nodes
#'
#' @return A vector of nodes indexes that are connected together with the ones
#'         provided in the `nodes` argument. The `nodes` themselves are not
#'         listed in this output
get_connected_nodes <- function(adj_list, nodes) {
  n_nodes <- length(adj_list)
  new_connections <- numeric(0)
  subnet <- nodes
  while (length(nodes) > 0 && length(subnet) < n_nodes) {
    nodes <- unlist(adj_list[nodes])
    nodes <- setdiff(nodes, subnet)
    subnet <- c(subnet, nodes)
    new_connections <- c(new_connections, nodes)
  }
  new_connections
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

#' Return an adjacency list of subnets
#'
#' @param adj_list A normal adjacency list
#'
#' @return An adjacency list where only the first node of a subnet contains the
#' subnet and all other contain only the first node
#'
#' @details
#' A graph with 4 components: 1, 2, 3, 4, and 5 and 6, 7, 8  would yield a list
#' like so:
#' 1: 2, 3, 4
#' 2: 1
#' 3: 1
#' 4: 1
#' 5: numeric(0)
#' 6: 7, 8
#' 7: 6,
#' 8: 6
#'
#' This format speeds up the construction of reachable sets on dense networks
get_subnet_adj_list <- function(adj_list) {
  for (i in seq_along(adj_list)) {
    # if first elt connected to before himself, already in a subnet
    if (length(adj_list[[i]]) == 0 || adj_list[[i]][1] < i)
      next
    # get current subnet of i
    subnet <- get_connected_nodes(adj_list, i)
    adj_list[[i]] <- subnet
    # all members of subnet point to i
    adj_list[subnet] <- i
  }
  adj_list
}

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
    dplyr::filter(.data$n == 1) |>
    dplyr::select(-"n")

  e_dup <- el_n |>
    dplyr::filter(.data$n > 1) |>
    dplyr::select(-"n") |>
    dplyr::arrange(.data$head, .data$tail, .data$start, .data$stop)

  e_dedup <- e_dup |>
    dplyr::group_by(.data$head, .data$tail) |>
    dplyr::mutate(
      lstart = dplyr::lag(.data$start),
      lstop = dplyr::lag(.data$stop),
      overlap = !is.na(.data$lstop) & !is.na(.data$lstart) &
                .data$start <= .data$lstop,
      stop = ifelse(.data$overlap, max(.data$stop, .data$lstop, na.rm = TRUE), .data$stop),
      start = ifelse(.data$overlap, min(.data$start, .data$lstart, na.rm = TRUE), .data$start)
    ) |>
    dplyr::select(-c("lstart", "lstop", "overlap")) |>
    dplyr::ungroup() |>
    unique()

  dplyr::bind_rows(e_unique, e_dedup)
}
