#' @title Get an Edgelist From the Specified Network
#'
#' @description This function outputs an edgelist from the specified network,
#'              selecting the method depending on the stored network type.

#' @inheritParams recovery.net
#' @param network Numerical index of the network from which the edgelist should
#'                be extracted. (May be > 1 for models with multiple overlapping
#'                networks.)
#'
#' @return
#' An edgelist in matrix form with two columns. Each column contains the
#' posit_ids (see \code{get_posit_ids}) of the nodes in each edge.
#'
#' @export
get_edgelist <- function(dat, network) {

  if (get_control(dat, "tergmLite")) {
    if (!network %in% seq_len(dat$num.nw)) {
      stop("There is no network '", network, "' to get an edgelist from")
    }
    el <- dat$run$el[[network]]
  } else {
    if (!network %in% seq_len(dat$num.nw)) {
      stop("There is no network '", network, "' to get an edgelist from")
    }
    at <- get_current_timestep(dat)
    el <- networkDynamic::get.dyads.active(dat$run$nw[[network]], at = at)
  }

  return(el)
}

## Discordant Edgelist ---------------------------------------------------------

#' @title Get Discordant Edgelist Based on Specified Status Variable
#'
#' @description This function returns a `data.frame` with a discordant
#'              edgelist, defined as the set of edges for which the status attribute
#'              of interest is discordant between the two partners.
#'
#' @inheritParams get_cumulative_edgelists_df
#' @param status.attr The name of the status attribute of interest.
#' @param head.status The value(s) of `status.attr` for which to look for the head of the edge.
#'        Can be a single value or a vector.
#' @param tail.status  The value(s) of `status.attr` for which to look for the tail of the edge.
#'        Can be a single value or a vector.
#'
#' @details
#' This is a generalized version of the `discord_edgelist` function.
#' It creates an edgelist of current partnerships in which the status attribute
#' of interest (as specified by the parameter `status.attr`) of one partner matches
#' the value (or one of the values) of the `head.status` parameter while the
#' corresponding status attribute of the other partner matches the value (or
#' one of the values) of the `tail.status` parameter.
#'
#' @return
#' A `data.frame` with the following columns:
#'  * `head`: Positional ID of the head node.
#'  * `tail`: Positional ID of the tail node.
#'  * `head_status`: Status of the head node.
#'  * `tail_status`: Status of the tail node.
#'  * `network`: The numerical index of the network on which the partnership is located.
#'
#' @seealso \code{\link{discord_edgelist}}
#'
#' @export
#' @keywords netMod internal
#'
get_discordant_edgelist <- function(dat, status.attr, head.status,
                                    tail.status, networks = NULL) {
  if (get_current_timestep(dat) == get_control(dat, "start") + 1 &&
        length(intersect(head.status, tail.status)) > 0) {
    warning("The head.status and tail.status arguments should be discordant.")
  }

  status <- get_attr(dat, status.attr)

  d_el <- get_edgelists_df(dat, networks)

  d_el_ordered <- dplyr::filter(
    d_el,
    status[d_el$head] %in% head.status &
      status[d_el$tail] %in% tail.status
  )

  d_el_rev <- dplyr::filter(
    d_el,
    status[d_el$tail] %in% head.status &
      status[d_el$head] %in% tail.status
  ) |>
    dplyr::select(head = "tail", tail = "head", "network")

  d_el <- dplyr::bind_rows(d_el_ordered, d_el_rev) |>
    dplyr::mutate(
      head_status = status[.data$head],
      tail_status = status[.data$tail]
    ) |>
    dplyr::select(
      "head", "tail", "head_status", "tail_status", "network"
    )

  return(d_el)

}

#' @title Get the Edgelist(s) from the Specified Network(s)
#'
#' @inheritParams get_cumulative_edgelists_df
#'
#' @return
#' A `data.frame` with the following columns:
#'  * `head`: Positional ID of the head node.
#'  * `tail`: Positional ID of the tail node.
#'  * `network`: The numerical index of the network on which the edge is located.
#'
#' @export
get_edgelists_df <- function(dat, networks = NULL) {

  networks <- if (is.null(networks)) seq_len(dat$num.nw) else networks
  el_list <- lapply(networks, get_edgelist, dat = dat)
  el_tibble <- lapply(el_list, as_tibble_edgelist)
  d_el <- dplyr::bind_rows(el_tibble)

  el_sizes <- vapply(el_list, nrow, numeric(1))
  d_el[["network"]] <- rep(networks, el_sizes)

  return(d_el)
}

#' @title Convert an Edgelist into a Tibble
#'
#' @param el An edgelist in matrix or data frame form.
#'
#' @return
#' The edgelist in tibble form with two columns named `head` and `tail`.
#'
#' @export
as_tibble_edgelist <- function(el) {

  if (nrow(el) > 0) {
    t_el <- tibble::tibble(head = el[, 1], tail = el[, 2])
  } else {
    t_el <- tibble::tibble(head = integer(0), tail = integer(0))
  }

  return(t_el)
}

## Cumulative Edgelists --------------------------------------------------------

#' @title Get a Cumulative Edgelist From a Specified Network
#'
#' @inheritParams recovery.net
#' @param network Numerical index of the network from which the cumulative
#'                edgelist should be extracted. (May be > 1 for models with
#'                multiple overlapping networks.)
#'
#' @return
#' A cumulative edgelist in \code{data.frame} form with 4 columns:
#' \itemize{
#'   \item \code{head}: the unique ID (see \code{get_unique_ids}) of the
#'         head node on the edge.
#'   \item \code{tail}: the unique ID (see \code{get_unique_ids}) of the
#'         tail node on the edge.
#'   \item \code{start}: the time step in which the edge started.
#'   \item \code{stop}: the time step in which the edge stopped; if ongoing,
#'         then \code{NA} is returned.
#' }
#'
#' @export
get_cumulative_edgelist <- function(dat, network) {
  if (!network %in% seq_len(dat$num.nw)) {
    stop("There is no network '", network,
         "' from which to get the cumulative edgelist.")
  }

  if (!get_control(dat, "cumulative.edgelist")) {
    stop("Failed to get the cumulative edgelist. It is likely not stored because the
         `cumulative.edgelist` control setting is set to `FALSE`.")
  }

  el_cuml <- dplyr::bind_rows(
    get_raw_elcuml(dat, network, active = FALSE),
    get_raw_elcuml(dat, network, active = TRUE)
  )

  return(el_cuml)
}

#' @title Update a Cumulative Edgelist of the Specified Network
#'
#' @inheritParams recovery.net
#' @param network Numerical index of the network for which the cumulative
#'                edgelist will be updated. (May be > 1 for models with
#'                multiple overlapping networks.)
#' @param truncate After how many time steps a partnership that is no longer
#'                 active should be removed from the output.
#'
#' @section Truncation:
#' To avoid storing a cumulative edgelist too long, the \code{truncate}
#' parameter defines a number of steps after which an edge that is no longer
#' active is truncated out of the cumulative edgelist.
#' When \code{truncate = Inf}, no edges are ever removed. When
#' \code{truncate = 0}, only the active edges are kept. You may want this
#' behavior to keep track of the active edges' start step.
#'
#' @inherit recovery.net return
#'
#' @export
update_cumulative_edgelist <- function(dat, network, truncate = 0) {
  if (!get_control(dat, "cumulative.edgelist")) {
    return(dat)
  }

  el <- get_edgelist(dat, network)
  el_cuml_cur <- get_raw_elcuml(dat, network, active = TRUE)
  el_cuml_hist <- get_raw_elcuml(dat, network, active = FALSE)

  at <- get_current_timestep(dat)

  # truncate el_cuml_hist
  if (truncate != Inf && truncate != 0) {
    rel.age <- at - el_cuml_hist$stop
    el_cuml_hist <- el_cuml_hist[rel.age <= truncate, ]
  }

  el <- tibble::tibble(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = TRUE
  )

  el_cuml_cur <- dplyr::full_join(el_cuml_cur, el, by = c("head", "tail"))

  # new edges
  new_edges_ids <- which(is.na(el_cuml_cur$start))
  if (length(new_edges_ids) > 0) {
    el_cuml_cur[new_edges_ids, ]$start <- at
  }

  # terminated edges
  terminated_edges_ids <- which(is.na(el_cuml_cur$current))
  el_cuml_cur$current <- NULL

  if (length(terminated_edges_ids) > 0) {
    el_cuml_term <- el_cuml_cur[terminated_edges_ids, ]

    # with truncate == 0, don't save any historic edges
    if (truncate != 0) {
      el_cuml_term$stop <- at - 1
      el_cuml_hist <- dplyr::bind_rows(el_cuml_hist, el_cuml_term)
    }

    el_cuml_cur <- el_cuml_cur[-terminated_edges_ids, ]
  }

  dat <- set_raw_elcuml(dat, network, el_cuml_cur, active = TRUE)
  dat <- set_raw_elcuml(dat, network, el_cuml_hist, active = FALSE)

  return(dat)
}

#' @title Get the Cumulative Edgelists of a Model
#'
#' @inheritParams recovery.net
#' @param networks Numerical indexes of the networks to extract the partnerships
#'                 from. (May be > 1 for models with multiple overlapping
#'                 networks.) If \code{NULL}, extract from all networks.
#'
#' @return
#' A \code{data.frame} with 5 columns:
#' \itemize{
#'   \item \code{index}: the unique ID (see \code{get_unique_ids}) of the
#'         indexes.
#'   \item \code{partner}: the unique ID (see \code{get_unique_ids}) of the
#'         partners/contacts.
#'   \item \code{start}: the time step in which the edge started.
#'   \item \code{stop}: the time step in which the edge stopped; if ongoing,
#'         then \code{NA} is returned.
#'   \item \code{network}: the numerical index for the network on which the
#'         partnership/contact is located.
#'  }
#'
#' @export
get_cumulative_edgelists_df <- function(dat, networks = NULL) {
  networks <- if (is.null(networks)) seq_len(dat$num.nw) else networks

  el_cuml_list <- lapply(networks, get_cumulative_edgelist, dat = dat)
  el_cuml_df <- dplyr::bind_rows(el_cuml_list)

  el_sizes <- vapply(el_cuml_list, nrow, numeric(1))
  el_cuml_df[["network"]] <- rep(networks, el_sizes)

  return(el_cuml_df)
}

#' @title Return the Historical Contacts (Partners) of a Set of Index Nodes
#'
#' @description
#' From a full cumulative edgelist that contains the history of contacts (both persistent and
#' one-time), this function returns a data frame containing details of the index (head) and partner
#' (tail) nodes, along with start and stop time steps for the partnership and the network location.
#'
#' @param index_posit_ids The positional IDs of the indexes of interest.
#' @param networks Numerical indexes of the networks to extract the partnerships from. (May be > 1
#'        for models with multi-layer networks.) If `NULL`, extract from all networks.
#' @param only.active.nodes If `TRUE`, then inactive (e.g., deceased) partners will be removed from
#'        the output.
#' @inheritParams update_cumulative_edgelist
#'
#' @return
#' A `data.frame` with 5 columns:
#'   * `index`: the unique IDs of the indexes.
#'   * `partner`: the unique IDs of the partners/contacts.
#'   * `start`: the time step at which the edge started.
#'   * `stop`: the time step in which the edge stopped; if ongoing, then `NA` is returned.
#'   * `network`: the numerical index for the network on which the partnership/contact is located.
#'
#' @details
#' Note that `get_partners` takes as input the positional IDs of the indexes of interest but returns
#' the unique IDs. That is by design, because while `get_partners` would be expected to be called
#' for active nodes, some partners (contacts) of nodes may be inactive in the network history.
#' Therefore, both index and partner IDs are returned as unique IDs for consistency. To convert
#' between a positional to a unique ID, you may use [`get_posit_ids`]; to convert between a
#' unique ID to a positional ID, you may use [`get_unique_ids`].
#'
#' @export
#'
get_partners <- function(dat, index_posit_ids, networks = NULL,
                         truncate = Inf, only.active.nodes = FALSE) {

  el_cuml_df <- get_cumulative_edgelists_df(dat, networks)
  index_unique_ids <- get_unique_ids(dat, index_posit_ids)

  partner_head_df <- el_cuml_df[el_cuml_df[["head"]] %in% index_unique_ids, ]
  partner_tail_df <- el_cuml_df[
    el_cuml_df[["tail"]] %in% index_unique_ids,
    c(2, 1, 3:5) # switch the head and tail columns
  ]

  colnames(partner_head_df) <- c("index", "partner", "start", "stop", "network")
  colnames(partner_tail_df) <- colnames(partner_head_df)

  partner_df <- dplyr::bind_rows(partner_head_df, partner_tail_df)

  if (only.active.nodes) {
    active_partners <- is_active_unique_ids(dat, partner_df[["partner"]])
    partner_df <- partner_df[active_partners, ]
  }

  if (truncate != Inf) {
    at <- get_current_timestep(dat)
    rel.age <- at - partner_df[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    partner_df <- partner_df[rel.age <= truncate, ]
  }

  return(partner_df)
}

#' @title Return the Cumulative Degree of a Set of Index Nodes
#'
#' @inheritParams get_partners
#'
#' @return
#' A \code{data.frame} with 2 columns:
#' \itemize{
#'   \item \code{index_pid}: the positional ID (see \code{get_posit_ids}) of the
#'         indexes.
#'   \item \code{degree}: the cumulative degree of the index.
#'  }
#'
#' @section Cumulative Degree:
#' The cumulative degree of a node is the number of edges connected to this
#' node at during the time window. The time window is by default all the steps
#' stored in the `cumulative_edgelist` or set by the `truncate` parameter.
#'
#' @export
get_cumulative_degree <- function(dat, index_posit_ids, networks = NULL,
                                  truncate = Inf, only.active.nodes = FALSE) {
  get_partners(
    dat, index_posit_ids, networks,
    truncate, only.active.nodes
  ) |>
    dplyr::summarize(degree = dplyr::n(), .by = "index") |>
    dplyr::mutate(index = get_posit_ids(dat, .data$index)) |>
    dplyr::select(index_pid = "index", "degree")
}

# Helper functions to get and set the cumulative edgelists
#
# Cumulative edgelists are split in to:
#   - el_cuml_cur: for the edges not yet disolved
#   - el_cuml_hist: for the edges where the start and stop time are known
#   (dissolved edges)
get_raw_elcuml <- function(dat, network, active) {
  loc <- if (active) "el_cuml_cur" else "el_cuml_hist"

  if (length(dat$run[[loc]]) >= network) {
    el_cuml <- dat$run[[loc]][[network]]
  } else {
    el_cuml <- NULL
  }

  if (is.null(el_cuml)) {
    el_cuml <- empty_el_cuml()
  }

  return(el_cuml)
}

set_raw_elcuml <- function(dat, network, el_cuml, active) {
  loc <- if (active) "el_cuml_cur" else "el_cuml_hist"
  dat$run[[loc]][[network]] <- el_cuml
  return(dat)
}

# template for the cumulative edgelists
empty_el_cuml <- function() {
  tibble::tibble(
    head  = numeric(0),
    tail  = numeric(0),
    start = numeric(0),
    stop  = numeric(0)
  )
}
