
#' @title Initialize Network Object
#'
#' @description Initialize an undirected network object for use in EpiModel
#'              workflows.
#'
#' @param n Network size.
#'
#' @details
#' This function is used in `EpiModel` workflows to initialize an empty
#' network object.  The network attributes `directed`, `bipartite`,
#' `hyper`, `loops`, and `multiple` are set to `FALSE`.
#'
#' @return
#' Returns an object of class `network`.
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
#' @description Sets a vertex attribute on an object of class `network`.
#'              This function simplifies the related function in the
#'              `network` package.
#'
#' @param x An object of class network.
#' @param attrname The name of the attribute to set.
#' @param value A vector of values of the attribute to be set.
#' @param v IDs for the vertices whose attributes are to be altered.
#'
#' @details
#' This function is used in `EpiModel` workflows to set vertex attributes
#' on an initialized empty network object (see [network_initialize()]).
#'
#' @return
#' Returns an object of class `network`.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' nw <- set_vertex_attribute(nw, "age", runif(100, 15, 65))
#' nw
#'
set_vertex_attribute <- function(x, attrname, value, v = NULL) {
  if (is.null(v)) {
    v <- seq_len(network.size(x))
  }
  g <- set.vertex.attribute(x, attrname, value, v)
  return(g)
}


#' @title Get Vertex Attribute on Network Object
#'
#' @description Gets a vertex attribute from an object of class `network`.
#'              This function simplifies the related function in the
#'              `network` package.
#'
#' @param x An object of class network.
#' @param attrname The name of the attribute to get.
#'
#' @details
#' This function is used in `EpiModel` workflows to query vertex
#' attributes on an initialized empty network object (see
#' [network_initialize()]).
#'
#' @return
#' Returns a vector of vertex attribute values for the attribute specified
#' by `attrname`.
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

#' @title Get Network Attributes from a Network Object
#'
#' @description Gets all network attributes except `"mnext"` from its
#'              network argument.
#'
#' @param x An object of class `network` or `networkLite`.
#'
#' @details
#' This function is used in `EpiModel` workflows to copy relevant network
#' attributes from the network object to the `netsim_dat` object when
#' initializing `netsim` runs.
#'
#' @return
#' Returns the named list of network attributes.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' get_network_attributes(nw)
#'
get_network_attributes <- function(x) {
  out <- list()
  for (name in setdiff(list.network.attributes(x), c("mnext"))) {
    new <- list(get.network.attribute(x, name))
    names(new) <- name
    out <- c(out, new)
  }
  out
}

make_networkDynamic <- function(sim, sim_num = 1, network = NULL) {
  if (!sim$control$save.cumulative.edgelist &&
        !sim$control$cumulative.edgelist) {
    stop("cumulative edgelist must be used and saved")
  }
  if (! "active" %in% sim$control$tracked.attributes) {
    stop("The `active` attribute must be tracked")
  }

  attr_hist <- get_attr_history(sim)
  n_nodes <- max(dplyr::filter(attr_hist$active, sim == sim_num)$uids)
  n_steps <- sim$control$nsteps

  nw <- network::network.initialize(n_nodes, directed = FALSE)
  el <- sim$cumulative.edgelist[[paste0("sim", sim_num)]] |>
    dplyr::mutate(stop = ifelse(is.na(stop), n_steps + 1, stop))

  if (!is.null(network))
    el <- dplyr::filter(el, network == .env$network)

  networkDynamic::add.edges.active(
    nw,
    head = el$head,
    tail = el$tail,
    onset = el$start,
    terminus = el$stop + 1L,
    names.eval = rep(list("network"), nrow(el)),
    vals.eval = lapply(el$network, \(x) list(network = x))
  )

  # Manage when nodes are active
  d_active <- dplyr::filter(attr_hist$active, sim == sim_num) |>
    dplyr::select(time, uids, values) |>
    dplyr::mutate(values = ifelse(values == 1, "onset", "terminus")) |>
    dplyr::group_by(uids, values) |>
    dplyr::filter(time == min(time)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = values, values_from = time) |>
    dplyr::mutate(terminus = ifelse(is.na(terminus), Inf, terminus))

  networkDynamic::activate.vertices(
    nw,
    v = d_active$uids,
    onset = d_active$onset,
    terminus = d_active$terminus + 1
  )

  for (item in names(attr_hist)) {
    if (item == "active") next
    d_item <- dplyr::filter(attr_hist[[item]], sim == sim_num) |>
      dplyr::select(time, uids, values) |>
      dplyr::arrange(time)
    for (d_t in split(d_item, d_item$time)) {
      networkDynamic::activate.vertex.attribute(
        nw,
        item,
        v = d_t$uids,
        value = d_t$values,
        onset = d_t$time[[1]],
        terminus = Inf
      )
    }
  }

  nw
}
