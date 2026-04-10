
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

#' @title Build a networkDynamic Object from Simulation Output
#'
#' @description
#' Reconstructs a `networkDynamic` object from a completed `netsim`
#' simulation, using the cumulative edgelist for edge spells and the
#' tracked attribute history for vertex activity and time-varying
#' attributes. The resulting object can be used with `ndtv` for
#' visualization or with `tsna` for temporal network analysis.
#'
#' @param sim An `EpiModel` object of class `netsim`.
#' @param sim_num The simulation number to use (default 1).
#' @param network Optional network number to filter edges. If `NULL`
#'   (default), edges from all networks are included.
#'
#' @return A `networkDynamic` object with active vertex spells and
#'   time-varying vertex attributes.
#'
#' @seealso [get_attr_at()], [get_attr_history()]
#' @export
make_networkDynamic <- function(sim, sim_num = 1, network = NULL) {
  if (!inherits(sim, "netsim"))
    stop("`sim` must be of class netsim")
  if (sim_num > sim$control$nsims || sim_num < 1)
    stop("Specify a single sim_num between 1 and ", sim$control$nsims)
  if (is.null(sim$cumulative.edgelist)) {
    stop("Cumulative edgelist not saved in netsim object. ",
         "Check control.net settings.")
  }
  if (length(sim$control$tracked.attributes) == 0) {
    stop("At least the `active` attribute must be tracked. ",
         "Check control.net `tracked.attributes` settings.")
  }

  attr_hist <- get_attr_history(sim)
  d_active <- attr_hist$active[attr_hist$active$sim == sim_num, ]
  n_nodes <- max(d_active$uids)
  n_steps <- sim$control$nsteps

  # Setup Nodes
  nw <- network::network.initialize(n_nodes, directed = FALSE)
  spells <- get_nodes_spell(d_active)
  networkDynamic::activate.vertices(
    nw,
    v = spells$uid,
    onset = spells$onset,
    terminus = spells$terminus
  )

  # Setup Edges
  el <- sim$cumulative.edgelist[[paste0("sim", sim_num)]] |>
    dplyr::mutate(stop = ifelse(is.na(stop), Inf, stop))

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
  networkDynamic::reconcile.edge.activity(nw, "reduce.to.vertices")


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
