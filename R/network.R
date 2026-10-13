
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
#' @details
#' The simulation must be run with `cumulative.edgelist = TRUE`,
#' `save.cumulative.edgelist = TRUE`, and `truncate.el.cuml = Inf` in
#' [`control.net`], so that the cumulative edgelist holds every edge, and with
#' at least one attribute in `tracked.attributes` or `tracked.attributes.once`,
#' so that the activity of the nodes is recorded. The vertex ids of the
#' `networkDynamic` object are the `unique_id`s of the nodes. The edges of all
#' the network layers are merged into one network, with the layer of each edge
#' in the `network` edge attribute.
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
  if (!isTRUE(sim$control$truncate.el.cuml == Inf)) {
    stop("`make_networkDynamic()` needs the full cumulative edgelist, but ",
         "the simulation ran with `truncate.el.cuml = ",
         format(NVL(sim$control$truncate.el.cuml, 0)), "`, which drops the ",
         "edges that dissolved earlier than that many time steps before. ",
         "Run it with `truncate.el.cuml = Inf` in `control.net()`.",
         call. = FALSE)
  }
  if (length(sim$control$tracked.attributes) == 0 &&
        length(sim$control$tracked.attributes.once) == 0) {
    stop("At least the `active` attribute must be tracked. ",
         "Check control.net `tracked.attributes` settings.")
  }

  attr_hist <- get_attr_history(sim)
  check_attr_history(attr_hist, sim_num)
  d_active <- attr_hist$active[attr_hist$active$sim == sim_num, ]
  n_nodes <- max(d_active$uids)

  # Setup Nodes
  nw <- network::network.initialize(n_nodes, directed = FALSE)
  spells <- get_nodes_spell(d_active)
  networkDynamic::activate.vertices(
    nw,
    v = spells$uid,
    onset = spells$onset,
    terminus = spells$terminus
  )
  # A vertex without any spell (e.g., a node that departed before a restart,
  # and so is not in the history) would otherwise count as always active
  never_active <- setdiff(seq_len(n_nodes), spells$uid)
  if (length(never_active) > 0) {
    networkDynamic::deactivate.vertices(
      nw,
      v = never_active,
      onset = -Inf,
      terminus = Inf,
      deactivate.edges = FALSE
    )
  }

  # Setup Edges
  el <- sim$cumulative.edgelist[[paste0("sim", sim_num)]]
  el$stop <-  ifelse(is.na(el$stop), Inf, el$stop)

  if (!is.null(network))
    el <- el[el$network == network, , drop = FALSE]

  # Clip edges if they disapear due one of the node exiting
  # In these cases, the cumulative edgelist and networkDynamic store the end
  # differently. This block reconciles that
  head_term <- get_spell_terminus(spells, el$head, el$start)
  tail_term <- get_spell_terminus(spells, el$tail, el$start)
  el$stop <- pmin(el$stop, head_term - 1L, tail_term - 1L)
  # onset is already correct, no need to update `el$start`
  # An edge formed on the step one of its nodes departed is never active
  el <- el[which(el$stop >= el$start), , drop = FALSE]

  networkDynamic::add.edges.active(
    nw,
    head = el$head,
    tail = el$tail,
    onset = el$start,
    terminus = el$stop + 1L,
    names.eval = rep(list("network"), nrow(el)),
    vals.eval = lapply(el$network, \(x) list(network = x))
  )

  for (item in names(attr_hist)) {
    if (item == "active") next
    d_item <- attr_hist[[item]]
    d_item <- d_item[d_item$sim == sim_num, c("time", "uids", "values")]
    nw <- set_vertex_tea(nw, item, d_item)
  }

  nw
}

# Set a time-varying vertex attribute from its history
#
# Writes the `<item>.active` vertex attribute that `networkDynamic` reads,
# with one call for all the vertices: for each vertex, the list of its
# recorded values and the matrix of their spells, each value lasting until
# the next one. This stores the same object as one
# `networkDynamic::activate.vertex.attribute()` call per time step, with
# `terminus = Inf`, at a fraction of the cost. When a vertex has several
# records at the same time step, the last one is used.
#
# @param nw A `network` object, with the vertex ids as `unique_id`s.
# @param item The name of the attribute.
# @param d_item A `data.frame` with columns `time`, `uids` and `values`.
# @return The updated `network` object.
# @noRd
set_vertex_tea <- function(nw, item, d_item) {
  if (nrow(d_item) == 0) return(nw)
  d_item <- d_item[order(d_item$uids, d_item$time), , drop = FALSE]
  d_item <- d_item[!duplicated(d_item[c("uids", "time")], fromLast = TRUE), ,
                   drop = FALSE]
  by_v <- split(seq_len(nrow(d_item)), d_item$uids)
  tea <- lapply(by_v, function(ix) {
    onset <- d_item$time[ix]
    list(
      as.list(d_item$values[ix]),
      cbind(onset, c(onset[-1], Inf), deparse.level = 0)
    )
  })
  v <- d_item$uids[vapply(by_v, `[`, 1L, 1L)]
  network::set.vertex.attribute(nw, paste0(item, ".active"), unname(tea),
                                v = v)
}
