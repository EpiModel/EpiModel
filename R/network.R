
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
#' `networkDynamic` object are the `unique_id`s of the nodes. Each pair of
#' nodes tied on a layer has one edge, with one activity spell per period the
#' tie was active, and the layer in the `network` edge attribute. With
#' `network = NULL` and more than one layer, the layers are merged into one
#' network that allows multiple edges (`multiple = TRUE`), so that a pair of
#' nodes tied on two layers has one edge per layer.
#'
#' The nodes are active from the first record of the attribute history: for a
#' new simulation from time 0, as in the `networkDynamic` objects that
#' [netsim()] saves without `tergmLite`, and for a restarted one from
#' `start - 1`, the snapshot its history starts from. An edge that formed
#' before the activity spell of one of its nodes starts is active from the
#' start of that spell, so that no edge is active between inactive vertices.
#'
#' @return A `networkDynamic` object with active vertex spells and
#'   time-varying vertex attributes.
#'
#' @seealso [get_attr_at()], [get_attr_history()]
#' @examples
#' \donttest{
#' nw <- network_initialize(n = 100)
#' est <- netest(nw, formation = ~edges, target.stats = 50,
#'               coef.diss = dissolution_coefs(~offset(edges), 10),
#'               verbose = FALSE)
#' control <- control.net(type = "SI", nsteps = 20, nsims = 1,
#'                        tergmLite = TRUE, resimulate.network = TRUE,
#'                        cumulative.edgelist = TRUE,
#'                        save.cumulative.edgelist = TRUE,
#'                        truncate.el.cuml = Inf,
#'                        tracked.attributes = "status", verbose = FALSE)
#' sim <- netsim(est, param.net(inf.prob = 0.3, act.rate = 1),
#'               init.net(i.num = 10), control)
#'
#' # The full network history of a tergmLite run, with each node's status
#' nd <- make_networkDynamic(sim)
#' nd
#' x10 <- network.collapse(nd, at = 10)
#' table(get_vertex_attribute(x10, "status"))
#' }
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

  el <- sim$cumulative.edgelist[[paste0("sim", sim_num)]]
  el$stop <-  ifelse(is.na(el$stop), Inf, el$stop)
  if (!is.null(network))
    el <- el[el$network == network, , drop = FALSE]

  # Setup Nodes. The edges of several layers may join the same pair of nodes,
  # which needs a network allowing multiple edges
  multi_layer <- length(unique(el$network)) > 1
  nw <- network::network.initialize(n_nodes, directed = FALSE,
                                    multiple = multi_layer)
  spells <- get_nodes_spell(d_active)
  # A new simulation records its first snapshot at time 1, after the network
  # of time 0 that its initial edges come from: the nodes of that snapshot are
  # active from time 0, as in netsim's own networkDynamic
  t_first <- min(d_active$time)
  if (isTRUE(sim$control$start == 1) && epi_start(sim) == 1 && t_first == 1) {
    spells$onset[spells$onset == t_first] <- 0
  }
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
  # An edge that formed before the activity spell of one of its nodes starts
  # (before the first snapshot of a restarted simulation's history) starts
  # with that spell
  head_on <- get_spell_onset(spells, el$head, el$start)
  tail_on <- get_spell_onset(spells, el$tail, el$start)
  el$start <- pmax(el$start, head_on, tail_on)

  # Clip edges if they disapear due one of the node exiting
  # In these cases, the cumulative edgelist and networkDynamic store the end
  # differently. This block reconciles that
  head_term <- get_spell_terminus(spells, el$head, el$start)
  tail_term <- get_spell_terminus(spells, el$tail, el$start)
  el$stop <- pmin(el$stop, head_term - 1L, tail_term - 1L)
  # An edge formed on the step one of its nodes departed is never active
  el <- el[which(el$stop >= el$start), , drop = FALSE]

  # One edge per pair of nodes and layer, with one spell per row of the
  # cumulative edgelist: a tie that dissolves and forms again is the same
  # edge, active again
  key <- paste(pmin(el$head, el$tail), pmax(el$head, el$tail), el$network)
  first_row <- !duplicated(key)
  edge_of_row <- match(key, key[first_row])
  if (any(first_row)) {
    network::add.edges(
      nw,
      tail = el$tail[first_row],
      head = el$head[first_row],
      names.eval = rep(list("network"), sum(first_row)),
      vals.eval = lapply(el$network[first_row], \(x) list(x))
    )
    networkDynamic::activate.edges(
      nw,
      e = edge_of_row,
      onset = el$start,
      terminus = el$stop + 1L
    )
  }
  nw <- networkDynamic::as.networkDynamic(nw)

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
