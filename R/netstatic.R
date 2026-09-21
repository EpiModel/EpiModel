
#' @title Static Network Layer for Network Epidemic Models
#'
#' @description Builds a fixed contact layer, such as households, classrooms,
#'              hospital wards, cabins, or an observed census edgelist, that
#'              [netsim()] accepts anywhere in its layer list alongside the
#'              dynamic layers estimated with [netest()]. The edges of a static
#'              layer never form or dissolve: there is no formation model, no
#'              dissolution model, and no resimulation.
#'
#' @param nw An object of class `network` holding the node set and any vertex
#'        attributes, as passed to [netest()] for the dynamic layers of the
#'        same model. The network must be undirected and not bipartite.
#' @param group.attr Name of a vertex attribute on `nw` whose values partition
#'        the nodes into groups. Every pair of nodes sharing a non-missing
#'        value becomes an edge, so each group is a clique. Nodes with a
#'        missing (`NA`) value are isolates on this layer.
#' @param edgelist A two-column matrix or `data.frame` of node indices giving
#'        the edges of the layer directly. Use this for a fixed contact
#'        structure that is not a set of cliques, such as an observed network.
#' @param arrivals Rule for wiring nodes that arrive during a simulation with
#'        vital dynamics into the layer: `"isolate"` (the default; the new node
#'        has no edges on this layer and its group attribute is `NA`), `"new"`
#'        (each arrival starts a group of size one, with a fresh group id), or
#'        `"join"` (each arrival joins an existing group and is connected to
#'        all of its current members). See Details.
#' @param arrivals.FUN Optional function used with `arrivals = "join"` to
#'        choose the group each arrival joins. See Details.
#'
#' @details
#' A `netest` object carries a formation model, a dissolution model, and a
#' starting network; `netsim` resimulates its layer each time step. Many
#' contact structures are better described as fixed: every pair of
#' co-residents is a household contact for the whole simulation, every pair
#' of pupils in a classroom is a classroom contact, and so on. Before this
#' function, such a layer had to be carried into `netsim` as a parameter and
#' walked by a custom infection module. `netstatic` instead builds a layer
#' object that fills a slot in the layer list, is skipped by the network
#' resimulation and the edges correction, and is read by the built-in
#' infection modules through [discord_edgelist()] like any other layer.
#'
#' The layer is specified in one of three ways. With `group.attr`, all pairs of
#' nodes sharing a value of that attribute are connected (the household case).
#' With `edgelist`, the edges are given directly. With neither, the edges
#' already on `nw` are used, which is the observed-network case. In the first
#' two cases, any edges already on `nw` are discarded.
#'
#' The group attribute must be an integer, numeric, or character vector; a
#' factor is refused because arrivals under `"new"` need to create values the
#' factor does not have. It should also be set on the network passed to
#' [netest()] for the other layers, since `netsim` reads nodal attributes from
#' the first layer in its list; when it is not, `netsim` copies it from the
#' static layer.
#'
#' @section Arrivals:
#' When a model has vital dynamics, each arriving node has to be placed on the
#' static layer. `netsim` applies the rule chosen with `arrivals` inside
#' [arrive_nodes()], once per static layer, after the arrivals module has
#' created the node and set its other attributes:
#'
#'  * `"isolate"`: the node gets no edges and its group attribute is `NA`. A
#'    custom module may place it later.
#'  * `"new"`: the node gets a fresh group id, one larger than the largest
#'    integer id in use (or a new unique string for character ids), and no
#'    edges. Groups then grow only through further `"join"` arrivals, so this
#'    rule suits models in which arrivals are single-person households.
#'  * `"join"`: the node joins an existing group and is connected to every
#'    active member of that group, including other nodes joining the same
#'    group in the same time step. By default the group is drawn with
#'    probability proportional to its current size, which is the distribution
#'    of a randomly chosen existing member. `arrivals.FUN` replaces that draw:
#'    it is called as `arrivals.FUN(dat, at, new_ids, network)` and must return
#'    one group id per element of `new_ids` (an `NA` leaves that node as an
#'    isolate). Because the function receives the full `dat` object it can
#'    read any nodal attribute; the example below places each newborn in a
#'    household that already has a young child, and a function that reads the
#'    group attribute back for `new_ids` lets a custom arrivals module set the
#'    group itself.
#'
#' The rule is applied to every static layer that has a `group.attr`. Layers
#' built from an `edgelist` or from the edges of `nw` support `"isolate"` only.
#'
#' @section Transmission over a static layer:
#' The built-in infection modules treat every layer alike, so a static layer
#' transmits with the model's `inf.prob` and `act.rate` unless those are given
#' per layer as [multilayer()] objects in [param.net()]:
#' `param.net(inf.prob = multilayer(0.45, 0.10), act.rate = multilayer(1, 2))`
#' with the household layer first and the community layer second. The
#' transmission matrix records the layer of each transmission in its
#' `network` column.
#'
#' Nodes departing the population are removed from a static layer with their
#' edges, as on every other layer. Network statistics for the layer are
#' recorded through `nwstats.formula` in [control.net()], with the default
#' `"formation"` meaning `~edges` for a static layer, and the cumulative
#' edgelist records its edges with a start time of 0. [netdx()] refuses a
#' `netstatic` object, since there is no model to diagnose; the group size
#' distribution and mean degree are shown by `print`.
#'
#' @return
#' An object of class `netstatic`, a list with elements:
#'
#'  * **newnetwork:** a `network` object holding the node set, the vertex
#'    attributes of `nw`, and the static edges.
#'  * **group.attr:** the grouping attribute name, or `NULL`.
#'  * **arrivals**, **arrivals.FUN:** the arrival rule.
#'  * **target.stats**, **target.stats.names:** the edge count of the layer,
#'    named `"edges"`, so that the network statistics table printed by
#'    [print.netsim()] and plotted by [plot.netsim()] shows the layer's own
#'    edge count as the target.
#'  * **edapprox:** `FALSE`.
#'  * **summary:** a list with the node count, edge count, mean degree,
#'    number of isolates, and (with `group.attr`) the group size distribution.
#'
#' @seealso [sample_groups()] builds a population from a table of group
#'   types, and [assign_groups()] assigns group ids to an existing population
#'   with composition constraints. [netsim()] runs the simulation, and
#'   [multilayer()] sets per-layer parameters and controls.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Households as cliques under a community TERGM layer
#' hh <- sample_groups(500, c("adult" = 0.3, "adult adult" = 0.4,
#'                            "adult adult child" = 0.2,
#'                            "adult adult child child" = 0.1),
#'                     attr.name = "age")
#' nw <- network_initialize(n = 500)
#' nw <- set_vertex_attribute(nw, "age", hh$age)
#' nw <- set_vertex_attribute(nw, "hh_id", hh$group)
#'
#' est_hh <- netstatic(nw, group.attr = "hh_id", arrivals = "join")
#' est_hh
#' print(est_hh, by = "age")
#'
#' est_com <- netest(nw, formation = ~edges + nodematch("age"),
#'                   target.stats = c(250, 175),
#'                   coef.diss = dissolution_coefs(~offset(edges), 20),
#'                   verbose = FALSE)
#'
#' param <- param.net(inf.prob = multilayer(0.3, 0.05), act.rate = 1)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 50, nsims = 1,
#'                        tergmLite = TRUE, resimulate.network = TRUE,
#'                        verbose = FALSE)
#' sim <- netsim(list(est_hh, est_com), param, init, control)
#'
#' # Share of transmissions on the household layer
#' tm <- get_transmat(sim)
#' mean(tm$network == 1)
#'
#' # Newborns join a household that already has a child
#' est_hh2 <- netstatic(nw, group.attr = "hh_id", arrivals = "join",
#'   arrivals.FUN = function(dat, at, new_ids, network) {
#'     hh <- get_attr(dat, "hh_id")
#'     age <- get_attr(dat, "age")
#'     pool <- hh[which(age == "child" & !is.na(hh))]
#'     pool[sample.int(length(pool), length(new_ids), replace = TRUE)]
#'   })
#' }
#'
netstatic <- function(nw, group.attr = NULL, edgelist = NULL,
                      arrivals = c("isolate", "new", "join"),
                      arrivals.FUN = NULL) {

  if (!inherits(nw, "network")) {
    stop("`nw` must be an object of class `network`.")
  }
  if (is.directed(nw)) {
    stop("`netstatic` supports undirected networks only.")
  }
  if (is.bipartite(nw)) {
    stop("`netstatic` does not support bipartite networks.")
  }
  arrivals <- match.arg(arrivals)
  if (!is.null(arrivals.FUN)) {
    if (!is.function(arrivals.FUN)) {
      stop("`arrivals.FUN` must be a function.")
    }
    if (arrivals != "join") {
      stop("`arrivals.FUN` applies only with `arrivals = \"join\"`.")
    }
  }
  if (!is.null(group.attr) && !is.null(edgelist)) {
    stop("Specify either `group.attr` or `edgelist`, not both.")
  }

  n <- network.size(nw)
  group.sizes <- NULL

  if (!is.null(group.attr)) {
    if (!is.character(group.attr) || length(group.attr) != 1) {
      stop("`group.attr` must be a single attribute name.")
    }
    if (!group.attr %in% list.vertex.attributes(nw)) {
      stop("There is no vertex attribute called `", group.attr, "` on `nw`.")
    }
    group <- get_vertex_attribute(nw, group.attr)
    if (is.factor(group)) {
      stop("The `", group.attr, "` attribute is a factor. Convert it with ",
           "`as.character()` or `as.integer()` before calling `netstatic`.")
    }
    if (!is.numeric(group) && !is.character(group)) {
      stop("The `", group.attr, "` attribute must be numeric or character.")
    }
    el <- group_edgelist(group)
    group.sizes <- table(tabulate(factor(group[!is.na(group)])),
                         dnn = NULL)
  } else if (!is.null(edgelist)) {
    if (arrivals != "isolate") {
      stop("`arrivals = \"", arrivals, "\"` requires `group.attr`; a layer ",
           "built from an edgelist supports `arrivals = \"isolate\"` only.")
    }
    el <- check_static_edgelist(edgelist, n)
  } else {
    if (network.edgecount(nw) == 0) {
      stop("Specify `group.attr` or `edgelist`, or pass a network that ",
           "already has edges.")
    }
    if (arrivals != "isolate") {
      stop("`arrivals = \"", arrivals, "\"` requires `group.attr`; a layer ",
           "built from the edges of `nw` supports `arrivals = \"isolate\"` ",
           "only.")
    }
    el <- check_static_edgelist(as.edgelist(nw)[, 1:2, drop = FALSE], n)
  }

  ## the layer network: node set and attributes of nw, static edges only
  newnetwork <- network_initialize(n)
  for (a in setdiff(list.vertex.attributes(nw), c("na", "vertex.names"))) {
    newnetwork <- set_vertex_attribute(newnetwork, a,
                                       get_vertex_attribute(nw, a))
  }
  network.vertex.names(newnetwork) <- network.vertex.names(nw)
  if (nrow(el) > 0) {
    newnetwork <- add.edges(newnetwork, tail = el[, 1], head = el[, 2])
  }

  deg <- tabulate(c(el[, 1], el[, 2]), nbins = n)

  out <- list()
  out$newnetwork <- newnetwork
  out$group.attr <- group.attr
  out$arrivals <- arrivals
  out$arrivals.FUN <- arrivals.FUN
  out$target.stats <- nrow(el)
  out$target.stats.names <- "edges"
  out$edapprox <- FALSE
  out$summary <- list(n = n,
                      edges = nrow(el),
                      mean.degree = mean(deg),
                      isolates = sum(deg == 0),
                      group.sizes = group.sizes)

  class(out) <- "netstatic"
  return(out)
}

#' @export
print.netstatic <- function(x, by = NULL, digits = 3, ...) {

  s <- x$summary

  cat("EpiModel Static Network Layer")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  if (!is.null(x$group.attr)) {
    cat("\nLayer type: cliques on `", x$group.attr, "`", sep = "")
  } else {
    cat("\nLayer type: fixed edgelist")
  }

  cat("\n\nLayer Summary")
  cat("\n-----------------------")
  cat("\nNodes:", s$n)
  cat("\nEdges:", s$edges)
  cat("\nMean degree:", round(s$mean.degree, digits))
  cat("\nIsolates:", s$isolates)

  if (!is.null(s$group.sizes)) {
    cat("\n\nGroup Size Distribution")
    cat("\n-----------------------\n")
    gs <- rbind(size = as.integer(names(s$group.sizes)),
                groups = as.integer(s$group.sizes))
    colnames(gs) <- rep("", ncol(gs))
    print(gs)
  }

  if (!is.null(by)) {
    if (is.null(x$newnetwork)) {
      cat("\n\nMean degree by `", by, "` is not available: the layer ",
          "network is not stored on this object.", sep = "")
    } else {
      if (!by %in% list.vertex.attributes(x$newnetwork)) {
        stop("There is no vertex attribute called `", by, "` on the layer.")
      }
      el <- as.edgelist(x$newnetwork)
      deg <- tabulate(c(el[, 1], el[, 2]), nbins = network.size(x$newnetwork))
      by.val <- get_vertex_attribute(x$newnetwork, by)
      cat("\n\nMean Degree by `", by, "`", sep = "")
      cat("\n-----------------------\n")
      print(round(tapply(deg, by.val, mean), digits))
    }
  }

  cat("\n\nArrivals:", x$arrivals)
  if (x$arrivals == "join") {
    if (is.null(x$arrivals.FUN)) {
      cat(" (group drawn in proportion to its size)")
    } else {
      cat(" (group chosen by arrivals.FUN)")
    }
  }
  cat("\n")

  invisible()
}

# Is a layer object, or its nwparam record on the dat object, a static layer?
is_static_layer <- function(x) {
  inherits(x, "netstatic")
}

# All within-group pairs for a grouping vector, as a sorted two-column integer
# matrix with tail < head. Nodes with NA are in no group. Groups of the same
# size are processed together, so the cost is one combn() per distinct size.
group_edgelist <- function(group) {
  ids <- which(!is.na(group))
  empty <- matrix(integer(0), ncol = 2)
  if (length(ids) == 0) {
    return(empty)
  }
  g <- group[ids]
  o <- order(g, ids)
  ids <- ids[o]
  g <- g[o]
  sizes <- tabulate(factor(g, levels = unique(g)))
  starts <- cumsum(c(1L, sizes))[seq_along(sizes)]

  pieces <- list()
  for (k in sort(unique(sizes[sizes > 1]))) {
    gk <- which(sizes == k)
    members <- matrix(ids[rep(starts[gk], each = k) + rep(0:(k - 1), length(gk))],
                      ncol = k, byrow = TRUE)
    pairs <- combn(k, 2)
    pieces[[length(pieces) + 1]] <- cbind(
      as.vector(members[, pairs[1, ], drop = FALSE]),
      as.vector(members[, pairs[2, ], drop = FALSE])
    )
  }
  if (length(pieces) == 0) {
    return(empty)
  }
  el <- do.call(rbind, pieces)
  el <- cbind(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]))
  el <- el[order(el[, 1], el[, 2]), , drop = FALSE]
  storage.mode(el) <- "integer"
  dimnames(el) <- NULL
  el
}

# Validate a user-supplied edgelist and put it in canonical form: integer,
# tail < head, no self-loops, no duplicates, sorted.
check_static_edgelist <- function(edgelist, n) {
  if (is.data.frame(edgelist)) {
    edgelist <- as.matrix(edgelist)
  }
  if (!is.matrix(edgelist) || ncol(edgelist) != 2 || !is.numeric(edgelist)) {
    stop("`edgelist` must be a two-column numeric matrix or data.frame of ",
         "node indices.")
  }
  if (nrow(edgelist) == 0) {
    return(matrix(integer(0), ncol = 2))
  }
  if (anyNA(edgelist) || any(edgelist != round(edgelist)) ||
        any(edgelist < 1) || any(edgelist > n)) {
    stop("`edgelist` must contain integer node indices between 1 and ", n, ".")
  }
  if (any(edgelist[, 1] == edgelist[, 2])) {
    stop("`edgelist` contains self-loops.")
  }
  el <- cbind(pmin(edgelist[, 1], edgelist[, 2]),
              pmax(edgelist[, 1], edgelist[, 2]))
  el <- unique(el)
  el <- el[order(el[, 1], el[, 2]), , drop = FALSE]
  storage.mode(el) <- "integer"
  dimnames(el) <- NULL
  el
}

# Append rows to a tergmLite edgelist, keeping its attributes (n and the
# as.edgelist class and metadata) and its sorted order.
add_edges_to_el <- function(el, new) {
  if (NROW(new) == 0) {
    return(el)
  }
  a <- attributes(el)
  out <- rbind(matrix(el, ncol = 2), matrix(new, ncol = 2))
  out <- out[order(out[, 1], out[, 2]), , drop = FALSE]
  a$dim <- dim(out)
  a$dimnames <- NULL
  attributes(out) <- a
  out
}

# Add edges to a static layer in either storage mode.
add_static_edges <- function(dat, network, el_new) {
  if (NROW(el_new) == 0) {
    return(dat)
  }
  at <- get_current_timestep(dat)
  el_new <- cbind(pmin(el_new[, 1], el_new[, 2]), pmax(el_new[, 1], el_new[, 2]))
  el_new <- unique(el_new)

  if (get_control(dat, "tergmLite") == TRUE) {
    dat$run$el[[network]] <- add_edges_to_el(dat$run$el[[network]], el_new)
    if (get_network_control(dat, network, "tergmLite.track.duration") == TRUE) {
      dat$run$net_attr[[network]][["lasttoggle"]] <- rbind(
        dat$run$net_attr[[network]][["lasttoggle"]],
        cbind(el_new, at)
      )
    }
  } else {
    nw <- get_network(dat, network = network)
    nw <- networkDynamic::add.edges.active(nw, tail = el_new[, 1],
                                           head = el_new[, 2],
                                           onset = at, terminus = Inf)
    dat <- set_network(dat, nw = nw, network = network)
  }
  return(dat)
}

# Fresh group ids for k arrivals, of the same type as the existing ids.
new_group_ids <- function(group, k) {
  existing <- group[!is.na(group)]
  if (is.numeric(group) || length(existing) == 0) {
    start <- if (length(existing) == 0) 0 else max(existing)
    return(start + seq_len(k))
  }
  ids <- make.unique(c(unique(existing), rep("arrival", k)), sep = "_")
  ids[length(ids) - k + seq_len(k)]
}

# Apply a static layer's arrival rule to the nodes new_ids, which arrive_nodes
# has already added to the layer as isolates. Sets the group attribute for
# the new nodes and, under "join", connects each of them to every active
# member of its group, including other nodes joining the same group in this
# time step.
static_layer_arrivals <- function(dat, network, new_ids) {
  nwparam <- get_nwparam(dat, network = network)
  group.attr <- nwparam$group.attr
  rule <- nwparam$arrivals
  if (is.null(group.attr) || rule == "isolate" || length(new_ids) == 0) {
    return(dat)
  }

  active <- get_attr(dat, "active")
  n <- length(active)
  group <- get_attr(dat, group.attr, override.null.error = TRUE)
  if (is.null(group)) {
    group <- rep(NA, n)
  } else if (length(group) < n) {
    group <- c(group, rep(NA, n - length(group)))
  }

  if (rule == "new") {
    group[new_ids] <- new_group_ids(group, length(new_ids))
    dat <- set_attr(dat, group.attr, group)
    return(dat)
  }

  ## the join rule: choose a group for each arrival, then connect it
  if (!is.null(nwparam$arrivals.FUN)) {
    at <- get_current_timestep(dat)
    new_group <- nwparam$arrivals.FUN(dat, at, new_ids, network)
    if (length(new_group) != length(new_ids)) {
      stop("`arrivals.FUN` for static layer ", network, " returned ",
           length(new_group), " group ids for ", length(new_ids),
           " arriving nodes.")
    }
  } else {
    old_ids <- setdiff(which(active == 1), new_ids)
    pool <- group[old_ids]
    pool <- pool[!is.na(pool)]
    if (length(pool) == 0) {
      new_group <- new_group_ids(group, length(new_ids))
    } else {
      new_group <- pool[sample.int(length(pool), length(new_ids),
                                   replace = TRUE)]
    }
  }
  group[new_ids] <- new_group
  dat <- set_attr(dat, group.attr, group)

  joining <- new_ids[!is.na(new_group)]
  if (length(joining) == 0) {
    return(dat)
  }
  in_group <- which(active == 1 & !is.na(group))
  members <- split(in_group, group[in_group])
  el_new <- lapply(joining, function(i) {
    m <- members[[as.character(group[i])]]
    m <- m[m != i]
    if (length(m) == 0) {
      return(NULL)
    }
    cbind(pmin(i, m), pmax(i, m))
  })
  el_new <- do.call(rbind, el_new)
  dat <- add_static_edges(dat, network, el_new)

  return(dat)
}


#' @title Sample a Population of Groups from a Table of Group Types
#'
#' @description Builds a population one group at a time from a table of group
#'              types, such as household compositions, returning each node's
#'              group id and the attributes its group type implies. The result
#'              supplies both the grouping attribute for [netstatic()] and the
#'              nodal attributes for the dynamic layers of the same model.
#'
#' @param n Number of nodes in the population.
#' @param types The group types. Either a character vector in which each
#'        element lists the members of one type separated by `sep`, such as
#'        `c("adult", "adult adult", "adult adult child")`; a named numeric
#'        vector whose names are those templates and whose values are the
#'        type weights, in which case `prob` is taken from the values; or a
#'        list of `data.frame`s, each with one row per member and any number of
#'        attribute columns, for types that carry more than one attribute per
#'        member.
#' @param prob Sampling weights for the types, one per type. They are
#'        normalized to sum to one. The default draws every type with equal
#'        probability.
#' @param attr.name Name of the attribute column in the result when `types` is
#'        a character vector.
#' @param sep Separator between member labels within a character template.
#'
#' @details
#' Types are drawn with replacement until the population reaches `n`. The
#' final group is drawn from the types whose size fits the remaining slots, so
#' that the population has exactly `n` nodes without cutting a group short.
#' When no type fits (no type has size one, and the remainder is smaller than
#' every type), the last group is truncated and a message says so.
#'
#' The realized distribution of group types, and so of group sizes and of the
#' attribute values, follows the `prob` weights up to sampling variation; the
#' expected mean group size is the weighted mean of the type sizes. The
#' expected share of nodes with a given attribute value is the weighted share
#' of that value among the members of all types. Tune the type table to hit
#' both.
#'
#' @return
#' A `data.frame` with one row per node and `n` rows, in group order:
#' `group` (integer group id, `1` to the number of groups), `type` (index of
#' the type the group was drawn from), and the attribute columns (one column
#' named `attr.name` for character templates, or the columns of the
#' `data.frame`s in `types`).
#'
#' @seealso [netstatic()] to turn the `group` column into a static clique
#'   layer, and [assign_groups()] for the reverse problem of assigning group
#'   ids to a population whose attributes already exist.
#'
#' @export
#'
#' @examples
#' # Household types by age group, weights chosen so that every child lives
#' # with at least one adult and about a quarter of older adults live alone
#' hh_types <- c("adult"                    = 0.15,
#'               "adult adult"              = 0.20,
#'               "elderly"                  = 0.12,
#'               "elderly elderly"          = 0.10,
#'               "adult elderly"            = 0.03,
#'               "adult child"              = 0.05,
#'               "adult adult child"        = 0.15,
#'               "adult adult child child"  = 0.15,
#'               "adult adult child elderly" = 0.05)
#' pop <- sample_groups(1000, hh_types, attr.name = "age")
#' head(pop, 10)
#' table(tabulate(pop$group))                 # household size distribution
#' prop.table(table(pop$age))                 # person-level age distribution
#'
#' # Every child shares a household with an adult
#' has_adult <- tapply(pop$age == "adult", pop$group, any)
#' all(has_adult[as.character(pop$group[pop$age == "child"])])
#'
#' # Types with more than one attribute per member
#' types <- list(
#'   data.frame(age = "adult", sex = "F"),
#'   data.frame(age = c("adult", "adult"), sex = c("F", "M")),
#'   data.frame(age = c("adult", "adult", "child"), sex = c("F", "M", "F"))
#' )
#' pop2 <- sample_groups(100, types, prob = c(0.3, 0.4, 0.3))
#' head(pop2)
#'
sample_groups <- function(n, types, prob = NULL, attr.name = "member",
                          sep = " ") {

  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n)) {
    stop("`n` must be a single positive integer.")
  }

  ## normalize types to a list of data.frames, one row per member
  if (is.character(types)) {
    if (!is.null(names(types))) {
      stop("`types` is a named character vector; give the weights in `prob` ",
           "or pass a named numeric vector.")
    }
    templates <- types
  } else if (is.numeric(types)) {
    if (is.null(names(types))) {
      stop("A numeric `types` must be named by the group-type templates.")
    }
    if (!is.null(prob)) {
      stop("Specify the type weights either as the values of `types` or in ",
           "`prob`, not both.")
    }
    prob <- unname(types)
    templates <- names(types)
  } else if (is.list(types)) {
    if (!all(vapply(types, is.data.frame, logical(1)))) {
      stop("A list `types` must contain one `data.frame` per group type.")
    }
    templates <- NULL
  } else {
    stop("`types` must be a character vector, a named numeric vector, or a ",
         "list of data.frames.")
  }
  if (!is.null(templates)) {
    members <- strsplit(trimws(templates), sep, fixed = TRUE)
    types <- lapply(members, function(m) {
      df <- data.frame(m, stringsAsFactors = FALSE)
      names(df) <- attr.name
      df
    })
  }
  K <- length(types)
  sizes <- vapply(types, nrow, integer(1))
  if (K == 0 || any(sizes == 0)) {
    stop("Every group type must have at least one member.")
  }
  if (any(c("group", "type") %in% names(types[[1]]))) {
    stop("The attribute columns of `types` cannot be named `group` or `type`.")
  }

  if (is.null(prob)) {
    prob <- rep(1, K)
  }
  if (length(prob) != K || any(is.na(prob)) || any(prob < 0) || sum(prob) == 0) {
    stop("`prob` must be one non-negative weight per group type.")
  }
  prob <- prob / sum(prob)

  ## draw whole groups until the population reaches n
  mean.size <- sum(prob * sizes)
  draw <- integer(0)
  total <- 0
  while (total < n) {
    n.draw <- ceiling(1.2 * (n - total) / mean.size) + 10
    batch <- sample.int(K, n.draw, replace = TRUE, prob = prob)
    cs <- total + cumsum(sizes[batch])
    keep <- which(cs <= n)
    draw <- c(draw, batch[keep])
    total <- if (length(keep) > 0) cs[max(keep)] else total
    if (length(keep) < n.draw) {
      break
    }
  }

  ## fill the remaining slots with types that fit; truncate as a last resort
  truncate <- 0L
  remaining <- n - total
  while (remaining > 0) {
    fit <- which(sizes <= remaining)
    if (length(fit) == 0) {
      t <- sample.int(K, 1, prob = prob)
      truncate <- as.integer(sizes[t] - remaining)
      draw <- c(draw, t)
      remaining <- 0
      message("No group type fits the ", remaining, " remaining slots; the ",
              "last group is truncated.")
    } else {
      t <- fit[sample.int(length(fit), 1, prob = prob[fit])]
      draw <- c(draw, t)
      remaining <- remaining - sizes[t]
    }
  }

  ## expand the drawn types into node rows, one block per type
  pieces <- lapply(seq_len(K), function(t) {
    gidx <- which(draw == t)
    if (length(gidx) == 0) {
      return(NULL)
    }
    k <- sizes[t]
    piece <- types[[t]][rep(seq_len(k), times = length(gidx)), , drop = FALSE]
    piece <- cbind(group = rep(gidx, each = k), type = t, piece)
    piece
  })
  out <- do.call(rbind, pieces)
  out <- out[order(out$group), , drop = FALSE]
  if (truncate > 0) {
    out <- out[seq_len(nrow(out) - truncate), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}


#' @title Assign Group Ids to an Existing Population
#'
#' @description Assigns nodes to groups of a target size distribution, subject
#'              to a composition rule that keeps dependent members (such as
#'              children) in groups that contain an anchor member (such as an
#'              adult). This is the reverse of [sample_groups()]: the nodal
#'              attributes already exist, from a census age distribution for
#'              example, and the group ids are built to match them.
#'
#' @param size.dist Distribution of group sizes: a numeric vector of weights
#'        named by group size, such as `c("1" = 0.28, "2" = 0.34, "3" = 0.15,
#'        "4" = 0.14, "5" = 0.09)`. An unnamed vector is taken over sizes
#'        `1, 2, ...`. Weights are normalized to sum to one.
#' @param role Optional vector of node roles, one per node, such as an age
#'        group. Its length sets the population size.
#' @param anchor Values of `role` that can anchor a group. Every group that
#'        contains a `dependent` member also contains an anchor, as long as
#'        there are enough anchors.
#' @param dependent Values of `role` that must share a group with an anchor.
#' @param n Population size when `role` is not given.
#'
#' @details
#' The algorithm has four steps:
#'
#'  1. Group sizes are drawn from `size.dist` until they sum to at least `n`;
#'     the last group is trimmed so that the sizes sum to exactly `n`.
#'  2. One anchor is placed in each group, in random order, until either the
#'     groups or the anchors run out. When there are fewer anchors than
#'     groups, the groups without an anchor cannot receive dependents.
#'  3. Dependents are placed in the open slots of anchored groups, each open
#'     slot being equally likely, so that a group receives dependents in
#'     proportion to its remaining size. When the dependents outnumber those
#'     slots, the surplus is added to random anchored groups, which then
#'     exceed their drawn size; this is the only case in which the realized
#'     size distribution departs from `size.dist`, and a message reports it.
#'  4. All remaining nodes (further anchors, and roles that are neither) fill
#'     the remaining open slots in the same way.
#'
#' Without `role`, or without `anchor`, nodes are assigned to the drawn group
#' sizes at random and no composition rule applies.
#'
#' The rule guarantees cross-generational contact for every dependent. It does
#' not otherwise control group composition: the number of dependents per
#' group, the pairing of anchors, and whether a third role (older adults, for
#' example) lives alone or with others all follow from the size distribution
#' and the population shares of the roles. When the composition itself is the
#' quantity to control, build the population from a table of group types with
#' [sample_groups()] instead.
#'
#' @return An integer vector of group ids, one per node, with ids running from
#'   `1` to the number of groups.
#'
#' @seealso [sample_groups()] and [netstatic()].
#'
#' @export
#'
#' @examples
#' # Ages from a census-like distribution, then households around them
#' set.seed(1)
#' n <- 1000
#' age <- sample(c("child", "adult", "elderly"), n, replace = TRUE,
#'               prob = c(0.22, 0.60, 0.18))
#' size.dist <- c("1" = 0.28, "2" = 0.34, "3" = 0.15, "4" = 0.14, "5" = 0.09)
#' hh_id <- assign_groups(size.dist, role = age,
#'                        anchor = c("adult", "elderly"), dependent = "child")
#'
#' table(tabulate(hh_id))                          # realized size distribution
#' has_anchor <- tapply(age %in% c("adult", "elderly"), hh_id, any)
#' all(has_anchor[as.character(hh_id[age == "child"])])  # every child has one
#'
#' nw <- network_initialize(n)
#' nw <- set_vertex_attribute(nw, "age", age)
#' nw <- set_vertex_attribute(nw, "hh_id", hh_id)
#' est_hh <- netstatic(nw, group.attr = "hh_id")
#' print(est_hh, by = "age")
#'
assign_groups <- function(size.dist, role = NULL, anchor = NULL,
                          dependent = NULL, n = NULL) {

  if (!is.numeric(size.dist) || length(size.dist) == 0 ||
        any(is.na(size.dist)) || any(size.dist < 0) || sum(size.dist) == 0) {
    stop("`size.dist` must be a vector of non-negative weights over group ",
         "sizes.")
  }
  sizes <- if (is.null(names(size.dist))) {
    seq_along(size.dist)
  } else {
    suppressWarnings(as.integer(names(size.dist)))
  }
  if (anyNA(sizes) || any(sizes < 1)) {
    stop("The names of `size.dist` must be positive integer group sizes.")
  }
  prob <- size.dist / sum(size.dist)

  if (!is.null(role)) {
    if (!is.null(n) && n != length(role)) {
      stop("`n` does not match the length of `role`.")
    }
    n <- length(role)
  }
  if (is.null(n) || !is.numeric(n) || length(n) != 1 || n < 1) {
    stop("Give the population size in `n` or through `role`.")
  }
  n <- as.integer(n)
  if (!is.null(dependent) && is.null(anchor)) {
    stop("`dependent` requires `anchor`.")
  }

  ## 1. group sizes summing to exactly n
  drawn <- integer(0)
  total <- 0L
  mean.size <- sum(prob * sizes)
  while (total < n) {
    k <- ceiling(1.2 * (n - total) / mean.size) + 10
    batch <- sizes[sample.int(length(sizes), k, replace = TRUE, prob = prob)]
    drawn <- c(drawn, batch)
    total <- total + sum(batch)
  }
  cs <- cumsum(drawn)
  last <- which(cs >= n)[1]
  drawn <- drawn[seq_len(last)]
  drawn[last] <- drawn[last] - (cs[last] - n)
  G <- length(drawn)

  ## no composition rule: random assignment to the drawn sizes
  if (is.null(role) || is.null(anchor)) {
    group <- integer(n)
    group[sample.int(n)] <- rep(seq_len(G), drawn)
    return(group)
  }

  is.anchor <- role %in% anchor
  is.dep <- if (is.null(dependent)) rep(FALSE, n) else role %in% dependent
  if (any(is.anchor & is.dep)) {
    stop("`anchor` and `dependent` must not share values.")
  }
  anchors <- sample(which(is.anchor))
  deps <- sample(which(is.dep))
  others <- sample(which(!is.anchor & !is.dep))

  group <- rep(NA_integer_, n)
  open <- drawn

  ## 2. one anchor per group while both last
  n.head <- min(G, length(anchors))
  heads <- anchors[seq_len(n.head)]
  group[heads] <- seq_len(n.head)
  open[seq_len(n.head)] <- open[seq_len(n.head)] - 1L
  anchors <- anchors[-seq_len(n.head)]
  anchored <- seq_len(n.head)

  ## 3. dependents into the open slots of anchored groups
  if (length(deps) > 0) {
    if (length(anchored) == 0) {
      message("No anchors in the population; dependents are placed at random.")
      others <- c(others, deps)
    } else {
      slots <- rep(anchored, open[anchored])
      n.fit <- min(length(deps), length(slots))
      if (n.fit > 0) {
        pick <- slots[sample.int(length(slots), n.fit)]
        group[deps[seq_len(n.fit)]] <- pick
        used <- tabulate(pick, nbins = G)
        open <- open - used
      }
      if (length(deps) > n.fit) {
        surplus <- deps[-seq_len(n.fit)]
        message(length(surplus), " dependents exceed the open slots of ",
                "anchored groups and are added to random anchored groups.")
        group[surplus] <- anchored[sample.int(length(anchored),
                                              length(surplus), replace = TRUE)]
      }
    }
  }

  ## 4. everyone else into the remaining open slots
  rest <- c(anchors, others)
  if (length(rest) > 0) {
    slots <- rep(seq_len(G), pmax(open, 0L))
    n.fit <- min(length(rest), length(slots))
    if (n.fit > 0) {
      group[rest[seq_len(n.fit)]] <- slots[sample.int(length(slots), n.fit)]
    }
    if (length(rest) > n.fit) {
      surplus <- rest[-seq_len(n.fit)]
      group[surplus] <- sample.int(G, length(surplus), replace = TRUE)
    }
  }

  ## drop ids of groups that ended up empty (possible only after surplus)
  group <- as.integer(factor(group, levels = sort(unique(group))))
  group
}
