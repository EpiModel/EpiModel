
#' @title Observed Network Layer for Network Epidemic Models
#'
#' @description Wraps a fully observed contact network, either a static
#'              `network` or a `networkDynamic` with edge spells, as a layer
#'              that [netsim()] accepts anywhere in its layer list alongside
#'              the layers estimated with [netest()]. The simulation reads the
#'              edges active at each time step from the observed object; there
#'              is no formation model, no dissolution model, and no
#'              resimulation. It is an addition to `netest` for the case in
#'              which the whole network has been observed, not a substitute
#'              for estimation from a sample (see Details).
#'
#' @param nw An object of class `network` (a static census: the edges are
#'        active at every time step) or `networkDynamic` (a dynamic census:
#'        the edges active at time step `at` are those whose spells cover
#'        `at`). The network must be undirected and not bipartite. Vertex
#'        attributes on `nw` are carried into the simulation as nodal
#'        attributes; temporally extended vertex attributes (`*.active`) are
#'        dropped with a message, since `netsim` keeps its own disease status.
#' @param window Observation window of a dynamic census, as `c(start, end)`.
#'        Defaults to the `net.obs.period` network attribute when `nw` has
#'        one, and otherwise to the range of the finite edge spell times. Used
#'        to warn when a simulation runs past the end of the observations.
#'
#' @details
#' [netest()] exists to turn partial, usually egocentric, network data into a
#' generative model that `netsim` can simulate from. That is the general
#' workflow, because samples are the norm and the fit is what lets a sample
#' stand in for a population. When the whole network has been observed, for
#' every node and every contact and, for a dynamic census, over every time
#' step, there is nothing to estimate: the observed object is what `netsim`
#' would otherwise have to generate. `netcensus` places it in the layer list
#' directly. Sensor, proximity-logger, contact-tracing, and animal-tracking
#' datasets are the usual sources. It is the wrong tool for a sample from
#' which one wants to generalize, and it is not a way to model how the
#' observed ties arise; a network whose ties are to be reproduced from their
#' predictors is estimated with `netest`, whatever their turnover.
#'
#' For a `networkDynamic` census, EpiModel time step `at` reads the edges
#' active at time `at` of the observed object, so the simulation clock is the
#' observation clock. Match `nsteps` in [control.net()] to the observation
#' window: by the `networkDynamic` convention, edges active at the last
#' observed time stay active indefinitely, so a simulation that runs past the
#' window sees a frozen edge set. [netsim()] warns when `nsteps` reaches the
#' end of `window`. Vertex activity spells are not used; every node is present
#' throughout.
#'
#' A census has a fixed node set, so vital dynamics are refused: [netsim()]
#' stops when the model has arrivals or departures and a `netcensus` layer,
#' and [arrive_nodes()] and [depart_nodes()] stop if a custom module tries to
#' add or remove nodes. Duration tracking under `tergmLite` is also refused
#' for a dynamic census, whose observed spells already carry the durations.
#'
#' The layer is skipped by the network resimulation and the edges correction,
#' and the built-in infection modules read it through [discord_edgelist()]
#' like any other layer, so per-layer `inf.prob` and `act.rate` may be given
#' with [multilayer()] in [param.net()]. Under `tergmLite`, the active
#' edgelist is extracted from the observed object at every time step; without
#' `tergmLite`, the `networkDynamic` object is used as stored. Network
#' statistics are recorded through `nwstats.formula` in [control.net()], with
#' the default `"formation"` meaning `~edges`. [netdx()] refuses a `netcensus`
#' object, since there is no model to diagnose; `print` summarizes the
#' observed edges over the window.
#'
#' @return
#' An object of class `netcensus`, a list with elements:
#'
#'  * **newnetwork:** the observed network, with temporally extended vertex
#'    attributes removed.
#'  * **dynamic:** `TRUE` for a `networkDynamic` census, `FALSE` for a static
#'    one.
#'  * **window:** the observation window of a dynamic census, or `NULL`.
#'  * **edapprox:** `FALSE`.
#'  * **summary:** a list with the node count, the vertex attribute names,
#'    and either the edge count, mean degree, and number of isolates (static)
#'    or the number of distinct edges ever observed and the number of edges
#'    active at each integer time in the window (dynamic).
#'
#' @seealso [netclique()] for a layer of cliques defined by a grouping
#'   attribute, which is the other model-free layer. [netsim()] runs the
#'   simulation, and [multilayer()] sets per-layer parameters and controls.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # A dynamic census: an observed networkDynamic with edge spells
#' library(networkDynamicData)
#' data(concurrencyComparisonNets)
#' obs <- netcensus(base)
#' obs
#'
#' param <- param.net(inf.prob = 0.5, act.rate = 1)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        resimulate.network = FALSE, verbose = FALSE)
#' sim <- netsim(obs, param, init, control)
#' plot(sim)
#'
#' # The same census under tergmLite, with the active edges read each step
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5,
#'                        tergmLite = TRUE, verbose = FALSE)
#' sim <- netsim(obs, param, init, control)
#'
#' # A static census next to an estimated layer
#' nw <- network_initialize(n = 100)
#' nw <- add.edges(nw, tail = 1:50, head = 51:100)
#' est <- netest(nw, ~edges, target.stats = 40,
#'               coef.diss = dissolution_coefs(~offset(edges), 10),
#'               verbose = FALSE)
#' sim <- netsim(list(netcensus(nw), est),
#'               param.net(inf.prob = multilayer(0.3, 0.1)),
#'               init.net(i.num = 5),
#'               control.net(type = "SI", nsteps = 20, nsims = 1,
#'                           tergmLite = TRUE, verbose = FALSE))
#' }
#'
netcensus <- function(nw, window = NULL) {

  if (!inherits(nw, "network")) {
    stop("`nw` must be an object of class `network` or `networkDynamic`.")
  }
  if (is.directed(nw)) {
    stop("`netcensus` supports undirected networks only.")
  }
  if (is.bipartite(nw)) {
    stop("`netcensus` does not support bipartite networks.")
  }
  dynamic <- networkDynamic::is.networkDynamic(nw)
  n <- network.size(nw)

  ## temporally extended vertex attributes are not nodal attributes netsim
  ## can use, and a stored disease status would collide with its own
  teas <- grep("\\.active$", list.vertex.attributes(nw), value = TRUE)
  if (length(teas) > 0) {
    message("Dropping temporally extended vertex attribute",
            if (length(teas) > 1) "s" else "", " ",
            paste0("`", teas, "`", collapse = ", "),
            " from the census network; netsim keeps its own disease status.")
    for (a in teas) {
      nw <- delete.vertex.attribute(nw, a)
    }
  }

  if (dynamic) {
    if (is.null(window)) {
      window <- census_window(nw)
    } else if (!is.numeric(window) || length(window) != 2 ||
                 anyNA(window) || window[1] >= window[2]) {
      stop("`window` must be a numeric vector `c(start, end)` with ",
           "`start < end`.")
    }
    times <- seq(ceiling(window[1]), floor(window[2]))
    if (length(times) > 1000) {
      times <- round(seq(times[1], times[length(times)], length.out = 1000))
    }
    edges.active <- vapply(times, function(at) {
      NROW(networkDynamic::get.dyads.active(nw, at = at))
    }, numeric(1))
    names(edges.active) <- times
    summary <- list(n = n,
                    dynamic = TRUE,
                    window = window,
                    edges.ever = network.edgecount(nw),
                    edges.active = edges.active,
                    attributes = census_attributes(nw))
  } else {
    if (!is.null(window)) {
      stop("`window` applies to a `networkDynamic` census only.")
    }
    if (network.edgecount(nw) == 0) {
      stop("`nw` has no edges.")
    }
    el <- as.edgelist(nw)
    deg <- tabulate(c(el[, 1], el[, 2]), nbins = n)
    summary <- list(n = n,
                    dynamic = FALSE,
                    edges = nrow(el),
                    mean.degree = mean(deg),
                    isolates = sum(deg == 0),
                    attributes = census_attributes(nw))
  }

  out <- list()
  out$newnetwork <- nw
  out$dynamic <- dynamic
  out$window <- window
  out$edapprox <- FALSE
  out$summary <- summary

  class(out) <- "netcensus"
  return(out)
}

#' @export
print.netcensus <- function(x, digits = 3, ...) {

  s <- x$summary

  cat("EpiModel Observed Network Layer")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  nwtype <- if (s$dynamic) "dynamic (networkDynamic)" else "static (network)"
  cat("\nNetwork type:", nwtype)

  cat("\n\nLayer Summary")
  cat("\n-----------------------")
  cat("\nNodes:", s$n)
  if (s$dynamic) {
    cat("\nObservation window:", s$window[1], "to", s$window[2])
    cat("\nEdges ever observed:", s$edges.ever)
    ea <- s$edges.active
    cat("\nEdges active per step: min", min(ea), ", median",
        round(stats::median(ea), digits), ", max", max(ea))
    cat("\nMean degree per step:", round(2 * mean(ea) / s$n, digits))
  } else {
    cat("\nEdges:", s$edges)
    cat("\nMean degree:", round(s$mean.degree, digits))
    cat("\nIsolates:", s$isolates)
  }
  cat("\nVertex attributes:",
      if (length(s$attributes) > 0) paste(s$attributes, collapse = ", ")
      else "none")
  cat("\n")

  invisible()
}

# A census layer specifically.
is_census_layer <- function(x) {
  inherits(x, "netcensus")
}

# Vertex attributes a census network carries into the simulation.
census_attributes <- function(nw) {
  setdiff(list.vertex.attributes(nw), c("na", "vertex.names", "active"))
}

# Observation window of a networkDynamic: its net.obs.period when it has one,
# otherwise the range of the finite spell times.
census_window <- function(nw) {
  obs <- nw %n% "net.obs.period"
  if (!is.null(obs) && length(obs$observations) > 0) {
    o <- do.call(rbind, obs$observations)
    return(c(min(o[, 1]), max(o[, 2])))
  }
  sp <- networkDynamic::get.edge.activity(nw, as.spellList = TRUE)
  onset <- sp$onset[is.finite(sp$onset)]
  terminus <- sp$terminus[is.finite(sp$terminus)]
  c(if (length(onset) > 0) min(onset) else 0,
    if (length(terminus) > 0) max(terminus) else Inf)
}

# The edges of a dynamic census active at time `at`, as a tergmLite edgelist:
# sorted, tail < head, with the network size in the `n` attribute.
census_edgelist_at <- function(nw, at, n = network.size(nw)) {
  el <- networkDynamic::get.dyads.active(nw, at = at)
  if (NROW(el) == 0) {
    el <- matrix(integer(0), ncol = 2)
  } else {
    el <- cbind(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2]))
    el <- unique(el)
    el <- el[order(el[, 1], el[, 2]), , drop = FALSE]
  }
  storage.mode(el) <- "integer"
  dimnames(el) <- NULL
  attr(el, "n") <- as.integer(n)
  el
}

# Stop when a census layer is present and a caller wants to change the node
# set, which the observed network cannot follow.
stop_if_census_layer <- function(dat, what) {
  census <- vapply(dat$nwparam, is_census_layer, logical(1))
  if (any(census)) {
    stop("Cannot ", what, ": network ", paste(which(census), collapse = ", "),
         " is a `netcensus` layer with a fixed node set.")
  }
  invisible(NULL)
}
