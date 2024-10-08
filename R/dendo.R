
#' @title Convert transmat Infection Tree into a phylo Object
#'
#' @description Converts a transmission matrix from the \code{get_transmat}
#'              function into a \code{phylo} class object.
#'
#' @param x An object of class \code{transmat}, the output from
#'        \code{\link{get_transmat}}.
#' @param vertex.exit.times Optional numeric vector providing the time of
#'        departure of vertices, to be used to scale the lengths of branches
#'        reaching to the tips. Index position on vector corresponds to network
#'        id. NA indicates no departure, so branch will extend to the end of the
#'        tree.
#' @param ... Further arguments (unused).
#'
#' @details
#' Converts a \code{\link{transmat}} object containing information about the
#' history of a simulated infection into a [`ape::phylo`] object
#' representation suitable for plotting as a tree with
#' \code{\link[ape]{plot.phylo}}. Each infection event becomes a 'node'
#' (horizontal branch) in the resulting \code{phylo} tree, and each network
#' vertex becomes a 'tip' of the tree. The infection events are labeled with the
#' vertex ID of the infector to make it possible to trace the path of infection.
#'
#' The infection timing information is included to position the phylo-nodes,
#' with the lines to the tips drawn to the max time value +1 (unless
#' \code{vertex.exit.times} are passed in it effectively assumes all vertices
#' are active until the end of the simulation).
#'
#' If the \code{transmat} contains multiple infection seeds (there are multiple
#' trees with separate root nodes), this function will return a list of class
#' \code{multiPhylo}, each element of which is a \code{phylo} object. See
#' \code{\link[ape]{read.tree}}.
#'
#' @return A \code{phylo} class object.
#'
#' @importFrom ape as.phylo
#' @export as.phylo.transmat
#'
#' @examples
#' set.seed(13)
#'
#' # Fit a random mixing TERGM with mean degree of 1 and mean edge
#' # duration of 20 time steps
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Parameterize the epidemic model as SI with one infected seed
#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 1)
#' control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose = FALSE)
#'
#' # Simulate the model
#' mod1 <- netsim(est, param, init, control)
#'
#' # Extract the transmission matrix
#' tm <- get_transmat(mod1)
#' head(tm, 15)
#'
#' # Convert to phylo object and plot
#' tmPhylo <- as.phylo.transmat(tm)
#' par(mar = c(1,1,1,1))
#' plot(tmPhylo, show.node.label = TRUE,
#'               root.edge = TRUE,
#'               cex = 0.75)
#'
as.phylo.transmat <- function(x,
                              vertex.exit.times,
                              ...) {

  # if not named properly, assume inf, sus at
  if (!all(c("inf", "sus", "at") %in% names(x))) {
    warning("input does not have appropriate column names for transmat,
            assuming first 3 should be 'inf','sus','at'")
    names(x) <- c("inf", "sus", "at")
  }
  tm <- x
  if (missing(vertex.exit.times)) {
    vertex.exit.times <- NULL
  }
  # find roots (infectors that never appear as sus)
  v <- setdiff(unique(tm$inf), unique(tm$sus))
  if (length(v) > 1) {
    message("found multiple trees, returning a list of ", length(v),
            "phylo objects")
    # need to extract the portions of the edgelist and call seperately
    sub_phylos <- lapply(v, function(v_sub) {
      # walk down the list to find elements below v_sub
      sub_rows <- which(tm$inf == v_sub)
      toFind <- v_sub
      while (length(toFind) > 0) {
        i <- toFind[1]
        sub_rows <- unique(c(sub_rows, which(tm$inf == i)))
        toFind <- c(toFind[-1], tm$sus[which(tm$inf == i)])
      }
      # call as.phylo on the subset of the edgelist
      as.phylo.transmat(tm[sub_rows, , drop = FALSE],
                        vertex.exit.times = vertex.exit.times)

    })
    names(sub_phylos) <- paste("seed", v, sep = "_")
    class(sub_phylos) <- c("multiPhylo", class(sub_phylos))
    return(sub_phylos)
  }

  el <- cbind(tm$inf, tm$sus)
  origNodes <- unique(as.vector(el))
  # if vertex.exit.times included check that it is consistant
  if (!is.null(vertex.exit.times)) {
    if (length(origNodes) > length(vertex.exit.times) ||
          any(origNodes > length(vertex.exit.times))) {
      stop("Vertex ids in edgelist imply a larger network size than
           vertex.exit.times")
    }
  }
  # translate ids in el to sequential integers starting from one
  el[, 1] <- match(el[, 1], origNodes)
  el[, 2] <- match(el[, 2], origNodes)

  maxTip <- max(el)  # need to know what phylo node ids will start

  maxTime <- max(x$at) + 1
  if (!is.null(vertex.exit.times)) {
    maxTime <- max(maxTime, vertex.exit.times, na.rm = TRUE)
  }
  # create new ids for phyloNodes
  phyloNodes <- seq(from = maxTip + 1, length.out = length(origNodes) - 1)
  Nnode <- length(phyloNodes)
  # create labels for each phylo node based on infector id
  phylo.label <- tm$inf

  # set default durations
  # since we don't know how long the graph vertices live, assume entire duration
  durations <- rep(NA, length(phyloNodes) * 2)
  tipExitTimes <- rep(maxTime, maxTip)
  if (!is.null(vertex.exit.times)) {
    # replace any NA values with max time
    vertex.exit.times[is.na(vertex.exit.times)] <- maxTime + 1
    # copy the vertex exit times into the appropriate positions in the
    # durations array
    durations[seq_len(maxTip)] <- vertex.exit.times[origNodes]
    # reorder the vertex.exit times to match new ids of tips
    tipExitTimes <- vertex.exit.times[origNodes]
  }

  # create a new edgelist by stepping through the existing edgelist
  # and creating the new links from phylo nodes to graph vertices (tips)
  # and from phylo node to phylo node
  # have to do this as progressive modifications

  # assume at least one xmit has occured
  # create the phylo node linking to the first
  # infector and infectee
  phyloEl <- rbind(cbind(phyloNodes[1], el[1, 1]),
                   cbind(phyloNodes[1], el[1, 2]))

  durations[1] <- tipExitTimes[el[1, 1]] - tm[["at"]][1]
  durations[2] <- tipExitTimes[el[1, 2]] - tm[["at"]][1]

  phyloN <- 1
  # loop over remaining rows
  if (nrow(el) > 1) {
    for (r in 2:nrow(el)) {
      # find id of infector
      infector <- el[r, 1]
      # find the phylo row of phylo node corresponding to the infector
      phyNRow <- which(phyloEl[, 2] == infector)
      # replace the infector with a new phylo node
      phyloEl[phyNRow, 2] <- phyloNodes[phyloN + 1]
      # link the new phylo node to the infector
      phyloEl <- rbind(phyloEl, cbind(phyloNodes[phyloN + 1], infector))
      # link the new phylo node to the infectee (tip)
      phyloEl <- rbind(phyloEl, cbind(phyloNodes[phyloN + 1], el[r, 2]))

      # update the timing on the replaced row that linked to tip
      durations[phyNRow] <-
        durations[phyNRow] - (tipExitTimes[infector] - tm[["at"]][r])
      # add timings for new rows equal to remaining time
      # infector
      durations[nrow(phyloEl) - 1] <- tipExitTimes[infector] - tm[["at"]][r]
      # infectee
      durations[nrow(phyloEl)] <- tipExitTimes[el[r, 2]] - tm[["at"]][r]


      # increment the phylo node counter
      phyloN <- phyloN + 1
    }
  }

  # format the output
  out <- list()
  out[["edge"]] <- phyloEl
  out[["Nnode"]] <- Nnode  # number of non-tip nodes
  out[["tip.label"]] <- origNodes
  out[["node.label"]] <- phylo.label
  out[["root.edge"]] <- x$at[1] # have to assume sim started at 0
  out[["edge.length"]] <- durations

  class(out) <- "phylo"
  return(out)
}

#' @title Convert transmat Infection Tree into a network Object
#'
#' @description Converts a transmission matrix from the \code{get_transmat}
#'              function into a [`network::network`] class object.
#'
#' @param x An object of class \code{transmat} to be converted into a
#'          \code{network} class object.
#' @param ... Unused.
#'
#' @details
#' When converting from a \code{transmat} to a \code{network} object, this
#' functions copies the edge attributes within the transmission matrix
#' (\code{'at'}, \code{'infDur'}, \code{'transProb'}, \code{'actRate'}, and
#' \code{'finalProb'}) into edge attributes on the network.
#'
#' @return A [`network::network`] object.
#'
#' @method as.network transmat
#' @export
#'
as.network.transmat <- function(x, ...) {
  tm <- x
  ids <- unique(c(tm$sus, tm$inf))

  # remap the ids to a new set
  tm$sus <- match(tm$sus, ids)
  tm$inf <- match(tm$inf, ids)
  net <- network(cbind(tm$inf, tm$sus), matrix.type = "edgelist")
  network.vertex.names(net) <- ids

  # add the other attributes to edges and vertices
  set.edge.attribute(net, "at", tm$at)
  set.edge.attribute(net, "infDur", tm$infDur)
  set.edge.attribute(net, "transProb", tm$transProb)
  set.edge.attribute(net, "actRate", tm$actRate)
  set.edge.attribute(net, "finalProb", tm$finalProb)

  net %v% "at" <- 0
  net <- set_vertex_attribute(net, "at", tm$at, v = tm$sus)

  return(net)
}

#' @title Plot transmat Infection Tree in Three Styles
#'
#' @description Plots the transmission matrix tree from from \code{get_transmat}
#'              in one of three styles: a phylogram, a directed network, or
#'              a transmission timeline.
#'
#' @param x A \code{\link{transmat}} object to be plotted.
#' @param style Character name of plot style. One of \code{"phylo"},
#'              \code{"network"}, or \code{"transmissionTimeline"}.
#' @param ... Additional plot arguments to be passed to lower-level plot
#'            functions (\code{plot.network}, \code{plot.phylo}, or
#'            \code{transmissionTimeline}).
#'
#' @details
#' The \code{phylo} plot requires the \code{ape} package. The
#' \code{transmissionTimeline} plot requires that the \code{ndtv} package.
#'
#' @export
#' @method plot transmat
#'
#' @seealso [`network::plot.network`], \code{\link[ape]{plot.phylo}},
#'          \code{\link[ndtv]{transmissionTimeline}}.
#'
plot.transmat <- function(x,
                          style = c("phylo", "network", "transmissionTimeline"),
                          ...) {

  style <- match.arg(style)

  switch(style,
    "transmissionTimeline" = tm_transsmissionTree_plot(x, ...),
    "network" = plot.network(as.network(x), ...),
    "phylo" = plot(as.phylo.transmat(x),
                   show.node.label = TRUE,
                   root.edge = TRUE,
                   label.offset = 0.1,
                   ...)
  )

}

#' @export
#' @aliases transmat
#' @rdname get_transmat
is.transmat <- function(x) {
  if (inherits(x, "transmat")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# this is a wrapper to load the namespace and call the transmissionTimeline plot
tm_transsmissionTree_plot <- function(x, ...) {
  requireNamespace("ndtv")
  ndtv::transmissionTimeline(x, ...)
}
