
#' @title Convert transmat infection tree into a phylo object
#'
#' @param x An object of class \code{"transmat"}, the output from
#'        \code{\link{get_transmat}}.
#' @param collapse.singles logical, (default TRUE) should
#'        \code{\link[ape]{collapse.singles}} be called on the phylo object
#'        before it is returned? (many infection trees contain intermediate
#'        nodes that must be removed to be a proper phylo tree)
#' @param ...  further arguments (unused)
#'
#' @details
#' Converts the edgelist matrix in the transmat object into a phylo object by
#' doing the required reordering and labeling.  Converts the infection timing into
#' elapsed time from parents' infections to be appropriate for the
#' \code{edge.length} component.  If the the tree does not have the appropriate
#' structure to be a phylogenetic tree (chains of multiple vertices with no
#' branches) the branches will be collapsed (depending on \code{collapse.singles})
#' and labeled with the latest vertex in the chain. Does not yet support infection
#' trees with multiple sources.
#'
#' @export as.phylo.transmat
#'
#' @examples
#' set.seed(10)
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 1)
#' control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose = FALSE,
#'                        use.pids = FALSE)
#'
#' mod1 <- netsim(est1, param, init, control)
#' tm <- get_transmat(mod1)
#' tmPhylo <- as.phylo.transmat(tm)
#' plot(tmPhylo, show.node.label = TRUE, cex = 0.7)
#'
as.phylo.transmat <- function(x, collapse.singles = TRUE, ...) {

  tm <- x
  el <- cbind(tm$inf,tm$sus)

  origNodes <- unique(el[,1])
  origTips <- setdiff(el[,2],origNodes)

  phyloTips <- 1:length(origTips)
  phyloNodes <- rep(NA,length(origNodes))
  phylo.label <- phyloNodes  # make a list to store the labels

  # store the duration of time between the infection and its parent's infection
  durations <- sapply(1:nrow(el), function(r) {

    # convert corresponding infection time to a duration rather than clock
    parentRow <- which(el[, 2] == el[r, 1]) # find edgelist row for v's parent

    #  parent will be missing if it is the root
    if (length(parentRow) > 0) {
      # set the time to the difference between v's infection time and v's
      # parent's infection time
      tm[["at"]][r] - tm[["at"]][parentRow]
    } else {
      tm[["at"]][r]  # assume infection started at time 0
    }
  })

  # find roots (infectors that never appear as sus)
  v <- setdiff(unique(el[, 1]), unique(el[, 2]))
  if (length(v) > 1) {
    warning("found multiple trees ", length(v), " not yet supported")
  }
  # figure out the ordering such that the root
  # node will be one larger than the tip nodes

  phyloN <- 1
  while (length(v) > 0) {
    origIndex <- which(origNodes == v[1])
    if (length(origIndex) > 0) {
      # add the element on the list of phylo nodes
      phyloNodes[origIndex] <- phyloN + length(origTips)
      # copy the old label to new position
      phylo.label[phyloN] <- origNodes[origIndex]
      phyloN <- phyloN + 1
    }
    kids <- el[el[, 1] == v[1], 2] # look up kids on the edgelist
    v <- c(v[-1], kids)
  }

  elTips <- el[, 2] %in% origTips
  # translate the edgelist to the new ids
  el[elTips, 2] <- phyloTips[match(el[elTips, 2], origTips)]
  el[!elTips, 2] <- phyloNodes[match(el[!elTips, 2], origNodes)]
  el[, 1] <- phyloNodes[match(el[, 1], origNodes)]

  out <- list()
  out[["edge"]] <- el
  out[["Nnode"]] <- length(phyloNodes)  # number of non-tip nodes
  out[["tip.label"]] <- origTips
  out[["node.label"]] <- phylo.label
  out[["root.edge"]] <- 0 # have to assume sim started at 0
  out[["edge.length"]] <- durations

  class(out) <- "phylo"
  if (collapse.singles == TRUE) {
    out <- ape::collapse.singles(out)
  }
  return(out)
}

#' @title Converts transmat infection tree into a network object
#'
#' @description Converts the edges of the infection tree described in the
#'              transmat object into a \code{\link{network}} object, copying in
#'              appropriate edge attributes for 'at', 'infDur', 'transProb',
#'              'actRate', and 'finalProb' and constructing a vertex attribute
#'              for 'at'.
#'
#' @param x an object of class \code{transmat} to be converted into a network
#'        object
#' @param ... unused
#'
#' @method as.network transmat
#' @export
#'
as.network.transmat <- function(x, ...){
  tm <- x
  ids <- unique(c(tm$sus, tm$inf))

  # remap the ids to a new set
  tm$sus <- match(tm$sus,ids)
  tm$inf <- match(tm$inf,ids)
  net <- network(cbind(tm$inf, tm$sus), matrix.type = "edgelist")
  network.vertex.names(net) <- ids

  # add the other attributes to edges and vertices
  #net%e%"at"<-tm$at
  set.edge.attribute(net, "at", tm$at)
  set.edge.attribute(net, "infDur", tm$infDur)
  set.edge.attribute(net, "transProb", tm$transProb)
  set.edge.attribute(net, "actRate", tm$actRate)
  set.edge.attribute(net, "finalProb", tm$finalProb)

  net %v% "at" <- 0
  set.vertex.attribute(net, "at", tm$at, v = tm$sus)

  return(net)
}

#' @title Plot transmat infection tree in one of several styles
#'
#' @description Plots the infection tree described in a \code{\link{transmat}}
#'              object in one of several styles: phylogentic tree, an un-rooted
#'              network, a hierarchical tree, or a transmissionTimeline.
#'
#' @param x a \code{\link{transmat}} object to be plotted
#' @param style character name of plot style. One of "phylo", "network", "gv_tree"
#'        or "transmissionTimeline"
#' @param ...  additional plot arguments to be passed to lower-level plot
#'        functions (plot.network, etc)
#'
#' @details The phylo plot requires the \code{ape} package. The gv_tree and
#' \code{ndtv::transmissionTimeline} require that the \code{ndtv} package
#' is installed, and the gv_tree requires a working Graphviz installation on the
#' system. \code{ndtv::install.graphviz}. All of the options are essentially
#' wrappers to other plot calls with some appropriate preset arguments.
#'
#' @export
#' @method plot transmat
#'
#' @seealso \code{\link{plot.network}},\code{\link[ape]{plot.phylo}}
#'
plot.transmat <- function(x,
                          style=c("phylo", "network", "transmissionTimeline"),
                          ...) {

  style <- match.arg(style)

  switch(style,
    "transmissionTimeline" = tm_transsmissionTree_plot(x, ...),
    "network" = plot.network(as.network(x), ...),
    # "gv_tree" = tm_gv_tree_plot(x, ...),
    "phylo" = plot(as.phylo.transmat(x), show.node.label = TRUE, cex = 0.7)
  )

}


# this is a wrapper that uses ndtv and graphviz to make a tree plot
tm_gv_tree_plot <- function(tm, ...) {

  # assumes graphviz is installed
  net <- as.network(tm)

  # calculate coords for transmission tree
  treeCoords <- ndtv::network.layout.animate.Graphviz(net,
                                              layout.par = list(gv.engine = "dot",
                                                                gv.args = "-Granksep=2"))
  treeCoords <- ndtv::layout.normalize(treeCoords, keep.aspect.ratio = FALSE)

  # peek at it
  plot(net, coord = treeCoords, displaylabels = TRUE, jitter = FALSE,
       label.pos = 2, label.cex = 0.7, ...)
}

#' @export
#' @aliases transmat
#' @rdname get_transmat
is.transmat <- function(x) {
  if ("transmat" %in% class(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# this is a wrapper to load the namespace and call the transmissionTimeline
tm_transsmissionTree_plot <- function(x, ...) {
  ndtv::transmissionTimeline(x, ...)
}

