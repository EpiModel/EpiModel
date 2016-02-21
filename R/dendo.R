
#' @title Convert transmat infection tree into a phylo object
#'
#' @description Converts the edgelist matrix in the transmat object into a phylo
#'              object by doing the required reordering and labeling.
#'
#' @param x An object of class \code{"transmat"}, the output from
#'        \code{\link{get_transmat}}.
#' @param collapse.singles logical, DEPRECATED
#' @param ...  further arguments (unused)
#'
#' @details
#' Converts a \code{\link{transmat}} object containing information about the history of a 
#' simulated infection into a \code{\link{phylo}} object representation suitable for plotting
#'  as a tree with \code{\link[ape]{plot.phylo}}.  Each infection 
#' event becomes a 'node' (horizontal branch) in the resulting phylo tree, and each 
#' network vertex becomes a 'tip' of the tree.  The infection events are labled with 
#' the vertex id of the infector to make it possible to trace the path of infection.  
#' 
#' The infection timing information is included to position the phylo-nodes, with the
#'  lines to the tips drawn to the max time value +1 (the transmat does not contain info
#'  on vertex activity, so it effectively assumes all vertices are active/alive until 
#'  the end of the simulation). The function does not yet support infection trees with multiple infection seeds
#' 
#' Note that in EpiModel versions <= 1.2.4, the phylo tree was constructed differently, translating network 
#' vertices to both phylo-nodes and tips and requiring 'collapse.singles' to prune it to an appropriate branching structure.
#'
#' @importFrom ape as.phylo
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
#' plot(tmPhylo, show.node.label = TRUE,
#'               root.edge=TRUE, 
#'               cex = 0.5)
#'
as.phylo.transmat <- function(x, collapse.singles, ...) {
  
  # warnings if somone tries to use old args that are no longer supported
  if(!missing(collapse.singles)){
    warning("the 'collapse.singles' argument to as.phylo.transmat is no longer supported and will be ignored")
  }
  
  tm <- x
  el <- cbind(tm$inf,tm$sus)
  origNodes <- unique(as.vector(el))
  # translate ids in el to sequential integers starting from one
  el[,1]<-match(el[,1],origNodes)
  el[,2]<-match(el[,2],origNodes)
  
  maxTip<-max(el)  # need to know what phylo node ids will start
  maxTime<-max(x$at)+1
  # create new ids for phyloNodes
  phyloNodes <- seq(from=maxTip+1,length.out=length(origNodes)-1)
  Nnode<-length(phyloNodes)
  # create labels for each phylo node based on infector id
  phylo.label <- tm$inf
  # this is an alternate label form like i_j
  #phylo.label <- sapply(1:length(phyloNodes),function(r){
  #  paste(tm[r,'inf'],tm[r,'sus'],sep='_')
  #})  
  
  # set default durations
  # since we don't know how long the graph vertices live, assume entire duration
  durations <-rep(maxTime+1,length(phyloNodes)*2)
  # find roots (infectors that never appear as sus)
  v <- setdiff(unique(el[, 1]), unique(el[, 2]))
  if (length(v) > 1) {
    warning("found multiple trees ", length(v), " not yet supported")
  }
  
  # create a new edgelist by stepping through the existing edgelist
  # and creating the new links from phylo nodes to graph vertices (tips)
  # and from phylo node to phylo node
  # have to do this as progressive modifications
  
  # assume at least one xmit has occured
  # create the phylo node linking to the first 
  # infector and infectee
  phyloEl <-rbind(cbind(phyloNodes[1],el[1,1]),
                  cbind(phyloNodes[1],el[1,2]))
  durations[1]<-maxTime-x$at[1]
  durations[2]<-maxTime-x$at[1]
  phyloN<-1
  # loop over remaining rows
  if(nrow(el)>1){
    for(r in 2:nrow(el)){
      # find id of infector
      infector<-el[r,1]
      # find the phylo row of phylo node corresponding to the infector
      phyNRow<-which(phyloEl[,2]==infector)
      # replace the infector with a new phylo node
      phyloEl[phyNRow,2]<-phyloNodes[phyloN+1]
      # link the new phylo node to the infector
      phyloEl<-rbind(phyloEl,cbind(phyloNodes[phyloN+1],infector))
      # link the new phylo node to the infectee
      phyloEl<-rbind(phyloEl,cbind(phyloNodes[phyloN+1],el[r,2]))
      
      # update the timing on the replaced row that linked to tip
      durations[phyNRow] <- durations[phyNRow] - (maxTime - tm[["at"]][r])
      # add timings for new rows equal to remaining time
      durations[nrow(phyloEl)-1]<- maxTime - tm[["at"]][r] 
      durations[nrow(phyloEl)]<- maxTime - tm[["at"]][r] 
      
      # increment the phylo node counter
      phyloN<- phyloN+1
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
#' @param style character name of plot style. One of "phylo", "network",
#'        or "transmissionTimeline"
#' @param ...  additional plot arguments to be passed to lower-level plot
#'        functions (plot.network, plot.phylo, etc)
#'
#' @details The phylo plot requires the \code{ape} package. The
#' \code{ndtv::transmissionTimeline} requires that the \code{ndtv} package
#' is installed. All of the options are essentially
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
    "phylo" = plot(as.phylo.transmat(x), 
                   show.node.label = TRUE, 
                   root.edge=TRUE,
                   label.offset=0.1,
                   ...
                   )
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
  requireNamespace('ndtv')
  ndtv::transmissionTimeline(x, ...)
}

