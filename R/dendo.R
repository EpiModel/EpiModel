# classes for working with dendograms from ape

# require(tsna)
# data(moodyContactSim)
# v10path<-tPath(moodyContactSim,v=10,start=0)
# 
# getChildren(10,17,v10path$previous)
# plot(v10path)
# as.phylo.tPath(moodyContactSim,v10path)
# this is a draft version (used with tsna code above)
as.phylo.tPath<-function(path){
  
  el<-cbind(path$previous,1:length(path$previous))
  # remove the seed row
  el<-el[path$previous!=0,]
  #phylo doesn't like "non-splitting" nodes that have no children
  singleKidNodes<-as.numeric(names(which(table(el[,1])==1)))
  if (length(singleKidNodes)>0){
    singleKidsRows<-which(el[,1]==singleKidNodes)
    warning('found vertices ',paste(singleKidNodes,collapse=','), ' that create "non-splitting" phylo nodes')
    #el<-el[-singleKidsRows,,drop=FALSE]
  }
  
  origNodes<-unique(el[,1])
  origTips<-setdiff(el[,2],origNodes)
  
  phyloTips<-1:length(origTips)
  phyloNodes<-rep(NA,length(origNodes))
  # figure out the ordering such that the root
  # node will be one larger than the tip nodes
  v<-which(path$previous==0)
  phyloN<-length(phyloTips)+1
  while(length(v)>0){
    origIndex<-which(origNodes==v[1])
    if(length(origIndex)>0){
      phyloNodes[origIndex]<-phyloN
      phyloN<-phyloN+1
    }
    kids<-which(path$previous==v[1])
    v<-c(v[-1],kids)
  }
  
  elTips<-el[,2]%in%origTips
  # translate the edgelist to the new ids
  el[elTips,2]<-phyloTips[match(el[elTips,2],origTips)]
  el[!elTips,2]<-phyloNodes[match(el[!elTips,2],origNodes)]
  el[,1]<-phyloNodes[match(el[,1],origNodes)]
  
  out<-list()
  out[['edge']]<-el
  out[['Nnode']]<-length(phyloNodes)  # number of non-tip nodes
  out[['tip.label']]<-origTips
  out[['node.label']]<-origNodes
  class(out)<-'phylo'
  return(out)
}

#' @title convert transmat infection tree into a phylo object
#' @method as.phylo transmat
#' @export as.phylo.transmat
#' @param x An object of class \code{'transmat'}, the output from \code{\link{get_transmat}}. 
#' @param collapse.singles logical, (default TRUE) should \code{\link[ape]{collapse.singles}} be called on the phylo object before it is returned? (many infection trees contain intermediate nodes that must be removed to be a proper phylo tree)
#' @param ...  further arguments (unused)
#' @details Converts the edgelist matrix in the transmat object into a phylo object by doing the required reordering and labeling.  Converts the infection timing into elapsed time from parents' infections to be appropriate for the \code{edge.length} component.  If the the tree does not have the appropriate structure to be a phylogenetic tree (chains of multiple vertices with no branches) the branches will be collapsed (depending on \code{collapse.singles}) and labeled with the latest vertex in the chain. Does not yet support infection trees with multiple sources.
#' @examples
#' library(EpiModel)
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

#' param <- param.net(inf.prob = 0.5)
#' init <- init.net(i.num = 1)
#' control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose.int = 0, use.pids = FALSE)
#' mod1 <- netsim(est1, param, init, control)
#' tm <- get_transmat(mod1)
#' tmPhylo<-as.phylo.transmat(tm)
#' plot(tmPhylo,show.node.label = TRUE,cex=0.7)
as.phylo.transmat<-function(x,collapse.singles=TRUE,...){
  requireNamespace('ape')
  tm<-x
  el<-cbind(tm$inf,tm$sus)

  origNodes<-unique(el[,1])
  origTips<-setdiff(el[,2],origNodes)
  
  phyloTips<-1:length(origTips)
  phyloNodes<-rep(NA,length(origNodes))
  phylo.label<-phyloNodes  # make a list to store the labels
  # store the duration of time between the infection and its parent's infection
  durations<-sapply(1:nrow(el),function(r){
    # convert corresponding infection time to a duration rather than clock
    parentRow<-which(el[,2]==el[r,1]) # find edgelist row for v's parent
    #  parent will be missing if it is the root
    if (length(parentRow)>0){
      # set the time to the difference between v's infection time and v's parent's infection time
      tm[['at']][r]-tm[['at']][parentRow]
    } else {
      tm[['at']][r]  # assume infection started at time 0
    }
  })

  # find roots (infectors that never appear as sus)
  v<-setdiff(unique(el[,1]),unique(el[,2]))
  if(length(v)>1){
    warning('found multiple trees ',length(v),' not yet supported')
  }
  # figure out the ordering such that the root
  # node will be one larger than the tip nodes
  
  phyloN<-1
  while(length(v)>0){
    origIndex<-which(origNodes==v[1])
    if(length(origIndex)>0){
      # add the element on the list of phylo nodes
      phyloNodes[origIndex]<-phyloN+length(origTips)
      # copy the old label to new position
      phylo.label[phyloN]<-origNodes[origIndex] 
      phyloN<-phyloN+1
    }
    kids<-el[el[,1]==v[1],2] # look up kids on the edgelist
    v<-c(v[-1],kids)
  }
  
  elTips<-el[,2]%in%origTips
  # translate the edgelist to the new ids
  el[elTips,2]<-phyloTips[match(el[elTips,2],origTips)]
  el[!elTips,2]<-phyloNodes[match(el[!elTips,2],origNodes)]
  el[,1]<-phyloNodes[match(el[,1],origNodes)]
  
  out<-list()
  out[['edge']]<-el
  out[['Nnode']]<-length(phyloNodes)  # number of non-tip nodes
  out[['tip.label']]<-origTips
  out[['node.label']]<-phylo.label
  out[['root.edge']]<-0 # have to assume sim started at 0
  out[['edge.length']]<-durations
  class(out)<-'phylo'
  if(collapse.singles){
    out<-ape::collapse.singles(out)
  }
  return(out)
}

#' @title converting transmat infection tree into a network object
#' @method as.network transmat
#' @export as.network.transmat
#' @param x an object of class \code{transmat} to be converted into a network object
#' @param ... unused
#' @description converts the edges of the infection tree described in the transmat object into a \code{\link{network}} object, copying in appropriate edge attributes for 'at', 'infDur', 'transProb', 'actRate', and 'finalProb' and constructing a vertex attribute for 'at'. 
as.network.transmat<-function(x,...){
  tm<-x
  ids<-unique(c(tm$sus,tm$inf))
  # remap the ids to a new set
  tm$sus<-match(tm$sus,ids)
  tm$inf<-match(tm$inf,ids)
  net<-network(cbind(tm$inf,tm$sus),matrix.type='edgelist')
  network.vertex.names(net)<-ids
  # add the other attributes to edges and vertices
  #net%e%'at'<-tm$at
  set.edge.attribute(net,'at',tm$at)
  set.edge.attribute(net,'infDur',tm$infDur)
  set.edge.attribute(net,'transProb',tm$transProb)
  set.edge.attribute(net,'actRate',tm$actRate)
  set.edge.attribute(net,'finalProb',tm$finalProb)

  
  net%v%'at'<-0
  set.vertex.attribute(net,'at',tm$at,v=tm$sus)
  return(net)
}

#' @title plot transmat infection tree in one of several styles
#' @method plot transmat
#' @export plot.transmat
#' @param x a \code{\link{transmat}} object to be plotted
#' @param style character name of plot style
#' @param ...  additional plot arguments to be passed to lower-level plot functions (plot.network, etc)
#' @description plots the infection tree described in a transmat object in one of several styles: phylogentic tree, a network, a hierarchical tree (gv_tree'), or a transmissionTimeline. The gv_tree and transmissionTimeline require that the ndtv package is installed, and the gv_tree requires a working Graphviz installation on the system. \code{\link[ndtv]{install.graphviz}}. All of the options are essentially wrappers to other plot calls with some appropriate preset arguments. 
plot.transmat<-function(x,style=c('phylo','network','gv_tree','transmissionTimeline'),...){
  style<-match.arg(style)
  switch (style,
    'transmissionTimeline' = tm_cascade_plot(x,...),
    'network' = plot.network(as.network(x),...),
    'gv_tree' = tm_gv_tree_plot(x,...),
    'phylo' = plot(as.phylo(x),show.node.label = TRUE,cex=0.7)
  )
}



# this is a wrapper that uses ndtv and graphviz to make a tree plot
tm_gv_tree_plot<-function(tm,...){
# assumes graphviz is installed
requireNamespace('ndtv')
net<-as.network(tm)
# calculate coords for transmission tree
treeCoords<-ndtv::network.layout.animate.Graphviz(net,
                                            layout.par=list(gv.engine='dot',
                                                            gv.args='-Granksep=2'))
treeCoords<-ndtv::layout.normalize(treeCoords,keep.aspect.ratio = FALSE)
# peek at it
plot(net,coord=treeCoords,displaylabels=TRUE,jitter=FALSE,label.pos=2,label.cex=0.7,...)
}

#' @export is.transmat
#' @aliases transmat
#' @rdname get_transmat
is.transmat<-function(x){
  if ('transmat'%in%class(x)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# TODO: this needs to be replaced to call the diffusionTimeline version in ndtv, once that is released to cran. 
tm_cascade_plot<-function(x,time.attr,
                          label = NULL,
                          displaylabels = TRUE,
                          label.cex = 0.7,
                          label.col = 1,
                          vertex.col = 2,
                          vertex.sides = 50,
                          edge.col = 'gray',
                          edge.lty = 1,
                          edge.lwd = 1,
                          xlab='time',
                          ylab='generation',
                          ...){
  net<-as.network(x)
  if(missing(time.attr)){
    if(is.transmat(x)){
      time.attr<-'at'
    } else if ('tPath'%in%class(x)){
      time.attr<-'tdist'
    }
  }
  el<-as.edgelist(net)
  times<-get.vertex.attribute(net,time.attr)
  coords<-matrix(0,nrow=network.size(net),ncol=2)
  yBin<-0
  # find roots
  v<-setdiff(unique(el[,1]),unique(el[,2]))
 
  visited<-integer(0)
  while(length(v)>0){
    visited<-c(visited,v[1])
    coords[v[1],]<-c(times[v[1]],yBin)
    kids<-el[el[,1]==v[1],2] # look up kids on the edgelist
    # in case network was not actually a tree, make sure we don't loop forever
    if(any(kids%in%visited)){
      stop('vertex was revisited: network does not appear to be a tree')
    }
    v<-c(v[-1],kids)
    yBin<-yBin+1
  }
  op <- par(no.readonly = TRUE)
  
  # set up the plotting window
  plot(coords,pch=NA,xlab=xlab,ylab=ylab,...)
  
  # expand the various plot parameters using network defaults
  if(is.null(label)){
    label<-network.vertex.names(net)
  }
  vertex.col<- plotArgs.network(net,'vertex.col',vertex.col)
  vertex.sides<-plotArgs.network(net,'vertex.sides',vertex.sides)
  # remap vertex.sides to a vertex pch approximation
  vertex.pch <- sapply(vertex.sides,function(sides){
    switch (as.character(sides),
    '3' = 24,
    '4' = 22,
    '50' = 21)})
  
  edge.col<-plotArgs.network(net,'edge.col',edge.col)
  edge.lty <- plotArgs.network(net,'edge.lty',edge.lty)
  edge.lwd <- plotArgs.network(net,'edge.lwd',edge.lwd)
  labels <-plotArgs.network(net,'labels',labels)
  label.col<-plotArgs.network(net,'label.col',label.col)
  label.pos<-4
  label.cex<-plotArgs.network(net,'label.cex',label.cex)
  
  # plot the 'edges'
  if(nrow(el)>0){
    lapply(1:nrow(el),function(e){
      lines(c(coords[el[e,1],1],
            coords[el[e,2],1]), # xcoords
            c(coords[el[e,1],2],
            coords[el[e,2],2]),# ycoords
            col=edge.col[e],
            lwd=edge.lwd[e],
            lty=edge.lty[e]
          )
    })
  }
  
  # plot the vertices
  if(network.size(net)>0){
    points(coords[,1],coords[,2],pch=21,bg=vertex.col)
      
    # plot labels
    if(displaylabels){
      text(coords[,1],coords[,2],labels = label,col=label.col,pos=label.pos,cex=label.cex)
    }
  }
  par(op)
}

