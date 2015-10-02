# classes for working with dendograms from ape
require(ape)

#' @title functions for converting to and from phylo class
#' @import ape
# require(tsna)
# data(moodyContactSim)
# v10path<-tPath(moodyContactSim,v=10,start=0)
# 
# getChildren(10,17,v10path$previous)
# plot(v10path)
# as.phylo.tPath(moodyContactSim,v10path)
 

 
# edgelist version 
# as.phylo.tPath<-function(net, path){
# 
#   el<-cbind(path$previous,1:length(path$previous))
#   # remove the seed row
#   el<-el[path$previous!=0,]
#   # 'nodes' are all the vertices with decendents
#   nodes<-unique(path$previous[path$previous!=0])
#   # tips are all the vertices that are not nodes
#   tips<-setdiff(1:length(path$previous),nodes)
#   
#   ids<-c(tips,nodes)
#   
#   # remap numbering
#   el[,1]<-match(el[,1],ids)
#   el[,2]<-match(el[,2],ids)
#   
#   # delete the singles
#   
#   
#   out<-list()
#   out[['edge']]<-el
#   out[['Nnode']]<-length(nodes)  # number of non-tip nodes
#   out[['tip.label']]<-network.vertex.names(net)[tips]
#   out[['node.label']]<-network.vertex.names(net)[nodes]
#   out[['edge.length']]<-path$tdist
#   class(out)<-'phylo'
#   return(out)
# }

# crawler version
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



as.phylo.transmat<-function(tm){
  el<-cbind(tm$inf,tm$sus)
  #phylo doesn't like "non-splitting" nodes that have no children
#   dummyTips<-numeric(0)
#   singleKidNodes<-as.numeric(names(which(table(el[,1])==1)))
#   if (length(singleKidNodes)>0){
#     singleKidsRows<-el[,1]%in%singleKidNodes
#     warning('found vertices ',paste(singleKidNodes,collapse = ','),' that create "non-splitting" phylo nodes, dummies marked with "*" added ')
#     #el<-el[!singleKidsRows,,drop=FALSE]
#     # add some dummy nodes
#     for(row in which(singleKidsRows)){
#       dummyTips<-c(dummyTips,max(el+1))
#       el<-rbind(el,c(el[row,1],max(el)+1))
#     }
#   }
  
  origNodes<-unique(el[,1])
  origTips<-setdiff(el[,2],origNodes)
  
  phyloTips<-1:length(origTips)
  phyloNodes<-rep(NA,length(origNodes))
  phylo.label<-phyloNodes  # make a list to store the labels

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

  # also translate the times
  times<-tm[['at']] # get the times 
  #times[elTips]<-phyloTips[match(el[elTips,2],origTips)]
  #times[!elTips]<-phyloNodes[match(el[!elTips,2],origNodes)]
  
  # if we have dummy nodes, change their labels
  tip.label<-origTips
#   if(length(dummyTips)>0){
#     tip.label[match(dummyTips,origTips)]<-'*'
#   }
  
  out<-list()
  out[['edge']]<-el
  out[['Nnode']]<-length(phyloNodes)  # number of non-tip nodes
  out[['tip.label']]<-tip.label
  out[['node.label']]<-phylo.label
  out[['root.edge']]<-0 # have to assume sim started at 0
  out[['edge.length']]<-times
  class(out)<-'phylo'
  return(out)
}

#' @title converting transmat infection tree into a network object
#' @method as.network transmat
#' @export as.network.transmat
as.network.transmat<-function(tm){
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
plot.transmat<-function(tm,style=c('cascade','network','gv_tree','phylo'),...){
  style<-match.arg(style)
  switch (style,
    'cascade' = tm_cascade_plot(tm,...),
    'network' = plot.network(as.network(tm),...),
    'gv_tree' = tm_gv_tree_plot(tm,...),
    'phylo' = plot(as.phylo(tm),show.node.label = TRUE,cex=0.7)
  )
}

if(FALSE){


phytest<-list(edge=matrix(c(5,4,
                            5,6,
                            6,3,
                            6,7,
                            7,2,
                            7,1),ncol=2,byrow=TRUE),
              Nnode=3,
              tip.label=1:4,
              node.label=5:7,
              edge.length=c(5,1,1,1,1,1))
class(phytest)<-'phylo'
plot(phytest,show.node.label=TRUE,use.edge.length = TRUE)



# this works
phytest<-list(edge=matrix(c(5,1,
                            5,6,
                            6,2,
                            6,7,
                            7,3,
                            7,4),ncol=2,byrow=TRUE),
              Nnode=3,
              tip.label=1:4,
              node.label=5:7,
              edge.length=c(5,1,1,1,1,1))
class(phytest)<-'phylo'
plot(phytest,show.node.label=TRUE,use.edge.length = TRUE)

}



### example code
transmat_example<-function(){
library(EpiModel)
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

param <- param.net(inf.prob = 0.5)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 40, nsims = 1, verbose.int = 0, use.pids = FALSE)
mod1 <- netsim(est1, param, init, control)

tm <- get_transmat(mod1)
tm
}

tm_gv_tree_plot<-function(tm,...){
# assumes graphviz is installed
requireNamespace('ndtv')
net<-as.network(tm)
# calculate coords for transmission tree
treeCoords<-network.layout.animate.Graphviz(net,
                                            layout.par=list(gv.engine='dot',
                                                            gv.args='-Granksep=2'))
treeCoords<-layout.normalize(treeCoords,keep.aspect.ratio = FALSE)
# peek at it
plot(net,coord=treeCoords,displaylabels=TRUE,jitter=FALSE,label.pos=2,label.cex=0.6,...)
}

is.transmat<-function(x){
  if ('transmat'%in%class(x)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

tm_cascade_plot<-function(x,time.attr,
                          label = network.vertex.names(x),
                          displaylabels = !missing(label),
                          label.cex = 1,
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
  if(is.tPath(x)){
    # assume tPath for now
    yBin<-net%v%'gsteps'
    coords<-cbind(times,yBin)
    yBin<-0
  } else {
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
  }
  op <- par(no.readonly = TRUE)
  
  # set up the plotting window
  plot(coords,pch=NA,xlab=xlab,ylab=ylab,...)
  
  # expand the various plot parameters using network defaults
  
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
  label.cex<-plotArgs.network(net,'label.col',label.col)
  
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
      text(coords[,1],coords[,2],labels = labels,col=label.col,pos=label.pos,cex=label.cex)
    }
  }
  par(op)
}

mat<-matrix(c(1,0,0,0,0,
         1,1,0,0,0,
         1,1,1,0,0,
         1,1,0,1,0,
         1,1,0,1,0),ncol=5)
rownames(mat)<-LETTERS[1:5]
plot(nj(dist.gene(mat)),show.node.label = TRUE)
plot(bionj(dist.gene(mat)),show.node.label = TRUE)
plot(fastme.bal(dist.gene(mat)),show.node.label = TRUE)
plot(njs(dist.gene(mat)),show.node.label = TRUE)
