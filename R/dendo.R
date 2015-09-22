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
as.phylo.tPath<-function(net, path){
  
  el<-cbind(path$previous,1:length(path$previous))
  # remove the seed row
  el<-el[path$previous!=0,]
  #phylo doesn't like "non-splitting" nodes that have no children
  singleKidNodes<-as.numeric(names(which(table(el[,1])==1)))
  if (length(singleKidNodes)>0){
    singleKidsRows<-which(el[,1]==singleKidNodes)
    warning('removed row(s) ',singleKidsRows, ' corresponding to single child')
    el<-el[-singleKidsRows,,drop=FALSE]
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
    phyloNodes[origIndex]<-phyloN
    phyloN<-phyloN+1
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
  out[['tip.label']]<-network.vertex.names(net)[origTips]
  out[['node.label']]<-network.vertex.names(net)[origNodes]
  class(out)<-'phylo'
  return(out)
}

#' @title convert transmat infection tree into a phylo tree object
#' @method as.phylo transmat
#' @export as.phylo.transmat
as.phylo.transmat<-function(tm){
  el<-cbind(tm$inf,tm$sus)
  #phylo doesn't like "non-splitting" nodes that have no children
  singleKidNodes<-as.numeric(names(which(table(el[,1])==1)))
  if (length(singleKidNodes)>0){
    singleKidsRows<-el[,1]%in%singleKidNodes
    warning('found vertices ',paste(singleKidNodes,collapse = ','),' that create "non-splitting" phlo nodes ')
    #el<-el[!singleKidsRows,,drop=FALSE]
  }
  
  origNodes<-unique(el[,1])
  origTips<-setdiff(el[,2],origNodes)
  
  phyloTips<-1:length(origTips)
  phyloNodes<-rep(NA,length(origNodes))

  # find roots (infectors that never appear as sus)
  v<-setdiff(unique(el[,1]),unique(el[,2]))
  if(length(v)>1){
    warning('found multiple trees ',length(v),' not yet supported')
  }
  # figure out the ordering such that the root
  # node will be one larger than the tip nodes
  phyloN<-length(origTips)+1
  while(length(v)>0){
    origIndex<-which(origNodes==v[1])
    if(length(origIndex)>0){
      # add the element on the list of phylo nodes
      phyloNodes[origIndex]<-phyloN
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
  out[['node.label']]<-origNodes
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
  return(net)
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
class(tm)<-c('transmat',class(tm))
tm
}