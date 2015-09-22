# classes for working with dendograms from ape
require(ape)

#' @title functions for converting to and from phylo class
#' @examples 
 require(tsna)
 data(moodyContactSim)
 v10path<-tPath(moodyContactSim,v=10,start=0)
 
 getChildren(10,17,v10path$previous)
 plot(v10path)
 as.phylo.tPath(moodyContactSim,v10path)
 

 
# edgelist version 
as.phylo.tPath<-function(net, path){

  el<-cbind(path$previous,1:length(path$previous))
  # remove the seed row
  el<-el[path$previous!=0,]
  # 'nodes' are all the vertices with decendents
  nodes<-unique(path$previous[path$previous!=0])
  # tips are all the vertices that are not nodes
  tips<-setdiff(1:length(path$previous),nodes)
  
  ids<-c(tips,nodes)
  
  # remap numbering
  el[,1]<-match(el[,1],ids)
  el[,2]<-match(el[,2],ids)
  
  # delete the singles
  
  
  out<-list()
  out[['edge']]<-el
  out[['Nnode']]<-length(nodes)  # number of non-tip nodes
  out[['tip.label']]<-network.vertex.names(net)[tips]
  out[['node.label']]<-network.vertex.names(net)[nodes]
  out[['edge.length']]<-path$tdist
  class(out)<-'phylo'
  return(out)
}

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

as.phylo.transmat<-function(tm){
  el<-cbind(tm$inf,tm$sus)
  #phylo doesn't like "non-splitting" nodes that have no children
  singleKidNodes<-as.numeric(names(which(table(el[,1])==1)))
  if (length(singleKidNodes)>0){
    singleKidsRows<-el[,1]%in%singleKidNodes
    warning('removed row(s) ',paste(which(singleKidsRows)), ' corresponding to single child')
    el<-el[!singleKidsRows,,drop=FALSE]
  }
  
  origNodes<-unique(el[,1])
  origTips<-setdiff(el[,2],origNodes)
  
  phyloTips<-1:length(origTips)
  phyloNodes<-rep(NA,length(origNodes))

  # find roots (infectors that never appear as sus)
  v<-setdiff(unique(el[,1]),unique(el[,2]))
  # figure out the ordering such that the root
  # node will be one larger than the tip nodes
  phyloN<-length(origTips)+1
  while(length(v)>0){
    origIndex<-which(origNodes==v[1])
    phyloNodes[origIndex]<-phyloN
    phyloN<-phyloN+1
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

as.network.transmat<-function(tm){
  ids<-unique(c(tm$sus,tm$inf))
  # remap the ids to a new set
  tm$sus<-match(tm$sus,ids)
  tm$inf<-match(tm$inf,ids)
  net<-network(cbind(tm$inf,tm$sus))
  network.vertex.names(net)<-ids
  return(net)
}


# function to recursively loop through infection tree
# recursive version
getChildren<-function(v,parent,previous){
  
  vKids<-previous[previous[,1]==v,2]
  if(length(vKids)>0){
    # if the child has children, call deepr
    vC<-lapply(vKids,function(vKid){
      getChildren(vKid,parent=parent+1,previous=previous)  # look up my kids, with a new parent
    })
    el<-do.call(rbind,vC)
    el<-rbind(el,cbind(parent,parent+1,-1)) # link my parent to my kids parent
    #TODO: can't figure out how to include a parent with kids in the correct group
    el<-rbind(el,cbind(parent,v,-2))
  } 
  else {
   # if I have no children, just link me to my parent
    el<-cbind(parent,v,-3)
  }
  return(el)
}

# sequential version
getChildren<-function(v,parent,previous){
  vKids<-previous[previous[,1]==v,2]
  newNodes<-max(previous[,2])+1
  # add a link to the root node?
  el<-matrix(c(newNodes[1],v,0),ncol=3,byrow=TRUE)
  while(length(vKids)>0){
    # find which children have children
    vC<-lapply(vKids,function(vKid){
      previous[previous[,1]==vKid,2]  # look up my kids, with a new parent
    })
    # for each kid found...
    for(kid in vC){
      switch (length(kid),
        1 =   el<-rbind(el,cbind(parent,v,-2)) # this is the end of branch, add with parent
      )
    }
    
    
    el<-rbind(el,cbind(parent,parent+1,-1)) # link my parent to my kids parent
    #TODO: can't figure out how to include a parent with kids in the correct group
    el<-rbind(el,cbind(parent,v,-2))
    el<-rbind(el,do.call(rbind,vC))
  } 

  return(el)
}


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

cat('(10,(7,(c,d)));',file='test.tre')





### example code

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