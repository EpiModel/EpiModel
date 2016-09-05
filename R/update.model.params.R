
# function to appropriately update model params based on ergm model terms
#p <- dat$p
#mf <- p$model.form
#md <- p$model.diss
#mhf <- p$MHproposal.form
#mhd <- p$MHproposal.diss


#TODO: check speed costs of executing the appropriate check.ergmTerm function in order to provide consistent checking of network types, argument names, etc. 

updateModelTermInputs<-function(dat){
  p <- dat$p
  mf <- p$model.form
  md <- p$model.diss
  mhf <- p$MHproposal.form
  mhd <- p$MHproposal.diss
  n <- attributes(dat$el)$n
  combindMaxDyads <- 0
  # we assume that model.form and model.diss have allready been set up in the appropriate structure
  # by ergm.getmodel using the network object, and this has also validated the terms
  # so for most terms we only need to setup the specific values of the input vectors
  # loop over formation model terms and update
  for (t in seq_along(mf$terms)){
    term<-mf$terms[[t]]
  
    if (term$name=='edges'){
# ---- EDGES -----------------------------  
      maxdyads <- choose(n, 2)
      mf$terms[[t]]$maxval <- maxdyads  # TODO: can we pull maxdyads from the term$maxval?
      combindMaxDyads <- combindMaxDyads + maxdyads
    } else if (term$name=='nodematch'){
 # ---- NODEMATCH -------------------
      # see ergm:::InitErgmTerm.nodematch
      # TODO: implement diff, don't match
      # get the name of the attribute to be used for nodecov
      attrname <- strsplit(term$coef.names,'.',fixed=TRUE)[[1]][2]
      # collect the values for the attribute
      nodecov <- dat$attr[[attrname]]
      u <- sort(unique(nodecov))
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      #ui <- seq(along = u)   <-- this is in original code, I don't understand why
      #inputs <- c(ui, nodecov)
      inputs <-nodecov
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    } else if (term$name=='concurrent'){
      # ---- CONCURRENT -------------------
      coef.names <- "concurrent"
      name <- "concurrent"
      inputs <- NULL
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
      
    } else {
      stop('fast_edgelist mode does not know how to update the term ',term$name,' in the formation model formula')
    }
    
  }
  # update combinded maxval
  mf$maxval[1] <- combindMaxDyads
  
  # loop over dissolution model terms and update
  for (t in seq_along(md$terms)){
    term<-md$terms[[t]]
    # ---- EDGES -----------------------------    
    if (term$name=='edges'){
      maxdyads <- choose(n, 2)
      md$terms[[t]]$maxval <- maxdyads  # TODO: can we pull maxdyads from the term$maxval?
      combindMaxDyads <- maxdyads
    } else {
      stop('fast_edgelist mode does not know how to update the term ',term$name,' in the dissolution model formula')
    }
    
  }
  md$maxval <- combindMaxDyads
  
  # update MHproposal.form
  
  # update MHproposal.dis
  # update the elements of the parameter list and return
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)
  dat$p <-p
  return(dat)
}