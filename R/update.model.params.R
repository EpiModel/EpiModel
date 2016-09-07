
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
      # need to get the formation formula to try to parse the params
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      # get the name of the attribute to be used for nodecov
      attrname <- args[[1]]
      # collect the values for the attribute
      nodecov <- dat$attr[[attrname]]
      u <- sort(unique(nodecov))
      # optionally remove values not indicated by 'keep'
      if (!is.null(args$keep)) {
        u <- u[args$keep]
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      ui <- seq(along = u)
      if (args$diff==TRUE) {
        inputs <- c(ui, nodecov)
      } else {
        inputs <- nodecov
      }
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
      
    } else if (term$name=='nodefactor'){
      
      # ---- NODEFACTOR -------------------
      # see ergm:::InitErgmTerm.nodefactor
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      # get the name of the attribute to be used for nodecov
      attrname <- args[[1]]      # collect the values for the attribute
      nodecov <- dat$attr[[attrname]]
      u <- sort(unique(nodecov))
      if (any(NVL(args$base, 0) != 0)) {
        u <- u[-args$base]
        if (length(u) == 0) {
          stop(" nodefactor term should be deleted because it contributes no statistics")
        }
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)
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
      # this is not one of the hardcoded terms, so stop
      stop("EpiModel's fast_edgelist mode does not know how to update the term '",term$name,"' in the formation model formula")
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

# returns a list of the arguments to the terms in the formula
# with offsets removed, evaluated in the formula calling environment
# returns a list where the first element is the term name and subsequent 
# (named) elements are the argument values named by the argument names.
get_formula_term_args_in_formula_env <-function(form,termIndex){
  # get the calling environment of the formula in case
  # there are substitutions
  formula.env<-environment(formula)
  args<-term.list.formula(form[[2]])[[termIndex]]
  # remove the offset term if it exists
  if(args[1]=='offset()'){
    args <- args[[-1]]
  }
  # hack to convert from a call to a list when evaluated
  args[[1]]<- as.name("list")
  # evaluate in formula's calling environment
  outlist <- eval(args,formula.env)
  return(outlist)
}