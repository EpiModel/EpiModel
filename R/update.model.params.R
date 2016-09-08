
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
      # concurrent doesn't actually accept any inputs
      inputs <- NULL
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
      
    } else if (term$name=='degree'){
      # ---- DEGREEE -------------------
      # see ergm:::InitErgmTerm.degree
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      
      d <- args[[1]]
      byarg <- args$byarg
      homophily <- args$homophily
      emptynwstats <- NULL
      if (!is.null(byarg)) {
        nodecov <- dat$attr[[byarg]]
        u <- sort(unique(nodecov))
        if (any(is.na(nodecov))) {
          u <- c(u, NA)
        }
        nodecov <- match(nodecov, u)
        if (length(u) == 1) {
          stop("Attribute given to degree() has only one value", 
               call. = FALSE)
        }
      } 
      if (!is.null(byarg) && !homophily) {
        lu <- length(u)
        du <- rbind(rep(d, lu), rep(1:lu, rep(length(d), lu)))
        if (any(du[1, ] == 0)) {
          emptynwstats <- rep(0, ncol(du))
          tmp <- du[2, du[1, ] == 0]
          for (i in 1:length(tmp)) tmp[i] <- sum(nodecov == 
                                                   tmp[i])
          emptynwstats[du[1, ] == 0] <- tmp
        }
      }  else {
        if (any(d == 0)) {
          emptynwstats <- rep(0, length(d))
          emptynwstats[d == 0] <- attr(dat$el,'n') # network size
        }
      } 
      if (is.null(byarg)) {
        if (length(d) == 0) {
          return(NULL)
        }
        inputs <- c(d)
      }  else if (homophily) {
        if (length(d) == 0) {
          return(NULL)
        }
        inputs <- c(d, nodecov)
      }  else {
        if (ncol(du) == 0) {
          return(NULL)
        }
        inputs <- c(as.vector(du), nodecov)
      }
      if (!is.null(emptynwstats)) {
        mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                  length(inputs), inputs)
        mf$terms[[t]]$emptynwstats <- emptynwstats
  
      } else {
        mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                  length(inputs), inputs)
        # belive it is also necessary to update the maxval for this statistic?
        mf$terms[[t]]$maxval<- attr(dat$el,'n') # network size
      }
     
      
    } 
    else if (term$name=='absdiff'){
      # ---- ABSDIFF -------------------
      # see ergm:::InitErgmTerm.absdiff
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      attrname <- args[[1]]
      # get the transformation function
      pow <- args$pow
      nodecov <- dat$attr[[attrname]]
      #TODO: check of pow passed in correctly
      mf$terms[[t]]$inputs <- c(pow, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    } else if (term$name=='nodecov'){
      # ---- NODECOV -------------------
      # see ergm:::InitErgmTerm.nodecov
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      attrname <- args[[1]]
      # get the transformation function
      f <- args$transform
      nodecov <- dat$attr[[attrname]]
      inputs <- f(nodecov)
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    } else if (term$name=='nodemix'){
      # ---- NODEMIX -------------------
      # see ergm:::InitErgmTerm.nodemix
      form <- dat$nwparam[[1]]$formation
      args<-get_formula_term_args_in_formula_env(form,t)
      attrname <- args[[1]]
      nodecov <- dat$attr[[attrname]]
      base <- args$base
      # ASSUMES NETWORK IS NOT BIPARTITE
      u <- sort(unique(nodecov))
      if (any(is.na(nodecov))) {
        u <- c(u, NA)
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      ucount <- sapply(ui, function(x) {
        sum(nodecov == x, na.rm = TRUE)
      })
      uui <- matrix(1:length(ui)^2, length(ui), length(ui))
      urm <- t(sapply(ui, rep, length(ui)))
      ucm <- sapply(ui, rep, length(ui))
      uun <- outer(u, u, paste, sep = ".")
      if (!is.directed(nw)) {
        uui <- uui[upper.tri(uui, diag = TRUE)]
        urm <- urm[upper.tri(urm, diag = TRUE)]
        ucm <- ucm[upper.tri(ucm, diag = TRUE)]
        uun <- uun[upper.tri(uun, diag = TRUE)]
      }
      if (any(NVL(a$base, 0) != 0)) {
        urm <- as.vector(urm)[-a$base]
        ucm <- as.vector(ucm)[-a$base]
        uun <- as.vector(uun)[-a$base]
      }
      inputs <- c(urm, ucm, nodecov)
      attr(inputs, "ParamsBeforeCov") <- 2 * length(uun)
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