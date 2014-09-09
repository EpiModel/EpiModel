
net.update.attribute <- function(all, attribute, type, params, include.inactive=F) {
  
  # TODO Check all
  # TODO Check attr
  
  switch(type, 
         
    ####### if type = uniform #######
    # Params should have a single numeric
    uniform = {
      if (!is.numeric(params) | length(params) != 1) stop('Params argument to function net.update.attributes 
              must be a numeric of length 1 when type=\'fixed\'.')
      switch(1+include.inactive,
        {all$attr[attribute][all$attr$active==1,] <- all$attr[attribute][all$attr$active==1,] + params},
        {all$attr[attribute] <- all$attr[attribute] + params}
      )  
    },

    
    ####### if type = by #######
    # Params should be a list; [[1]] should be an attribute name, [[2]] should be a set of values for that attr;
    #   [[3]] should be a vector of numerics of the same length as [[2]]
    by = {
      if (!is.list(params) | length(params) != 3) stop('Params argument to function net.update.attributes 
              must be a list of length 3 when type=\'by\'.')
      # TODO More error checking.
      switch(1+include.inactive,
         {all$attr[attribute][all$attr$active==1,] <- all$attr[attribute][all$attr$active==1,] + 
            params[[3]][match(all$attr[params[[1]]][all$attr$active==1,],params[[2]])]},
         {all$attr[attribute] <- all$attr[attribute] + 
            params[[3]][match(all$attr[params[[1]]][,],params[[2]])]}
      )  
    },
    
    ####### if type = transmat #######
    # Params should be a list; [[1]] should be a set of values for the main attribute;
    #   [[2]] should be a square matrix whose dimension lengths are the same as the length of [[1]], and
    #   whose rows each sum to 1.
    transmat = {
      if (!is.list(params) | length(params) != 2) stop('Params argument to function net.update.attributes 
              must be a list of length 2 when type=\'by\'.')
      # TODO More error checking.      
      new.vals <- sample(1:3,1,T,my.transmat[1,])
      switch(1+include.inactive,
        {all$attr[attribute][all$attr$active==1,] <- all$attr[attribute][all$attr$active==1,] + 
           },
        {all$attr[attribute] <- all$attr[attribute] + 
           }
      )
    },         
      
    # if type does not match any of the options
    stop("Argument type to function net.update.attributes must be one of:
            \'uniform\', \'by\', \'transmat\'.")    
  )
  
  return(all)

}
