
net.update.attributes <- function(all, attribute, type, params, include.inactive=F) {
  
  # Check all
  # Check attr
  
  switch(type, 
         
    ####### if type = uniform #######
    # Params should have a single numeric
    uniform = {
      if (!is.numeric(params) | length(params) != 1) stop('Params argument to function net.update.attributes 
              must be a numeric of length 1 when type=\'fixed\'.')
      switch(include.inactive,
        F= all$attr[attribute][all$attr$active==1] <- all$attr[attribute][all$attr$active==1] + params,
        T= all$attr[attribute] <- all$attr[attribute] + params
      )  
    },

    
    ####### if type = by #######
    # Params should be a list; [[1]] should be an attribute name, [[2]] should be a set of values for that attr;
    #   [[3]] should be a vector of numerics of the same length as [[2]]
    by = {
      if (!is.list(params) | length(params) != 1) stop('Params argument to function net.update.attributes 
              must be a list of length 3 when type=\'by\'.')
      # More error checking would probably be good.
      #switch(include.inactive,
      #       F= all$attr[attribute][all$attr$active==1] <- all$attr[attribute][all$attr$active==1] + params,
      #       T= all$attr[attribute] <- all$attr[attribute] + params[]
      #)  
    },
    
    ####### if type = transmat #######
    transmat = {
      
    },         
      
    # if type does not match any of the options
    stop('Argument type to function net.update.attributes was not of an acceptable option.')    
  )
  
  return(all)

}
