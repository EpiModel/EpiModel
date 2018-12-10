
#' @title Get Epidemic Output from netsim Model
#'
#' @description Provides all active model state sizes from the network at the
#'              specified time step, output to a list of vectors.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @details
#' This network utility is used during the \code{\link{netsim}} simulation
#' process to efficiently query the current size of each state or compartment
#' in the model at any given timestep. For a bipartite network, the current state
#' size for each mode, and overall is provided.
#'
#' @export
#' @keywords netUtils internal
#'
get_prev.net <- function(dat, at) {

  active <- dat$attr$active

  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL

  status <- l$status

  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(status == "s")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "s" &
                                                     get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num <- sum(status == "i")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == "i" &
                                                     get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]] <- sum(status == "r" &
                                                       get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num <- length(status)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]] <- sum(get(ebn) == ebv[i])
      }
    }
  } else {
    # at > 1
    dat$epi$s.num[at] <- sum(status == "s")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "s" &
                                                         get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num[at] <- sum(status == "i")
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == "i" &
                                                         get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(status == "r")
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num", ebun[i])]][at] <- sum(status == "r" &
                                                           get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num[at] <- length(status)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]][at] <- sum(get(ebn) == ebv[i])
      }
    }
  }
  return(dat)
}


#' @title Get Epidemic Output from netsim Model
#'
#' @description Provides all active model state sizes from the network at the
#'              specified time step, output to a list of vectors.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @details
#' This network utility is used during the \code{\link{netsim}} simulation
#' process to efficiently query the current  state size for each mode and overall.
#'
#' @export
#' @keywords netUtils internal
#'
get_prev.net.bip <- function(dat, at) {

  active <- dat$attr$active

  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL

  status <- l$status
  mode <- idmode(dat$nw)[active == 1]


  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(status == "s" & mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "s" &
                                                     mode == 1 &
                                                     get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num <- sum(status == "i" & mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == "i" &
                                                     mode == 1 &
                                                     get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(status == "r" & mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "r" &
                                                       mode == 1 &
                                                       get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num <- sum(mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]] <- sum(mode == 1 &
                                                   get(ebn) == ebv[i])
      }
    }
    dat$epi$s.num.m2 <- sum(status == "s" & mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num.m2", ebun[i])]] <- sum(status == "s" &
                                                        mode == 2 &
                                                        get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num.m2 <- sum(status == "i" & mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num.m2", ebun[i])]] <- sum(status == "i" &
                                                        mode == 2 &
                                                        get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num.m2 <- sum(status == "r" & mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num.m2", ebun[i])]] <- sum(status == "r" &
                                                          mode == 2 &
                                                          get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num.m2 <- sum(mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num.m2", ebun[i])]] <- sum(mode == 2 &
                                                      get(ebn) == ebv[i])
      }
    }
  } else {
    # at > 1
    dat$epi$s.num[at] <- sum(status == "s" & mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "s" &
                                                         mode == 1 &
                                                         get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num[at] <- sum(status == "i" & mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == "i" &
                                                         mode == 1 &
                                                         get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(status == "r" & mode == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "r" &
                                                           mode == 1 &
                                                           get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num[at] <- sum(mode == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]][at] <- sum(mode == 1 &
                                                       get(ebn) == ebv[i])
      }
    }
    dat$epi$s.num.m2[at] <- sum(status == "s" & mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num.m2", ebun[i])]][at] <- sum(status == "s" &
                                                            mode == 2 &
                                                            get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num.m2[at] <- sum(status == "i" & mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num.m2", ebun[i])]][at] <- sum(status == "i" &
                                                            mode == 2 &
                                                            get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num.m2[at] <- sum(status == "r" & mode == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num.m2", ebun[i])]][at] <- sum(status == "r" &
                                                              mode == 2 &
                                                              get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num.m2[at] <- sum(mode == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num.m2", ebun[i])]][at] <- sum(mode == 2 &
                                                          get(ebn) == ebv[i])
      }
    }
  }
  return(dat)
}

