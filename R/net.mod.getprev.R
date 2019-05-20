
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
#' in the model at any given timestep. For a two-group network, the current state
#' size for each group and overall is provided.
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
get_prev.net.grp <- function(dat, at) {

  active <- dat$attr$active

  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL

  status <- l$status
  group <- idgroup(dat$nw)[active == 1]

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
    dat$epi$s.num <- sum(status == "s" & group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "s" &
                                                     group == 1 &
                                                     get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num <- sum(status == "i" & group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]] <- sum(status == "i" &
                                                     group == 1 &
                                                     get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(status == "r" & group == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]] <- sum(status == "r" &
                                                       group == 1 &
                                                       get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num <- sum(group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]] <- sum(group == 1 &
                                                   get(ebn) == ebv[i])
      }
    }
    dat$epi$s.num.g2 <- sum(status == "s" & group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num.g2", ebun[i])]] <- sum(status == "s" &
                                                        group == 2 &
                                                        get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num.g2 <- sum(status == "i" & group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num.g2", ebun[i])]] <- sum(status == "i" &
                                                        group == 2 &
                                                        get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num.g2 <- sum(status == "r" & group == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num.g2", ebun[i])]] <- sum(status == "r" &
                                                          group == 2 &
                                                          get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num.g2 <- sum(group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num.g2", ebun[i])]] <- sum(group == 2 &
                                                      get(ebn) == ebv[i])
      }
    }
  } else {
    # at > 1
    dat$epi$s.num[at] <- sum(status == "s" & group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "s" &
                                                         group == 1 &
                                                         get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num[at] <- sum(status == "i" & group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num", ebun[i])]][at] <- sum(status == "i" &
                                                         group == 1 &
                                                         get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(status == "r" & group == 1)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("s.num", ebun[i])]][at] <- sum(status == "r" &
                                                           group == 1 &
                                                           get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num[at] <- sum(group == 1)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num", ebun[i])]][at] <- sum(group == 1 &
                                                       get(ebn) == ebv[i])
      }
    }
    dat$epi$s.num.g2[at] <- sum(status == "s" & group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("s.num.g2", ebun[i])]][at] <- sum(status == "s" &
                                                            group == 2 &
                                                            get(ebn) == ebv[i])
      }
    }
    dat$epi$i.num.g2[at] <- sum(status == "i" & group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("i.num.g2", ebun[i])]][at] <- sum(status == "i" &
                                                            group == 2 &
                                                            get(ebn) == ebv[i])
      }
    }
    if (dat$control$type == "SIR") {
      dat$epi$r.num.g2[at] <- sum(status == "r" & group == 2)
      if (eb == TRUE) {
        for (i in 1:length(ebun)) {
          dat$epi[[paste0("r.num.g2", ebun[i])]][at] <- sum(status == "r" &
                                                              group == 2 &
                                                              get(ebn) == ebv[i])
        }
      }
    }
    dat$epi$num.g2[at] <- sum(group == 2)
    if (eb == TRUE) {
      for (i in 1:length(ebun)) {
        dat$epi[[paste0("num.g2", ebun[i])]][at] <- sum(group == 2 &
                                                          get(ebn) == ebv[i])
      }
    }
  }
  return(dat)
}

