
#' @title Simulate Dynamic Network at Time 1
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps, to be used in \code{\link{netsim}} models.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param nw A \code{networkDynamic} object.
#' @param nsteps Number of time steps to simulate the network over.
#' @param nsims Number of independent network simulations.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets <- function(x, nw, nsteps, control) {

  suppressWarnings(
    sim <- simulate(nw,
                    formation = x$formation,
                    dissolution = x$coef.diss$dissolution,
                    coef.form = x$coef.form,
                    coef.diss = x$coef.diss$coef.crude,
                    time.slices = nsteps,
                    time.start = 1,
                    time.offset = 0,
                    constraints = x$constraints,
                    monitor = control$nwstats.formula,
                    nsim = 1,
                    control = control$set.control.stergm))

  return(sim)
}


#' @title Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @param x A master object passed through \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at) {

  idsActive <- which(dat$attr$active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)

  if (dat$param$groups == 2) {
    groupids.1 <- which(get.vertex.attribute(dat$nw, "group") == 1)
    groupids.2 <- which(get.vertex.attribute(dat$nw, "group") == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  # Pull network model parameters
  nwparam <- get_nwparam(dat)

  # Serosorting model check
  statOnNw <- ("status" %in% dat$temp$fterms)
  status <- dat$attr$status
  if (statOnNw == TRUE && length(unique(status)) == 1) {
    stop("Stopping simulation because status in formation formula and ",
         "no longer any discordant nodes",
         call. = TRUE)
  }

  # Set up nwstats df
  if (dat$control$save.nwstats == TRUE) {
    if (at == 2) {
      nwstats <- attributes(dat$nw)$stats
      dat$stats$nwstats <- as.data.frame(nwstats)
    }
  }

  # Network simulation
  if (anyActive > 0 & dat$control$depend == TRUE) {
    suppressWarnings(
      dat$nw <- simulate(dat$nw,
                         formation = nwparam$formation,
                         dissolution = nwparam$coef.diss$dissolution,
                         coef.form = nwparam$coef.form,
                         coef.diss = nwparam$coef.diss$coef.adj,
                         constraints = nwparam$constraints,
                         time.start = at,
                         time.slices = 1,
                         time.offset = 0,
                         monitor = dat$control$nwstats.formula,
                         control = dat$control$set.control.stergm))

    # Set up nwstats df
    if (dat$control$save.nwstats == TRUE) {
      dat$stats$nwstats <- rbind(dat$stats$nwstats,
                                 tail(attributes(dat$nw)$stats, 1)[,])
    }

    if (dat$control$delete.nodes == TRUE) {
      dat$nw <- network.extract(dat$nw, at = at)
      inactive <- which(dat$attr$active == 0) #Merge
      dat$nw.update$resim$inactive <- inactive
    }

    #Edges Correction
    if (dat$param$groups == 1) {
      old.num <- dat$epi$num[at - 1]
      new.num <- sum(dat$attr$active == 1)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(old.num) -
        log(new.num)
    }
    if (dat$param$groups == 2) {
      group <- idgroup(dat$nw)
      old.num.g1 <- dat$epi$num[at - 1]
      old.num.g2 <- dat$epi$num.g2[at - 1]
      new.num.g1 <- sum(dat$attr$active == 1 & group == 1)
      new.num.g2 <- sum(dat$attr$active == 1 & group == 2)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
        log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
    }
  }
  return(dat)
}

#' @title Simulate Dynamic Network at Time 1
#'
#' @description This function simulates a dynamic network over one or multiple
#'              time steps, to be used in \code{\link{netsim}} models.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param nw A \code{networkDynamic} object.
#' @param nsteps Number of time steps to simulate the network over.
#' @param nsims Number of independent network simulations.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#'
#' @export
#' @keywords netUtils internal
#'
sim_nets <- function(x, nw, nsteps, control) {

  suppressWarnings(
    sim <- simulate(nw,
                    formation = x$formation,
                    dissolution = x$coef.diss$dissolution,
                    coef.form = x$coef.form,
                    coef.diss = x$coef.diss$coef.crude,
                    time.slices = nsteps,
                    time.start = 1,
                    time.offset = 0,
                    constraints = x$constraints,
                    monitor = control$nwstats.formula,
                    nsim = 1,
                    control = control$set.control.stergm))

  return(sim)
}


#' @title TergmLite: Resimulate Dynamic Network at Time 2+
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @param x A master object passed through \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets.tgl <- function(dat, at) {

  idsActive <- which(dat$attr$active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
  if (dat$param$groups == 2) {
    groupids.1 <- which(dat$attr$groups == 1)
    groupids.2 <- which(dat$attr$groups == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  nwparam <- get_nwparam(dat)

  dat$nw$el[[1]] <- tergmLite::simulate_network(p = dat$nw$p[[1]],
                                                el = dat$nw$el[[1]],
                                                coef.form = nwparam$coef.form,
                                                coef.diss = nwparam$coef.diss$coef.adj,
                                                save.changes = TRUE)

  #Edges Correction
  if (dat$param$groups == 1) {
    old.num <- dat$epi$num[at - 1]
    new.num <- sum(dat$attr$active == 1)
    dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
      log(old.num) -
      log(new.num)
  }
  if (dat$param$groups == 2) {
    group <- idgroup(dat$nw)
    old.num.g1 <- dat$epi$num[at - 1]
    old.num.g2 <- dat$epi$num.g2[at - 1]
    new.num.g1 <- sum(dat$attr$active == 1 & group == 1)
    new.num.g2 <- sum(dat$attr$active == 1 & group == 2)
    dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
      log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
      log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
  }


  return(dat)
}


#' @title EpiModel Network Writes
#'
#' @description This function handles all call to network object contained in
#' \code{dat$nw} during simulation, whether a direct or indirect manipulation.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#'
nw.update.net <- function(dat, at) {

  tea.status <- dat$control$tea.status

  #Resimulate Network----

  #Deactive inactive nodes
  inactive <- dat$nw.update$resim$inactive
  if (length(inactive) > 0) {
    dat$attr <- deleteAttr(dat$attr, inactive)
  }

  if (dat$param$vital != FALSE) {

  #Departures----

    idsDpt <- unlist(dat$nw.update$idsDpt)
    idsDpt <- as.vector(idsDpt)
    #Deactive all departures on the network -

    if (length(idsDpt) > 0) {
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDpt, deactivate.edges = TRUE)
    }

  #Arrivals----
    nArrivals <- dat$nw.update$arr$nArrivals
    if (sum(nArrivals) > 0) {
      nCurr <- network.size(dat$nw)
      #New Arrivals
      dat$nw <- add.vertices(dat$nw, nv = sum(nArrivals))
      newNodes <- (nCurr + 1):(nCurr + sum(nArrivals))
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
      if (length(nArrivals) > 1) {
        dat$nw <- set.vertex.attribute(dat$nw, "group",
                                       rep(1:2, c(nArrivals[1], nArrivals[2])),
                                       newNodes)
      }

      # Set attributes on nw
      fterms <- dat$temp$fterms
      curr.tab <- get_attr_prop(dat$nw, fterms)
      if (length(curr.tab) > 0) {
        dat$nw <- update_nwattr(dat$nw, newNodes, dat$control$attr.rules,
                                curr.tab, dat$temp$t1.tab)
      }

      # Save any val on attr
      dat <- copy_toall_attr(dat, at, fterms)

      if (tea.status == TRUE) {
        if ("status" %in% fterms) {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = dat$attr$status[newNodes],
                                              onset = at, terminus = Inf,
                                              v = newNodes)
        } else {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = "s", onset = at, terminus = Inf,
                                              v = newNodes)
        }
      }
      if (!("status" %in% fterms)) {
        dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
      }
      dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
      dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
      dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))

      ## Handles infTime when incoming nodes are infected
      newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
      dat$attr$infTime[newNodesInf] <- at

      if (length(unique(sapply(dat$attr, length))) != 1) {
        stop("Attribute list of unequal length. Check arrivals.net module.")
      }
    }
  }

  #Recovery----

  idsRecov <- dat$nw.update$rec$idsRecov
  recovState <- dat$nw.update$rec$recovState
  status <- dat$attr$status

  if (length(idsRecov) > 0) {
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = recovState, onset = at,
                                          terminus = Inf, v = idsRecov)
    }
  }

  if ("status" %in% dat$temp$fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }

  #Infection----

  #Active and set vertex attribute of infected
  idsNewInf <- dat$nw.update$inf$idsNewInf
  tea.status <- dat$control$tea.status
  if (length(idsNewInf) > 0) {
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
    }

    if ("status" %in% dat$temp$fterms) {
      dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
    }
  }

  # Discordant Edgelist Calculattion----
  dat$temp$del <- discord_edgelist(dat, at)

  #Output-----
  return(dat)
}

#' @title TergmLite: EpiModel Network Writes
#'
#' @description This function handles all call to network object contained in
#' \code{dat$nw} during simulation, whether a direct or indirect manipulation.
#'
#' @param dat Master list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @export
#'
nw.update.net.tgl <- function(dat, at) {

  tea.status <- dat$control$tea.status

  #Resimulate Network----

  #Deactive inactive nodes
  inactive <- dat$nw.update$resim$inactive
  if (length(inactive) > 0) {
    dat$attr <- deleteAttr(dat$attr, inactive)
  }

  if (dat$param$vital != FALSE) {

    #Departures----

    idsDpt <- unlist(dat$nw.update$idsDpt)
    idsDpt <- as.vector(idsDpt)
    #Deactive all departures on the network -

    if (length(idsDpt) > 0) {
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                    v = idsDpt, deactivate.edges = TRUE)
    }

    #Arrivals----
    nArrivals <- dat$nw.update$arr$nArrivals
    if (sum(nArrivals) > 0) {
      nCurr <- network.size(dat$nw)
      #New Arrivals
      dat$nw <- add.vertices(dat$nw, nv = sum(nArrivals))
      newNodes <- (nCurr + 1):(nCurr + sum(nArrivals))
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
      if (length(nArrivals) > 1) {
        dat$nw <- set.vertex.attribute(dat$nw, "group",
                                       rep(1:2, c(nArrivals[1], nArrivals[2])),
                                       newNodes)
      }

      # Set attributes on nw
      fterms <- dat$temp$fterms
      curr.tab <- get_attr_prop(dat$nw, fterms)
      if (length(curr.tab) > 0) {
        dat$nw <- update_nwattr(dat$nw, newNodes, dat$control$attr.rules,
                                curr.tab, dat$temp$t1.tab)
      }

      # Save any val on attr
      dat <- copy_toall_attr(dat, at, fterms)

      if (tea.status == TRUE) {
        if ("status" %in% fterms) {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = dat$attr$status[newNodes],
                                              onset = at, terminus = Inf,
                                              v = newNodes)
        } else {
          dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                              value = "s", onset = at, terminus = Inf,
                                              v = newNodes)
        }
      }
      if (!("status" %in% fterms)) {
        dat$attr$status <- c(dat$attr$status, rep("s", length(newNodes)))
      }
      dat$attr$active <- c(dat$attr$active, rep(1, length(newNodes)))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, length(newNodes)))
      dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, length(newNodes)))
      dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, length(newNodes)))

      ## Handles infTime when incoming nodes are infected
      newNodesInf <- intersect(newNodes, which(dat$attr$status == "i"))
      dat$attr$infTime[newNodesInf] <- at

      if (length(unique(sapply(dat$attr, length))) != 1) {
        stop("Attribute list of unequal length. Check arrivals.net module.")
      }
    }
  }

  #Recovery----

  idsRecov <- dat$nw.update$rec$idsRecov
  recovState <- dat$nw.update$rec$recovState
  status <- dat$attr$status

  if (length(idsRecov) > 0) {
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                          value = recovState, onset = at,
                                          terminus = Inf, v = idsRecov)
    }
  }

  if ("status" %in% dat$temp$fterms) {
    dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  }

  #Infection----

  #Active and set vertex attribute of infected
  idsNewInf <- dat$nw.update$inf$idsNewInf
  tea.status <- dat$control$tea.status
  if (length(idsNewInf) > 0) {
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw,
                                          prefix = "testatus",
                                          value = "i",
                                          onset = at,
                                          terminus = Inf,
                                          v = idsNewInf)
    }

    if ("status" %in% dat$temp$fterms) {
      dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
    }
  }

  # Discordant Edgelist Calculattion----
  dat$temp$del <- discord_edgelist.tgl(dat, at)

  #Output-----
  return(dat)
}
