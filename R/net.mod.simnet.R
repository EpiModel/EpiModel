
#' @title Resimulate Dynamic Network at Time `at`
#'
#' @description This function resimulates the dynamic network in stochastic
#'              network models simulated in \code{\link{netsim}} with dependence
#'              between the epidemic and demographic processes and the network
#'              structure.
#'
#' @inheritParams recovery.net
#' @param nsteps number of time steps to simulate
#'
#' @inherit recovery.net return
#'
#' @export
#' @keywords netUtils internal
#'
resim_nets <- function(dat, at, nsteps = 1) {

  # Control settings for resimulation
  tergmLite <- get_control(dat, "tergmLite")
  isTERGM <- get_control(dat, "isTERGM")
  save.nwstats <- get_control(dat, "save.nwstats")
  resimulate.network <- get_control(dat, "resimulate.network")
  nwstats.formula <- get_control(dat, "nwstats.formula")
  set.control.tergm <- get_control(dat, "set.control.tergm")
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")

  # Edges Correction
  dat <- edges_correct(dat, at)
    
  # active attribute (all models)
  active <- get_attr(dat, "active")
  idsActive <- which(active == 1)
  anyActive <- ifelse(length(idsActive) > 0, TRUE, FALSE)
    
  # group attribute (built-in models)
  if (dat$param$groups == 2) {
    group <- get_attr(dat, "group")
    groupids.1 <- which(group == 1)
    groupids.2 <- which(group == 2)
    nActiveG1 <- length(intersect(groupids.1, idsActive))
    nActiveG2 <- length(intersect(groupids.2, idsActive))
    anyActive <- ifelse(nActiveG1 > 0 & nActiveG2 > 0, TRUE, FALSE)
  }

  # Pull network model parameters
  nwparam <- get_nwparam(dat)

  # Network resimulation
  if (anyActive == TRUE & resimulate.network == TRUE) {

    if (tergmLite == FALSE) {
      # Full tergm/network Method
      nw <- dat$nw[[1]]
      output <- "networkDynamic"
    } else {
      # tergmLite/networkLite Method
      nw <- networkLite(dat$el[[1]], dat$attr)
      output <- "final"
      if (tergmLite.track.duration == TRUE) {
        nw %n% "time" <- dat$nw[[1]] %n% "time"
        nw %n% "lasttoggle" <- dat$nw[[1]] %n% "lasttoggle"
      }
    }

    # TERGM simulation
    if (isTERGM == TRUE) {
      dat$nw[[1]] <- simulate(nw ~
                                Form(nwparam$formation) +
                                Persist(nwparam$coef.diss$dissolution),
                              coef = c(nwparam$coef.form,
                                       nwparam$coef.diss$coef.adj),
                              constraints = nwparam$constraints,
                              time.start = at - 1,
                              time.slices = nsteps,
                              output = output,
                              monitor = nwstats.formula,
                              control = set.control.tergm,
                              dynamic = TRUE)
      
    } else {
      dat$nw[[1]] <- simulate(nwparam$formation,
                              basis = nw,
                              coef = c(nwparam$coef.form),
                              constraints = nwparam$constraints,
                              time.start = at - 1,
                              time.slices = nsteps,
                              output = output,
                              monitor = nwstats.formula,
                              control = set.control.tergm,
                              dynamic = TRUE)
    }

    # Update nwstats data frame
    if (save.nwstats == TRUE) {
      new.nwstats <- attributes(dat$nw[[1]])$stats
      keep.cols <- which(!duplicated(colnames(new.nwstats)))
      new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
      dat$stats$nwstats[[1]] <- rbind(dat$stats$nwstats[[1]], new.nwstats)
    }
    
    if (tergmLite == TRUE) {
      dat$el[[1]] <- as.edgelist(dat$nw[[1]])
    }
  }

  return(dat)
}


#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficient in a dynamic network model
#'              simulated in \code{\link{netsim}} to preserve the mean
#'              degree of nodes in the network.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @keywords internal
#' @export
#'
edges_correct <- function(dat, at) {

  resimulate.network <- get_control(dat, "resimulate.network")
  groups <- get_param(dat, "groups")
  active <- get_attr(dat, "active")

  if (at > 1 && resimulate.network == TRUE) {

    if (groups == 1) {
      index <- at - 1
      old.num <- get_epi(dat, "num", index)
      new.num <- sum(active == 1)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(old.num) -
        log(new.num)
    }
    if (groups == 2) {
      index <- at - 1
      group <- get_attr(dat, "group")
      old.num.g1 <- get_epi(dat, "num", index)
      old.num.g2 <- get_epi(dat, "num.g2", index)
      new.num.g1 <- sum(active == 1 & group == 1)
      new.num.g2 <- sum(active == 1 & group == 2)
      dat$nwparam[[1]]$coef.form[1] <- dat$nwparam[[1]]$coef.form[1] +
        log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
        log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
    }
  }
  return(dat)
}
