# Exported Functions ------------------------------------------------------

#' @title Vertex Attributes for Two-Group Network
#'
#' @description Outputs static vertex attributes for a bipartite network for one
#'              specified mode.
#'
#' @param nw An object of class \code{network} or \code{networkDynamic}.
#' @param mode Mode number.
#' @param val Nodal attribute to return.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(n = 10, bipartite = 5)
#' nw <- set.vertex.attribute(nw, "male", rep(0:1, each = 5))
#' bipvals(nw, mode = 1, "male")
#'
bipvals <- function(nw, mode, val) {
  if (!is.numeric(nw$gal$bipartite)) {
    stop("nw must be a bipartite network", call. = FALSE)
  }
  if (missing(mode)) {
    stop("Specify mode=1 or mode=2", call. = FALSE)
  }
  nw %s% modeids(nw, mode) %v% val
}


# Exported Functions ------------------------------------------------------

#' @title Vertex Attributes for Two-Group Network
#'
#' @description Outputs static vertex attributes for a two-group network for one
#'              specified group.
#'
#' @param nw An object of class \code{network} or \code{networkDynamic}.
#' @param group group number.
#' @param val Nodal attribute to return.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(n = 10)
#' nw <- set.vertex.attribute(nw, "male", rep(0:1, each = 5))
#' grpvals(nw, group = 1, "male")
#'
grpvals <- function(nw, group, val) {
  flag <- "group" %in% names(nw$val[[1]])
  if (flag == FALSE) {
    stop("nw must be a two-group network")
  }
  if (missing(group)) {
    stop("Specify group=1 or group=2", call. = FALSE)
  }
  nw %s% groupids(nw, group) %v% val
}


#' @title Calculate Equilibrium for Infection Prevalence
#'
#' @description Calculates the relative change in infection prevalence across a
#'              time series of an epidemic model to assess equilibrium.
#'
#' @param x An \code{EpiModel} object of class \code{dcm}, \code{icm}, or
#'        \code{netsim}.
#' @param numer Numerator for prevalence calculation.
#' @param denom Denominator for prevalence calculation.
#' @param nsteps Number of time steps at end of model simulation to calculate
#'        equilibrium as the absolute value of the difference between the
#'        minimum prevalence value and the maximum prevalence value over that
#'        time range.
#' @param threshold Threshold value for determining equilibrium.
#' @param digits Number of digits to round for table output.
#' @param invisible If \code{TRUE}, function will suppress output to console and
#'        return summary statistics invisibly.
#'
#' @details
#' This function calculates whether equilibrium in disease prevalence, or any
#' other fraction of two compartments contained in an epidemic model, have
#' reached equilibrium over a time series. Equilibrium is calculated as the
#' absolute value of the difference of the maximum prevalence and minimum
#' prevalence over a specified time series. That time series is specified as the
#' final \code{nsteps} time steps of an epidemic model. A larger \code{nsteps}
#' specification will therefore calculate differences over a longer time series.
#'
#' @export
#'
#' @examples
#' # Calculate equilibrium for a DCM
#' param <- param.dcm(inf.prob = 0.2, inf.prob.g2 = 0.1, act.rate = 0.5,
#'                    balance = "g1", rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
#'                    a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
#'                    ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
#'                    di.rate.g2 = 1 / 90)
#' init <- init.dcm(s.num = 500, i.num = 1,
#'                  s.num.g2 = 500, i.num.g2 = 1)
#' control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
#' x <- dcm(param, init, control)
#' plot(x)
#'
#' # Different calculation options
#' calc_eql(x, nsteps = 100)
#' calc_eql(x, nsteps = 250)
#' calc_eql(x, nsteps = 100, numer = "i.num.g2", denom = "num.g2")
#' calc_eql(x, nsteps = 100, numer = "i.num.g2", denom = "num.g2",
#'          threshold = 0.00001)
#'
calc_eql <- function(x, numer = "i.num", denom = "num",
                     nsteps, threshold = 0.001, digits = 4, invisible = FALSE) {
  if (!(class(x) %in% c("dcm", "icm", "netsim"))) {
    stop("x must an object of class dcm, icm, or netsim", call. = FALSE)
  }
  # Change scipen based on digits
  old.scipen <- options()$scipen
  options(scipen = digits + 1)
  # Convert model to df and calculate prevalence
  df <- as.data.frame(x)
  if (!(numer %in% names(df))) {
    stop("numer must be an output compartment on x", call. = FALSE)
  }
  if (!(denom %in% names(df))) {
    stop("denom must be an output compartment on x", call. = FALSE)
  }
  prev <- df[[numer]]/df[[denom]]

  # Truncate vector and calculate difference
  tprev <- tail(prev, nsteps)
  diff <- abs(max(tprev) - min(tprev))
  if (invisible == FALSE) {
    cat("Equilibrium Results")
    cat("\n=================================")
    cat("\nStart Prev.:  ", round(head(tprev, 1), digits))
    cat("\nEnd Prev.:    ", round(tail(tprev, 1), digits))
    cat("\nMax Prev.:    ", round(max(tprev), digits))
    cat("\nMin Prev.:    ", round(min(tprev), digits))
    cat("\nRel. Diff.:   ", round(diff, digits))
    cat("\n<= Threshold: ", diff <= threshold)
  }
  on.exit(options(scipen = old.scipen))
  out <- list(strprev = round(head(tprev, 1), digits),
              endprev = round(tail(tprev, 1), digits),
              maxprev = round(max(tprev), digits),
              minprev = round(min(tprev), digits),
              reldiff = round(diff, digits),
              thresh = diff <= threshold)
  invisible(out)
}


#' @title Check Degree Distribution for Bipartite Target Statistics
#'
#' @description Checks for consistency in the implied network statistics
#'              of a bipartite network in which the mode size and mode-specific
#'              degree distributions are specified.
#'
#' @param num Number of nodes in mode 1.
#' @param num.g2 Number of nodes in mode 2.
#' @param deg.dist.g1 Vector with fractional degree distribution for mode 1.
#' @param deg.dist.g2 Vector with fractional degree distribution for mode 2.
#'
#' @details
#' This function outputs the number of nodes of degree 0 to m, where m is the
#' length of a fractional degree distribution vector, given that vector and the
#' size of the mode. This utility is used to check for balance in implied degree
#' given that fractional distribution within bipartite network simulations, in
#' which the degree-constrained counts must be equal across modes.
#'
#' @seealso
#' For a detailed explanation of this function, see the tutorial:
#' \href{http://statnet.github.io/tut/NetUtils.html}{EpiModel Network
#' Utility Functions}.
#'
#' @export
#' @keywords netUtils
#'
#' @examples
#' # An imbalanced distribution
#' check_bip_degdist(num.g1 = 500, num.g2 = 500,
#'                   deg.dist.g2 = c(0.40, 0.55, 0.03, 0.02),
#'                   deg.dist.g1 = c(0.48, 0.41, 0.08, 0.03))
#'
#' # A balanced distribution
#' check_bip_degdist(num.g1 = 500, num.g2 = 500,
#'                   deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
#'                   deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03))
#'
check_bip_degdist <- function(num.g1, num.g2,
                              deg.dist.g1, deg.dist.g2) {
  deg.counts.g1 <- deg.dist.g1 * num.g1
  deg.counts.g2 <- deg.dist.g2 * num.g2
  tot.deg.g1 <- sum(deg.counts.g1 * (1:length(deg.dist.g1) - 1))
  tot.deg.g2 <- sum(deg.counts.g2 * (1:length(deg.dist.g2) - 1))
  mat <- matrix(c(deg.dist.g1, deg.counts.g1,
                  deg.dist.g2, deg.counts.g2), ncol = 4)
  mat <- rbind(mat, c(sum(deg.dist.g1), tot.deg.g1, sum(deg.dist.g2), tot.deg.g2))
  colnames(mat) <- c("m1.dist", "m1.cnt", "m2.dist", "m2.cnt")
  rownames(mat) <- c(paste0("Deg", 0:(length(deg.dist.g1) - 1)), "Edges")
  cat("Bipartite Degree Distribution Check\n")
  cat("=============================================\n")
  print(mat, print.gap = 3)
  cat("=============================================\n")
  reldiff <- (tot.deg.g1 - tot.deg.g2) / tot.deg.g2
  absdiff <- abs(tot.deg.g1 - tot.deg.g2)
  if (sum(deg.dist.g1) <= 0.999 | sum(deg.dist.g1) >= 1.001 |
      sum(deg.dist.g2) <= 0.999 | sum(deg.dist.g2) >= 1.001 | absdiff > 1) {
    if (sum(deg.dist.g1) <= 0.999 | sum(deg.dist.g1) >= 1.001) {
      cat("** deg.dist.g1 TOTAL != 1 \n")
    }
    if (sum(deg.dist.g2) <= 0.999 | sum(deg.dist.g2) >= 1.001) {
      cat("** deg.dist.g2 TOTAL != 1 \n")
    }
    if (absdiff > 1) {
      if (tot.deg.g1 > tot.deg.g2) {
        msg <- "Mode 1 Edges > Mode 2 Edges:"
      } else {
        msg <- "Mode 1 Edges < Mode 2 Edges:"
      }
      cat("**", msg, round(reldiff, 3), "Rel Diff \n")
    }
  } else {
    cat("** Edges balanced ** \n")
  }
  invisible(c(tot.deg.g1, deg.counts.g1, deg.counts.g2))
}


#' @title Creates a TEA Variable for Infection Status for \code{ndtv} Animations
#'
#' @description Creates a new color-named temporally-extended attribute (TEA)
#'              variable in a \code{networkDynamic} object containing a disease
#'              status TEA in numeric format.
#'
#' @param nd An object of class \code{networkDynamic}.
#' @param old.var Old TEA variable name.
#' @param old.sus Status value for susceptible in old TEA variable.
#' @param old.inf Status value for infected in old TEA variable.
#' @param old.rec Status value for recovered in old TEA variable.
#' @param new.var New TEA variable name to be stored in \code{networkDynamic}
#'        object.
#' @param new.sus Status value for susceptible in new TEA variable.
#' @param new.inf Status value for infected in new TEA variable.
#' @param new.rec Status value for recovered in new TEA variable.
#' @param verbose Print progress to console.
#'
#' @details
#' The \code{ndtv} package (\url{https://cran.r-project.org/package=ndtv}) produces
#' animated visuals for dynamic networks with evolving edge structures and nodal
#' attributes. Nodal attribute dynamics in \code{ndtv} movies require a temporally
#' extended attribute (TEA) containing a standard R color for each node at each
#' time step. By default, the \code{EpiModel} package uses TEAs to store disease
#' status history in network model simulations run in \code{\link{netsim}}. But,
#' that status TEA is in numeric format (0, 1, 2). The \code{color_tea} function
#' transforms those numeric values of that disease status TEA into a TEA with
#' color values in order to visualize status changes in \code{ndtv}.
#'
#' The convention in \code{\link{plot.netsim}} is to color the susceptible
#' nodes as blue, infected nodes as red, and recovered nodes as green. Alternate
#' colors may be specified using the \code{new.sus}, \code{new.inf}, and
#' \code{new.rec} parameters, respectively.
#'
#' Using the \code{color_tea} function with a \code{netsim} object requires that
#' TEAs for disease status be used and that the \code{networkDynamic} object be
#' saved in the output: both \code{tea.status} and \code{save.network} must be
#' set to \code{TRUE} in \code{\link{control.net}}.
#'
#' @seealso \code{\link{netsim}} and the \code{ndtv} package documentation.
#' @keywords colorUtils
#' @export
#'
color_tea <- function(nd, old.var = "testatus", old.sus = "s", old.inf = "i",
                      old.rec = "r", new.var = "ndtvcol", new.sus, new.inf,
                      new.rec, verbose = TRUE) {
  if (missing(new.inf)) {
    new.inf <- transco("firebrick", 0.75)
  }
  if (missing(new.sus)) {
    new.sus <- transco("steelblue", 0.75)
  }
  if (missing(new.rec)) {
    new.rec <- transco("seagreen", 0.75)
  }
  times <- 1:max(get.change.times(nd))
  for (at in times) {
    stat <- get.vertex.attribute.active(nd, old.var, at = at)
    infected <- which(stat == old.inf)
    uninfected <- which(stat == old.sus)
    recovered <- which(stat == old.rec)
    nd <- activate.vertex.attribute(nd, prefix = new.var, value = new.inf,
                                    onset = at, terminus = Inf, v = infected)
    nd <- activate.vertex.attribute(nd, prefix = new.var, value = new.sus,
                                    onset = at, terminus = Inf, v = uninfected)
    nd <- activate.vertex.attribute(nd, prefix = new.var, value = new.rec,
                                    onset = at, terminus = Inf, v = recovered)
    if (verbose == TRUE) {
      cat("\n", at, "/", max(times), "\t", sep = "")
    }
  }
  return(nd)
}


#' @title Copies Vertex Attributes in Formation Formula to attr List
#'
#' @description Copies the vertex attributes stored on the network object to the
#'              master attr list in the dat data object.
#'
#' @param dat Master data object passed through \code{netsim} simulations.
#' @param at Current time step.
#' @param fterms Vector of attributes used in formation formula, usually as
#'        output of \code{\link{get_formula_term_attr}}.
#'
#' @seealso \code{\link{get_formula_term_attr}}, \code{\link{get_attr_prop}},
#'          \code{\link{update_nwattr}}.
#' @keywords netUtils internal
#' @export
#'
copy_toall_attr <- function(dat, at, fterms) {
  otha <- names(dat$nw$val[[1]])
  otha <- otha[which(otha %in% fterms)]
  if (length(otha) > 0) {
    for (i in seq_along(otha)) {
      va <- get.vertex.attribute(dat$nw, otha[i])
      dat$attr[[otha[i]]] <- va
      if (at == 1) {
        if (!is.null(dat$control$epi.by) && dat$control$epi.by == otha[i]) {
          dat$temp$epi.by.vals <- unique(va)
        }
      }
    }
  }
  return(dat)
}


#' @title Dissolution Coefficients for Stochastic Network Models
#'
#' @description Calculates dissolution coefficients, given a dissolution model
#'              and average edge duration, to pass as offsets to an ERGM/STERGM
#'              model fit in \code{netest}.
#'
#' @param dissolution Right-hand sided STERGM dissolution formula
#'        (see \code{\link{netest}}). See below for list of supported dissolution
#'        models.
#' @param duration A vector of mean edge durations in arbitrary time units.
#' @param d.rate Departure or exit rate from the population, as a single homogenous
#'        rate that applies to the entire population.
#'
#' @details
#' This function performs two calculations for dissolution coefficients
#' used in a network model estimated with \code{\link{netest}}:
#' \enumerate{
#'  \item \strong{Transformation:} the mean duration of edges in a network are
#'        mathematically transformed to logit coefficients.
#'  \item \strong{Adjustment:} in a dynamic network simulation in an open
#'        population (in which there are departures), it is further necessary to
#'        adjust these coefficients for dynamic simulations; this upward adjustment
#'        accounts for departure as a competing risk to edge dissolution.
#' }
#'
#' The current dissolution models supported by this function and in network
#' model estimation in \code{\link{netest}} are as follows:
#' \itemize{
#'  \item \code{~offset(edges)}: a homogeneous dissolution model in which the
#'         edge duration is the same for all partnerships. This requires
#'         specifying one duration value.
#'  \item \code{~offset(edges) + offset(nodematch("<attr>"))}: a heterogeneous
#'         model in which the edge duration varies by whether the nodes in the
#'         dyad have similar values of a specified attribute. The duration vector
#'         should now contain two values: the first is the mean edge duration of
#'         non-matched dyads, and the second is the duration of the matched dyads.
#'  \item \code{~offset(edges) + offset(nodemix("<attr>"))}: a heterogenous model
#'         that extends the nodematch model to include non-binary attributes for
#'         homophily. The duration vector should first contain the base value,
#'         then the values for every other possible combination in the term.
#'  \item \code{~offset(edges) + offset(nodefactor("<attr>"))}: a heterogenous
#'         model in which the edge duration varies by a specified attribute. The
#'         duration vector should first contain the base value, then the values
#'         for every other value of that attribute in the term.
#' }
#'
#' @return
#' A list of class \code{disscoef} with the following elements:
#' \itemize{
#'  \item \strong{dissolution:} right-hand sided STERGM dissolution formula
#'         passed in the function call.
#'  \item \strong{duration:} mean edge durations passed into the function.
#'  \item \strong{coef.crude:} mean durations transformed into logit
#'        coefficients.
#'  \item \strong{coef.adj:} crude coefficients adjusted for the risk of
#'        departure on edge persistence, if the \code{d.rate} argument is supplied.
#'  \item \strong{d.rate:} the departure rate.
#' }
#'
#' @seealso
#' The theory and details of this function are explained in detail in the
#' \href{http://statnet.github.io/tut/NetUtils.html}{EpiModel Network
#' Utility Functions} tutorial.
#'
#' @export
#' @keywords netUtils
#'
#' @examples
#' # Homogeneous dissolution model with no departures
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#'
#' # Homogeneous dissolution model with departures
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25,
#'                   d.rate = 0.001)
#'
#' # Heterogeneous dissolution model in which same-race edges have
#' # shorter duration compared to mixed-race edges, with no departures
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10))
#'
#' # Heterogeneous dissolution model in which same-race edges have
#' # shorter duration compared to mixed-race edges, with departures
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10), d.rate = 0.001)
#'
#'
dissolution_coefs <- function(dissolution, duration, d.rate = 0) {
  # Error check for duration < 1
  if (any(duration < 1)) {
    stop("All values in duration must be >= 1", call. = FALSE)
  }
  # Check form of dissolution formula
  form.length <- length(strsplit(as.character(dissolution)[2], "[+]")[[1]])
  t1.edges <- grepl("offset[(]edges",
                    strsplit(as.character(dissolution)[2], "[+]")[[1]][1])
  if (form.length == 2) {
    t2 <- strsplit(as.character(dissolution)[2], "[+]")[[1]][2]
    t2.term <- NULL
    if (grepl("offset[(]nodematch", t2)) {
      t2.term <- "nodematch"
    } else if (grepl("offset[(]nodefactor", t2)) {
      t2.term <- "nodefactor"
    } else if (grepl("offset[(]nodemix", t2)) {
      t2.term <- "nodemix"
    }
  }
  model.type <- NA
  if (form.length == 1 && t1.edges == TRUE) {
    model.type <- "homog"
  } else if (form.length == 2 && t1.edges == TRUE &&
             t2.term %in% c("nodematch", "nodefactor", "nodemix")) {
    model.type <- "hetero"
  } else {
    model.type <- "invalid"
  }
  if (length(d.rate) > 1) {
    stop("Length of d.rate must be 1", call. = FALSE)
  }
  # Log transformation of duration to coefficent
  if (t1.edges == FALSE) {
    stop("Dissolution models must start with offset(edges)", call. = FALSE)
  }
  if (form.length == 1) {
    if (length(duration) > 1) {
      stop("Dissolution model length is 1, but number of duration was ",
           length(duration), call. = FALSE)
    }
    pg <- (duration[1] - 1) / duration[1]
    ps2 <- (1 - d.rate) ^ 2
    coef.crude <- log(pg / (1 - pg))
    if (ps2 <= pg) {
      d.rate_ <- round(1-sqrt(pg),5)
      str <- paste("The competing risk of departure is too high for the given",
                   " duration of ", duration[1], "; specify a d.rate lower than ",
                   d.rate_,".",sep="")
      stop(str, call. = FALSE)
    }
    coef.adj <- log(pg / (ps2 - pg))
  }
  if (form.length == 2) {
    if (t2.term %in% c("nodematch", "nodefactor", "nodemix")) {
      coef.crude <- coef.adj <- NA
      for (i in 1:length(duration)) {
        pg.thetaX <- (duration[i] - 1) / duration[i]
        ps2.thetaX <- (1 - d.rate) ^ 2
        if (sqrt(ps2.thetaX) <= pg.thetaX) {
          stop("The competing risk of departure is too high for the given the ",
               "duration in place ", i, ". Specify a lower d.rate", call. = FALSE)
        }
        if (i == 1) {
          coef.crude[i] <- log(pg.thetaX / (1 - pg.thetaX))
          coef.adj[i] <- log(pg.thetaX / (ps2.thetaX - pg.thetaX))
        } else {
          coef.crude[i] <- log(pg.thetaX / (1 - pg.thetaX)) - coef.crude[1]
          coef.adj[i] <- log(pg.thetaX / (ps2.thetaX - pg.thetaX)) - coef.adj[1]
        }
      }
    } else {
      stop("Supported heterogeneous dissolution model terms are nodematch, ",
           "nodefactor, or nodemix", call. = FALSE)
    }
  }
  out <- list()
  out$dissolution <- dissolution
  out$duration <- duration
  out$coef.crude <- coef.crude
  out$coef.adj <- coef.adj
  out$d.rate <- d.rate
  out$model.type <- model.type
  class(out) <- "disscoef"
  return(out)
}


#' @title Table of Edge Censoring
#'
#' @description Outputs a table of the number and percent of edges that are
#'              left-censored, right-censored, both-censored, or uncensored for a
#'              \code{networkDynamic} object.
#'
#' @param el Timed edgelist with start and end times extracted from a
#'        \code{networkDynamic} object using the \code{as.data.frame.networkDynamic}
#'        function.
#'
#' @export
#' @keywords netUtils
#'
#' @details
#' Given a STERGM simulation over a specified number of time steps, the edges
#' within that simulation may be left-censored (started before the first step),
#' right-censored (continued after the last step), right and left-censored, or
#' uncensored. The amount of censoring will increase when the average edge
#' duration approaches the length of the simulation.
#'
#' @examples
#' # Initialize and parameterize network model
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the network and extract a timed edgelist
#' dx <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE, verbose = FALSE)
#' el <- as.data.frame(dx)
#'
#' # Calculate censoring
#' edgelist_censor(el)
#'
edgelist_censor <- function(el) {
  # left censored
  leftcens <- el$onset.censored
  leftcens.num <- sum(leftcens)
  leftcens.pct <- leftcens.num / nrow(el)
  # right censored
  rightcens <- el$terminus.censored
  rightcens.num <- sum(rightcens)
  rightcens.pct <- rightcens.num / nrow(el)
  # partnership lasts for entire window (left and right censored)
  lrcens <- el$onset.censored & el$terminus.censored
  lrcens.num <- sum(lrcens)
  lrcens.pct <- lrcens.num / nrow(el)
  # fully observed
  nocens <- el$onset.censored == FALSE & el$terminus.censored == FALSE
  nocens.num <- sum(nocens)
  nocens.pct <- nocens.num / nrow(el)
  ## Table
  nums <- rbind(leftcens.num, rightcens.num, lrcens.num, nocens.num)
  pcts <- rbind(leftcens.pct, rightcens.pct, lrcens.pct, nocens.pct)
  out <- cbind(nums, pcts)
  rownames(out) <- c("Left Cens.", "Right Cens.", "Both Cens.", "No Cens.")
  colnames(out) <- c("num", "pct")
  return(out)
}


#' @title Mean Age of Partnerships over Time
#'
#' @description Outputs a vector of mean ages of edges at a series of timesteps
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param el If not passing \code{x}, a timed edgelist from a \code{networkDynamic}
#'        object extracted with the \code{as.data.frame.networkDynamic} function.
#'
#' @details
#' This function calculates the mean partnership age at each time step over
#' a dynamic network simulation from \code{\link{netest}}. These objects
#' contain the network, edgelist, and dissolution objects needed for the
#' calculation. Alternatively, one may pass in these objects separately if
#' \code{netest} was not used, or statistics were not run requested after
#' the estimation.
#'
#' Currently, the calculations are limited to those dissolution formulas with a single
#' homogenous dissolution (\code{~offset(edges)}). This functionality will be
#' expanded in future releases.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' # Initialize and parameterize the network model
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the network and extract a timed edgelist
#' dx <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE, verbose = FALSE)
#' el <- as.data.frame(dx)
#'
#' # Calculate ages directly from edgelist
#' mean_ages <- edgelist_meanage(el = el)
#' mean_ages
#'
#' # Alternatively, netdx calculates these
#' dx$pages
#' identical(dx$pages[[1]], mean_ages)
#'
edgelist_meanage <- function(x, el) {
  # If passing a netest object directly
  if (!(missing(x))) {
    el <- x$edgelist
  }
  terminus <- el$terminus
  onset <- el$onset
  minterm <- min(terminus)
  maxterm <- max(terminus)
  meanpage <- rep(NA, maxterm)
  for (at in minterm:maxterm) {
    actp <- (onset <= at & terminus > at) |
      (onset == at & terminus == at);
    page <- at - onset[actp] + 1
    meanpage[at] <- mean(page)
  }
  meanpage <- meanpage[1:(length(meanpage) - 1)]
  return(meanpage)
}


#' @title Proportional Table of Vertex Attributes
#'
#' @description Calculates the proportional distribution of each vertex attribute
#'              contained on the network, with a possible limitation to those
#'              attributes contained in the formation formula only.
#'
#' @param nw The \code{networkDynamic} object contained in the \code{netsim}
#'        simulation.
#' @param fterms Vector of attributes used in formation formula, usually as
#'        output of \code{\link{get_formula_term_attr}}.
#' @param only.formula Limit the tables to those terms only in \code{fterms},
#'        otherwise output proportions for all attributes on the network object.
#'
#' @seealso \code{\link{get_formula_term_attr}}, \code{\link{copy_toall_attr}},
#'          \code{\link{update_nwattr}}.
#' @keywords netUtils internal
#' @export
#'
get_attr_prop <- function(nw, fterms, only.formula = TRUE) {
  if (is.null(fterms)) {
    return(NULL)
  }
  nwVal <- names(nw$val[[1]])
  if (only.formula == TRUE) {
    nwVal <- nwVal[which(nwVal %in% fterms)]
  }
  out <- list()
  for (i in 1:length(nwVal)) {
    tab <- prop.table(table(nw %v% nwVal[i]))
    out[[i]] <- tab
  }
  names(out) <- nwVal
  return(out)
}


#' @title Outputs ERGM Formula Attributes into a Character Vector
#'
#' @description Given a formation formula for a network model, outputs it into
#'              a character vector of vertex attributes to be used in \code{netsim}
#'              simulations.
#'
#' @param form an ergm model formula
#' @param nw a network object
#'
#' @export
#'
get_formula_term_attr <- function(form, nw) {

  nw_attr <- names(nw$val[[1]])
  nw_attr <- setdiff(nw_attr, c("active", "vertex.names", "na"))

  if (length(nw_attr) == 0) {
    return(NULL)
  }
  matches <- sapply(nw_attr, function(x) grepl(x, form))
  matches <- colSums(matches)

  out <- names(matches)[which(matches == 1)]
  if (length(out) == 0) {
    return(NULL)
  } else {
    return(out)
  }

}

#' @title Mode Numbers for Bipartite Network
#'
#' @description Outputs mode numbers give ID numbers for a bipartite network.
#'
#' @param nw Object of class \code{network} or \code{networkDynamic}.
#' @param ids Vector of ID numbers for which the mode number
#'        should be returned.
#'
#' @seealso \code{\link{modeids}} provides the reverse functionality.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(10, bipartite = 5)
#' idmode(nw)
#' idmode(nw, ids = c(3, 6))
#'
idmode <- function(nw, ids) {
  n <- network.size(nw)
  if (missing(ids)) {
    ids <- seq_len(n)
  }
  if (any(ids > n)) {
    stop("Specify ids between 1 and ", n)
  }
  if (!is.bipartite(nw)) {
    out <- rep(1, n)
  } else {
    m1size <- nw$gal$bipartite
    modes <- c(rep(1, m1size),
               rep(2, n - m1size))
    out <- modes[ids]
  }
  return(out)
}


#' @title Mode Numbers for Two-Group Network
#'
#' @description Outputs group numbers give ID numbers for a two-group network.
#'
#' @param nw Object of class \code{network} or \code{networkDynamic}.
#' @param ids Vector of ID numbers for which the group number
#'        should be returned.
#'
#' @seealso \code{\link{groupids}} provides the reverse functionality.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(10)
#' nw <- set.vertex.attribute(nw, "group", rep(c(1,2), each = 5))
#' idgroup(nw)
#' idgroup(nw, ids = c(3, 6))
#'
idgroup <- function(nw, ids) {
  n <- network.size(nw)
  if (missing(ids)) {
    ids <- seq_len(n)
  }
  if (any(ids > n)) {
    stop("Specify ids between 1 and ", n)
  }
  flag <- "group" %in% names(nw$val[[1]])
  if (!flag) {
    out <- rep(1, n)
  } else {
    g1size <- length(which(get.vertex.attribute(nw, "group") == 1))
    groups <- c(rep(1, g1size),
                rep(2, n - g1size))
    out <- groups[ids]
  }
  return(out)
}


#' @title ID Numbers for Bipartite Network
#'
#' @description Outputs ID numbers for a mode number for a bipartite network.
#'
#' @param nw Object of class \code{network} or \code{networkDynamic}.
#' @param mode Mode number to return ID numbers for.
#'
#' @seealso \code{\link{idmode}} provides the reverse functionality.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(10, bipartite = 5)
#' modeids(nw, mode = 2)
#'
modeids <- function(nw, mode) {
  if (!is.numeric(nw$gal$bipartite)) {
    stop("nw must be a bipartite network")
  }
  if (missing(mode)) {
    stop("Specify mode=1 or mode=2")
  }
  n <- network.size(nw)
  m1size <- nw$gal$bipartite
  if (mode == 1) {
    out <- 1:m1size
  }
  if (mode == 2) {
    out <- (m1size + 1):n
  }
  return(out)
}


#' @title ID Numbers for Two-Group Network
#'
#' @description Outputs ID numbers for a group number for a two-group network.
#'
#' @param nw Object of class \code{network} or \code{networkDynamic}.
#' @param group Group number to return ID numbers for.
#'
#' @seealso \code{\link{idmode}} provides the reverse functionality.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(10)
#' nw <- set.vertex.attribute(nw, "group", rep(c(1,2), each = 5))
#' groupids(nw, group = 2)
#'
groupids <- function(nw, group) {
  flag <- "group" %in% names(nw$val[[1]])
  if (!flag) {
    stop("nw must be a two-group network")
  }
  if (missing(group)) {
    stop("Specify group=1 or group=2")
  }

  if (group == 1) {
    out <- which(get.vertex.attribute(nw, "group") == 1)
  }
  if (group == 2) {
    out <- which(get.vertex.attribute(nw, "group") == 2)
  }
  return(out)
}


#' @title Update Attribute Values for a Two-Group Network
#'
#' @description Adds new values for attributes in a two-group network in which
#'              there may be arrivals/entries in the first mode, which requires
#'              splitting the attribute vector into two, adding the new values,
#'              and re-concatenating the two updated vectors.
#'
#' @param dat Master data object passed through \code{netsim} simulations.
#' @param var Variable to update.
#' @param val Fixed value to set for all incoming nodes.
#' @param nCurrG1 Number currently in group 1.
#' @param nCurrG2 Number currently in group 2.
#' @param narrivals Number of arrivals/entries in group 1.
#' @param narrivalsG2 Number of arrivals/entries in group 2.
#'
#' @export
#' @keywords netUtils internal
#'
split_bip <- function(dat, var, val, nCurrG1, nCurrG2, nArrivals, nArrivalsG2) {
  VarG1  <- dat$attr[[var]][groupids(dat$nw, group = 1)]
  VarG2  <- dat$attr[[var]][groupids(dat$nw, group = 2)]
  newVarG1 <- as.vector(sapply(VarG1, function(x) ifelse(is.na(x) == TRUE, val, x)))
  newVarG2 <- as.vector(sapply(VarG2, function(x) ifelse(is.na(x) == TRUE, val, x)))
  #FLAG 4/26
  #oldVarG1 <- dat$attr[[var]][1:nCurrG1]
  #oldVarG2 <- dat$attr[[var]][(nCurrG1 + 1):(nCurrG1 + nCurrG2)]
  #newVarG1 <- c(oldVarG1, rep(val, nArrivals))
  #newVarG2 <- c(oldVarG2, rep(val, nArrivalsG2))
  newVar <- c(newVarG1, newVarG2)
  dat$attr[[var]] <- newVar
  return(dat)
}


#' @title Updates Vertex Attributes for Incoming Vertices
#'
#' @description Updates the vertex attributes on a network for new nodes incoming
#'              into that network, based on a set of rules for each attribute
#'              that the user specifies in \code{control.net}.
#'
#' @param nw The \code{networkDynamic} object used in \code{netsim} simulations.
#' @param newNodes Vector of nodal IDs for incoming nodes at the current time
#'        step.
#' @param rules List of rules, one per attribute to be set, governing how to set
#'        the values of each attribute.
#' @param curr.tab Current proportional distribution of all vertex attributes.
#' @param t1.tab Proportional distribution of all vertex attributes at the outset
#'        of the simulation.
#'
#' @seealso \code{\link{copy_toall_attr}}, \code{\link{get_attr_prop}},
#'          \code{\link{update_nwattr}}.
#' @keywords netUtils internal
#' @export
#'
update_nwattr <- function(nw, newNodes, rules, curr.tab, t1.tab) {
  for (i in 1:length(curr.tab)) {
    vname <- names(curr.tab)[i]
    rule <- rules[[vname]]
    if (is.null(rule)) {
      rule <- "current"
    }
    if (rule == "current") {
      vclass <- class(nw %v% vname)
      if (vclass == "character") {
        nattr <- sample(names(curr.tab[[vname]]),
                        size = length(newNodes),
                        replace = TRUE,
                        prob = curr.tab[[vname]])
      } else {
        nattr <- sample(as.numeric(names(curr.tab[[i]])),
                        size = length(newNodes),
                        replace = TRUE,
                        prob = curr.tab[[i]])
      }
    } else if (rule == "t1") {
      vclass <- class(nw %v% vname)
      if (vclass == "character") {
        nattr <- sample(names(t1.tab[[vname]]),
                        size = length(newNodes),
                        replace = TRUE,
                        prob = t1.tab[[vname]])
      } else {
        nattr <- sample(as.numeric(names(t1.tab[[i]])),
                        size = length(newNodes),
                        replace = TRUE,
                        prob = t1.tab[[i]])
      }
    } else {
      nattr <- rep(rules[[vname]], length(newNodes))
    }
    nw <- set.vertex.attribute(nw, attrname = vname,
                               value = nattr, v = newNodes)
  }
  return(nw)
}


#' @title Get Individual Degree from Network or Edgelist
#'
#' @description A fast method for querying the current degree of all individuals
#'              within a network.
#'
#' @param x Either an object of class \code{network} or \code{edgelist} generated
#'        from a network. If \code{x} is an edgelist, then it must contain
#'        an attribute for the total network size, \code{n}.
#'
#' @details
#' Individual-level data on the current degree of nodes within a network is
#' often useful for summary statistics and modeling complex interactions between
#' degree. Given a \code{network} class object, \code{net}, one way to look
#' up the current degree is to get a summary of the ERGM term, \code{sociality},
#' as in: \code{summary(net ~ sociality(base = 0))}. But that is computionally
#' inefficient for a number of reasons. This function provide a fast method
#' for generating the vector of degree using a query of the edgelist. It is
#' even faster if the parameter \code{x} is already transformed as an edgelist.
#'
#' @export
#'
#' @examples
#' nw <- network.initialize(500, directed = FALSE)
#'
#' set.seed(1)
#' fit <- ergm(nw ~ edges, target.stats = 250)
#' sim <- simulate(fit)
#'
#' # Slow ERGM-based method
#' ergm.method <- unname(summary(sim ~ sociality(base = 0)))
#' ergm.method
#'
#' # Fast tabulate method with network object
#' deg.net <- get_degree(sim)
#' deg.net
#'
#' # Even faster if network already transformed into an edgelist
#' el <- as.edgelist(sim)
#' deg.el <- get_degree(el)
#' deg.el
#'
#' identical(as.integer(ergm.method), deg.net, deg.el)
#'
get_degree <- function(x) {
  if (inherits(x, "network")) {
    x <- as.edgelist(x)
  }
  if (is.null(attr(x, "n"))) {
    stop("x missing an n attribute")
  }
  n <- attr(x, "n")
  out <- tabulate(x, nbins = n)
  return(out)
}


#' @title Truncate Simulation Time Series
#'
#' @description Left-truncates a simulation epidemiological summary statistics and
#'              network statistics at a specified time step.
#'
#' @param x Object of class \code{netsim} or \code{icm}.
#' @param at Time step at which to left-truncate the time series.
#'
#' @details
#' This function would be used when running a follow-up simulation from time steps
#' \code{b} to \code{c} after a burnin period from time \code{a} to \code{b},
#' where the final time window of interest for data analysis is \code{b} to \code{c}
#' only.
#'
#' @export
#'
#' @examples
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SI", nsteps = 200, nsims = 1)
#' mod1 <- icm(param, init, control)
#' df <- as.data.frame(mod1)
#' print(df)
#' plot(mod1)
#' mod1$control$nsteps
#'
#' mod2 <- truncate_sim(mod1, at = 150)
#' df2 <- as.data.frame(mod2)
#' print(df2)
#' plot(mod2)
#' mod2$control$nsteps
#'
truncate_sim <- function(x, at) {
  if (class(x) != "icm" && class(x) != "netsim") {
    stop("x must be either an object of class icm or class netsim",
         call. = FALSE)
  }
  rows <- at:(x$control$nsteps)
  # epi
  x$epi <- lapply(x$epi, function(r) r[rows, ])
  # control settings
  x$control$start <- 1
  x$control$nsteps <- max(seq_along(rows))
  return(x)
}
