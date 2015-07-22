
# Exported Functions ------------------------------------------------------

#' @title Vertex Attributes for Bipartite Network
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
#' reach equilibrium over a time series. Equilibrium is calculated as the
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
#'                    b.rate = 1 / 100, b.rate.g2 = NA, ds.rate = 1 / 100,
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
#' @param num.m1 Number of nodes in mode 1.
#' @param num.m2 Number of nodes in mode 2.
#' @param deg.dist.m1 Vector with fractional degree distribution for mode 1.
#' @param deg.dist.m2 Vector with fractional degree distribution for mode 2.
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
#' check_bip_degdist(num.m1 = 500, num.m2 = 500,
#'                   deg.dist.m2 = c(0.40, 0.55, 0.03, 0.02),
#'                   deg.dist.m1 = c(0.48, 0.41, 0.08, 0.03))
#'
#' # A balanced distribution
#' check_bip_degdist(num.m1 = 500, num.m2 = 500,
#'                   deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
#'                   deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))
#'
check_bip_degdist <- function(num.m1, num.m2,
                              deg.dist.m1, deg.dist.m2) {

  deg.counts.m1 <- deg.dist.m1 * num.m1
  deg.counts.m2 <- deg.dist.m2 * num.m2

  tot.deg.m1 <- sum(deg.counts.m1 * (1:length(deg.dist.m1) - 1))
  tot.deg.m2 <- sum(deg.counts.m2 * (1:length(deg.dist.m2) - 1))

  mat <- matrix(c(deg.dist.m1, deg.counts.m1,
                  deg.dist.m2, deg.counts.m2), ncol = 4)
  mat <- rbind(mat, c(sum(deg.dist.m1), tot.deg.m1, sum(deg.dist.m2), tot.deg.m2))

  colnames(mat) <- c("m1.dist", "m1.cnt", "m2.dist", "m2.cnt")
  rownames(mat) <- c(paste0("Deg", 0:(length(deg.dist.m1) - 1)), "Edges")

  cat("Bipartite Degree Distribution Check\n")
  cat("=============================================\n")
  print(mat, print.gap = 3)
  cat("=============================================\n")

  reldiff <- (tot.deg.m1 - tot.deg.m2) / tot.deg.m2
  absdiff <- abs(tot.deg.m1 - tot.deg.m2)

  if (sum(deg.dist.m1) <= 0.999 | sum(deg.dist.m2) <= 0.999 | absdiff > 1) {
    if (sum(deg.dist.m1) <= 0.999) {
      cat("** deg.dist.m1 TOTAL != 1 \n")
    }
    if (sum(deg.dist.m2) <= 0.999) {
      cat("** deg.dist.m2 TOTAL != 1 \n")
    }

    if (absdiff > 1) {
      if (tot.deg.m1 > tot.deg.m2) {
        msg <- "Mode 1 Edges > Mode 2 Edges:"
      } else {
        msg <- "Mode 1 Edges < Mode 2 Edges:"
      }
      cat("**", msg, round(reldiff, 3), "Rel Diff \n")
    }
  } else {
    cat("** Edges balanced ** \n")
  }
  invisible(c(tot.deg.m1, deg.counts.m1, deg.counts.m2))
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
#' The \code{ndtv} package (\url{http://cran.r-project.org/package=ndtv}) produces
#' animated visuals for dynamic networks with evolving edge structures and nodal
#' attributes. Nodal attribute dynamics in \code{ndtv} movies require a temporally
#' extended attribute (TEA) containing a standard R color for each node at each
#' time step. By default, the \code{EpiModel} package uses TEAs to store disease
#' status history in network model simulations run in \code{\link{netsim}}. But
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
#'        output of \code{\link{get_formula_terms}}.
#'
#' @seealso \code{\link{get_formula_terms}}, \code{\link{get_attr_prop}},
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
#' @param d.rate Death or exit rate from the population, as a single homogenous
#'        rate that applies to the entire population.
#'
#' @details
#' This function performs two calculations for dissolution coefficients
#' used in a network model estimated with \code{\link{netest}}:
#' \enumerate{
#'  \item \strong{Transformation:} the mean duration of edges in a network are
#'        mathematically transformed to logit coefficients.
#'  \item \strong{Adjustment:} In a dynamic network simulation in an open
#'        population (in which there are deaths), it is further necessary to
#'        adjust these coefficients for dynamic simulations; this upward adjustment
#'        accounts for death as a competing risk to edge dissolution.
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
#'  \item \strong{coef.crude:} mean durations transformed into a logit
#'        coefficients.
#'  \item \strong{coef.adj:} crude coefficients adjusted for the risk of
#'        death on edge persistence, if the \code{d.rate} argument is supplied.
#'  \item \strong{d.rate:} the death rate.
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
#' # Homogeneous dissolution model with no deaths
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#'
#' # Homogeneous dissolution model with deaths
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25,
#'                   d.rate = 0.001)
#'
#' # Heterogeneous dissolution model in which same-race edges have
#' # shorter duration compared to mixed-race edges, with no deaths
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10))
#'
#' # Heterogeneous dissolution model in which same-race edges have
#' # shorter duration compared to mixed-race edges, with deaths
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10), d.rate = 0.001)
#'
#'
dissolution_coefs <- function(dissolution, duration, d.rate = 0) {

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
    if (sqrt(ps2) <= pg) {
      stop("The competing risk of mortality is too high given the duration. Specify a lower d.rate",
           call. = FALSE)
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
         stop("The competing risk of mortality is too high for the given the ",
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
#' sim <- netdx(est, nsims = 1, nsteps = 100, verbose = FALSE)
#' el <- sim$edgelist[[1]]
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
#' sim <- netdx(est, nsims = 1, nsteps = 100, verbose = FALSE)
#' el <- sim$edgelist[[1]]
#'
#' # Calculate ages directly from edgelist
#' ( ma <- edgelist_meanage(el = el) )
#'
#' # Alternatively, netdx calculates these
#' sim$pages
#' identical(sim$pages[[1]], ma)
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
#'        output of \code{\link{get_formula_terms}}.
#' @param only.formula Limit the tables to those terms only in \code{fterms},
#'        otherwise output proportions for all attributes on the network object.
#'
#' @seealso \code{\link{get_formula_terms}}, \code{\link{copy_toall_attr}},
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


#' @title Outputs Formula Terms into a Character Vector
#'
#' @description Given a formation formula for a network model, outputs it into
#'              a character vector of terms to be used in \code{netsim}
#'              simulations.
#'
#' @param formula Right-hand sided formation formula.
#'
#' @seealso \code{\link{copy_toall_attr}}, \code{\link{get_attr_prop}},
#'          \code{\link{update_nwattr}}.
#' @keywords netUtils internal
#' @export
#'
get_formula_terms <- function(formula) {

  fterms <- attributes(terms.formula(formula))$term.labels
  fterms <- strsplit(fterms, split = "[\"]")
  tl <- sapply(fterms, length)
  if (all(tl == 1)) {
    fterms <- NULL
  } else {
    fterms <- fterms[tl > 1]
    fterms <- unique(sapply(fterms, function(x) x[2]))
  }

  return(fterms)
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


#' @title Update Attribute Values for a Bipartite Network
#'
#' @description Adds new values for attributes in a bipartite network in which
#'              there may be births/entries in the first mode, which requires
#'              splitting the attribute vector into two, adding the new values,
#'              and re-concatenating the two updated vectors.
#'
#' @param dat Master data object passed through \code{netsim} simulations.
#' @param var Variable to update.
#' @param val Fixed value to set for all incoming nodes.
#' @param nCurrM1 Number currently in mode 1.
#' @param nCurrM2 Number currently in mode 2.
#' @param nBirths Number of births/entries in mode 1.
#' @param nBirthsM2 Number of births/entries in mode 2.
#'
#' @export
#' @keywords netUtils internal
#'
split_bip <- function(dat, var, val, nCurrM1, nCurrM2, nBirths, nBirthsM2) {

  oldVarM1 <- dat$attr[[var]][1:nCurrM1]
  oldVarM2 <- dat$attr[[var]][(nCurrM1 + 1):(nCurrM1 + nCurrM2)]

  newVarM1 <- c(oldVarM1, rep(val, nBirths))
  newVarM2 <- c(oldVarM2, rep(val, nBirthsM2))

  newVar <- c(newVarM1, newVarM2)

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
