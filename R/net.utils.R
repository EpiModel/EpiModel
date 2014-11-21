
# Exported Functions ------------------------------------------------------

#' @title Vertex Attributes for Bipartite Network
#'
#' @description Outputs static vertex attributes for a bipartite network for one
#'              specified mode.
#'
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param mode mode number to extract values for.
#' @param val static attribute values to return.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network.initialize(10, bipartite = 5)
#' nw %v% "sex" <- rep(c("M", "F"), each = 5)
#' bipvals(nw, mode = 1, "sex")
#'
bipvals <- function(nw,
                    mode,
                    val
                    ) {

  if (!is.numeric(nw$gal$bipartite))
    stop("nw must be a bipartite network")
  if (missing(mode))
    stop("Specify mode=1 or mode=2")

  nw %s% modeids(nw, mode) %v% val
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
#' \href{http://statnet.org/EpiModel/vignette/NetUtils.html}{EpiModel Network
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

  num <- num.m1 + num.m2

  deg.counts.m1 <- deg.dist.m1 * num.m1
  deg.counts.m2 <- deg.dist.m2 * num.m2

  tot.deg.m1 <- sum(deg.counts.m1 * (1:length(deg.dist.m1) - 1))
  tot.deg.m2 <- sum(deg.counts.m2 * (1:length(deg.dist.m2) - 1))

  mat <- matrix(c(deg.dist.m1, deg.counts.m1,
                  deg.dist.m2, deg.counts.m2), ncol=4)
  mat <- rbind(mat, c(sum(deg.dist.m1), tot.deg.m1, sum(deg.dist.m2), tot.deg.m2))

  colnames(mat) <- c("m1.dist", "m1.cnt", "m2.dist", "m2.cnt")
  rownames(mat) <- c(paste0("Deg", 0:(length(deg.dist.m1)-1)), "Edges")

  cat("Bipartite Degree Distribution Check\n")
  cat("=============================================\n")
  print(mat, print.gap = 3)
  cat("=============================================\n")

  if (sum(deg.dist.m1) != 1 | sum(deg.dist.m2) != 1 | round(tot.deg.m1) != round(tot.deg.m2)) {
    if (sum(deg.dist.m1) != 1) {
      cat("** deg.dist.m1 TOTAL != 1 \n")
    }
    if (sum(deg.dist.m2) != 1) {
      cat("** deg.dist.m2 TOTAL != 1 \n")
    }

    if (round(tot.deg.m1) != round(tot.deg.m2)) {
      reldiff <- (tot.deg.m1-tot.deg.m2)/tot.deg.m2
      if (reldiff > 0) {
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
#' @param nd an object of class \code{networkDynamic}.
#' @param old.var old TEA variable name.
#' @param old.sus status value for susceptible in old TEA variable.
#' @param old.inf status value for infected in old TEA variable.
#' @param old.rec status value for recovered in old TEA variable.
#' @param new.var new TEA variable name to be stored in \code{networkDynamic}
#'        object.
#' @param new.sus status value for susceptible in new TEA variable.
#' @param new.inf status value for infected in new TEA variable.
#' @param new.rec status value for recovered in new TEA variable.
#' @param verbose print progress to console.
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
color_tea <- function(nd,
                      old.var = "testatus",
                      old.sus = "s",
                      old.inf = "i",
                      old.rec = "r",
                      new.var = "ndtvcol",
                      new.sus,
                      new.inf,
                      new.rec,
                      verbose = TRUE) {

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
      cat("\n", at, "/", max(times), "\t", sep="")
    }
  }

  return(nd)
}



#' @title Copies Vertex Attributes in Formation Formula to attr List
#'
#' @description Copies the vertex attributes stored on the network object to the
#'              master attr list in the dat data object.
#'
#' @param dat master data object passed through \code{netsim} simulations.
#' @param at current time step.
#' @param fterms vector of attributes used in formation formula, usually as
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
#' @param dissolution a right-hand sided STERGM dissolution formula
#'        (see \code{\link{netest}}); currently limited to a \code{~offset(edges)}
#'        dissolution model.
#' @param duration an average edge duration in arbitrary time units.
#' @param d.rate death/exit rate in the absence of disease.
#'
#' @details
#' This function performs two calculations to obtain a dissolution coefficient
#' used in a network model estimated with \code{\link{netest}}:
#' \enumerate{
#'  \item \strong{Transformation:} the average duration of edges in a network are
#'        mathematically transformed to logit coefficients.
#'  \item \strong{Adjustment:} In a dynamic network simulation in an open
#'        population (in which there are births and deaths), it is necessary to
#'        adjust the dissolution coefficient for the STERGM simulations to account
#'        for the death as a competing risk to edge dissolution.
#' }
#'
#' Future releases of \code{EpiModel} will allow for more flexibility in the
#' possible dissolution models that may be calculated, including models with
#' heterogenous dissolution probabilities conditional on nodal or edge attributes.
#'
#' @return
#' A list of class \code{disscoef} with the following elements:
#' \itemize{
#'  \item \strong{dissolution:} the right-hand sided STERGM dissolution formula
#'         passed in the function call.
#'  \item \strong{duration:} the average edge duration passed in the function
#'        call.
#'  \item \strong{coef.crude:} the average duration transformed into a logit
#'        coefficient for use in STERGM simulations.
#'  \item \strong{coef.adj:} the crude coefficient adjusted for the impact of
#'        death on edge persistence, if the \code{d.rate} argument is supplied.
#'  \item \strong{d.rate:} the death rate, supplied via the \code{d.rate} argument.
#' }
#'
#' @seealso
#' The theory and details of this function are explained in detail in the
#' \href{http://statnet.org/EpiModel/vignette/NetUtils.html}{EpiModel Network
#' Utility Functions} tutorial.
#'
#' @export
#' @keywords netUtils
#'
#' @examples
#' dissolution <- ~offset(edges)
#' duration <- 25
#' dissolution_coefs(dissolution, duration)
#' dissolution_coefs(dissolution, duration, d.rate = 0.001)
#'
dissolution_coefs <- function(dissolution,
                              duration,
                              d.rate = 0
                              ) {


  # Check form of dissolution formula
  form.length <- length(strsplit(as.character(dissolution)[2], "[+]")[[1]])
  t1.edges <- grepl("offset[(]edges",
                    strsplit(as.character(dissolution)[2], "[+]")[[1]][1])

  # Log transformation of duration to coefficent
  if (t1.edges == TRUE && form.length == 1) {
    coef.diss <- log(duration[1] - 1)
  } else {
    stop("Only ~offset(edges) dissolution models currently supported",
         call. = FALSE)
  }

  if (d.rate > 0) {
    # Exogenous death correction to coefficient
    exp.dur <- 1 + exp(coef.diss)
    prob.diss <- 1 / exp.dur

    prob.neither.dying <- (1 - d.rate)^2
    prob.either.dying <- 2*d.rate - d.rate^2

    prob <- 1 - ((prob.diss - prob.either.dying) / prob.neither.dying)
    if (prob >= 1) {
      stop("The competing risk of mortality is too high for the given duration. Specify a lower d.rate",
           call. = FALSE)
    }
    coef.diss.adj <- logit(prob)
  } else {
    coef.diss.adj <- coef.diss
  }

  out <- list()
  out$dissolution <- dissolution
  out$duration <- duration
  out$coef.adj <- coef.diss.adj
  out$coef.crude <- coef.diss
  out$d.rate <- d.rate

  class(out) <- "disscoef"
  return(out)
}



#' @title Table of Edge Censoring
#'
#' @description Outputs a table of the number and percent of edges that are
#'              left-censored, right-censored, both-censored, or uncensored for a
#'              \code{networkDynamic} object.
#'
#' @param el a timed edgelist with start and end times extracted from a
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
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution_coefs(dissolution, duration = 20)
#'
#' # Model estimation
#' est <- netest(nw,
#'               formation,
#'               dissolution,
#'               target.stats,
#'               coef.diss,
#'               verbose = FALSE)
#'
#' # Simulate the network and extract a timed edgelist
#' sim <- netdx(est, nsims = 1, nsteps = 100, verbose = FALSE)
#' el <- sim$edgelist[[1]]
#'
#' # Calculate censoring
#' edgelist_censor(el)
#'
edgelist_censor <- function(el) {

  time.steps <- max(el$terminus)
  min.step <- min(el$onset)

  # left censored
  leftcens <- el$onset.censored
  leftcens.num <- sum(leftcens)
  leftcens.pct <- leftcens.num/nrow(el)

  # right censored
  rightcens <- el$terminus.censored
  rightcens.num <- sum(rightcens)
  rightcens.pct <- rightcens.num/nrow(el)

  # partnership lasts for entire window (left and right censored)
  lrcens <- el$onset.censored & el$terminus.censored
  lrcens.num <- sum(lrcens)
  lrcens.pct <- lrcens.num/nrow(el)

  # fully observed
  nocens <- el$onset.censored == FALSE & el$terminus.censored == FALSE
  nocens.num <- sum(nocens)
  nocens.pct <- nocens.num/nrow(el)

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
#' @param x an \code{EpiModel} object of class \code{\link{netest}}.
#' @param el if not passing \code{x}, a timed edgelist from a \code{networkDynamic}
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
#' formation <- ~ edges
#' target.stats <- 50
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution_coefs(dissolution, duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, dissolution,
#'               target.stats, coef.diss, verbose = FALSE)
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

  meanpage <- meanpage[1:(length(meanpage)-1)]

  return(meanpage)
}


#' @title Proportional Table of Vertex Attributes
#'
#' @description Calculates the proportional distribution of each vertex attribute
#'              contained on the network, with a possible limitation to those
#'              attributes contained in the formation formula only.
#'
#' @param nw the \code{networkDynamic} object contained in the \code{netsim}
#'        simulation.
#' @param fterms vector of attributes used in formation formula, usually as
#'        output of \code{\link{get_formula_terms}}.
#' @param only.formula limit the tables to those terms only in \code{fterms},
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
#' @param formula a right-hand sided formation formula.
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



#' @title Get Epidemic Output from netsim Model
#'
#' @description Provides all active model state sizes from the network at the
#'              specified time step, output to a list of vectors.
#'
#' @param dat a list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{\link{netsim}}.
#' @param at current time step.
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
  modes <- dat$param$modes

  # Subset attr to active == 1
  l <- lapply(1:length(dat$attr), function(x) dat$attr[[x]][active == 1])
  names(l) <- names(dat$attr)
  l$active <- l$infTime <- NULL

  status <- l$status

  if (modes == 2) {
    mode <- idmode(dat$nw)[active == 1]
  }

  ## Subsetting for epi.by control
  eb <- !is.null(dat$control$epi.by)
  if (eb == TRUE) {
    ebn <- dat$control$epi.by
    ebv <- dat$temp$epi.by.vals
    ebun <- paste0(".", ebn, ebv)
    assign(ebn, l[[ebn]])
  }

  ## One mode networks
  if (modes == 1) {
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

  } else {
    # Bipartite networks
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
  }

  return(dat)
}



#' @title Mode Numbers for Bipartite Network
#'
#' @description Outputs mode numbers give ID numbers for a bipartite network.
#'
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param ids a vector of ID numbers for which the mode number
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
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param mode mode number to return ID numbers for.
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
    out <- (m1size+1):n
  }

  return(out)
}


#' @title Query Active Nodes in NetworkDynamic Object
#'
#' @description Outputs information on the active nodes in a \code{networkDynamic}
#'              object.
#'
#' @param nw an object of class \code{networkDynamic}.
#' @param at current time step.
#' @param out function output, with options of \code{out="vec"} for
#'        a T/F vector of whether the node is active, \code{out="ids"} for
#'        a vector of IDs active, \code{out="prev"} for the number of
#'        nodes that are active, and \code{out="all"} to return a list of
#'        the prior three elements.
#' @param mode if \code{nw} is bipartite, the mode number for status (may
#'        be ignored if requesting output for both modes).
#' @param active.default if \code{TRUE}, elements without an activity attribute
#'        will be regarded as active.
#'
#' @details
#' This is a specialized version of \code{\link{is.active}} from the
#' \code{networkDynamic} package that allows for key output to be efficiently
#' generated for use in \code{\link{netsim}} simulations.
#'
#' @seealso \code{\link{is.active}}, \code{\link{get_prev.net}}
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' # Initialize NW and activate vertices
#' nw <- network.initialize(20)
#'
#' # Activate all vertices, then deactive half at time 5
#' activate.vertices(nw, onset = 1, terminus = 10)
#' deactivate.vertices(nw, onset = 5, terminus = 10, v = 1:10)
#'
#' # Output all information for vertices at time 1 and time 5
#' node_active(nw, at = 1, out = "all")
#' node_active(nw, at = 5, out = "all")
#'
node_active <- function(nw,
                        at,
                        out,
                        mode,
                        active.default = FALSE
                        ) {

  if (!(missing(mode)) && !is.numeric(nw$gal$bipartite))
    stop("nw must be bipartite if mode argument is used")

  if (out %in% c("vec", "ids", "prev")) {
    if (missing(mode)) {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.vec <- node.active
      out.ids <- which(node.active)
      out.prev <- sum(node.active)
    } else {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      if (mode == 1) {
        out.vec <- node.active[ids.m1]
        out.ids <- intersect(which(node.active), ids.m1)
        out.prev <- sum(node.active[ids.m1])
      }
      if (mode == 2) {
        out.vec <- node.active[ids.m2]
        out.ids <- intersect(which(node.active), ids.m2)
        out.prev <- sum(node.active[ids.m2])
      }
    }
  }
  if (out == "all") {
    if (!is.numeric(nw$gal$bipartite)) {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.all <- list()
      out.all$vec$all <- node.active
      out.all$ids$all <- which(node.active)
      out.all$prev$all <- sum(node.active)
    } else {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.all <- list()
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      out.all$vec$m1 <- node.active[ids.m1]
      out.all$vec$m2 <- node.active[ids.m2]
      out.all$vec$all <- node.active
      out.all$ids$m1 <- intersect(which(node.active), ids.m1)
      out.all$ids$m2 <- intersect(which(node.active), ids.m2)
      out.all$ids$all <- which(node.active)
      out.all$prev$m1 <- sum(node.active[ids.m1])
      out.all$prev$m2 <- sum(node.active[ids.m2])
      out.all$prev$all <- sum(node.active)
    }
  }

  if (out == "vec") {
    return(out.vec)
  }
  if (out == "ids") {
    return(out.ids)
  }
  if (out == "prev") {
    return(out.prev)
  }
  if (out == "all") {
    return(out.all)
  }

}


#' @title Update Attribute Values for a Bipartite Network
#'
#' @description Adds new values for attributes in a bipartite network in which
#'              there may be births/entries in the first mode, which requires
#'              splitting the attribute vector into two, adding the new values,
#'              and re-concatenating the two updated vectors.
#'
#' @param dat master data object passed through \code{netsim} simulations.
#' @param var variable to update.
#' @param val fixed value to set for all incoming nodes.
#' @param nCurrM1 number currently in mode 1.
#' @param nCurrM2 number currently in mode 2.
#' @param nBirths number of births/entries in mode 1.
#' @param nBirthsM2 number of births/entries in mode2.
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
#' @param nw the \code{networkDynamic} object used in \code{netsim} simulations.
#' @param newNodes vector of nodal IDs for incoming nodes at the current time
#'        step.
#' @param rules list of rules, one per attribute to be set, governing how to set
#'        the values of each attribute.
#' @param curr.tab current proportional distribution of all vertex attributes.
#' @param t1.tab proportional distribution of all vertex attributes at the outset
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


# Unexported Functions ----------------------------------------------------

# logit transformation of a probability
logit <- function(x) {
  log(x / (1 - x))
}

