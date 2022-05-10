# Network-related Utility Functions -----------------------------------------

#' @title Check Degree Distribution for Balance in Target Statistics
#'
#' @description Checks for consistency in the implied network statistics
#'              of a two-group network in which the group size and
#'              group-specific degree distributions are specified.
#'
#' @param num.g1 Number of nodes in group 1.
#' @param num.g2 Number of nodes in group 2.
#' @param deg.dist.g1 Vector with fractional degree distribution for group 1.
#' @param deg.dist.g2 Vector with fractional degree distribution for group 2.
#'
#' @details
#' This function outputs the number of nodes of degree 0 to g, where g is the
#' length of a fractional degree distribution vector, given that vector and the
#' size of the group. This utility is used to check for balance in implied
#' degree given that fractional distribution within two-group network
#' simulations, in which the degree-constrained counts must be equal across
#' groups.
#'
#' @export
#' @keywords netUtils
#'
#' @examples
#' # An unbalanced distribution
#' check_degdist_bal(num.g1 = 500, num.g2 = 500,
#'                   deg.dist.g2 = c(0.40, 0.55, 0.03, 0.02),
#'                   deg.dist.g1 = c(0.48, 0.41, 0.08, 0.03))
#'
#' # A balanced distribution
#' check_degdist_bal(num.g1 = 500, num.g2 = 500,
#'                   deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
#'                   deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03))
#'
check_degdist_bal <- function(num.g1, num.g2,
                              deg.dist.g1, deg.dist.g2) {
  deg.counts.g1 <- deg.dist.g1 * num.g1
  deg.counts.g2 <- deg.dist.g2 * num.g2
  tot.deg.g1 <- sum(deg.counts.g1 * (seq_along(deg.dist.g1) - 1))
  tot.deg.g2 <- sum(deg.counts.g2 * (seq_along(deg.dist.g2) - 1))
  mat <- matrix(c(deg.dist.g1, deg.counts.g1,
                  deg.dist.g2, deg.counts.g2), ncol = 4)
  mat <- rbind(mat, c(sum(deg.dist.g1), tot.deg.g1, sum(deg.dist.g2),
                      tot.deg.g2))
  colnames(mat) <- c("g1.dist", "g1.cnt", "g2.dist", "g2.cnt")
  rownames(mat) <- c(paste0("Deg", 0:(length(deg.dist.g1) - 1)), "Edges")
  cat("Degree Distribution Check\n")
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
        msg <- "Group 1 Edges > Group 2 Edges:"
      } else {
        msg <- "Group 1 Edges < Group 2 Edges:"
      }
      cat("**", msg, round(reldiff, 3), "Rel Diff \n")
    }
  } else {
    cat("** Edges balanced ** \n")
  }
  invisible(c(tot.deg.g1, deg.counts.g1, deg.counts.g2))
}


#' @title Create a TEA Variable for Infection Status for \code{ndtv} Animations
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
#' @param verbose If \code{TRUE}, print progress to console.
#'
#' @details
#' The \code{ndtv} package (\url{https://cran.r-project.org/package=ndtv})
#' produces animated visuals for dynamic networks with evolving edge structures
#' and nodal attributes. Nodal attribute dynamics in \code{ndtv} movies require
#' a temporally extended attribute (TEA) containing a standard R color for each
#' node at each time step. By default, the \code{EpiModel} package uses TEAs to
#' store disease status history in network model simulations run in
#' \code{\link{netsim}}. But that status TEA is in numeric format (0, 1, 2).
#' The \code{color_tea} function transforms those numeric values of that disease
#' status TEA into a TEA with color values in order to visualize status changes
#' in \code{ndtv}.
#'
#' The convention in \code{\link{plot.netsim}} is to color the susceptible
#' nodes as blue, infected nodes as red, and recovered nodes as green. Alternate
#' colors may be specified using the \code{new.sus}, \code{new.inf}, and
#' \code{new.rec} parameters, respectively.
#'
#' Using the \code{color_tea} function with a \code{netsim} object requires that
#' TEAs for disease status be used and that the \code{networkDynamic} object be
#' saved in the output: \code{tergmListe} must be  set to \code{FALSE} in
#' \code{\link{control.net}}.
#'
#' @return The updated object of class \code{networkDynamic}.
#'
#' @seealso \code{\link{netsim}} and the \code{ndtv} package documentation.
#' @keywords colorUtils
#' @export
#'
color_tea <- function(nd, old.var = "testatus", old.sus = "s", old.inf = "i",
                      old.rec = "r", new.var = "ndtvcol", new.sus, new.inf,
                      new.rec, verbose = TRUE) {
  if (missing(new.inf)) {
    new.inf <- adjustcolor(2, 0.75)
  }
  if (missing(new.sus)) {
    new.sus <- adjustcolor(4, 0.75)
  }
  if (missing(new.rec)) {
    new.rec <- adjustcolor(3, 0.75)
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


#' @title Copy Vertex Attributes From Network to \code{dat} List
#'
#' @description Copies the vertex attributes stored on the network object to the
#'              main \code{attr} list in the \code{dat} data object.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @seealso \code{\link{get_formula_term_attr}}, \code{\link{get_attr_prop}},
#'          \code{\link{auto_update_attr}}, and
#'          \code{\link{copy_datattr_to_nwattr}}.
#' @keywords netUtils internal
#' @export
#'
copy_nwattr_to_datattr <- function(dat) {
  otha <- list.vertex.attributes(dat$nw[[1]])
  otha <- setdiff(otha, c("na", "vertex.names", "active",
                          "testatus.active", "tergm_pid"))
  if (length(otha) > 0) {
    for (i in seq_along(otha)) {
      va <- get_vertex_attribute(dat$nw[[1]], otha[i])
      dat$attr[[otha[i]]] <- va
      if (!is.null(dat$control$epi.by) && dat$control$epi.by == otha[i]) {
        dat$temp$epi.by.vals <- unique(va)
      }
    }
  }
  return(dat)
}


#' @title Copy Vertex Attributes from the \code{dat} List to the Network Object
#'
#' @description Copies the vertex attributes stored on the main \code{attr} list
#'              on \code{dat} to the network object on \code{dat}.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @seealso \code{\link{get_formula_term_attr}}, \code{\link{get_attr_prop}},
#'          \code{\link{auto_update_attr}}, and
#'          \code{\link{copy_nwattr_to_datattr}}.
#' @keywords netUtils internal
#' @export
#'
copy_datattr_to_nwattr <- function(dat) {
  nwterms <- dat$temp$nwterms
  special.attr <- "status"
  if (dat$param$groups == 2) {
    special.attr <- c(special.attr, "group")
  }
  nwterms <- union(nwterms, special.attr)
  attr.to.copy <- union(nwterms, special.attr)
  attr <- dat$attr[attr.to.copy]
  if (length(attr.to.copy) > 0) {
    if (length(attr.to.copy) == 1) {
      dat$nw[[1]] <- set_vertex_attribute(dat$nw[[1]], names(attr), attr[[1]])
    } else {
      dat$nw[[1]] <- set_vertex_attribute(dat$nw[[1]], names(attr), attr)
    }
  }
  return(dat)
}

#' @title Dissolution Coefficients for Stochastic Network Models
#'
#' @description Calculates dissolution coefficients, given a dissolution model
#'              and average edge duration, to pass as offsets to an ERGM/TERGM
#'              model fit in \code{netest}.
#'
#' @param dissolution Right-hand sided STERGM dissolution formula
#'        (see \code{\link{netest}}). See below for list of supported
#'        dissolution models.
#' @param duration A vector of mean edge durations in arbitrary time units.
#' @param d.rate Departure or exit rate from the population, as a single
#'        homogeneous rate that applies to the entire population.
#'
#' @details
#' This function performs two calculations for dissolution coefficients
#' used in a network model estimated with \code{\link{netest}}:
#' \enumerate{
#'  \item \strong{Transformation:} the mean durations of edges in a network are
#'        mathematically transformed to logit coefficients.
#'  \item \strong{Adjustment:} in a dynamic network simulation in an open
#'        population (in which there are departures), it is further necessary to
#'        adjust these coefficients; this upward adjustment accounts for
#'        departure as a competing risk to edge dissolution.
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
#'         dyad have similar values of a specified attribute. The duration
#'         vector should now contain two values: the first is the mean edge
#'         duration of non-matched dyads, and the second is the duration of the
#'         matched dyads.
#'  \item \code{~offset(edges) + offset(nodemix("<attr>"))}: a heterogeneous
#'         model that extends the nodematch model to include non-binary
#'         attributes for homophily. The duration vector should first contain
#'         the base value, then the values for every other possible combination
#'         in the term.
#'  \item \code{~offset(edges) + offset(nodefactor("<attr>"))}: a heterogeneous
#'         model in which the edge duration varies by a specified attribute. The
#'         duration vector should first contain the base value, then the values
#'         for every other value of that attribute in the term. This option is
#'         deprecated.
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
#'        departure on edge persistence, if the \code{d.rate} argument is
#'        supplied.
#'  \item \strong{coef.form.corr:} corrections to be subtracted from formation
#'        coefficients.
#'  \item \strong{d.rate:} the departure rate.
#'  \item \strong{diss.model.type:} the form of the dissolution model; options
#'        include \code{edgesonly}, \code{nodematch}, \code{nodemix}, and
#'        \code{nodefactor}.
#' }
#'
#' @export
#' @keywords netUtils
#'
#' @examples
#' ## Homogeneous dissolution model with no departures
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25)
#'
#' ## Homogeneous dissolution model with departures
#' dissolution_coefs(dissolution = ~offset(edges), duration = 25,
#'                   d.rate = 0.001)
#'
#' ## Heterogeneous dissolution model in which same-race edges have
#' ## shorter duration compared to mixed-race edges, with no departures
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10))
#'
#' ## Heterogeneous dissolution model in which same-race edges have
#' ## shorter duration compared to mixed-race edges, with departures
#' dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
#'                   duration = c(20, 10), d.rate = 0.001)
#'
#' \dontrun{
#' ## Extended example for differential homophily by age group
#' # Set up the network with nodes categorized into 5 age groups
#' nw <- network_initialize(n = 1000)
#' age.grp <- sample(1:5, 1000, TRUE)
#' nw <- set_vertex_attribute(nw, "age.grp", age.grp)
#'
#' # durations = non-matched, age.grp1 & age.grp1, age.grp2 & age.grp2, ...
#' # TERGM will include differential homophily by age group with nodematch term
#' # Target stats for the formation model are overall edges, and then the number
#' # matched within age.grp 1, age.grp 2, ..., age.grp 5
#' form <- ~edges + nodematch("age.grp", diff = TRUE)
#' target.stats <- c(450, 100, 125, 40, 80, 100)
#'
#' # Target stats for the dissolution model are duration of non-matched edges,
#' # then duration of edges matched within age.grp 1, age.grp 2, ..., age.grp 5
#' durs <- c(60, 30, 80, 100, 125, 160)
#' diss <- dissolution_coefs(~offset(edges) +
#'                             offset(nodematch("age.grp", diff = TRUE)),
#'                           duration = durs)
#'
#' # Fit the TERGM
#' fit <- netest(nw, form, target.stats, diss)
#'
#' # Full diagnostics to evaluate model fit
#' dx <- netdx(fit, nsims = 10, ncores = 4, nsteps = 300)
#' print(dx)
#'
#' # Simulate one long time series to examine timed edgelist
#' dx <- netdx(fit, nsims = 1, nsteps = 5000, keep.tedgelist = TRUE)
#'
#' # Extract timed-edgelist
#' te <- as.data.frame(dx)
#' head(te)
#'
#' # Limit to non-censored edges
#' te <- te[which(te$onset.censored == FALSE & te$terminus.censored == FALSE),
#'          c("head", "tail", "duration")]
#' head(te)
#'
#' # Look up the age group of head and tail nodes
#' te$ag.head <- age.grp[te$head]
#' te$ag.tail <- age.grp[te$tail]
#' head(te)
#'
#' # Recover average edge durations for age-group pairing
#' mean(te$duration[te$ag.head != te$ag.tail])
#' mean(te$duration[te$ag.head == 1 & te$ag.tail == 1])
#' mean(te$duration[te$ag.head == 2 & te$ag.tail == 2])
#' mean(te$duration[te$ag.head == 3 & te$ag.tail == 3])
#' mean(te$duration[te$ag.head == 4 & te$ag.tail == 4])
#' mean(te$duration[te$ag.head == 5 & te$ag.tail == 5])
#' durs
#' }
#'
dissolution_coefs <- function(dissolution, duration, d.rate = 0) {
  # Error check for duration < 1
  if (any(duration < 1)) {
    stop("All values in duration must be >= 1", call. = FALSE)
  }
  # Check form of dissolution formula
  diss.model.type <- NA
  form.length <- length(strsplit(as.character(dissolution)[2], "[+]")[[1]])
  t1.edges <- grepl("offset[(]edges",
                    strsplit(as.character(dissolution)[2], "[+]")[[1]][1])
  if (form.length == 1 && t1.edges == TRUE) {
    diss.model.type <- "edgesonly"
  } else {
    if (form.length == 2 && t1.edges == TRUE) {
      t2 <- strsplit(as.character(dissolution)[2], "[+]")[[1]][2]
      t2.term <- NULL
      if (grepl("offset[(]nodematch", t2)) {
        t2.term <- diss.model.type <- "nodematch"
      } else {
        if (grepl("offset[(]nodefactor", t2)) {
          t2.term <- diss.model.type <- "nodefactor"
          warning("Support for dissolution models containing a nodefactor term
                  is deprecated, and will be removed in a future release.")
          # TODO: remove functionality and deprecation message in future release
        } else {
          if (grepl("offset[(]nodemix", t2)) {
          t2.term <- diss.model.type <- "nodemix"
          } else stop("The form of the dissolution argument is invalid. Type
                      help(\'dissolution_coefs\') to see the set of options
                      allowed.")
        }
      }
    }
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
    if (ps2 <= pg) {
      d.rate_ <- round(1 - sqrt(pg), 5)
      str <- paste("The competing risk of departure is too high for the given",
                   " duration of ", duration[1],
                   "; specify a d.rate lower than ", d.rate_, ".", sep = "")
      stop(str, call. = FALSE)
    }

    coef.crude <- log(pg / (1 - pg))
    coef.adj <- log(pg / (ps2 - pg))
    coef.form.corr <- log(1 + pg / (1 - pg))
  }
  if (form.length == 2) {
    if (t2.term %in% c("nodematch", "nodefactor", "nodemix")) {
      coef.crude <- coef.adj <- coef.form.corr <- NA
      for (i in seq_along(duration)) {
        pg <- (duration[i] - 1) / duration[i]
        ps2 <- (1 - d.rate) ^ 2

        if (ps2 <= pg) {
          d.rate_ <- round(1 - sqrt(pg), 5)
          stop("The competing risk of departure is too high for the given",
               " edge duration of ", duration[i], " in place ", i, ". ",
               "Specify a d.rate lower than ", d.rate_, ".", sep = "")
        }
        if (i == 1) {
          coef.crude[i] <- log(pg / (1 - pg))
          coef.adj[i] <- log(pg / (ps2 - pg))
          coef.form.corr[i] <- log(1 + pg / (1 - pg))
        } else {
          coef.crude[i] <- log(pg / (1 - pg)) - coef.crude[1]
          coef.adj[i] <- log(pg / (ps2 - pg)) - coef.adj[1]
          coef.form.corr[i] <- log(1 + pg / (1 - pg)) - coef.form.corr[1]
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
  out$coef.form.corr <- coef.form.corr
  out$d.rate <- d.rate
  out$diss.model.type <- diss.model.type
  class(out) <- "disscoef"
  return(out)
}


#' @title Table of Edge Censoring
#'
#' @description Outputs a table of the number and percent of edges that are
#'              left-censored, right-censored, both-censored, or uncensored for
#'              a \code{networkDynamic} object.
#'
#' @param el A timed edgelist with start and end times extracted from a
#'        \code{networkDynamic} object using the
#'        \code{as.data.frame.networkDynamic} function.
#'
#' @return A 4 x 2 table containing the number and percent of edges in \code{el}
#'         that are left-censored, right-censored, both-censored, or uncensored.
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
#' nw <- network_initialize(n = 100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#'
#' # Model estimation
#' est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' # Simulate the network and extract a timed edgelist
#' dx <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE,
#'       verbose = FALSE)
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
#' @description Outputs a matrix of mean ages of edges at a series of timesteps.
#'
#' @param el A timed edgelist with start and end times extracted from a
#'        \code{networkDynamic} object using the
#'        \code{as.data.frame.networkDynamic} function.
#' @param diss_term A string indicating the form of heterogeneity present
#'        in the dissolution model (options \code{nodematch}, \code{nodemix},
#'        and \code{nodefactor})
#'        or \code{NULL} for homogeneous (edges-only) dissolution model
#' @param diss_attr A vector containing values of the nodal attribute
#'        associated with the heterogeneity in the dissolution model,
#'        or \code{NULL} for homogeneous (edges-only) dissolution model
#' @param diss_arg A vector containing any additional argument(s) to the
#'        term associated with the heterogeneity in the dissolution model,
#'        or \code{NULL} if there is no relevant argument. At present, this is
#'        only used when the dissolution model contains a \code{nodematch} term,
#'        in which case \code{diss_arg} should equal either \code{TRUE} or
#'        \code{FALSE} to indicate the value of the \code{diff} argument to
#'        that term. Additional uses may be added in the future.
#'
#' @details
#' This function calculates the mean partnership age at each time step over
#' a dynamic network simulation expressed in the form of a timed edgelist.
#' Means may be calculated for all edges, or disaggregated by nodal attribute
#' combinations.
#'
#' @return A matrix containing the mean edge age at each time step (rows),
#' with either one column (for homogeneous models, i.e. when
#' \code{diss_term = NULL}) or one column per attribute value combination
#' (for heterogeneous models).
#'
#' @keywords netUtils internal
#'
edgelist_meanage <- function(el, diss_term = NULL,
                             diss_attr = NULL, diss_arg = NULL) {
# TODO: remove nodefactor from documentation in later version
  terminus <- el$terminus
  onset <- el$onset
  minterm <- 1
  maxterm <- max(terminus)

  if (is.null(diss_term) || diss_term=="nodefactor") {
    # TODO: remove nodefactor in future release
    meanpage <- matrix(NA, maxterm, 1)
  } else {
    attrvalues <- sort(unique(diss_attr))
    n.attrvalues <- length(attrvalues)
    attr1 <- diss_attr[el$head]
    attr2 <- diss_attr[el$tail]
    if(diss_term == "nodematch" && diss_arg == FALSE) {
      meanpage <- matrix(NA, maxterm, 2)
    }
    if(diss_term == "nodematch" && diss_arg == TRUE) {
      meanpage <- matrix(NA, maxterm, n.attrvalues+1)
    }
    if(diss_term == "nodemix") {
      n.attrcombos <- n.attrvalues * (n.attrvalues + 1) / 2
      meanpage <- matrix(NA, maxterm, n.attrcombos)
      indices2.grid <- expand.grid(row = 1:n.attrvalues, col = 1:n.attrvalues)
      rowleqcol <- indices2.grid$row <= indices2.grid$col #assumes undirected
      indices2.grid <- indices2.grid[rowleqcol, ]
  }}

  for (at in minterm:maxterm) {
    actp <- (onset <= at & terminus > at) |
      (onset == at & terminus == at);
    page <- at - onset[actp] + 1

    if (is.null(diss_term) || diss_term == "nodefactor") {
      # TODO: remove nodefactor in future release
      meanpage[at,1] <- mean(page)
    } else {
        attr1a <- attr1[actp]
        attr2a <- attr2[actp]
        if(diss_term == "nodematch") {
          if(diss_arg == TRUE) {
            su_diss_attr <- sort(unique(diss_attr))
            meanpage[at,1] <- mean(page[attr1a != attr2a])
            for (k in seq_along(su_diss_attr)) {
              meanpage[at,k+1] <-
                mean(page[attr1a==su_diss_attr[k] & attr2a==su_diss_attr[k]])
            }
          } else {
            meanpage[at, 1] <- mean(page[attr1a!=attr2a])
            meanpage[at, 2] <- mean(page[attr1a==attr2a])
          }
        }
        if(diss_term == "nodemix") {
          for(i in 1:nrow(indices2.grid)) {
            if(indices2.grid$row[i]==indices2.grid$col[i]) {
              meanpage[at,i] <- mean(
                page[attr1a==attrvalues[indices2.grid$row[i]] &
                         attr2a==attrvalues[indices2.grid$col[i]]]
              )
            } else {
            meanpage[at, i] <- mean(
              c(page[attr1a==attrvalues[indices2.grid$row[i]] &
                     attr2a==attrvalues[indices2.grid$col[i]]],
                page[attr1a==attrvalues[indices2.grid$col[i]] &
                     attr2a==attrvalues[indices2.grid$row[i]]]
                ))
            }
          }
        }
    }
  }
  meanpage <- head(meanpage, -1)  # remove last row
  return(meanpage)
}


#' @title Proportional Table of Vertex Attributes
#'
#' @description Calculates the proportional distribution of each vertex
#'              attribute contained in a network.
#'
#' @inheritParams recovery.net
#' @param nwterms Vector of attributes on the network object, usually as
#'        output of \code{\link{get_formula_term_attr}}.
#'
#' @return A table containing the proportional distribution of each attribute in
#'         \code{nwterms}.
#'
#' @seealso \code{\link{get_formula_term_attr}},
#'          \code{\link{copy_nwattr_to_datattr}},
#'          \code{\link{auto_update_attr}}.
#' @keywords netUtils internal
#' @export
#'
get_attr_prop <- function(dat, nwterms) {

  if (is.null(nwterms)) {
    return(NULL)
  }

  nwVal <- names(dat$attr)
  nwVal <- setdiff(nwVal, c("na", "vertex.names", "active", "entrTime",
                            "exitTime", "infTime", "group", "status"))
  out <- list()
  if (length(nwVal) > 0) {
    for (i in seq_along(nwVal)) {
      tab <- prop.table(table(dat$attr[[nwVal[i]]]))
      out[[i]] <- tab
    }
    names(out) <- nwVal
  }

  return(out)
}


#' @title Output ERGM Formula Attributes into a Character Vector
#'
#' @description Given a formation formula for a network model, outputs a
#'              character vector of vertex attributes to be used in
#'              \code{\link{netsim}} simulations.
#'
#' @param form An ERGM model formula.
#' @param nw A network object.
#'
#' @return A character vector of vertex attributes.
#'
#' @export
#'
get_formula_term_attr <- function(form, nw) {

  nw_attr <- list.vertex.attributes(nw)
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

#' @title Output Network Attributes into a Character Vector
#'
#' @description Given a simulated network, outputs a character vector of vertex
#'              attributes to be used in \code{\link{netsim}} simulations.
#'
#' @param nw A network object.
#'
#' @return A character vector of vertex attributes.
#'
#' @export
#'
get_network_term_attr <- function(nw) {

  nw_attr <- list.vertex.attributes(nw)
  nw_attr <- setdiff(nw_attr, c("active", "vertex.names", "na",
                                "testatus.active", "tergm_pid"))

  if (length(nw_attr) == 0) {
    return(NULL)
  }

  out <- nw_attr
  if (length(out) == 0) {
    return(NULL)
  } else {
    return(out)
  }

}

#' @title Group Numbers for Two-Group Network
#'
#' @description Outputs group numbers given ID numbers for a two-group network.
#'
#' @param nw Object of class \code{network} or \code{networkDynamic}.
#' @param ids Vector of ID numbers for which the group number
#'        should be returned.
#'
#' @return A vector containing the group number for each of the specified nodes.
#'
#' @export
#' @keywords netUtils internal
#'
#' @examples
#' nw <- network_initialize(n = 10)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 5))
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

  flag <- "group" %in% list.vertex.attributes(nw)
  if (!flag) {
    out <- rep(1, n)
  } else {
    groups <- get_vertex_attribute(nw, "group")
    out <- groups[ids]
  }

  return(out)
}

#' @title Update Vertex Attributes for Incoming Vertices
#'
#' @description Updates the vertex attributes on a network for new nodes
#'              incoming into that network, based on a set of rules for each
#'              attribute that the user specifies in \code{\link{control.net}}.
#'
#' @inheritParams recovery.net
#' @param newNodes Vector of nodal IDs for incoming nodes at the current time
#'        step.
#' @param curr.tab Current proportional distribution of all vertex attributes.
#'
#' @inherit recovery.net return
#'
#' @seealso \code{\link{copy_nwattr_to_datattr}}, \code{\link{get_attr_prop}},
#'          \code{\link{auto_update_attr}}.
#' @keywords netUtils internal
#' @export
#'
auto_update_attr <- function(dat, newNodes, curr.tab) {

  rules <- get_control(dat, "attr.rules")
  active <- get_attr(dat, "active")
  t1.tab <- dat$temp$t1.tab

  for (i in seq_along(curr.tab)) {
    vname <- names(curr.tab)[i]
    needs.updating <- ifelse(length(get_attr(dat, vname)) < length(active),
                             TRUE, FALSE)
    if (length(vname) > 0 & needs.updating == TRUE) {
      rule <- rules[[vname]]

      if (is.null(rule)) {
        rule <- "current"
      }
      if (rule == "current") {
        vclass <- class(get_attr(dat, vname))
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
        vclass <- class(get_attr(dat, vname))
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
      dat$attr[[vname]] <- c(dat$attr[[vname]], nattr)
    }
  }

  return(dat)
}


#' @title Get Individual Degree from Network or Edgelist
#'
#' @description A fast method for querying the current degree of all individuals
#'              within a network.
#'
#' @param x Either an object of class \code{network} or \code{edgelist}
#'        generated from a network. If \code{x} is an edgelist, then it must
#'        contain an attribute for the total network size, \code{n}.
#'
#' @details
#' Individual-level data on the current degree of nodes within a network is
#' often useful for summary statistics. Given a \code{network} class object,
#' \code{net}, one way to look up the current degree is to get a summary of the
#' ERGM term, \code{sociality}, as in:
#' \code{summary(net ~ sociality(nodes = NULL))}. But that is computationally
#' inefficient for a number of reasons. This function provides a fast method for
#' generating the vector of degrees using a query of the edgelist. It is even
#' faster if the parameter \code{x} is already transformed into an edgelist.
#'
#' @return A vector of length equal to the total network size, containing the
#'         current degree of each node in the network.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(n = 500)
#'
#' set.seed(1)
#' fit <- ergm(nw ~ edges, target.stats = 250)
#' sim <- simulate(fit)
#'
#' # Slow ERGM-based method
#' ergm.method <- unname(summary(sim ~ sociality(nodes = NULL)))
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
#' @description Left-truncates simulation epidemiological summary statistics
#'              and network statistics at a specified time step.
#'
#' @param x Object of class \code{netsim} or \code{icm}.
#' @param at Time step at which to left-truncate the time series.
#'
#' @details
#' This function would be used when running a follow-up simulation from time
#' steps \code{b} to \code{c} after a burn-in period from time \code{a} to
#' \code{b}, where the final time window of interest for data analysis is
#' \code{b} to \code{c} only.
#'
#' @return The updated object of class \code{netsim} or \code{icm}.
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
  if (!inherits(x, c("icm", "netsim"))) {
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
