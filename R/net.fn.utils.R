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
  if (sum(deg.dist.g1) <= 0.999 || sum(deg.dist.g1) >= 1.001 ||
        sum(deg.dist.g2) <= 0.999 || sum(deg.dist.g2) >= 1.001 || absdiff > 1) {
    if (sum(deg.dist.g1) <= 0.999 || sum(deg.dist.g1) >= 1.001) {
      cat("** deg.dist.g1 TOTAL != 1 \n")
    }
    if (sum(deg.dist.g2) <= 0.999 || sum(deg.dist.g2) >= 1.001) {
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


#' @title Create a TEA Variable for Infection Status for `ndtv` Animations
#'
#' @description Creates a new color-named temporally-extended attribute (TEA)
#'              variable in a `networkDynamic` object containing a disease
#'              status TEA in numeric format.
#'
#' @param nd An object of class `networkDynamic`.
#' @param old.var Old TEA variable name.
#' @param old.sus Status value for susceptible in old TEA variable.
#' @param old.inf Status value for infected in old TEA variable.
#' @param old.rec Status value for recovered in old TEA variable.
#' @param new.var New TEA variable name to be stored in `networkDynamic`
#'        object.
#' @param new.sus Status value for susceptible in new TEA variable.
#' @param new.inf Status value for infected in new TEA variable.
#' @param new.rec Status value for recovered in new TEA variable.
#' @param verbose If `TRUE`, print progress to console.
#'
#' @details
#' The `ndtv` package (<https://cran.r-project.org/package=ndtv>)
#' produces animated visuals for dynamic networks with evolving edge structures
#' and nodal attributes. Nodal attribute dynamics in `ndtv` movies require
#' a temporally extended attribute (TEA) containing a standard R color for each
#' node at each time step. By default, the `EpiModel` package uses TEAs to
#' store disease status history in network model simulations run in
#' [netsim()]. But that status TEA is in numeric format (0, 1, 2).
#' The `color_tea` function transforms those numeric values of that disease
#' status TEA into a TEA with color values in order to visualize status changes
#' in `ndtv`.
#'
#' The convention in [plot.netsim()] is to color the susceptible
#' nodes as blue, infected nodes as red, and recovered nodes as green. Alternate
#' colors may be specified using the `new.sus`, `new.inf`, and
#' `new.rec` parameters, respectively.
#'
#' Using the `color_tea` function with a `netsim` object requires that
#' TEAs for disease status be used and that the `networkDynamic` object be
#' saved in the output: `tergmListe` must be  set to `FALSE` in
#' [control.net()].
#'
#' @return The updated object of class `networkDynamic`.
#'
#' @seealso [netsim()] and the `ndtv` package documentation.
#' @keywords colorUtils
#' @export
#'
color_tea <- function(nd, old.var = "testatus", old.sus = "s", old.inf = "i",
                      old.rec = "r", new.var = "ndtvcol", new.sus = NULL,
                      new.inf = NULL, new.rec = NULL, verbose = TRUE) {
  if (is.null(new.inf)) {
    new.inf <- adjustcolor(2, 0.75)
  }
  if (is.null(new.sus)) {
    new.sus <- adjustcolor(4, 0.75)
  }
  if (is.null(new.rec)) {
    new.rec <- adjustcolor(3, 0.75)
  }
  times <- seq_len(max(get.change.times(nd)))
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


#' @title Copy Vertex Attributes From Network to `netsim_dat` List
#'
#' @description Copies the vertex attributes stored on the network object to the
#'              main `attr` list in the `netsim_dat` data object.
#'
#' @inheritParams recovery.net
#' @param nw Network from which to copy vertex attributes.
#'
#' @inherit recovery.net return
#'
#' @seealso [get_formula_term_attr()], [get_attr_prop()],
#'          [auto_update_attr()], and
#'          [copy_datattr_to_nwattr()].
#' @keywords netUtils internal
#' @export
#'
copy_nwattr_to_datattr <- function(dat, nw) {
  otha <- list.vertex.attributes(nw)
  otha <- setdiff(otha, c("na", "vertex.names", "active",
                          "testatus.active", "tergm_pid"))
  if (length(otha) > 0) {
    epi_by <- get_control(dat, "epi.by", override.null.error = TRUE)
    for (i in seq_along(otha)) {
      va <- get_vertex_attribute(nw, otha[i])
      dat <- set_attr(dat, otha[i], va)
      if (!is.null(epi_by) && epi_by == otha[i]) {
        dat$run$epi.by.vals <- unique(va)
      }
    }
  }
  return(dat)
}


#' @title Copy Vertex Attributes from the `netsim_dat` List to the Network
#'        Objects
#'
#' @description Copies the vertex attributes stored on the main `attr` list
#'              of the `netsim_dat` object to each of the network objects
#'              stored on the `netsim_dat` object.
#'
#' @inheritParams recovery.net
#'
#' @inherit recovery.net return
#'
#' @seealso [get_formula_term_attr()], [get_attr_prop()],
#'          [auto_update_attr()], and
#'          [copy_nwattr_to_datattr()].
#' @keywords netUtils internal
#' @export
#'
copy_datattr_to_nwattr <- function(dat) {
  nwterms <- dat$run$nwterms
  special.attr <- "status"
  if (dat$param$groups == 2) {
    special.attr <- c(special.attr, "group")
  }
  nwterms <- union(nwterms, special.attr)
  attr.to.copy <- union(nwterms, special.attr)
  attr <- get_attr_list(dat, attr.to.copy)
  if (length(attr.to.copy) > 0) {
    if (length(attr.to.copy) == 1) {
      for (network in seq_len(dat$num.nw)) {
        dat$run$nw[[network]] <- set_vertex_attribute(dat$run$nw[[network]],
                                                      names(attr),
                                                      attr[[1]])
      }
    } else {
      for (network in seq_len(dat$num.nw)) {
        dat$run$nw[[network]] <- set_vertex_attribute(dat$run$nw[[network]],
                                                      names(attr),
                                                      attr)
      }
    }
  }
  return(dat)
}

#' @title Dissolution Coefficients for Stochastic Network Models
#'
#' @description Calculates dissolution coefficients, given a dissolution model
#'              and average edge duration, to pass as offsets to an ERGM/TERGM
#'              model fit in `netest`.
#'
#' @param dissolution Right-hand sided STERGM dissolution formula
#'        (see [netest()]). See below for list of supported
#'        dissolution models.
#' @param duration A vector of mean edge durations in arbitrary time units.
#' @param d.rate Departure or exit rate from the population, as a single
#'        homogeneous rate that applies to the entire population.
#'
#' @details
#' This function performs two calculations for dissolution coefficients
#' used in a network model estimated with [netest()]:
#'
#'  1. **Transformation:** the mean durations of edges in a network are
#'     mathematically transformed to logit coefficients.
#'  2. **Adjustment:** in a dynamic network simulation in an open
#'     population (in which there are departures), it is further necessary to
#'     adjust these coefficients; this upward adjustment accounts for
#'     departure as a competing risk to edge dissolution.
#'
#'
#' The current dissolution models supported by this function and in network
#' model estimation in [netest()] are as follows:
#'
#'  * `~offset(edges)`: a homogeneous dissolution model in which the
#'         edge duration is the same for all partnerships. This requires
#'         specifying one duration value.
#'  * `~offset(edges) + offset(nodematch("<attr>"))`: a heterogeneous
#'         model in which the edge duration varies by whether the nodes in the
#'         dyad have similar values of a specified attribute. The duration
#'         vector should now contain two values: the first is the mean edge
#'         duration of non-matched dyads, and the second is the duration of the
#'         matched dyads.
#'  * `~offset(edges) + offset(nodemix("<attr>"))`: a heterogeneous
#'         model that extends the nodematch model to include non-binary
#'         attributes for homophily. The duration vector should first contain
#'         the base value, then the values for every other possible combination
#'         in the term.
#'
#'
#' @return
#' A list of class `disscoef` with the following elements:
#'
#'  * **dissolution:** right-hand sided STERGM dissolution formula
#'         passed in the function call.
#'  * **duration:** mean edge durations passed into the function.
#'  * **coef.crude:** mean durations transformed into logit
#'        coefficients.
#'  * **coef.adj:** crude coefficients adjusted for the risk of
#'        departure on edge persistence, if the `d.rate` argument is
#'        supplied.
#'  * **coef.form.corr:** corrections to be subtracted from formation
#'        coefficients.
#'  * **d.rate:** the departure rate.
#'  * **diss.model.type:** the form of the dissolution model; options
#'        include `edgesonly`, `nodematch`, and `nodemix`.
#'
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
      } else if (grepl("offset[(]nodemix", t2)) {
        t2.term <- diss.model.type <- "nodemix"
      } else {
        stop("The form of the dissolution argument is invalid. Type
              help(\'dissolution_coefs\') to see the set of options
              allowed.")
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
    if (t2.term %in% c("nodematch", "nodemix")) {
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
      stop("Supported heterogeneous dissolution model terms are nodematch ",
           "or nodemix", call. = FALSE)
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
#'              a `networkDynamic` object.
#'
#' @param el A timed edgelist with start and end times extracted from a
#'        `networkDynamic` object using the
#'        `as.data.frame.networkDynamic` function.
#'
#' @return A 4 x 2 table containing the number and percent of edges in `el`
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


#' @title Proportional Table of Vertex Attributes
#'
#' @description Calculates the proportional distribution of each vertex
#'              attribute contained in a network.
#'
#' @inheritParams recovery.net
#' @param nwterms Vector of attributes on the network object, usually as
#'        output of [get_formula_term_attr()].
#'
#' @return
#' A table containing the proportional distribution of each attribute in
#' `nwterms`.
#'
#' @seealso [get_formula_term_attr()],
#'          [copy_nwattr_to_datattr()],
#'          [auto_update_attr()].
#' @keywords netUtils internal
#' @export
#'
get_attr_prop <- function(dat, nwterms) {

  if (is.null(nwterms)) {
    return(NULL)
  }
  attr_list <- get_attr_list(dat)
  nwVal <- names(attr_list)
  nwVal <- setdiff(nwVal, c("na", "vertex.names", "active", "entrTime",
                            "exitTime", "infTime", "group", "status"))
  out <- list()
  if (length(nwVal) > 0) {
    for (i in seq_along(nwVal)) {
      tab <- prop.table(table(attr_list[[nwVal[i]]]))
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
#'              [netsim()] simulations.
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
#'              attributes to be used in [netsim()] simulations.
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
#' @param nw Object of class `network` or `networkDynamic`.
#' @param ids Vector of ID numbers for which the group number
#'        should be returned. If `NULL` (default), return all IDs.
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
idgroup <- function(nw, ids = NULL) {
  n <- network.size(nw)
  if (is.null(ids)) {
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
#'              attribute that the user specifies in [control.net()].
#'
#' @inheritParams recovery.net
#' @param newNodes Vector of nodal IDs for incoming nodes at the current time
#'        step.
#' @param curr.tab Current proportional distribution of all vertex attributes.
#'
#' @inherit recovery.net return
#'
#' @seealso [copy_nwattr_to_datattr()], [get_attr_prop()],
#'          [auto_update_attr()].
#' @keywords netUtils internal
#' @export
#'
auto_update_attr <- function(dat, newNodes, curr.tab) {

  rules <- get_control(dat, "attr.rules")
  active <- get_attr(dat, "active")
  t1.tab <- dat$run$t1.tab

  for (i in seq_along(curr.tab)) {
    vname <- names(curr.tab)[i]
    needs.updating <- ifelse(length(get_attr(dat, vname)) < length(active),
                             TRUE, FALSE)
    if (length(vname) > 0 && needs.updating == TRUE) {
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
      dat <- append_attr(dat, vname, nattr, length(nattr))
    }
  }

  return(dat)
}


#' @title Get Individual Degree from Network or Edgelist
#'
#' @description A fast method for querying the current degree of all individuals
#'              within a network.
#'
#' @param x Either an object of class `network` or `edgelist`
#'        generated from a network. If `x` is an edgelist, then it must
#'        contain an attribute for the total network size, `n`.
#'
#' @details
#' Individual-level data on the current degree of nodes within a network is
#' often useful for summary statistics. Given a `network` class object,
#' `net`, one way to look up the current degree is to get a summary of the
#' ERGM term, `sociality`, as in:
#' `summary(net ~ sociality(nodes = NULL))`. But that is computationally
#' inefficient for a number of reasons. This function provides a fast method for
#' generating the vector of degrees using a query of the edgelist. It is even
#' faster if the parameter `x` is already transformed into an edgelist.
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
#' @param x Object of class `dcm`, `netsim`, or `icm`.
#' @param at Time step at which to left-truncate the time series.
#'
#' @details
#' This function would be used when running a follow-up simulation from time
#' steps `b` to `c` after a burn-in period from time `a` to
#' `b`, where the final time window of interest for data analysis is
#' `b` to `c` only.
#'
#' @return The updated object of class `dcm`, `netsim`, or `icm`.
#'
#' @export
#'
#' @examples
#' # DCM example
#' param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.dcm(s.num = 500, i.num = 1)
#' control <- control.dcm(type = "SI", nsteps = 200)
#' mod1 <- dcm(param, init, control)
#' mod1$control$nsteps
#'
#' mod2 <- truncate_sim(mod1, at = 150)
#' mod2$control$nsteps
#'
#' # ICM example
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
#' init <- init.icm(s.num = 500, i.num = 1)
#' control <- control.icm(type = "SI", nsteps = 200, nsims = 1)
#' mod1 <- icm(param, init, control)
#' mod1$control$nsteps
#'
#' mod2 <- truncate_sim(mod1, at = 150)
#' mod2$control$nsteps
#'
truncate_sim <- function(x, at) {
  if (!inherits(x, c("dcm", "icm", "netsim"))) {
    stop("x must be an object of class dcm, icm, or netsim",
         call. = FALSE)
  }

  if (inherits(x, "dcm")) {
    row_start <- which(x$control$timesteps == at)
    if (length(row_start) == 0) {
      stop("Specified value of at is not in the control$timesteps vector",
           call. = FALSE)
    }
    rows <- row_start:nrow(x$epi[[1]])
    # epi
    x$epi <- lapply(x$epi, function(r) r[rows, , drop = FALSE])
    # control settings
    x$control$timesteps <- x$control$timesteps[rows]
    x$control$nsteps <- max(x$control$timesteps)
  } else {
    rows <- at:(x$control$nsteps)
    # epi
    x$epi <- lapply(x$epi, function(r) r[rows, , drop = FALSE])
    # control settings
    x$control$start <- 1
    x$control$nsteps <- max(seq_along(rows))
  }

  return(x)
}

#' Make a Lightweight Restart Point From a `netsim` Object with tergmLite
#'
#' Extract the elements required for re-initializing a `netsim` simulation from
#' a completed simulation. This function also resets the Unique IDs and Time
#' values to reduce the size of the simulation. This function only works for
#' simulations where `control$tergmLite = TRUE`
#'
#' @param sim_obj a `netsim` object from an ended `netsim` call.
#' @param sim_num the number of the simulation to extract from the `netsim`
#'        object (default = 1).
#' @param keep_steps The number of simulation steps to keep from the previous
#'        run. By default only keep one but more is possible if some
#'        back-history is wanted.
#' @param time_attrs a `character` vector containing the names of the attributes
#'        that are expressed in time-steps. These will be offsetted so the last
#'        step in the original simulation become the step 1 (default) in the new
#'        ones. If no such attributes exist, pass `c()`.
#'
#' @details
#' The restart point created always contains a single simulation and drops the
#' `attr.history`, the `raw.records` and `stats` from the initial simulation.
#'
#' The epi trackers, cumulative edgelists, transmission matrix and `nwstats` are
#' truncated to only contain the last `keep_steps` entries.
#'
#' Warning: the `time_attrs` argument is mandatory. Almost all simulation worth
#' restarting have such attributes (e.g. time.of.hiv.infection). If no such
#' argument exists, passing `c()` will allow the function to run while ensuring
#' that this was done on purpose.
#'
#' When restarting from the output of this function, it is suggested to express
#' the time steps in a relative maner in `control.net`:
#' ```
#' control.net(
#'   start = restart_point$control$nsteps + 1),
#'   nsteps = restart_point$control$nsteps + 1 + 104)
#' )
#' ```
#'
#' @examples
#' \dontrun{
#' # With  pre-existing `sim`, `param` and `init` object (see `netsim`)
#'
#' # List all attributes that store a time step
#' time_attrs <- c(
#'   "inf.time",
#'   "stage.time",
#'   "aids.time",
#'   "prep.start.last"
#' )
#' # Make a restart point a re-run for 10 more timesteps
#' x <- make_restart_point(sim, time_attrs, sim_num = 1, keep_steps = 1)
#' control <- control_msm(
#'   start = x$control$nsteps + 1,
#'   nsteps = x$control$nsteps + 1 + 10
#' )
#' sim <- netsim(x, param, init, control)
#' }
#'
#' @return a trimed `netsim` object with only one simulation that is ready to be
#'         used as a restart point.
#'
#' @export
make_restart_point <- function(sim_obj, time_attrs,
                               sim_num = 1, keep_steps = 1) {
  if (!inherits(sim_obj, c("netsim"))) {
    stop("`sim_obj` must be  an object of class `netsim`")
  }
  required_names <- c(
    "control", "param", "nwparam", "epi", "run", "coef.form", "num.nw"
  )
  missing_names <- setdiff(required_names, names(sim_obj))
  if (length(missing_names) > 0) {
    stop(
      "`sim_obj` is missing the following elements required for",
      " re-initialization: ", paste.and(missing_names)
    )
  }
  if (sim_num < 1 || sim_num > sim_obj$control$nsims) {
    stop("`sim_num` must be be >= 1 and <= `sim_obj$control$nsims`")
  }
  if (!sim_obj$control$tergmLite) {
    stop("Only `netsim` object with `tergmLite == TRUE` are supported")
  }

  # Select  the simulation of interest, that renames the selected sim: `sim1`
  x <- get_sims(sim_obj, sims = sim_num)
  n_steps <- x$control$nsteps
  run_ls <- x$run$sim1

  # Keep only the last `keep_steps` rows of each epi
  if (keep_steps < 1 || keep_steps > n_steps) {
    stop("`keep_steps` must be >= 1 and <= `sim_obj$control$nsteps`")
  }
  keep_rows <- (n_steps - keep_steps + 1):n_steps
  x$epi <- lapply(x$epi, function(r) r[keep_rows, , drop = FALSE])

  # If `nwstats` are saved, keep only the last rows
  if (x$control$save.nwstats) {
    x$stats$nwstats$sim1 <- lapply(
      x$stats$nwstats$sim1,
      function(d) d[keep_rows, , drop = FALSE]
    )
  }

  # Fix UIDs
  n_nodes <- length(run_ls$attr$active)
  uid_offset <- min(run_ls$attr$unique_id) - 1
  run_ls$attr$unique_id <- run_ls$attr$unique_id - uid_offset
  run_ls$last_unique_id <- run_ls$last_unique_id - uid_offset

  # Time correction
  time_offset <- n_steps - keep_steps
  x$control$start <- 1
  x$control$nsteps <- keep_steps

  # Time attributes - offset so last step is now `keep_steps`
  time_attrs <- union(c("entrTime", "exitTime"), time_attrs)
  missing_attrs <- setdiff(time_attrs, names(run_ls$attr))
  if (length(missing_attrs) > 0) {
    stop(
      "Some time attributes are not present in the attributes list:",
      paste.and(missing_names)
    )
  }
  run_ls$attr[time_attrs] <- lapply(
    run_ls$attr[time_attrs],
    function(v) v - time_offset
  )

  # Cumulative Edgelist - fix time and UIDs
  run_ls$el_cuml_cur <- lapply(
    run_ls$el_cuml_cur,
    function(el) {
      el$head <- el$head - uid_offset
      el$tail <- el$tail - uid_offset
      el$start <- el$start - time_offset
      el
    }
  )
  # For Historical one - truncate to 1 (only edges in the kept history)
  run_ls$el_cuml_hist <- lapply(
    run_ls$el_cuml_hist,
    function(el) {
      el$head <- el$head - uid_offset
      el$tail <- el$tail - uid_offset
      el$start <- el$start - time_offset
      el$stop <- el$stop - time_offset
      el[el$stop >= 1, , drop = FALSE]
    }
  )

  # the edgelist stores the name of the vertices. We don't use it with
  # `tergmLite` and it takes a lot of space
  run_ls$el <- lapply(run_ls$el, function(x) {
    attr(x, "vnames") <- NULL
    x
  })

  # If transmat was saved, trim it and offset the `at` column
  if (x$control$save.transmat) {
    tsmt <- x$stats$transmat$sim1
    tsmt$at <- tsmt$at - time_offset
    x$stats$transmat$sim1 <- tsmt[tsmt$at > 0, , drop = FALSE]
  }

  x$run$sim1 <- run_ls
  x$attr.history <- list()
  x$raw.records <- list()
  return(x)
}

#' @title Function to Reduce the Size of a `netest` Object
#'
#' @description Trims formula environments from the `netest` object.
#'              Optionally converts the `newnetwork` element of the
#'              `netest` object to a `networkLite` class, and removes
#'              the `fit` element (if present) from the `netest`
#'              object.
#'
#' @param object A `netest` class object.
#' @param as.networkLite If `TRUE`, converts `object$newnetwork`
#'        to a `networkLite`.
#' @param keep.fit If `FALSE`, removes the `object$fit` (if present)
#'        on the `netest` object.
#' @param keep Character vector of object names to keep in formula environments.
#'        By default, all objects are removed.
#'
#' @details
#' With larger, more complex network structures with epidemic models, it is
#' generally useful to reduce the memory footprint of the fitted TERGM model
#' object (estimated with [netest()]). This utility function removes
#' all but the bare essentials needed for simulating a network model with
#' [netsim()].
#'
#' The function always trims the environments of `object$constraints` and
#' `object$coef.diss$dissolution`.
#'
#' When both `edapprox = TRUE` and `nested.edapprox = TRUE` in the
#' `netest` call, also trims the environments of `object$formula`
#' and `object$formation`.
#'
#' When both `edapprox = TRUE` and `nested.edapprox = FALSE` in the
#' `netest` call, also trims the environments of `object$formula`,
#' `environment(object$formation)$formation`, and
#' `environment(object$formation)$dissolution`.
#'
#' When `edapprox = FALSE` in the `netest` call, also trims the
#' environments of `object$formation`,
#' `environment(object$formula)$formation` and
#' `environment(object$formula)$dissolution`.
#'
#' By default all objects are removed from these trimmed environments. Specific
#' objects may be retained by passing their names as the `keep` argument.
#' For the output of `trim_netest` to be usable in [netsim()]
#' simulation, any objects referenced in the formulas should be included in the
#' `keep` argument.
#'
#' If `as.networkLite = TRUE`, converts `object$newnetwork` to a
#' `networkLite` object. If `keep.fit = FALSE`, removes `fit` (if
#' present) from `object`.
#'
#' @return
#' A `netest` object with formula environments trimmed, optionally with the
#' `newnetwork` element converted to a `networkLite` and the
#' `fit` element removed.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(n = 100)
#' formation <- ~edges + concurrent
#' target.stats <- c(50, 25)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
#' est <- netest(nw, formation, target.stats, coef.diss,
#'               set.control.ergm = control.ergm(MCMC.burnin = 1e5,
#'                                               MCMC.interval = 1000))
#' print(object.size(est), units = "KB")
#'
#' est.small <- trim_netest(est)
#' print(object.size(est.small), units = "KB")
#'
trim_netest <- function(object, as.networkLite = TRUE, keep.fit = FALSE,
                        keep = character(0)) {
  if (object$edapprox == TRUE) {
    object$formula <- trim_env(object$formula, keep = keep)
    if (object$nested.edapprox == TRUE) {
      object$formation <- trim_env(object$formation, keep = keep)
    } else {
      # trim environments for formation and dissolution inside formation
      environment(object$formation)$formation <-
        trim_env(environment(object$formation)$formation, keep = keep)
      environment(object$formation)$dissolution <-
        trim_env(environment(object$formation)$dissolution, keep = keep)
    }
  } else {
    object$formation <- trim_env(object$formation, keep = keep)
    # trim environments for formation and dissolution inside formula
    environment(object$formula)$formation <-
      trim_env(environment(object$formula)$formation, keep = keep)
    environment(object$formula)$dissolution <-
      trim_env(environment(object$formula)$dissolution, keep = keep)
  }

  object$coef.diss$dissolution <- trim_env(object$coef.diss$dissolution, keep = keep)
  object$constraints <- trim_env(object$constraints, keep = keep)

  if (keep.fit == FALSE) {
    object$fit <- NULL
  }

  if (as.networkLite == TRUE) {
    object$newnetwork <- as.networkLite(object$newnetwork)
  }

  return(object)
}