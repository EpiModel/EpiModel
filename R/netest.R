
#' @title Dynamic Network Model Estimation
#'
#' @description Estimates statistical network models using the exponential
#'              random graph modeling (ERGM) framework with extensions for
#'              dynamic/temporal models (STERGM).
#'
#' @param nw An object of class \code{network}.
#' @param formation Right-hand sided STERGM formation formula in the form
#'        \code{~edges + ...}, where \code{...} are additional network
#'        statistics.
#' @param target.stats Vector of target statistics for the formation model, with
#'        one number for each network statistic in the model.
#' @param coef.diss An object of class \code{disscoef} output from the
#'        \code{\link{dissolution_coefs}} function.
#' @param constraints Right-hand sided formula specifying constraints for the
#'        modeled network, in the form \code{~...}, where \code{...} are
#'        constraint terms. By default, no constraints are set.
#' @param coef.form Vector of coefficients for the offset terms in the formation
#'        formula.
#' @param edapprox If \code{TRUE}, use the indirect edges dissolution
#'        approximation  method for the dynamic model fit, otherwise use the
#'        more time-intensive full STERGM estimation (see details).
#' @param set.control.ergm Control arguments passed to \code{simulate.ergm} (see
#'        details).
#' @param set.control.stergm Control arguments passed to \code{simulate.stergm}
#'        (see details).
#' @param verbose Print model fitting progress to console.
#'
#' @details
#' \code{netest} is a wrapper function for the \code{ergm} and \code{stergm}
#' functions that estimate static and dynamic network models, respectively.
#' Network model estimation is the first step in simulating a stochastic network
#' epidemic model in \code{EpiModel}. The output from \code{netest} is a
#' necessary input for running the epidemic simulations in \code{\link{netsim}}.
#' With a fitted network model, one should always first proceed to model
#' diagnostics, available through the \code{\link{netdx}} function, to check
#' model fit. A detailed description of fitting these models, along with
#' examples, may be found in the
#' \href{http://www.epimodel.org/tut.html}{Basic Network Models} tutorials.
#'
#' @section Edges Dissolution Approximation:
#' The edges dissolution approximation method is described in Carnegie et al.
#' This approximation requires that the dissolution coefficients are known, that
#' the formation model is being fit to cross-sectional data conditional on those
#' dissolution coefficients, and that the terms in the dissolution model are a
#' subset of those in the formation model. Under certain additional conditions,
#' the formation coefficients of a STERGM model are approximately equal to the
#' coefficients of that same model fit to the observed cross-sectional data as
#' an ERGM, minus the corresponding coefficients in the dissolution model. The
#' approximation thus estimates this ERGM (which is typically much faster than
#' estimating a STERGM) and subtracts the dissolution coefficients.
#'
#' The conditions under which this approximation best hold are when there are
#' few relational changes from one time step to another; i.e. when either
#' average relational durations are long, or density is low, or both.
#' Conveniently, these are the same conditions under which STERGM estimation is
#' slowest.  Note that the same approximation is also used to obtain starting
#' values for the STERGM estimate when the latter is being conducted.  The
#' estimation does not allow for calculation of standard errors, p-values, or
#' likelihood for the formation model; thus, this approach is of most use when
#' the main goal of estimation is to drive dynamic network simulations rather
#' than to conduct inference on the formation model. The user is strongly
#' encouraged to examine the behavior of the resulting simulations to confirm
#' that the approximation is adequate for their purposes. For an example, see
#' the vignette for the package \code{tergm}.
#' 
#' It has recently been found that subtracting a modified version of the 
#' dissolution coefficients from the formation coefficients provides a more 
#' principled approximation, and this is now the form of the approximation 
#' applied by \code{netest}.  (The modified values subtracted from the formation
#' coefficients are equivalent to the (crude) dissolution coefficients with 
#' their target durations increased by 1.)
#' 
#' @section Control Arguments:
#' The \code{ergm} and \code{stergm} functions allow control settings for the
#' model fitting process. When fitting a STERGM directly (setting
#' \code{edapprox} to \code{FALSE}), control parameters may be passed to the
#' \code{stergm} function with the \code{set.control.stergm} argument in
#' \code{netest}. The controls should be input through the
#' \code{control.stergm()} function, with the available parameters listed in the
#' \code{\link{control.stergm}} help page in the \code{tergm} package.
#'
#' When fitting a STERGM indirectly (setting \code{edapprox} to \code{TRUE}),
#' control settings may be passed to the \code{ergm} function using
#' \code{set.control.ergm} in \code{netest}. The controls should be input
#' through the \code{control.ergm()} function, with the available parameters
#' listed in the \code{\link[ergm:control.ergm]{control.ergm}} help page in the
#' \code{ergm} package. An example is below.
#'
#' @references
#' Krivitsky PN, Handcock MS. "A separable model for dynamic networks." JRSS(B).
#' 2014; 76.1:29-46.
#'
#' Carnegie NB, Krivitsky PN, Hunter DR, Goodreau SM. An approximation method
#' for improving dynamic network model fitting. Journal of Computational and
#' Graphical Statistics. 2014; 24(2): 502-519.
#'
#' Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical
#' Modeling of Infectious Disease over Networks. Journal of Statistical
#' Software. 2018; 84(8): 1-47.
#'
#' @keywords model
#' @seealso Use \code{\link{netdx}} to diagnose the fitted network model, and
#'          \code{\link{netsim}} to simulate epidemic spread over a simulated
#'          dynamic network consistent with the model fit.
#'
#' @export
#'
#' @examples
#' # Initialize a network of 100 nodes
#' nw <- network_initialize(n = 100)
#'
#' # Set formation formula
#' formation <- ~edges + concurrent
#'
#' # Set target statistics for formation
#' target.stats <- c(50, 25)
#'
#' # Obtain the offset coefficients
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
#'
#' # Estimate the STERGM using the edges dissolution approximation
#' est <- netest(nw, formation, target.stats, coef.diss,
#'               set.control.ergm = control.ergm(MCMC.burnin = 1e5,
#'                                               MCMC.interval = 1000))
#' est
#'
#' # To estimate the STERGM directly, use edapprox = FALSE
#' # est2 <- netest(nw, formation, target.stats, coef.diss, edapprox = FALSE)
#'
netest <- function(nw, formation, target.stats, coef.diss, constraints,
                   coef.form = NULL, edapprox = TRUE,
                   set.control.ergm, set.control.stergm,
                   verbose = FALSE) {

  if (missing(constraints)) {
    constraints	<- trim_env(~.)
  }

  if (class(coef.diss) != "disscoef") {
    stop("dissolution must be of input through dissolution_coefs function",
         call. = FALSE)
  }
  dissolution <- coef.diss$dissolution
  if (coef.diss$coef.crude[1] == -Inf) {
    is.tergm <- FALSE
  } else {
    is.tergm <- TRUE
  }

  if (is.tergm == TRUE) {
    diss_check(formation, dissolution)
  }

  if (edapprox == FALSE) {

    if (missing(set.control.stergm)) {
      set.control.stergm <- control.stergm()
    }

    fit <- stergm(nw,
                  formation = formation,
                  dissolution = dissolution,
                  targets = "formation",
                  target.stats = target.stats,
                  offset.coef.form = coef.form,
                  offset.coef.diss = coef.diss$coef.crude,
                  constraints = constraints,
                  estimate = "EGMME",
                  eval.loglik = FALSE,
                  control = set.control.stergm,
                  verbose = verbose)

    coef.form <- fit$formation.fit

    out <- list()
    out$fit <- fit
    out$formation <- formation
    out$target.stats <- target.stats
    out$target.stats.names <-
      names(fit$formation.fit$coef)[!fit$formation.fit$offset]
    out$coef.form <- coef.form$coef
    out$dissolution <- dissolution
    out$coef.diss <- coef.diss
    out$constraints <- constraints
    out$edapprox <- edapprox

  } else {

    if (missing(set.control.ergm)) {
      set.control.ergm <- control.ergm()
    }

    formation.nw <- nonsimp_update.formula(formation, nw ~ ., from.new="nw")

    fit <- ergm(formation.nw,
                target.stats = target.stats,
                constraints = constraints,
                offset.coef = coef.form,
                eval.loglik = FALSE,
                control = set.control.ergm,
                verbose = verbose)

    coef.form <- fit$coef
    coef.form.crude <- coef.form
    if (is.tergm == TRUE) {
      l.cfc <- length(coef.diss$coef.form.corr)
      coef.form[1:l.cfc] <- coef.form[1:l.cfc] - coef.diss$coef.form.corr
    }

    # Reduce size of output object
    fit$initialfit <- NULL
    fit$constrained <- NULL
    environment(fit$sample.obs) <- NULL
    environment(fit$reference) <- NULL

    out <- list()
    out$fit <- fit
    out$formation <- formation
    out$target.stats <- target.stats
    if (length(names(fit$coef)) == length(target.stats)) {
      out$target.stats.names <- names(fit$coef)
    } else {
      out$target.stats.names <- names(fit$coef)[!fit$offset]
    }
    out$coef.form <- coef.form
    out$coef.form.crude <- coef.form.crude
    out$coef.diss <- coef.diss
    out$constraints <- constraints
    out$edapprox <- edapprox

  }

  class(out) <- "netest"
  return(out)
}


diss_check <- function(formation, dissolution) {

  # Split formulas into separate terms
  form.terms <- strsplit(as.character(formation)[2], "[+]")[[1]]
  diss.terms <- strsplit(as.character(dissolution)[2], "[+]")[[1]]

  # Remove whitespace
  form.terms <- gsub("\\s", "", form.terms)
  diss.terms <- gsub("\\s", "", diss.terms)

  offpos.f <- grep("offset(", form.terms, fixed = TRUE)
  form.terms[offpos.f] <- substr(form.terms[offpos.f], nchar("offset(") + 1,
                                 nchar(form.terms[offpos.f]) - 1)
  offpos.d <- grep("offset(", diss.terms, fixed = TRUE)
  diss.terms[offpos.d] <- substr(diss.terms[offpos.d], nchar("offset(") + 1,
                                 nchar(diss.terms[offpos.d]) - 1)

  argpos.f <- regexpr("\\(", form.terms)
  argpos.d <- regexpr("\\(", diss.terms)

  # Matrix with terms in row 1, args in row 2
  form.terms <- vapply(regmatches(form.terms, argpos.f, invert = TRUE),
                       function(x) {
                         if (length(x) < 2) {
                           x <- c(x, "")
                         } else {
                           x[2] <- substr(x[2], 1, nchar(x[2]) - 1)
                         }
                         x
                       },
                       c(term = "", args = ""))
  diss.terms <- vapply(regmatches(diss.terms, argpos.d, invert = TRUE),
                       function(x) {
                         if (length(x) < 2) {
                           x <- c(x, "")
                         } else {
                           x[2] <- substr(x[2], 1, nchar(x[2]) - 1)
                         }
                         x
                       },
                       c(term = "", args = ""))


  matchpos <- match(diss.terms[1, ], form.terms[1, ])

  if (any(is.na(matchpos))) {
    stop("Dissolution model is not a subset of formation model.", call. = FALSE)
  }
  if (!all(diss.terms[1, ] %in% c("edges", "nodemix",
                                  "nodematch", "nodefactor"))) {
    stop("The only allowed dissolution terms are edges, nodemix,
         nodematch and ", "nodefactor", call. = FALSE)
  }
  if (any(matchpos != 1:ncol(diss.terms))) {
    stop("Order of terms in the dissolution model does not correspond to the ",
         "formation model.", call. = FALSE)
  }
  if (any(diss.terms[2, ] != form.terms[2, 1:ncol(diss.terms)])) {
    stop("Term options for one or more terms in dissolution model do not ",
         "match the options in the formation model.", call. = FALSE)
  }

}


#' @title Adjust Dissolution Component of Network Model Fit
#'
#' @description Adjusts the dissolution component of an dynamic ERGM fit using
#'              the \code{netest} function with the edges dissolution
#'              approximation method.
#'
#' @param old.netest An object of class \code{netest}, from the
#'        \code{\link{netest}} function.
#' @param new.coef.diss An object of class \code{disscoef}, from the
#'        \code{\link{dissolution_coefs}} function.
#'
#' @details
#' Fitting an ERGM is a computationally intensive process when the model
#' includes dyadic dependent terms. With the edges dissolution approximation
#' method of Carnegie et al, the coefficients for a temporal ERGM are
#' approximated by fitting a static ERGM and adjusting the formation
#' coefficients to account for edge dissolution. This function provides a very
#' efficient method to adjust the coefficients of that model when one wants to
#' use a different dissolution model; a typical use case may be to fit several
#' different models with different average edge durations as targets. The
#' example below exhibits that case.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' nw <- network_initialize(n = 1000)
#'
#' # Two dissolutions: an average duration of 300 versus 200
#' diss.300 <- dissolution_coefs(~offset(edges), 300, 0.001)
#' diss.200 <- dissolution_coefs(~offset(edges), 200, 0.001)
#'
#' # Fit the two reference models
#' est300 <- netest(nw = nw,
#'                 formation = ~edges,
#'                 target.stats = c(500),
#'                 coef.diss = diss.300)
#'
#' est200 <- netest(nw = nw,
#'                 formation = ~edges,
#'                 target.stats = c(500),
#'                 coef.diss = diss.200)
#'
#' # Alternatively, update the 300 model with the 200 coefficients
#' est200.compare <- update_dissolution(est300, diss.200)
#'
#' identical(est200$coef.form, est200.compare$coef.form)
#'}
#'
update_dissolution <- function(old.netest, new.coef.diss) {

  if (class(old.netest) != "netest") {
    stop("old.netest must be an object of class netest", call. = FALSE)
  }
  if (class(new.coef.diss) != "disscoef") {
    stop("new.coef.diss must be an object of class disscoef", call. = FALSE)
  }
  if (old.netest$edapprox != TRUE) {
    stop("Edges dissolution approximation must be used for this adjustment",
         call. = FALSE)
  }

  out <- old.netest

  ## remove old correction
  l.cd.o <- length(out$coef.diss$coef.form.corr)
  out$coef.form[1:l.cd.o] <- out$coef.form[1:l.cd.o] + out$coef.diss$coef.form.corr
  ## apply new correction (may differ in length from the old correction)
  l.cd.n <- length(new.coef.diss$coef.form.corr)
  out$coef.form[1:l.cd.n] <- out$coef.form[1:l.cd.n] - new.coef.diss$coef.form.corr

  out$coef.diss <- new.coef.diss

  return(out)
}

