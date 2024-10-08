
#' @title Dynamic Network Model Estimation
#'
#' @description Estimates statistical network models using the exponential
#'              random graph modeling (ERGM) framework with extensions for
#'              dynamic/temporal models (STERGM).
#'
#' @param nw An object of class `network` or `egor`, with the latter
#'        indicating an `ergm.ego` fit.
#' @param formation Right-hand sided STERGM formation formula in the form
#'        `~edges + ...`, where `...` are additional network statistics.
#' @param target.stats Vector of target statistics for the formation model, with
#'        one number for each network statistic in the model. Ignored if
#'        fitting via `ergm.ego`.
#' @param coef.diss An object of class `disscoef` output from the
#'        [`dissolution_coefs`] function.
#' @param constraints Right-hand sided formula specifying constraints for the
#'        modeled network, in the form `~...`, where `...` are
#'        constraint terms. By default, no constraints are set.
#' @param coef.form Vector of coefficients for the offset terms in the formation
#'        formula.
#' @param edapprox If `TRUE`, use the indirect edges dissolution
#'        approximation  method for the dynamic model fit, otherwise use the
#'        more time-intensive full STERGM estimation (see details).  For
#'        `nw` of class `egor`, only `edapprox = TRUE` is supported.
#' @param set.control.ergm Control arguments passed to `ergm` (see details).
#' @param set.control.ergm.ego Control arguments passed to `ergm.ego` (see
#'        details).
#' @param set.control.tergm Control arguments passed to `tergm` (see details).
#' @param verbose If `TRUE`, print model fitting progress to console.
#' @param nested.edapprox Logical. If `edapprox = TRUE` the dissolution
#'        model is an initial segment of the formation model (see details).
#' @param ... Additional arguments passed to other functions.
#'
#' @details
#' `netest` is a wrapper function for the `ergm`, `ergm.ego`, and `tergm`
#' functions that estimate static and dynamic network models. Network model
#' estimation is the first step in simulating a stochastic network epidemic model
#' in `EpiModel`. The output from `netest` is a necessary input for running the
#' epidemic simulations in [`netsim`]. With a fitted network model, one should
#' always first proceed to model diagnostics, available through the [`netdx`]
#' function, to check model fit. A detailed description of fitting these
#' models, along with examples, may be found in the
#' [Network Modeling for Epidemics](https://epimodel.github.io/sismid/)
#' tutorials.
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
#' slowest. Note that the same approximation is also used to obtain starting
#' values for the STERGM estimate when the latter is being conducted. The
#' estimation does not allow for calculation of standard errors, p-values, or
#' likelihood for the formation model; thus, this approach is of most use when
#' the main goal of estimation is to drive dynamic network simulations rather
#' than to conduct inference on the formation model. The user is strongly
#' encouraged to examine the behavior of the resulting simulations to confirm
#' that the approximation is adequate for their purposes. For an example, see
#' the vignette for the package `tergm`.
#'
#' It has recently been found that subtracting a modified version of the
#' dissolution coefficients from the formation coefficients provides a more
#' principled approximation, and this is now the form of the approximation
#' applied by `netest`. The modified values subtracted from the formation
#' coefficients are equivalent to the (crude) dissolution coefficients with
#' their target durations increased by 1. The `nested.edapprox` argument
#' toggles whether to implement this modified version by appending the
#' dissolution terms to the formation model and appending the relevant values to
#' the vector of formation model coefficients (value = `FALSE`), whereas
#' the standard version subtracts the relevant values from the initial formation
#' model coefficients (value = `TRUE`).
#'
#' @section Control Arguments:
#' The `ergm`, `ergm.ego`, and `tergm` functions allow control settings for the
#' model fitting process. When fitting a STERGM directly (setting
#' `edapprox` to `FALSE`), control parameters may be passed to the
#' `tergm` function with the `set.control.tergm` argument in `netest`.
#' The controls should be input through the `control.tergm()` function,
#' with the available parameters listed in the [`tergm::control.tergm`] help
#' page in the `tergm` package.
#'
#' When fitting a STERGM indirectly (setting `edapprox` to `TRUE`), control
#' settings may be passed to the `ergm` function using `set.control.ergm`,
#' or to the `ergm.ego` function using `set.control.ergm.ego`. The controls
#' should be input through the `control.ergm()` and `control.ergm.ego()`
#' functions, respectively, with the available parameters listed in the
#' [`ergm::control.ergm`] help page in the `ergm` package and the
#' [`ergm.ego::control.ergm.ego`] help page in the `ergm.ego`
#' package. An example is below.
#'
#' @return A fitted network model object of class `netest`.
#'
#' @references
#' Krivitsky PN, Handcock MS. "A Separable Model for Dynamic Networks." JRSS(B).
#' 2014; 76.1: 29-46.
#'
#' Carnegie NB, Krivitsky PN, Hunter DR, Goodreau SM. An Approximation Method
#' for Improving Dynamic Network Model Fitting. Journal of Computational and
#' Graphical Statistics. 2014; 24(2): 502-519.
#'
#' Jenness SM, Goodreau SM and Morris M. EpiModel: An R Package for Mathematical
#' Modeling of Infectious Disease over Networks. Journal of Statistical
#' Software. 2018; 84(8): 1-47.
#'
#' @keywords model
#' @seealso Use [`netdx`] to diagnose the fitted network model, and [`netsim`]
#'          to simulate epidemic spread over a simulated dynamic network
#'          consistent with the model fit.
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
                   set.control.ergm = control.ergm(),
                   set.control.tergm = control.tergm(),
                   set.control.ergm.ego = NULL,
                   verbose = FALSE, nested.edapprox = TRUE, ...) {

  if (missing(constraints)) {
    constraints	<- trim_env(~.)
  }

  if (!inherits(coef.diss, "disscoef")) {
    stop("dissolution must be of input through dissolution_coefs function",
         call. = FALSE)
  }
  dissolution <- coef.diss$dissolution
  if (coef.diss$duration[1] == 1) {
    is.tergm <- FALSE
  } else {
    is.tergm <- TRUE
  }

  if (is.tergm == TRUE && nested.edapprox == TRUE) {
    diss_check(formation, dissolution)
  }

  if (edapprox == FALSE) {

    fit <- tergm(~Form(formation) + Persist(dissolution),
                 basis = nw,
                 targets = "formation",
                 target.stats = target.stats,
                 offset.coef = c(coef.form, coef.diss$coef.crude),
                 constraints = constraints,
                 estimate = "EGMME",
                 eval.loglik = FALSE,
                 control = set.control.tergm,
                 verbose = verbose)

    coef.form <- fit # there is no longer a separate formation fit
    which_form <- which(grepl("^Form~", names(coef(fit))) |
                          grepl("^offset\\(Form~", names(coef(fit))))
    form_names <- names(coef(fit)[which_form])[!fit$offset[which_form]]

    out <- list()
    out$formation <- formation
    out$target.stats <- target.stats
    out$target.stats.names <-
      substr(form_names, 6, nchar(form_names))
    out$coef.form <- coef(coef.form)[which_form]
    out$dissolution <- dissolution
    out$coef.diss <- coef.diss
    out$constraints <- constraints
    out$edapprox <- edapprox
    # convert ergm_state to network
    out$newnetwork <- as.network(fit$newnetwork)
    delete.network.attribute(out$newnetwork, "time")
    delete.network.attribute(out$newnetwork, "lasttoggle")
    out$formula <- trim_env(fit$formula, keep = c("formation", "dissolution"))

  } else {

    if (inherits(nw, "egor")) {
      # ergm.ego case
      if (system.file(package = "ergm.ego") == "")
        stop("The `ergm.ego` package is required to estimate from `egor` objects.\n ")
      if (is.null(set.control.ergm.ego))
        set.control.ergm.ego <- ergm.ego::control.ergm.ego()

      fit <- ergm.ego::ergm.ego(
        formation,
        basis = nw,
        popsize = 0,
        constraints = constraints,
        offset.coef = coef.form,
        control = set.control.ergm.ego,
        verbose = verbose
      )
      target.stats <- fit$m

    } else {
      # ergm case
      fit <- ergm(formation,
                  basis = nw,
                  target.stats = target.stats,
                  constraints = constraints,
                  offset.coef = coef.form,
                  eval.loglik = FALSE,
                  control = set.control.ergm,
                  verbose = verbose)
    }

    coef.form <- coef(fit)
    coef.form.crude <- coef.form
    if (is.tergm == TRUE) {
      if (nested.edapprox == TRUE) {
        l.cfc <- length(coef.diss$coef.form.corr)
        coef.form[1:l.cfc] <- coef.form[1:l.cfc] - coef.diss$coef.form.corr
      } else {
        ## implement the edapprox by appending the dissolution model to the
        ## formation model and appending the relevant values to the vector of
        ## formation model coefficients
        formation <- trim_env(~Passthrough(formation) + Passthrough(dissolution),
                              keep = c("formation", "dissolution"))
        coef.form <- c(coef.form, -coef.diss$coef.form.corr)
      }
    }

    out <- list()
    out$formation <- formation
    out$target.stats <- target.stats
    ## subselect coef names for targeted statistics, including extremal targets
    out$target.stats.names <- names(coef(fit))[!fit$offset | (fit$drop != 0)]
    out$coef.form <- coef.form
    out$coef.form.crude <- coef.form.crude
    out$coef.diss <- coef.diss
    out$constraints <- constraints
    out$edapprox <- edapprox
    out$nested.edapprox <- nested.edapprox
    out$newnetwork <- NVL(fit$newnetwork, fit$network)
    out$formula <- fit$formula
  }

  out$summary <- summary(fit, ...)
  out$fit <- fit

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
                                  "nodematch"))) {
    stop("The only allowed dissolution terms are edges, nodemix,
         and nodematch", call. = FALSE)
  }
  if (any(matchpos != seq_len(ncol(diss.terms)))) {
    stop("Order of terms in the dissolution model does not correspond to the ",
         "formation model.", call. = FALSE)
  }
  if (any(diss.terms[2, ] != form.terms[2, seq_len(ncol(diss.terms))])) {
    stop("Term options for one or more terms in dissolution model do not ",
         "match the options in the formation model.", call. = FALSE)
  }

}


#' @title Adjust Dissolution Component of Network Model Fit
#'
#' @description Adjusts the dissolution component of a dynamic ERGM fit using
#'              the \code{\link{netest}} function with the edges dissolution
#'              approximation method.
#'
#' @param old.netest An object of class \code{netest}, from the
#'        \code{\link{netest}} function.
#' @param new.coef.diss An object of class \code{disscoef}, from the
#'        \code{\link{dissolution_coefs}} function.
#' @param nested.edapprox Logical. If \code{edapprox = TRUE} the dissolution
#'        model is an initial segment of the formation model (see details in
#'        \code{\link{netest}}).
#'
#' @details
#' Fitting an ERGM is a computationally intensive process when the model
#' includes dyad dependent terms. With the edges dissolution approximation
#' method of Carnegie et al, the coefficients for a temporal ERGM are
#' approximated by fitting a static ERGM and adjusting the formation
#' coefficients to account for edge dissolution. This function provides a very
#' efficient method to adjust the coefficients of that model when one wants to
#' use a different dissolution model; a typical use case may be to fit several
#' different models with different average edge durations as targets. The
#' example below exhibits that case.
#'
#' @return An updated network model object of class \code{netest}.
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
update_dissolution <- function(old.netest, new.coef.diss,
                               nested.edapprox = TRUE) {

  if (!inherits(old.netest, "netest")) {
    stop("old.netest must be an object of class netest", call. = FALSE)
  }
  if (!inherits(new.coef.diss, "disscoef")) {
    stop("new.coef.diss must be an object of class disscoef", call. = FALSE)
  }
  if (old.netest$edapprox != TRUE) {
    stop("Edges dissolution approximation must be used for this adjustment",
         call. = FALSE)
  }

  out <- old.netest

  ## if an old correction was applied...
  if (old.netest$coef.diss$coef.crude[1] != -Inf) {
    ## remove the old correction
    if (out$nested.edapprox == TRUE) {
      ## adjust the formation model coefficients to remove the old edapprox
      l.cd.o <- length(out$coef.diss$coef.form.corr)
      out$coef.form[1:l.cd.o] <- out$coef.form[1:l.cd.o] +
        out$coef.diss$coef.form.corr
    } else {
      ## remove the part of the formation model and coefficient vector
      ## corresponding to the old edapprox
      out$formation <- eval(quote(formation),
                            envir = environment(out$formation))
      out$coef.form <-
        out$coef.form[seq_len(length(out$coef.form) -
                                length(out$coef.diss$coef.form.corr))]
    }
  }

  ## if a new correction should be applied...
  if (new.coef.diss$coef.crude[1] != -Inf) {
    ## apply the new correction
    if (nested.edapprox == TRUE) {
      ## check that the new dissolution model is an initial segment of the
      ## formation model
      diss_check(out$formation, new.coef.diss$dissolution)

      ## implement new edapprox by adjusting the formation model coefficients
      l.cd.n <- length(new.coef.diss$coef.form.corr)
      out$coef.form[1:l.cd.n] <- out$coef.form[1:l.cd.n] -
        new.coef.diss$coef.form.corr
    } else {
      ## implement new edapprox by appending the new dissolution model to the
      ## formation model and appending the relevant values to the vector of
      ## formation model coefficients
      # nolint start
      formation <- out$formation
      dissolution <- new.coef.diss$dissolution
      # nolint end
      out$formation <- trim_env(~Passthrough(formation) + Passthrough(dissolution),
                                keep = c("formation", "dissolution"))
      out$coef.form <- c(out$coef.form, -new.coef.diss$coef.form.corr)
    }
  }

  out$coef.diss <- new.coef.diss
  out$nested.edapprox <- nested.edapprox

  return(out)
}
