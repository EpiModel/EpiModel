
#' @title Merge Data across Stochastic Individual Contact Model Simulations
#'
#' @description Merges epidemiological data from two independent simulations of
#'              stochastic individual contact models from [icm()].
#'
#' @param x An `EpiModel` object of class [icm()].
#' @param y Another `EpiModel` object of class [icm()], with the
#'        identical model parameterization as `x`.
#' @param ...  Additional merge arguments (not used).
#'
#' @details
#' This merge function combines the results of two independent simulations of
#' [icm()] class models, simulated under separate function calls. The
#' model parameterization between the two calls must be exactly the same, except
#' for the number of simulations in each call. This allows for manual
#' parallelization of model simulations.
#'
#' This merge function does not work the same as the default merge, which allows
#' for a combined object where the structure differs between the input elements.
#' Instead, the function checks that objects are identical in model
#' parameterization in every respect (except number of simulations) and binds
#' the results.
#'
#' @return An `EpiModel` object of class [icm()] containing the
#'         data from both `x` and `y`.
#'
#' @method merge icm
#' @keywords extract
#' @export
#'
#' @examples
#' param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
#' init <- init.icm(s.num = 1000, i.num = 100)
#' control <- control.icm(type = "SI", nsteps = 10,
#'                        nsims = 3, verbose = FALSE)
#' x <- icm(param, init, control)
#'
#' control <- control.icm(type = "SI", nsteps = 10,
#'                        nsims = 1, verbose = FALSE)
#' y <- icm(param, init, control)
#'
#' z <- merge(x, y)
#'
#' # Examine separate and merged data
#' as.data.frame(x)
#' as.data.frame(y)
#' as.data.frame(z)
#'
merge.icm <- function(x, y, ...) {

  ## Check structure
  if (length(x) != length(y) || !identical(names(x), names(y))) {
    stop("x and y have different structure")
  }
  if (x$control$nsims > 1 && y$control$nsims > 1 &&
        !identical(sapply(x, class), sapply(y, class))) {
    stop("x and y have different structure")
  }

  ## Check params
  check1 <- identical(x$param, y$param)
  check2 <- identical(x$control[-which(names(x$control) == "nsims")],
                      y$control[-which(names(y$control) == "nsims")])

  if (check1 == FALSE) {
    stop("x and y have different parameters")
  }
  if (check2 == FALSE) {
    stop("x and y have different controls")
  }


  z <- x
  new.range <- (x$control$nsims + 1):(x$control$nsims + y$control$nsims)

  # Merge data
  for (i in seq_along(x$epi)) {
    if (x$control$nsims == 1) {
      x$epi[[i]] <- data.frame(x$epi[[i]])
    }
    if (y$control$nsims == 1) {
      y$epi[[i]] <- data.frame(y$epi[[i]])
    }
    z$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(z$epi[[i]])[new.range] <- paste0("sim", new.range)
  }

  z$control$nsims <- max(new.range)

  return(z)
}


#' @title Merge Model Simulations across netsim Objects
#'
#' @description Merges epidemiological data from two independent simulations of
#'              stochastic network models from `netsim`.
#'
#' @param x An `EpiModel` object of class [netsim()].
#' @param y Another `EpiModel` object of class [netsim()],
#'        with the identical model parameterization as `x`.
#' @param keep.transmat If `TRUE`, keep the transmission matrices from the
#'        original `x` and `y` elements. Note: transmission matrices
#'        only saved when (`save.transmat == TRUE`).
#' @param keep.network If `TRUE`, keep the `networkDynamic` objects
#'        from the original `x` and `y` elements. Note: network
#'        only saved when (`tergmLite == FALSE`).
#' @param keep.nwstats If `TRUE`, keep the network statistics (as set by
#'        the `nwstats.formula` parameter in `control.net`) from
#'        the original `x` and `y` elements.
#' @param keep.other If `TRUE`, keep the other simulation elements (as set
#'        by the `save.other` parameter in `control.net`) from the
#'        original `x` and `y` elements.
#' @param param.error If `TRUE`, if `x` and `y` have different
#'        params (in [param.net()]) or controls (passed in
#'        [control.net()]) an error will prevent the merge. Use
#'        `FALSE` to override that check.
#' @param keep.diss.stats If `TRUE`, keep `diss.stats` from the
#'        original `x` and `y` objects.
#' @param ...  Additional merge arguments (not currently used).
#'
#' @details
#' This merge function combines the results of two independent simulations of
#' [netsim()] class models, simulated under separate function calls.
#' The model parameterization between the two calls must be exactly the same,
#' except for the number of simulations in each call. This allows for manual
#' parallelization of model simulations.
#'
#' This merge function does not work the same as the default merge, which allows
#' for a combined object where the structure differs between the input elements.
#' Instead, the function checks that objects are identical in model
#' parameterization in every respect (except number of simulations) and binds
#' the results.
#'
#' @return An `EpiModel` object of class [netsim()] containing
#'         the data from both `x` and `y`.
#'
#' @method merge netsim
#' @keywords extract
#' @export
#'
#' @examples
#' \donttest{
#' # Network model
#' nw <- network_initialize(n = 100)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
#' est <- netest(nw, formation = ~edges, target.stats = 25,
#'               coef.diss = coef.diss, verbose = FALSE)
#'
#' # Epidemic models
#' param <- param.net(inf.prob = 1)
#' init <- init.net(i.num = 1)
#' control <- control.net(type = "SI", nsteps = 20, nsims = 2,
#'                        save.nwstats = TRUE,
#'                        nwstats.formula = ~edges + degree(0),
#'                        verbose = FALSE)
#' x <- netsim(est, param, init, control)
#' y <- netsim(est, param, init, control)
#'
#' # Merging
#' z <- merge(x, y)
#'
#' # Examine separate and merged data
#' as.data.frame(x)
#' as.data.frame(y)
#' as.data.frame(z)
#' }
#'
merge.netsim <- function(x, y, keep.transmat = TRUE, keep.network = TRUE,
                         keep.nwstats = TRUE, keep.other = TRUE,
                         param.error = TRUE, keep.diss.stats = TRUE, ...) {

  ## Check structure
  if (length(x) != length(y) || !identical(names(x), names(y))) {
    stop("x and y have different structure")
  }
  x$control$nsims <- as.integer(x$control$nsims)
  y$control$nsims <- as.integer(y$control$nsims)
  if (x$control$nsims > 1 && y$control$nsims > 1 &&
        !all(sapply(x, function(i) class(i)[1]) ==
               sapply(y, function(i) class(i)[1]))) {
    stop("x and y have different structure")
  }

  ## Check params
  check1 <- identical(param_without_random_values(x$param),
                      param_without_random_values(y$param))
  ## these carry formulas, whose environments differ between two calls to
  ## `control.net`, so they are compared by value below rather than by identity
  env.controls <- c("monitors", "nwstats.formula", "set.control.ergm",
                    "set.control.tergm")

  check2 <- identical(x$control[-which(names(x$control) %in%
                                         c("nsims", env.controls))],
                      y$control[-which(names(y$control) %in%
                                         c("nsims", env.controls))])

  check2 <- check2 &&
    all(vapply(env.controls,
               function(nm) isTRUE(all.equal(x$control[[nm]], y$control[[nm]])),
               logical(1)))

  if (check1 == FALSE && param.error == TRUE) {
    stop("x and y have different parameters")
  }
  if (check2 == FALSE && param.error == TRUE) {
    stop("x and y have different controls")
  }

  z <- x
  z$control$nsims <- x$control$nsims + y$control$nsims
  newnames <- paste0("sim", seq_len(z$control$nsims))
  z$param$random.params.values <- merge_random_params_values(
    x$param$random.params.values, y$param$random.params.values,
    x$control$nsims, y$control$nsims
  )

  # Merge epi data
  for (i in seq_along(x$epi)) {
    z$epi[[i]] <- cbind(x$epi[[i]], y$epi[[i]])
    names(z$epi[[i]]) <- newnames
  }

  ## `save.other` must name the same elements in both objects for them to be
  ## combined
  if (keep.other == TRUE && !is.null(x$control$save.other) &&
        !is.null(y$control$save.other) &&
        !identical(x$control$save.other, y$control$save.other)) {
    stop("Elements in save.other differ between x and y")
  }

  ## slots the caller asked to drop rather than merge
  dropped <- c(
    if (keep.transmat == FALSE) "stats$transmat",
    if (keep.network == FALSE) "network",
    if (keep.nwstats == FALSE) "stats$nwstats",
    if (keep.diss.stats == FALSE) "diss.stats",
    if (keep.other == FALSE) x$control$save.other
  )

  ## every other slot holding one element per simulation is concatenated. A
  ## slot that only one of the objects carries is dropped rather than left
  ## holding one object's simulations alone
  for (path in netsim_sim_slots(x)) {
    if (paste(path, collapse = "$") %in% dropped || is.null(y[[path]])) {
      z[[path]] <- NULL
    } else {
      merged <- setNames(c(x[[path]], y[[path]]), newnames)
      ## `c` drops attributes, so any class on the list itself is restored
      class(merged) <- oldClass(x[[path]])
      z[[path]] <- merged
    }
  }

  return(z)
}
