#' @title Parameters That Vary With the Duration of Infection
#'
#' @description Marks a vector of transmission probabilities, act rates, or
#'              recovery rates as varying with the time since infection, for
#'              the `inf.prob`, `inf.prob.g2`, `act.rate`, `rec.rate`, and
#'              `rec.rate.g2` arguments of [param.net()]. The built-in
#'              [netsim()] modules apply the value at each infected node's
#'              current duration of infection.
#'
#' @param x A numeric vector with one value per time step since infection: the
#'        first element applies in the first time step of an infection, the
#'        second in the second, and so on. The last element carries forward
#'        for the rest of the infection.
#'
#' @details
#' A duration of infection is the number of time steps since the infected
#' node was infected, counted from 1: a node infected at time step 10 is in
#' its first step of infection at time steps 10 and 11, its second at time
#' step 12, and so on. For transmission, the duration is that of the infected
#' partner; for recovery, that of the recovering node.
#'
#' `by_infection_duration(c(0.5, 0.5, 0.1))` gives a probability of 0.5 in the
#' first two time steps of an infection and 0.1 from the third time step on,
#' which can represent an acute stage of high infectiousness. This describes
#' variation over the course of each infection, not over calendar time: to
#' change a parameter at a given time step of the simulation, use scenarios
#' or parameter updaters (see `vignette("model-parameters", package =
#' "EpiModel")`).
#'
#' In multi-layer models, an entry of a [multilayer()] parameter may itself be
#' a `by_infection_duration` object:
#' `inf.prob = multilayer(by_infection_duration(c(0.5, 0.1)), 0.05)`.
#'
#' With the built-in modules, a plain vector passed to these arguments stops
#' [netsim()] with an error, since a plain vector carries no record of what its
#' positions mean; custom modules keep their own reading of plain vectors, such
#' as values by group or by layer. Parameter tables cannot represent these
#' objects yet, so [param.net_to_table()] stops on them.
#'
#' @return An object of class `by_infection_duration`: the numeric vector `x`
#'         with that class.
#'
#' @seealso [param.net()], [multilayer()].
#'
#' @export
#'
#' @examples
#' # Transmission probability of 0.05 in the first five time steps of an
#' # infection, then 0.15
#' param <- param.net(inf.prob = by_infection_duration(c(rep(0.05, 5), 0.15)),
#'                    act.rate = 1)
#' param
#'
#' \dontrun{
#' nw <- network_initialize(n = 100)
#' est <- netest(nw, formation = ~edges, target.stats = 50,
#'               coef.diss = dissolution_coefs(~offset(edges), 10),
#'               verbose = FALSE)
#' sim <- netsim(est, param, init.net(i.num = 10),
#'               control.net(type = "SI", nsteps = 25, nsims = 1,
#'                           verbose = FALSE))
#' tm <- get_transmat(sim)
#' table(tm$infDur, tm$transProb)
#' }
#'
by_infection_duration <- function(x) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x)) {
    stop("`x` must be a numeric vector of length one or more, without ",
         "missing values.", call. = FALSE)
  }
  structure(as.numeric(unclass(x)), class = "by_infection_duration")
}

#' @export
print.by_infection_duration <- function(x, ...) {
  cat("Values by duration of infection (time steps since infection),",
      "the last carried forward:\n")
  print(unclass(x), ...)
  invisible(x)
}

#' @export
format.by_infection_duration <- function(x, ...) {
  paste0("by_infection_duration(", paste(deparse(unclass(x)), collapse = ""),
         ")")
}

# The value of a parameter for each infected node given its duration of
# infection: a scalar applies at every duration, and a by_infection_duration()
# vector is indexed by duration, its last value carried forward. `name` and
# `layer` label the error for a plain vector.
infection_duration_value <- function(x, infDur, name, layer = NULL) {
  if (inherits(x, "by_infection_duration")) {
    x <- unclass(x)
    return(x[pmin(infDur, length(x))])
  }
  if (length(x) != 1) {
    stop_plain_duration_vector(name, x, layer)
  }
  rep(x, length(infDur))
}

# As infection_duration_value(), for a parameter that may also vary by layer
# through multilayer(): each discordant edge gets the entry for its layer.
layer_duration_value <- function(x, network, infDur, name) {
  if (!inherits(x, "multilayer")) {
    return(infection_duration_value(x, infDur, name))
  }
  out <- numeric(length(network))
  for (k in unique(network)) {
    idx <- which(network == k)
    out[idx] <- infection_duration_value(x[[k]], infDur[idx], name, layer = k)
  }
  out
}

stop_plain_duration_vector <- function(name, x, layer = NULL) {
  stop("`", name, "`", if (!is.null(layer)) paste0(" (layer ", layer, ")"),
       " is a plain vector of length ", length(x), ". The built-in modules ",
       "no longer read a plain vector as varying with the duration of ",
       "infection; wrap it in by_infection_duration() to give one value per ",
       "time step since infection.", call. = FALSE)
}

# Check the parameters that the built-in netsim modules read by duration of
# infection: each must be a single value or a by_infection_duration() object,
# and inf.prob, inf.prob.g2, and act.rate may also be multilayer() objects of
# those. Called by crosscheck.net() only when a built-in model type is used,
# since custom modules give plain vectors their own meaning.
check_duration_params <- function(param) {
  for (name in c("inf.prob", "inf.prob.g2", "act.rate", "rec.rate",
                 "rec.rate.g2")) {
    x <- param[[name]]
    if (is.null(x)) {
      next
    }
    if (inherits(x, "multilayer")) {
      if (name %in% c("rec.rate", "rec.rate.g2")) {
        stop("`", name, "` cannot be a multilayer() object: recovery does ",
             "not depend on the network layer.", call. = FALSE)
      }
      for (k in seq_along(x)) {
        check_duration_value(x[[k]], name, layer = k)
      }
    } else {
      check_duration_value(x, name)
    }
  }
  invisible(TRUE)
}

check_duration_value <- function(x, name, layer = NULL) {
  if (inherits(x, "by_infection_duration")) {
    return(invisible(TRUE))
  }
  if (!is.numeric(x)) {
    stop("`", name, "`", if (!is.null(layer)) paste0(" (layer ", layer, ")"),
         " must be numeric.", call. = FALSE)
  }
  if (length(x) != 1) {
    stop_plain_duration_vector(name, x, layer)
  }
  invisible(TRUE)
}

# by_infection_duration() objects apply to netsim only; DCM and ICM read
# vector parameters differently (sensitivity runs for DCMs).
stop_if_by_infection_duration <- function(param, model) {
  bid <- names(param)[vapply(param, inherits, logical(1),
                             "by_infection_duration")]
  if (length(bid) > 0) {
    stop("by_infection_duration() parameters apply to network models ",
         "simulated with netsim() only, not to ", model, " models: `",
         paste(bid, collapse = "`, `"), "`.", call. = FALSE)
  }
  invisible(TRUE)
}
