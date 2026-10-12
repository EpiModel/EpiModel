#' Make a list of EpiModel scenarios from a data.frame of scenarios
#'
#' An EpiModel scenario allows one or multiple set of parameters to be applied
#' to a model a predefined timesteps. They are usually used by a researcher who
#' wants to model counterfactuals using a pre calibrated model.
#'
#' @param scenarios.df a `data.frame`
#'
#' @return a list of EpiModel scenarios
#'
#' @section scenarios.df:
#' The `scenarios.df` is a `data.frame` of values to be used as
#' parameters.
#'
#' It must contain a ".at" column, specifying when the changes should occur.
#' It requires the "updater" module of EpiModel. *See, vignette*. If the ".at"
#' value of a row is less than two, the changes will be applied to the
#' parameter list iteself. The second mandatory column is ".scenario.id". It
#' is used to distinguish the different scenarios. If multiple rows share the
#' same ".scenario.id", the resulting scenario will contain one updater per row.
#' This permits modifying parameters at multiple points in time. (e.g. an
#' intervention limited in time).
#'
#' The other column names must correspond either to:
#' the name of one parameter if this parameter is of size 1 or the name of the
#' parameter with "_1", "_2", "_N" with the second part being the position of
#' the value for a parameter of size > 1. This means that the parameter names
#' cannot contain any underscore "_". (e.g "a.rate", "d.rate_1", "d.rate_2")
#'
#' @section Vector parameters:
#' A column with a position suffix sets that position of the current value of
#' the parameter and leaves the other positions as they are. For example, with
#' `prep.start.rate = c(0.0055, 0.0042, 0.0066)` in `param`, a scenario with
#' only a `prep.start.rate_1` column of `0.5` gives
#' `prep.start.rate = c(0.5, 0.0042, 0.0066)`, both when it is applied before
#' the simulation by [use_scenario()] and when it is applied by an updater
#' during the simulation. Before EpiModel 2.7.0, the parameter was replaced by
#' a vector built from the scenario columns alone (`0.5` here), so every
#' position had to be given.
#'
#' A column without a suffix (`prep.start.rate`) still replaces the whole
#' parameter with its single value. Positions past the end of the current
#' vector extend it when they leave no gap: with a vector of length 3,
#' `prep.start.rate_4` alone extends it to length 4, while
#' `prep.start.rate_5` alone is an error. A scenario cannot shorten a vector.
#' Only plain vectors can be set by position: a parameter given as a
#' [multilayer()] or [by_infection_duration()] object, or as a list, is an
#' error. Use a parameter updater that sets the whole value instead.
#'
#' @export
create_scenario_list <- function(scenarios.df) {
  check_scenarios_df(scenarios.df)
  scenarios.names <- unique(scenarios.df[[".scenario.id"]])

  scenarios <- lapply(
    scenarios.names,
    function(id) {
      make_scenario(
                    scenarios.df[scenarios.df[[".scenario.id"]] == id, ])
    }
  )

  names(scenarios) <- scenarios.names

  return(scenarios)
}

#' a scenario is a list with the following elements:
#'   - `id`: the identifier of the scenario
#'   - `.param.updater.list`: as in the updater module
#'
#' when running "use_scenario", if any updater has an `at` value < 2, it is
#' applied immediatly (so before netsim)
#'
#' An updater made from columns with a position suffix (`x_2`) also holds a
#' `positions` list, giving for each such parameter the positions the row
#' sets (`list(x = 2)`). Only these positions of the current value are
#' changed when the updater is applied (see `merge_param_positions`).
#'
#' @noRd
make_scenario <- function(scenario.rows) {
  scenario <- list(
    id = scenario.rows[[".scenario.id"]][1],
    .param.updater.list = vector(mode = "list", length = nrow(scenario.rows))
  )

  elements.at <- scenario.rows[[".at"]]
  # ensures that `scenario.rows` is a `tibble`
  scenario.rows <- dplyr::select(scenario.rows, - c(".at", ".scenario.id"))

  check_params_names(names(scenario.rows))
  positions <- flat_params_positions(names(scenario.rows))

  for (i in seq_along(elements.at)) {
    updater <- list(
      at = elements.at[[i]],
      param = unflatten_params(scenario.rows[i, , drop = FALSE])
    )
    if (length(positions) > 0) {
      updater[["positions"]] <- positions
    }
    scenario[[".param.updater.list"]][[i]] <- updater
  }

  return(scenario)
}

#' Apply a scenario object to a param.net object
#'
#' @param scenario a scenario object usually created from a `data.frame` of
#' scenarios using the `create_scenario_list` function. See the vignette
#' "network-model-scenarios".
#'
#' @section scenario:
#' A scenario is a list containing an "id" field, the name of the scenario and
#' a ".param.updater.list" containing a list of updaters that modifies the
#' parameters of the model at given time steps. If a scenario contains a
#' parameter not defined in the `param` object, an error will be produced.
#' A scenario column with a position suffix (`x_2`) sets only that position of
#' a vector parameter; see the "Vector parameters" section of
#' [create_scenario_list()].
#' See the vignette "model-parameters" for the technical detail of their
#' implementation.
#'
#' @inheritParams update_params
#' @inherit update_params return
#'
#' @export
use_scenario <- function(param, scenario) {
  scenario.params <- unique(unlist(lapply(
    scenario[[".param.updater.list"]],
    function(element) names(element[["param"]])
  )))

  undef.params <- setdiff(scenario.params, names(param))
  if (length(undef.params) > 0) {
    stop("Some of the scenario parameters are not defined in `param`: \n'",
         paste0(undef.params, collapse = "', '"), "'")
  }

  elements.at <- vapply(
    scenario[[".param.updater.list"]],
    function(element) element[["at"]],
    numeric(1)
  )

  for (i in which(elements.at < 2)) {
    updater <- scenario[[".param.updater.list"]][[i]]
    new.param <- merge_param_positions(param, updater[["param"]],
                                       updater[["positions"]])
    param <- update_params(param, new.param)
  }

  param[[".param.updater.list"]] <- c(
    param[[".param.updater.list"]],
    scenario[[".param.updater.list"]][elements.at >= 2]
  )

  param[[".scenario.id"]] <- scenario[["id"]]

  return(param)
}

#' Helper function validating the format of a `scenarios.df`
#' @noRd
check_scenarios_df <- function(scenarios.df) {
  checks <-
    all(c(".scenario.id", ".at") %in% names(scenarios.df)) &&
    is.numeric(scenarios.df[[".at"]]) &&
    all(as.integer(scenarios.df[[".at"]]) == scenarios.df[[".at"]])

  if (!checks) {
    stop(
      "A `data.frame` of scenarios must have a '.scenario.id' column \n",
      "and a '.at' column containing integers."
    )
  }
}

#' helper function to make a ragged param list into a flat one
#' @noRd
flatten_params <- function(params) {
  params <- remove_special_params(params)
  params.length <- vapply(params, length, 0)
  params.n <- sum(params.length)
  params.flat <- vector(mode = "list", length = params.n)

  i <- 1
  n <- 1
  while (n <= params.n) {
    l <- params.length[i]
    cur.param <- as.list(params[[i]])
    cur.name <- names(params)[i]
    if (l > 1)
      cur.name <- paste0(cur.name, "_", seq_len(l))
    params.flat[n:(n + l - 1)] <- cur.param
    names(params.flat)[n:(n + l - 1)] <- cur.name

    n <- n + l
    i <- i + 1
  }
  return(params.flat)
}

#' list the "special parameters" from a param list: all parameters starting
#' with "."
#' @noRd
list_special_params <- function(params.names) {
  return(params.names[grep("^\\.", params.names)])
}

#' helper function to remove the "special parameters" from a param list.
#' @noRd
remove_special_params <- function(params) {
  special.params.names <- list_special_params(names(params))
  params[!names(params) %in% special.params.names]
}

#' helper function to split flat parameter names into the parameter name, the
#' position (1 when the name has no position suffix), and whether the name has
#' a position suffix
#' @noRd
split_flat_names <- function(flat.names) {
  elts <- strsplit(flat.names, "_")
  list(
    name = vapply(elts, function(x) x[1], ""),
    pos = vapply(elts, function(x) {
      as.numeric(if (length(x) == 2) x[2] else 1)
    }, numeric(1)),
    suffixed = lengths(elts) == 2
  )
}

#' helper function listing, for each parameter set by position in a flat
#' parameter list, the positions set. A name without a suffix next to suffixed
#' names of the same parameter counts as position 1, as in `unflatten_params`.
#' Returns an empty list when no name has a position suffix.
#' @noRd
flat_params_positions <- function(flat.names) {
  flat <- split_flat_names(flat.names)
  by.position <- unique(flat$name[flat$suffixed])
  positions <- lapply(by.position, function(nm) flat$pos[flat$name == nm])
  names(positions) <- by.position
  positions
}

#' helper function merging the values of parameters set by position into
#' their current values
#'
#' `current` and `new` are named lists of parameters, and `positions` a named
#' list giving for some parameters of `new` the positions it sets (see
#' `make_scenario`). For each of these parameters, the positions set are copied
#' from `new` into the current value, leaving its other positions as they
#' are, and the result replaces the value in `new`. A parameter missing from
#' `current` keeps its value in `new`, as before positions were merged.
#' @noRd
merge_param_positions <- function(current, new, positions) {
  for (nm in intersect(names(positions), names(new))) {
    if (is.null(current[[nm]])) next
    new[[nm]] <- merge_positions(current[[nm]], new[[nm]], positions[[nm]], nm)
  }
  new
}

#' helper function setting positions `pos` of a plain vector `current` to
#' the same positions of `value`. Positions past the end of `current` extend
#' it when they leave no gap.
#' @noRd
merge_positions <- function(current, value, pos, name) {
  if (!is.atomic(current) || is.object(current)) {
    stop("A scenario sets `", name, "` by position (`", name, "_",
         min(pos), "`), but `", name, "` is a ", class(current)[1],
         " object rather than a plain vector. Only plain vectors can be set ",
         "by position. Use a parameter updater that sets the whole value ",
         "instead.", call. = FALSE)
  }
  gaps <- setdiff(seq_len(max(pos)), c(seq_along(current), pos))
  if (length(gaps) > 0) {
    stop("A scenario sets `", name, "_", max(pos), "`, but `", name,
         "` has length ", length(current), " and the scenario gives no value ",
         "for position ", paste(gaps, collapse = ", "), ". Positions past the ",
         "end of a vector parameter must leave no gap.", call. = FALSE)
  }
  current[pos] <- value[pos]
  current
}

#' helper function to make a flat param list into a ragged one
#' @noRd
unflatten_params <- function(params.flat) {
  check_params_names(names(params.flat))
  check_params_flat(params.flat)

  # get the param name and position for each element.
  # With pos == 1 if none provided
  flat <- split_flat_names(names(params.flat))
  set.elements.names <- flat$name
  set.elements.pos <- flat$pos

  params.names <- unique(set.elements.names)
  params <- vector(mode = "list", length = length(params.names))
  names(params) <- params.names

  for (i in seq_along(set.elements.names)) {
    nme <- set.elements.names[[i]]
    params[[nme]][set.elements.pos[i]] <- params.flat[[i]]
  }

  return(params)
}

#' helper function to check if a list of flat parameters is actually flat
#' @noRd
check_params_flat <- function(params.flat) {
  params.length <- vapply(params.flat, length, 0)
  params.list <- vapply(params.flat, is.list, TRUE)
  if (any(params.length != 1) || any(params.list)) {
    stop("A flat parameter list should contain only length 1 non list elements")
  }
  invisible(TRUE)
}

#' helper function to check the correctness of the flat parameters names
#' @noRd
check_params_names <- function(params.names) {
  # The position suffix is a decimal integer with no leading zero. The bracket
  # expression must allow 0 after the first digit, otherwise any position
  # containing a 0 (`_10`, `_20`, `_100`) is rejected while its neighbors pass.
  params.pattern <- "^[[:alpha:]][[:alnum:].]*(_[1-9][0-9]*)?$"
  correct.format <- grepl(params.pattern, params.names)

  if (!all(correct.format)) {
    stop("The following flat parameter names are malformed: \n`",
      paste0(params.names[!correct.format], collapse = "`, `"), "`\n\n",
      "you can check the names with ",
      '`grepl("', params.pattern, '", your.names)` \n',
      "Example: 'unique.param', 'param.set_1', 'param.set_2'"
    )
  }

  invisible(return(TRUE))
}
