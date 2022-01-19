flatten_params <- function(params) {
  params <- remove_special_params(params)
  params.length <- vapply(params, length, 0)
  params.n <- sum(params.length)
  params.flat <- vector(mode = "list", length = params.n)

  i <- n <- 1
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

list_special_params <- function(params) {
  builtin.special.params <- c(
    "param.updater.list",
    "random.params",
    "random.params.values"
  )

  builtin.special.params <- intersect(builtin.special.params, names(params))
  dot.special.params <- names(params)[grep("^\\.", names(params))]

  return(unique(c(builtin.special.params, dot.special.params)))
}

remove_special_params <- function(params) {
  special.params.names <- list_special_params(params)
  params[!names(params) %in% special.params.names]
}

unflatten_params <- function(params.flat) {
  check_params_names(names(params.flat))
  check_params_flat(params.flat)

  set.elements.names <- vapply(
    names(params.flat),
    function(x) sub("_.*$", "", x),
    ""
  )

  params.names <- unique(set.elements.names)
  params <- vector(mode = "list", length = length(params.names))
  names(params) <- params.names

  for (i in seq_along(set.elements.names)) {
    nme <- set.elements.names[[i]]
    params[[nme]] <- c(params[[nme]], params.flat[[i]])
  }

  return(params)
}

check_params_flat <- function(params.flat) {
  params.length <- vapply(params.flat, length, 0)
  params.list <- vapply(params.flat, is.list, TRUE)
  if (any(params.length != 1) || any(params.list)) {
    stop("A flat parameter list should contain only length 1 non list elements")
  }
  invisible(TRUE)
}

check_params_names <- function(params.names) {
  params.pattern <- "^[[:alpha:]][[:alnum:].]*(_[1-9]+)?$"
  correct.format <- grepl(params.pattern, params.names)

  if (!all(correct.format)) {
    stop("The following flat parameter names are malformed: \n`",
      paste0(params.names[!correct.format], collapse = "`, `"), "`\n\n",
      "you can check the names with ",
      '`grepl("', params.pattern, '", your.names)` \n',
      "Example: 'unique.param', 'param.set_1', 'param.set_2'"
    )
  }

  special.params.names <- list_special_params(params.names)
  if (length(special.params.names) != 0) {
    stop("The following special parameter names are not allowed: \n`",
         paste0(special.params.names, collapse = "`, `"), "`\n\n"
    )
  }

  invisible(return(TRUE))
}

make_scenarios_list <- function(scenarios.df) {
  if (!is.null(scenarios.df[[".scenario.id"]])) {
    scenarios.names <- scenarios.df[[".scenario.id"]]
    scenarios.df[[".scenario.id"]] <- NULL
  } else {
    scenarios.names <- seq_len(nrow(scenarios.df))
  }

  scenarios.list <- lapply(
    seq_len(nrow(scenarios.df)),
    function(row) unflatten_params(scenarios.df[row, ])
  )
  names(scenarios.list) <- scenarios.names

  return(scenarios.list)
}

make_scenarios_updaters <- function(scenarios.df, at) {
  scenarios.list <- make_scenarios_list(scenarios.df)
  updaters.list <- lapply(
    scenarios.list,
    function(scenario) list(list(at = at, param = scenario))
  )
  names(updaters.list) <- names(scenarios.list)
  return(updaters.list)
}

make_scenarios_df <- function(scenarios.list) {
  scenarios.rows <- lapply(scenarios.list, flatten_params)
  scenarios.df <- dplyr::bind_rows(scenarios.rows)
  dplyr::mutate(scenarios.df, .scenario.id = names(scenarios.list))
}
# update params to create scenarios
#
# params_orig <- function() {
#   list(
#     aasd.ad = sample(seq_len(9), 1),
#     d.adk = sample(seq_len(9), 3, replace = TRUE),
#     eas.d = sample(letters, 1, replace = TRUE),
#     fd.dd = sample(LETTERS, 4, replace = TRUE)
#   )
# }
#
# scenarios.df <- dplyr::bind_rows(
#   lapply(1:5, function(x) tibble::as_tibble(flatten_params(params_orig())))
# )
# scenarios.df[[".scenario.id"]] <- replicate(
#   5,
#   paste0(sample(LETTERS, 10, replace = TRUE), collapse = "")
# )
#
# scenarios.list <- make_scenarios_list(scenarios.df)
# scenarios.updaters <- make_scenarios_updaters(scenarios.df, 100)
#
# param_df <- read.csv("/home/adrien/Documents/Projects/BigNets/data/input/calibration.csv")
# make_scenarios_list(param_df)
