
#' @title Functions to Access and Edit the Master List Object in Network Models
#'
#' @description These `get_`, `set_`, `append_` and `add` functions allow a safe
#'              and efficient way to retrieve and mutate the Master list object
#'              of network models (`dat`).
#'
#' @param dat a Master list object of network models
#' @param item a character vector conaining the name of the element to access.
#'        Can be of length > 1 for `get_*_list` functions
#' @param indexes for `get_epi` and `get_attr`, a numeric vector of indexes or
#'        a logical vector to subset the desired `item`
#' @param value new value to be attributed in the `set_` and `append_` functions
#' @param override.null.error if TRUE, `get_` return NULL if the `item` does not
#'        exist instead of throwing an error. (default = FALSE)
#' @param override.length.check if TRUE, `set_attr` allows the modification of
#'        the `item` size. (default = FALSE)
#' @return a vector or a list of vector for `get_` functions. And the Master
#'         list object for `set_` and `add_` functions
#'
#' @section Mutability:
#' The `set_`, `append_` and `add_` functions DO NOT mutate the dat object in
#' place. The result must be assigned back to `dat` in order to be registered
#' `dat <- set_*(dat, item, value)`
#'
#' @section `set_` vs `add_`:
#' The `set_` functions edit a pre-existing element or create a new one if it
#' does not exist already by calling the `add_` functions internally.
#'
#' @examples
#' dat <- list(
#'   attr = list(
#'     active = rbinom(100, 1, 0.9)
#'   ),
#'   epi = list(),
#'   param = list(),
#'   init = list(),
#'   control = list(
#'     nsteps = 150
#'   )
#' )
#'
#' dat <- add_attr(dat, "age")
#' dat <- set_attr(dat, "age", runif(100))
#' dat <- set_attr(dat, "status", rbinom(100, 1, 0.9))
#' dat <- set_attr(dat, "status", rep(1, 150), override.length.check = TRUE)
#' dat <- append_attr(dat, "status", 1, 10)
#' dat <- append_attr(dat, "age", NA, 10)
#' get_attr_list(dat)
#' get_attr_list(dat, c("age", "active"))
#' get_attr(dat, "status")
#' get_attr(dat, "status", c(1, 4))
#'
#' dat <- add_epi(dat, "i.num")
#' dat <- set_epi(dat, "i.num", 150, 10)
#' dat <- set_epi(dat, "s.num", 150, 90)
#' get_epi_list(dat)
#' get_epi_list(dat, c("i.num", "s.num"))
#' get_epi(dat, "i.num")
#' get_epi(dat, "i.num", c(1, 4))
#' get_epi(dat, "i.num", rbinom(150, 1, 0.2) == 1)
#'
#' dat <- add_param(dat, "x")
#' dat <- set_param(dat, "x", 0.4)
#' dat <- set_param(dat, "y", 0.8)
#' get_param_list(dat)
#' get_param_list(dat, c("x", "y"))
#' get_param(dat, "x")
#'
#' dat <- add_init(dat, "x")
#' dat <- set_init(dat, "x", 0.4)
#' dat <- set_init(dat, "y", 0.8)
#' get_init_list(dat)
#' get_init_list(dat, c("x", "y"))
#' get_init(dat, "x")
#'
#' dat <- add_control(dat, "x")
#' dat <- set_control(dat, "x", 0.4)
#' dat <- set_control(dat, "y", 0.8)
#' get_control_list(dat)
#' get_control_list(dat, c("x", "y"))
#' get_control(dat, "x")
#'
#' @name net-accessor
NULL

#' @rdname net-accessor
#' @export
get_attr_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["attr"]]

  } else {
    missing_item <- setdiff(item, names(dat[["attr"]]))
    if (length(missing_item) > 0) {
      stop(paste("There is no attributes called",
                 paste(missing_item, collapse = ", "),
                 "in the attributes list of the Master list object (dat)"))
    }

    out <- dat[["attr"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_attr <- function(dat, item, indexes = NULL, override.null.error = FALSE) {
  if (!item %in% names(dat[["attr"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop(paste("There is no attribute called", item,
                 "in the attributes list of the Master list object (dat)"))
    }
  } else {
    if (is.null(indexes)) {
      out <- dat[["attr"]][[item]]
    } else {
      if (is.logical(indexes)) {
        if (length(indexes) != length(dat[["attr"]][[item]])) {
          stop("(logical) `indexes` has to have a length equal to the number of
              nodes in the network")
        }
      } else if (is.numeric(indexes)) {
        if (any(indexes > length(dat[["attr"]][[item]]))) {
          stop("Some (numeric) `indexes` are larger than the number of nodes in
              the network")
        }
      } else {
        stop("`indexes` must be logical, numeric, or NULL")
      }

      out <- dat[["attr"]][[item]][indexes]
    }
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_attr <- function(dat, item) {
  if (item %in% names(dat[["attr"]])) {
    stop(paste0("Cannot create the attribute '", item,
               "': exists already"))
  }

  dat[["attr"]][[item]] <- rep(NA, length(dat[["attr"]][["active"]]))

  return(dat)
}

#' @rdname net-accessor
#' @export
set_attr <- function(dat, item, value, override.length.check = FALSE) {
  if (!item %in% names(dat[["attr"]])) {
    dat <- add_attr(dat, item)
  }

  if (!override.length.check &&
      length(value) != length(dat[["attr"]][["active"]])) {
    stop(
      "When trying to edit the ", `item`, " nodal attribute: ",
      "The size of the `value` vector is not equal to the number of node in",
      "the network. \n",
      "Expected: ", length(dat[["attr"]][["active"]]), "\n" ,
      "Given: ", length(value)
    )
  }

  dat[["attr"]][[item]] <- value

  return(dat)
}

#' @param n.new the number of new elements to append at the end of `item`
#' @rdname net-accessor
#' @export
append_attr <- function(dat, item, value, n.new) {
  if (!is.numeric(n.new) || n.new < 0) {
    stop("`n_new` must be numeric and greater than or equal to zero.")
  }

  if (length(value) == 1) {
    new_vals <- rep(value, n.new)
  } else if (length(value) == n.new) {
    new_vals <- value
  } else {
    stop("`value` must be of length one or `n.new`.")
  }

  old_vals <- get_attr(dat, item, override.null.error = TRUE)
  dat <- set_attr(dat, item, c(old_vals, new_vals),
                  override.length.check = TRUE)

  return(dat)
}

#' @rdname net-accessor
#' @export
get_epi_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["epi"]]

  } else {
    missing_item <- setdiff(item, names(dat[["epi"]]))
    if (length(missing_item) > 0) {
      stop(paste("There is no epi output called",
                 paste(missing_item, collapse = ", "),
                 "in the epi output list of the Master list object (dat)"))
    }

    out <- dat[["epi"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_epi <- function(dat, item, indexes = NULL, override.null.error = FALSE) {
  if (!item %in% names(dat[["epi"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop(paste("There is no epi out called", item,
                 "in the epi out list of the Master list object (dat)"))
    }
  } else {
    if (is.null(indexes)) {
      out <- dat[["epi"]][[item]]
    } else {
      if (is.logical(indexes)) {
        if (length(indexes) != dat[["control"]][["nsteps"]]) {
          stop("(logical) `indexes` has to have a length equal to the number of
              steps planned for for the simulation (control[['nsteps']])")
        }
      } else if (is.numeric(indexes)) {
        if (any(indexes > dat[["control"]][["nsteps"]])) {
          stop("Some (numeric) `indexes` are larger than the number of
              steps planned for for the simulation (control[['nsteps']])")
        }
      } else {
        stop("`indexes` must be logical, numeric, or NULL")
      }

      out <- dat[["epi"]][[item]][indexes]
    }
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_epi <- function(dat, item) {
  if (item %in% names(dat[["epi"]])) {
    stop(paste("Cannot create the epi output, ", item,
               ": exists already"))
  }

  dat[["epi"]][[item]] <- rep(NA, dat[["control"]][["nsteps"]])

  return(dat)
}

#' @param at timestep where to add the new value for the epi output `item`
#' @rdname net-accessor
#' @export
set_epi <- function(dat, item, at,  value) {
  if (length(at) != 1 || !is.numeric(at)) {
    stop("`at` must be numeric and of length one")
  }

  if (!item %in% names(dat[["epi"]])) {
    dat <- add_epi(dat, item)
  }

  if (at > length(dat[["epi"]][[item]])) {

      dat[["epi"]][[item]] <- c(
        dat[["epi"]][[item]],
        rep(NA, dat[["control"]][["nsteps"]] - length(dat[["epi"]][[item]]))
      )
  }

  dat[["epi"]][[item]][at] <- value

  return(dat)
}

#' @rdname net-accessor
#' @export
get_param_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["param"]]

  } else {
    missing_item <- setdiff(item, names(dat[["param"]]))
    if (length(missing_item) > 0) {
      stop(paste("There is no parameters called",
                 paste(missing_item, collapse = ", "),
                 "in the parameter list of the Master list object (dat)"))
    }

    out <- dat[["param"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_param <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat[["param"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop(paste("There is no parameter called", item,
                 "in the parameter list of the Master list object (dat)"))
    }
  } else {
    out <- dat[["param"]][[item]]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_param <- function(dat, item) {
  if (item %in% names(dat[["param"]])) {
    stop(paste("Cannot create the parameter, ", item,
               ": exists already"))
  }

  dat[["param"]][[item]] <- NA

  return(dat)
}

#' @rdname net-accessor
#' @export
set_param <- function(dat, item, value) {
  if (!item %in% names(dat[["param"]])) {
    dat <- add_param(dat, item)
  }

  dat[["param"]][[item]] <- value

  return(dat)
}

#' @rdname net-accessor
#' @export
get_control_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["control"]]

  } else {
    missing_item <- setdiff(item, names(dat[["control"]]))
    if (length(missing_item) > 0) {
      stop(paste("There is no control value called",
                 paste(missing_item, collapse = ", "),
                 "in the control list of the Master list object (dat)"))
    }

    out <- dat[["control"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_control <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat[["control"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop(paste("There is no control value called", item,
                 "in the control list of the Master list object (dat)"))
    }
  } else {
    out <- dat[["control"]][[item]]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_control <- function(dat, item) {
  if (item %in% names(dat[["control"]])) {
    stop(paste("Cannot create the control value, ", item,
               ": exists already"))
  }

  dat[["control"]][[item]] <- NA

  return(dat)
}

#' @rdname net-accessor
#' @export
set_control <- function(dat, item, value) {
  if (!item %in% names(dat[["control"]])) {
    dat <- add_control(dat, item)
  }

  dat[["control"]][[item]] <- value

  return(dat)
}

#' @rdname net-accessor
#' @export
get_init_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["init"]]

  } else {
    missing_item <- setdiff(item, names(dat[["init"]]))
    if (length(missing_item) > 0) {
      stop(paste("There is no init value called",
                 paste(missing_item, collapse = ", "),
                 "in the init list of the Master list object (dat)"))
    }

    out <- dat[["init"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_init <- function(dat, item, override.null.error = FALSE) {
  if (!item %in% names(dat[["init"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop(paste("There is no init value called", item,
                 "in the init list of the Master list object (dat)"))
    }
  } else {
    out <- dat[["init"]][[item]]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_init <- function(dat, item) {
  if (item %in% names(dat[["init"]])) {
    stop(paste("Cannot create the init value, ", item,
               ": exists already"))
  }

  dat[["init"]][[item]] <- NA

  return(dat)
}

#' @rdname net-accessor
#' @export
set_init <- function(dat, item, value) {
  if (!item %in% names(dat[["init"]])) {
    dat <- add_init(dat, item)
  }

  dat[["init"]][[item]] <- value

  return(dat)
}
