#' @title Functions to Access and Edit the Main List Object in Network Models
#'
#' @description These \code{get_}, \code{set_}, \code{append_}, and \code{add_}
#'              functions allow a safe and efficient way to retrieve and mutate
#'              the main list object of network models (\code{dat}).
#'
#' @inheritParams recovery.net
#' @param item A character vector containing the name of the element to access
#'        (for \code{get_} functions), create (for \code{add_} functions), or
#'        edit (for \code{set_} and \code{append_} functions). Can be of length
#'        > 1 for \code{get_*_list} functions.
#' @param posit_ids For \code{set_attr} and \code{get_attr}, a numeric vector of
#'        posit_ids or a logical vector to subset the desired \code{item}.
#' @param value New value to be attributed in the \code{set_} and \code{append_}
#'        functions.
#' @param override.null.error If TRUE, \code{get_} will return NULL if the
#'         \code{item} does not exist instead of throwing an error.
#'         (default = FALSE).
#' @param override.length.check If TRUE, \code{set_attr} allows the modification
#'        of the \code{item} size. (default = FALSE).
#' @param n.new For \code{append_core_attr}, the number of new nodes to initiate
#'        with core attributes; for \code{append_attr}, the number of new
#'        elements to append at the end of \code{item}.
#' @param at For \code{get_epi}, the timestep at which to access the specified
#'        \code{item}; for \code{set_epi}, the timestep at which to add the new
#'        value for the epi output \code{item}; for \code{append_core_attr}, the
#'        current time step.
#' @return A vector or a list of vectors for \code{get_} functions; the main
#'         list object for \code{set_}, \code{append_}, and \code{add_}
#'         functions.
#'
#' @section Core Attribute:
#' The \code{append_core_attr} function initializes the attributes necessary for
#' EpiModel to work (the four core attributes are: "active", "unique_id",
#' "entrTime", and "exitTime"). These attributes are used in the initialization
#' phase of the simulation, to create the nodes (see
#' \code{\link{initialize.net}}); and also used when adding nodes during the
#' simulation (see \code{\link{arrivals.net}}).
#'
#' @section Mutability:
#' The \code{set_}, \code{append_}, and \code{add_} functions DO NOT modify the
#' \code{dat} object in place. The result must be assigned back to \code{dat} in
#' order to be registered: \code{dat <- set_*(dat, item, value)}.
#'
#' @section \code{set_} and \code{append_} vs \code{add_}:
#' The \code{set_} and \code{append_} functions edit a pre-existing element or
#' create a new one if it does not exist already by calling the \code{add_}
#' functions internally.
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
      stop("There is no attributes called `",
           paste(missing_item, collapse = ", "),
           "` in the attributes list of the main list object (dat)")
    }

    out <- dat[["attr"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_attr <- function(dat, item, posit_ids = NULL, override.null.error = FALSE) {
  if (!item %in% names(dat[["attr"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop("There is no attribute called `", item,
           "` in the attributes list of the main list object (dat)")
    }
  } else {
    if (is.null(posit_ids)) {
      out <- dat[["attr"]][[item]]
    } else {
      if (is.logical(posit_ids)) {
        if (length(posit_ids) != length(dat[["attr"]][[item]])) {
          stop("(logical) `posit_ids` has to have a length equal to the ",
               "number of nodes in the network")
        }
      } else if (is.numeric(posit_ids)) {
        if (length(posit_ids > 0) &&
            any(posit_ids > length(dat[["attr"]][[item]]))) {
          stop("Some (numeric) `posit_ids` are larger than the number of ",
               "nodes in the network")
        }
      } else {
        stop("`posit_ids` must be logical, numeric, or NULL")
      }

      out <- dat[["attr"]][[item]][posit_ids]
    }
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_attr <- function(dat, item) {
  if (item %in% names(dat[["attr"]])) {
    stop("Cannot create the attribute '", item, "': exists already")
  }

  dat[["attr"]][[item]] <- rep(NA, length(dat[["attr"]][["active"]]))

  return(dat)
}

#' @rdname net-accessor
#' @export
set_attr <- function(dat, item, value, posit_ids = NULL,
  override.length.check = FALSE) {
  if (!item %in% names(dat[["attr"]])) {
    dat <- add_attr(dat, item)
  }

  if (is.null(posit_ids)) {
    if (!override.length.check &&
      length(value) != length(dat[["attr"]][["active"]])) {
      stop(
        "When trying to edit the ", `item`, " nodal attribute: ",
        "The size of the `value` vector is not equal to the number of nodes in",
        " the network. \n",
        "Expected: ", length(dat[["attr"]][["active"]]), "\n",
        "Given: ", length(value)
      )
    }
    dat[["attr"]][[item]] <- value
  } else {
    if (is.logical(posit_ids)) {
      if (length(posit_ids) != length(dat[["attr"]][[item]])) {
        stop("(logical) `posit_ids` has to have a length equal to the number ",
          "of nodes in the network")
      }
    } else if (is.numeric(posit_ids)) {
      if (length(posit_ids) == 0) {
        return(dat)
      } else if (any(posit_ids > length(dat[["attr"]][[item]]))) {
        stop("Some (numeric) `posit_ids` are larger than the number of nodes ",
          " in the network")
      }
    } else {
      stop("`posit_ids` must be logical, numeric, or NULL")
    }

    if (!override.length.check &&
      length(value) != 1 &&
      length(value) != length(dat[["attr"]][["active"]][posit_ids])) {
      stop(
        "When trying to edit the `", item, "` nodal attribute: ",
        "The size of the `value` vector is not equal to the number of nodes ",
        "selected by the `posit_ids` vector nor of length 1. \n",
        "Expected: ", length(dat[["attr"]][["active"]][posit_ids]), " or 1 \n",
        "Given: ", length(value)
      )
    }

    dat[["attr"]][[item]][posit_ids] <- value
  }

  return(dat)
}

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
      stop("There is no epi output called `",
           paste(missing_item, collapse = ", "),
           "` in the epi output list of the main list object (dat)")
    }

    out <- dat[["epi"]][item]
  }

  return(out)
}

#' @rdname net-accessor
#' @export
get_epi <- function(dat, item, at = NULL, override.null.error = FALSE) {
  if (!item %in% names(dat[["epi"]])) {
    if (override.null.error) {
      out <- NULL
    } else {
      stop("There is no epi out called `", item,
           "` in the epi out list of the main list object (dat)")
    }
  } else {
    if (is.null(at)) {
      out <- dat[["epi"]][[item]]
    } else {
      if (is.logical(at)) {
        if (length(at) != dat[["control"]][["nsteps"]]) {
          stop("(logical) `at` has to have a length equal to the number of
              steps planned for for the simulation (control[['nsteps']])")
        }
      } else if (is.numeric(at)) {
        if (any(at > dat[["control"]][["nsteps"]])) {
          stop("Some (numeric) `at` are larger than the number of
              steps planned for for the simulation (control[['nsteps']])")
        }
      } else {
        stop("`at` must be logical, numeric, or NULL")
      }

      out <- dat[["epi"]][[item]][at]
    }
  }

  return(out)
}

#' @rdname net-accessor
#' @export
add_epi <- function(dat, item) {
  if (item %in% names(dat[["epi"]])) {
    stop("Cannot create the epi output, ", item, ": exists already")
  }

  dat[["epi"]][[item]] <- rep(NA, dat[["control"]][["nsteps"]])

  return(dat)
}

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
      stop("There is no parameters called `",
           paste(missing_item, collapse = ", "),
           "` in the parameter list of the main list object (dat)")
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
      stop("There is no parameter called `", item,
           "` in the parameter list of the main list object (dat)")
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
    stop("Cannot create the parameter, ", item, ": exists already")
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
      stop("There is no control value called `",
           paste(missing_item, collapse = ", "),
           "` in the control list of the main list object (dat)")
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
      stop("There is no control value called `", item,
           "` in the control list of the main list object (dat)")
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
    stop("Cannot create the control value, ", item,
         ": exists already")
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
      stop("There is no init value called `",
           paste(missing_item, collapse = ", "),
           "` in the init list of the main list object (dat)")
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
      stop("There is no init value called `", item,
           "` in the init list of the main list object (dat)")
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
    stop("Cannot create the init value, ", item,
         ": exists already")
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

# Core Attributes --------------------------------------------------------------

#' @rdname net-accessor
#' @export
append_core_attr <- function(dat, at, n.new) {
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- append_attr(dat, "entrTime", at, n.new)
  dat <- append_attr(dat, "exitTime", NA, n.new)

  dat <- update_unique_ids(dat, n.new)

  return(dat)
}

#' @title Create the Unique Identifiers for New Nodes
#'
#' @description This function is called by \code{\link{append_core_attr}} and
#' appends new \code{unique_ids} to the created nodes. It also keeps track of
#' the already used \code{unique_ids} with the \code{dat[["_last_unique_id"]]}
#' variable.
#'
#' @inheritParams recovery.net
#' @param n.new The number of new nodes to give \code{unique_ids} to.
#'
#' @inherit recovery.net return
#'
#' @keywords internal
update_unique_ids <- function(dat, n.new) {
  last_unique_id <- if (is.null(dat[["_last_unique_id"]])) 0L
                    else dat[["_last_unique_id"]]
  next_unique_ids <- seq_len(n.new) + last_unique_id
  dat[["_last_unique_id"]] <- last_unique_id + as.integer(n.new)
  dat <- append_attr(dat, "unique_id", next_unique_ids, n.new)

  return(dat)
}

#' @title Check that All Attributes in the Main Object are of Equal
#'        Length
#'
#' @inheritParams recovery.net
#'
#' @return invisible(TRUE) if everything is correct; an error if not.
#'
#' @keywords internal not_used
check_attr_lengths <- function(dat) {
  attr_lengths <- vapply(dat[["attr"]], length, numeric(1))
  expected_length <- attr_lengths["active"]
  wrong_lengths <- which(attr_lengths != expected_length)

  if (length(wrong_lengths > 0)) {
    msg <- c(
      "Some attribute are not of the correct length \n",
      "Expected length: ", expected_length, "\n",
      "Wrong length attributes: \n"
    )

    for (i in seq_along(wrong_lengths)) {
      msg <- c(msg, "`", names(wrong_lengths)[i], "`: ", wrong_lengths[i], "\n")
    }

    stop(msg)
  }

  return(invisible(TRUE))
}

# Unique / Positional identifier converters ------------------------------------

#' @title Convert Unique Identifiers to/from Positional Identifiers
#'
#' @description EpiModel refers to its nodes either by positional identifiers
#'              (\code{posit_ids}), which describe the position of a node in the
#'              \code{attr} vector, or by unique identifiers
#'              (\code{unique_ids}), which allow references to nodes even after
#'              they are deactivated.
#'
#' @section All elements:
#'   When \code{unique_ids} or \code{posit_ids} is NULL (default)
#'   the full list of positional IDs or unique IDs is returned.
#'
#' @section Deactivated nodes:
#'   When providing \code{unique_ids} of deactivated nodes to
#'   \code{get_posit_ids}, \code{NA}s are returned instead and a warning is
#'   produced.
#'
#' @inheritParams recovery.net
#' @param unique_ids A vector of node unique identifiers (default = NULL).
#' @param posit_ids A vector of node positional identifiers (default = NULL).
#' @return A vector of unique or positional identifiers.
#'
#' @name unique_id-tools
NULL

#' @rdname unique_id-tools
#' @export
get_unique_ids <- function(dat, posit_ids = NULL) {
  if (is.null(posit_ids)) {
    return(get_attr(dat, "unique_id"))
  }

  unique_ids <- get_attr(dat, "unique_id", posit_ids = posit_ids)
  return(unique_ids)
}

#' @rdname unique_id-tools
#' @export
get_posit_ids <- function(dat, unique_ids = NULL) {
  if (is.null(unique_ids)) {
    return(seq_along(get_attr(dat, "active")))
  }
  posit_ids <- base::match(unique_ids, get_attr(dat, "unique_id"))

  if (any(is.na(posit_ids))) {
    warning(
      "While converting `unique_ids` to `posit_ids`, some `unique_ids`",
      " correspond to deactivated nodes and NAs were produced"
    )
  }

  return(posit_ids)
}

#' @title Are These Nodes Active (Unique IDs)
#'
#' @inheritParams recovery.net
#' @param unique_ids A vector of node unique identifiers.
#'
#' @return A logical vector with TRUE if the node is still active and FALSE
#' otherwise.
#'
#' @export
is_active_unique_ids <- function(dat, unique_ids) {
  suppressWarnings({
    posit_ids <- get_posit_ids(dat, unique_ids)
  })

  return(is_active_posit_ids(dat, posit_ids))
}

#' @title Are These Nodes Active (Positional IDs)
#'
#' @inheritParams recovery.net
#' @param posit_ids A vector of node positional identifiers.
#'
#' @return A logical vector with TRUE if the node is still active and FALSE
#' otherwise.
#'
#' @export
is_active_posit_ids <- function(dat, posit_ids) {
  active <- get_attr(dat, "active")
  return(active[posit_ids] %in% 1)
}
