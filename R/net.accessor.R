#' @title Functions to Access and Edit the Master List Object in Network Models
#'
#' @description These \code{get_}, \code{set_}, \code{append_} and \code{add} functions allow a safe
#'              and efficient way to retrieve and mutate the Master list object
#'              of network models (\code{dat}).
#'
#' @param dat a Master list object of network models
#' @param item a character vector conaining the name of the element to access.
#'        Can be of length > 1 for \code{get_*_list} functions
#' @param indexes for \code{get_epi} and \code{get_attr}, a numeric vector of
#'        indexes or a logical vector to subset the desired \code{item}
#' @param value new value to be attributed in the \code{set_} and \code{append_}
#'        functions
#' @param override.null.error if TRUE, \code{get_} return NULL if the
#'         \code{item} does not exist instead of throwing an error.
#'         (default = FALSE)
#' @param override.length.check if TRUE, \code{set_attr} allows the modification
#'        of the \code{item} size. (default = FALSE)
#' @return a vector or a list of
#' vector for \code{get_} functions. And the Master list object for \code{set_}
#' and \code{add_} functions
#'
#' @section Core Attribute:
#' The \code{append_core_attr} function initialize the attributes necessary for
#' EpiModel to work (Currently "active" and "uid"). It is used in the
#' inilization phase of the simulation, to create the nodes (see
#' \code{initialize.net}). It is also used when adding nodes during the simulation
#' (see \code{arrival.net})
#'
#' @section Mutability:
#' The \code{set_}, \code{append_} and \code{add_} functions DO NOT modify the
#' dat object in place. The result must be assigned back to \code{dat} in order
#' to be registered \code{dat <- set_*(dat, item, value)}
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
           "` in the attributes list of the Master list object (dat)")
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
      stop("There is no attribute called `", item,
           "` in the attributes list of the Master list object (dat)")
    }
  } else {
    if (is.null(indexes)) {
      out <- dat[["attr"]][[item]]
    } else {
      if (is.logical(indexes)) {
        if (length(indexes) != length(dat[["attr"]][[item]])) {
          stop("(logical) `indexes` has to have a length equal to the number ",
               "of nodes in the network")
        }
      } else if (is.numeric(indexes)) {
        if (any(indexes > length(dat[["attr"]][[item]]))) {
          stop("Some (numeric) `indexes` are larger than the number of nodes ",
               " in the network")
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
    stop("Cannot create the attribute '", item, "': exists already")
  }

  dat[["attr"]][[item]] <- rep(NA, length(dat[["attr"]][["active"]]))

  return(dat)
}

#' @rdname net-accessor
#' @export
set_attr <- function(dat, item, value, indexes = NULL,
                     override.length.check = FALSE) {
  if (!item %in% names(dat[["attr"]])) {
    dat <- add_attr(dat, item)
  }

  if (is.null(indexes)) {
    dat <- set_attr_full(dat, item, value, override.length.check)
  } else {
    dat <- set_attr_indexes(dat, item, value, indexes, override.length.check)
  }

  return(dat)
}

#' Version of \code{set_attr} without indexes
#' keywords internal
set_attr_full <- function(dat, item, value, override.length.check) {
  if (!override.length.check &&
    length(value) != length(dat[["attr"]][["active"]])) {
    stop(
      "When trying to edit the ", `item`, " nodal attribute: ",
      "The size of the `value` vector is not equal to the number of node in",
      " the network. \n",
      "Expected: ", length(dat[["attr"]][["active"]]), "\n",
      "Given: ", length(value)
    )
  }

  dat[["attr"]][[item]] <- value

  return(dat)
}

#' Version of \code{set_attr} with indexes
#' keywords internal
set_attr_indexes <- function(dat, item, value, indexes, override.length.check) {
  if (is.logical(indexes)) {
    if (length(indexes) != length(dat[["attr"]][[item]])) {
      stop("(logical) `indexes` has to have a length equal to the number ",
        "of nodes in the network")
    }
  } else if (is.numeric(indexes)) {
    if (any(indexes > length(dat[["attr"]][[item]]))) {
      stop("Some (numeric) `indexes` are larger than the number of nodes ",
        " in the network")
    }
  } else {
    stop("`indexes` must be logical, numeric, or NULL")
  }

  if (!override.length.check &&
      length(value) != 1 &&
      length(value) != length(dat[["attr"]][["active"]][indexes])) {
    stop(
      "When trying to edit the ", `item`, " nodal attribute: ",
      "The size of the `value` vector is not equal to the number of node ",
      "selected by the `indexes` vector nor of length 1. \n",
      "Expected: ", length(dat[["attr"]][["active"]][indexes]), " or 1 \n",
      "Given: ", length(value)
    )
  }

  dat[["attr"]][[item]][indexes] <- value

  return(dat)
}

#' @param n.new the number of new elements to append at the end of \code{item}
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
           "` in the epi output list of the Master list object (dat)")
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
      stop("There is no epi out called `", item,
           "` in the epi out list of the Master list object (dat)")
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
    stop("Cannot create the epi output, ", item, ": exists already")
  }

  dat[["epi"]][[item]] <- rep(NA, dat[["control"]][["nsteps"]])

  return(dat)
}

#' @param at timestep where to add the new value for the epi output \code{item}
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
           "` in the parameter list of the Master list object (dat)")
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
           "` in the parameter list of the Master list object (dat)")
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
           "` in the control list of the Master list object (dat)")
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
           "` in the control list of the Master list object (dat)")
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
           "` in the init list of the Master list object (dat)")
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
           "` in the init list of the Master list object (dat)")
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

#' @param n.new the number of new nodes to initiate with core attributes
#' @param at current time step
#' @rdname net-accessor
#' @export
append_core_attr <- function(dat, at, n.new) {
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- append_attr(dat, "entrTime", at, n.new)
  dat <- append_attr(dat, "exitTime", NA, n.new)

  dat <- update_uids(dat, n.new)

  return(dat)
}

#' @title Create the uids for the new nodes
#'
#' @description This function is called by `append_core_attr` and append new
#' uids to the created nodes. It also keeps track of the already used uids with
#' the /code{dat[["_last_uid"]]} variable
#'
#' @param dat a Master list object of network models
#' @param n.new the number of new nodes to give \code{uid} to
#'
#' @return the Master list object of network models (\code{dat})
#'
#' @keywords internal
update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0L else dat[["_last_uid"]]
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + as.integer(n.new)
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}

#' @title Check that all \code{attr}ibutes in the master object are of equal
#'        length
#'
#' @param dat a Master list object of network models
#'
#' @return invisible(TRUE) if everythin is correct. It throws an error otherwise
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
