#' Helper functions to access and edit the Master list object of network models
#'
#' These `get_`, `set_` and `add` functions allow a safe and efficient way to
#' retrieve and mutate the Master list object of network models (`dat`).
#'
#' @section mutability:
#' The `set_` and `add_` functions DO NOT mutate the dat object in place.
#' The result must be assigned back to `dat` in order to be registered
#' `dat <- set_*(dat, item, value)`
#'
#' @section `set_` vs `add_`:
#' The `set_` functions edit a pre-existing element or create a new one if it
#' does not exist already by calling the `add_` functions internally.
#'
#' @param dat a Master list object of network models
#' @param item a character vector conaining the name of the element to access.
#'        Can be of length > 1 for `get_*_list` functions
#' @param indexes for `get_epi` and `get_attr`, a numeric vector of indexes or
#'        a logical vector to subset the desired `item`
#' @return a vector or a list of vector for `get_` functions. And the Master
#'         list object for `set_` and `add_` functions
#'
#' @name dat_get_set
NULL

#' @rdname dat_get_set
#' @export
get_attr_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["attr"]]

  } else {
    missing_item <- setdiff(item, names(dat[["attr"]]))
    if (length(missing_item) > 0)
      stop(paste("There is no attributes called",
                 paste(missing_item, collapse = ", "),
                 "in the attributes list of the Master list object (dat)"))

    out <- dat[["attr"]][item]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
get_attr <- function(dat, item, indexes = NULL) {
  if (!item %in% names(dat[["attr"]]))
      stop(paste("There is no attribute called", item,
                 "in the attributes list of the Master list object (dat)"))

  if (is.null(indexes)) {
    out <- dat[["attr"]][[item]]

  } else {
    if (is.logical(indexes)) {
      if (length(indexes) != length(dat[["attr"]][[item]]))
        stop("(logical) `indexes` has to have a length equal to the number of
              nodes in the network")
    } else if(is.numeric(indexes)) {
      if (any(indexes > length(dat[["attr"]][[item]])))
        stop("Some (numeric) `indexes` are larger than the number of nodes in
              the network")
    } else {
      stop("`indexes` must be logical, numeric, or NULL")
    }

    out <- dat[["attr"]][[item]][indexes]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
add_attr <- function(dat, item) {
  if (item %in% names(dat[["attr"]]))
    stop(paste0("Cannot create the attribute '", item,
               "': exists already"))

  dat[["attr"]][[item]] <- rep(NA, length(dat$attr$active))

  return(dat)
}

#' @rdname dat_get_set
#' @export
set_attr <- function(dat, item, value) {
  if (!item %in% names(dat[["attr"]])) {
    message(paste0("Creating attribute: '", item,
               "' in the attributes list of the Master list object (dat)"))


    dat <- add_attr(dat, item)
  }

  if (length(value) != length(dat$attr$active))
    stop(paste0(
      "When trying to edit the ", `item`, " nodale attribute: The size",
       " of the `value` vector is not equal to the number of node in
       the network. Expected: ", length(dat$attr$active), ", given : ",
       length(value)))

  dat[["attr"]][[item]] <- value

  return(dat)
}

#' @rdname dat_get_set
#' @export
get_epi_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["epi"]]

  } else {
    missing_item <- setdiff(item,names(dat[["epi"]]))
    if (length(missing_item) > 0)
      stop(paste("There is no Epidemic output called",
                 paste(missing_item, collapse = ", "),
                 "in the Epidemic output list of the Master list object (dat)"))

    out <- dat[["epi"]][item]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
get_epi <- function(dat, item, indexes = NULL) {
  if (!item %in% names(dat[["epi"]]))
      stop(paste("There is no Epidemic output called", item,
                 "in the Epidemic output list of the Master list object (dat)"))

  if (is.null(indexes)) {
    out <- dat[["epi"]][[item]]

  } else {
    if (is.logical(indexes)) {
      if (length(indexes) != dat$control$nsteps)
        stop("(logical) `indexes` has to have a length equal to the number of
              steps planned for for the simulation (control$nsteps)")
    } else if(is.numeric(indexes)) {
      if (any(indexes > dat$control$nsteps))
        stop("Some (numeric) `indexes` are larger than the number of
              steps planned for for the simulation (control$nsteps)")
    } else {
      stop("`indexes` must be logical, numeric, or NULL")
    }

    out <- dat[["epi"]][[item]][indexes]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
add_epi <- function(dat, item) {
  if (item %in% names(dat[["epi"]]))
    stop(paste("Cannot create the Epidemic outpout, ", item,
               ": exists already"))

  dat[["epi"]][[item]] <- rep(NA, dat$control$nsteps)

  return(dat)
}

#' @rdname dat_get_set
#' @export
set_epi <- function(dat, item, value) {
  if (!item %in% names(dat[["epi"]])) {
    message(paste("Creating Epidemic output: '", item,
               "'in the Epidemic output list of the Master list object (dat)"))

    dat <- add_epi(dat, item)
  }

  dat[["epi"]][[item]] <- value

  return(dat)
}

#' @param at timestep where to add the new value for the Epidemic Outuput `item`
#' @rdname dat_get_set
#' @export
set_epi_at <- function(dat, item, at,  value) {
  if (length(at) != 1 || !is.numeric(at))
    stop("`at` must be numeric and of length one")

  if (!item %in% names(dat[["epi"]])) {
    message(paste("Creating Epidemic output: '", item,
               "'in the Epidemic output list of the Master list object (dat)"))

    dat <- add_epi(dat, item)
  }

  if (at > length(dat[["epi"]][[item]])) {

      dat[["epi"]][[item]] <- c(
        dat[["epi"]][[item]],
        rep(NA, dat$control$nsteps - length(dat[["epi"]][[item]]))
      )
  }

  dat[["epi"]][[item]][at] <- value

  return(dat)
}

#' @rdname dat_get_set
#' @export
get_param_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["param"]]

  } else {
    missing_item <- setdiff(item,names(dat[["param"]]))
    if (length(missing_item) > 0)
      stop(paste("There is no parameters called",
                 paste(missing_item, collapse = ", "),
                 "in the parameter list of the Master list object (dat)"))

    out <- dat[["param"]][item]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
get_param <- function(dat, item) {
  if (!item %in% names(dat[["param"]]))
      stop(paste("There is no parameter called", item,
                 "in the parameter list of the Master list object (dat)"))

  out <- dat[["param"]][[item]]

  return(out)
}

#' @rdname dat_get_set
#' @export
add_param <- function(dat, item) {
  if (item %in% names(dat[["param"]]))
    stop(paste("Cannot create the parameter, ", item,
               ": exists already"))

  dat[["param"]][[item]] <- NA

  return(dat)
}

#' @rdname dat_get_set
#' @export
set_param <- function(dat, item, value) {
  if (!item %in% names(dat[["param"]])) {
    message(paste("Creating parameter: '", item,
               "'in the parameter list of the Master list object (dat)"))

    dat <- add_param(dat, item)
  }

  dat[["param"]][[item]] <- value

  return(dat)
}

#' @rdname dat_get_set
#' @export
get_control_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["control"]]

  } else {
    missing_item <- setdiff(item,names(dat[["control"]]))
    if (length(missing_item) > 0)
      stop(paste("There is no control value called",
                 paste(missing_item, collapse = ", "),
                 "in the control list of the Master list object (dat)"))

    out <- dat[["control"]][item]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
get_control <- function(dat, item) {
  if (!item %in% names(dat[["control"]]))
      stop(paste("There is no control value called", item,
                 "in the control list of the Master list object (dat)"))

  out <- dat[["control"]][[item]]

  return(out)
}

#' @rdname dat_get_set
#' @export
add_control <- function(dat, item) {
  if (item %in% names(dat[["control"]]))
    stop(paste("Cannot create the control value, ", item,
               ": exists already"))

  dat[["control"]][[item]] <- NA

  return(dat)
}

#' @rdname dat_get_set
#' @export
set_control <- function(dat, item, value) {
  if (!item %in% names(dat[["control"]])) {
    message(paste("Creating control value: '", item,
               "'in the control list of the Master list object (dat)"))

    dat <- add_control(dat, item)
  }

  dat[["control"]][[item]] <- value

  return(dat)
}

#' @rdname dat_get_set
#' @export
get_init_list <- function(dat, item = NULL) {
  if (is.null(item)) {
    out <- dat[["init"]]

  } else {
    missing_item <- setdiff(item,names(dat[["init"]]))
    if (length(missing_item) > 0)
      stop(paste("There is no init value called",
                 paste(missing_item, collapse = ", "),
                 "in the init list of the Master list object (dat)"))

    out <- dat[["init"]][item]
  }

  return(out)
}

#' @rdname dat_get_set
#' @export
get_init <- function(dat, item) {
  if (!item %in% names(dat[["init"]]))
      stop(paste("There is no init value called", item,
                 "in the init list of the Master list object (dat)"))

  out <- dat[["init"]][[item]]

  return(out)
}

#' @rdname dat_get_set
#' @export
add_init <- function(dat, item) {
  if (item %in% names(dat[["init"]]))
    stop(paste("Cannot create the init value, ", item,
               ": exists already"))

  dat[["init"]][[item]] <- NA

  return(dat)
}

#' @rdname dat_get_set
#' @export
set_init <- function(dat, item, value) {
  if (!item %in% names(dat[["init"]])) {
    message(paste("Creating init value: '", item,
               "'in the init list of the Master list object (dat)"))

    dat <- add_init(dat, item)
  }

  dat[["init"]][[item]] <- value

  return(dat)
}

