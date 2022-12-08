#' @title Function to Modify the List of Modules During the Simulation
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @details
#' This function removes modules from the module list by reading the
#' ".modules.to.remove" parameter. If the parameter is empty, the function
#' immediatly returns. Otherwise, the modules specified are removed. After the
#' removal, the ".modules.to.remove" parameter is set back to NULL.
#' Therefore, this functionnality is to be used with the "scenarios" API or the
#' ".param.updater" parameter as module removal should happen at specified
#' timestep.
#'
#' @keywords internal
modules_updater <- function(dat) {
  to_remove <- get_param(dat, ".modules.to.remove", override.null.error = TRUE)
  if (!is.null(to_remove)) {
    dat <- remove_modules(dat, to_remove)
    dat <- set_param(dat, ".modules.to.remove", NULL)
 }
  return(dat)
}

#' @title Get the List of Modules
#'
#' @inheritParams recovery.net
#'
#' @return A named list of modules to be run by the model
#'
#' @keywords internal
get_modules <- function(dat) {
  dat[["_modules_list"]]
}


#' @title Set the List of Modules
#'
#' @inheritParams recovery.net
#' @param modules A named list of modules to be run by the model
#'
#' @inherit recovery.net return
#'
#' @keywords internal
set_modules <- function(dat, modules) {
  dat[["_modules_list"]] <- modules
  return(dat)
}

#' @title Populate the Module List After the Initialization Module is run
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @keywords internal
make_module_list <- function(dat) {
  ## Module order
  morder <- get_control(dat, "module.order", override.null.error = TRUE)
  if (is.null(morder)) {
    bi.mods <- get_control(dat, "bi.mods")
    user.mods <- get_control(dat, "user.mods")
    lim.bi.mods <- bi.mods[
      -which(bi.mods %in% c("initialize.FUN", "verbose.FUN"))
      ]
    morder <- c(user.mods, lim.bi.mods)
  }

  ## Make the `modules` list
  modules <- vector(mode = "list", length = length(morder))
  for (i in seq_along(morder)) {
    modules [[i]] <- get_control(dat, morder[i])
  }
  names(modules) <- morder
  dat <- set_modules(dat, modules)

  return(dat)
}


#' @title Remove a Set of Modules From the Module List
#'
#' @inheritParams recovery.net
#' @param names.to.remove a character vector containing the name of the modules
#'   to remove.
#' @inherit recovery.net return
#'
#' @keywords internal
remove_modules <- function(dat, names.to.remove) {
  modules <- get_modules(dat)
  at <- get_current_timestep(dat)

  # Trying to remove modules not present => warning
  missing_mods <- setdiff(names.to.remove, names(modules))
  if (length(missing) > 0) {
    stop(
      "\n\nAt timestep = ", at, ":\n",
      "    an attempt was made to remove the following modules: \n'", missing_mods,
      paste0(missing_mods, collapse = "', '"), "'\n",
      "But they were not in the modules list. Maybe a typo was made"
    )

    names.to.remove <- intersect(names.to.remove, names(modules))
  }

  # In any case => message
  modules <- modules[!names(modules) %in% names.to.remove]
  message(
      "\n\nAt timestep = ", at, ":\n",
      "    the following modules: \n'", missing_mods,
      paste0(missing_mods, collapse = "', '"), "'\n",
      "were disabled."
  )

  dat <- set_modules(dat, modules)
  return(dat)
}
