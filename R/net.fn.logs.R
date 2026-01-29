#' @title Message to Find in Which Module a `condition` Occurred
#'
#' @description This function returns a formatted string describing when, where,
#'              and why an error, message, or warning occurred.
#'
#' @param cond The type of `condition` handled (message, warning, error).
#' @param module The name of the module where the `condition` occurred.
#' @param at The time step the `condition` occurred.
#' @param msg The `condition`'s message.
#'
#' @return A formatted string describing where and when the `condition`
#'         occurred as well as the `condition`'s message.
#'
#' @keywords internal
netsim_cond_msg <- function(cond, module, at, msg) {
  paste0("\n\tA ", cond, " occured in module '", module, "' at step ", at)
}

#'  Handle the Logging of Traceback and Dumping of Frames on Error
#'
#'  If `control$.traceback.on.error == TRUE`, this function prints the traceback
#'  of the current simulation to STDIN. This is useful when `ncores > 1` or in
#'  HPC settings.
#'  If `control$.dump.frames.on.error == TRUE`, this function saves a debugging
#'  dump for "postmortem debugging". The dumps are named
#'  "dump_%Y%m%d_%H%M%S_s.rda" and stored at the root of the working directory.
#'
#' @inheritParams recovery.net
#' @param s The number of the simulation that failed
#'
#' @return Nothing, after logging and dumping frames, the function gives the
#'   control back to the general error handler
#'
#' @keywords internal
netsim_error_logger <- function(dat, s) {
  if (get_control(dat, ".traceback.on.error")) {
    message("\n",
      "***************\n",
      "** TRACEBACK **\n",
      "***************"
    )
    traceback(0)
  }

  if (get_control(dat, ".dump.frame.on.error")) {
    dump_name <- format(Sys.time(), format = "dump_%Y%m%d_%H%M%S")
    dump_name <- paste0(dump_name, "_", s, ".rda")
    star_header <- paste0(rep("*", nchar(dump_name)))
    message("\n",
      star_header, "\n",
      "DUMP FILE:\n",
      dump_name, "\n",
      star_header, "\n"
    )
    utils::dump.frames()
    save.image(file = dump_name)
  }
}