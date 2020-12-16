#' @param n.new the number of new nodes to initiate with core attributes
#' @rdname net-accessor
#' @export
append_core_attr <- function(dat, n.new) {
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- update_uids(dat, n.new)

	return(dat)
}

#' @title Create the uids for the new nodes
#'
#' @description This function is called by `append_core_attr` and append new
#' uids to the created nodes. It also keeps track of the already used uids with
#' the `dat[["_last_uid"]]` variable
#'
#' @param dat a Master list object of network models
#' @param n.new the number of new nodes to give `uid` to
#'
#' @return the Master list object of network models (`dat`)
#'
#' @keywords internal
update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0L else dat[["_last_uid"]]
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + n.new
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}

#' @title Check that all `attr`ibutes in the master object are of equal length
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

