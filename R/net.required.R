#' @keywords internal
append_required_attr <- function(dat, n.new) {
	required_attr <- c("active", "uid")

	for (item in required_attr) {
		if (!item %in% names(dat[["attr"]])) {
			dat <- add_attr(dat, item)
		}
	}

  dat <- append_attr(dat, "active", 1, n.new)
  dat <- update_uids(dat, n.new)

	return(dat)
}

#' @keywords internal
update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0 else dat[["_last_uid"]]
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + n.new
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}


