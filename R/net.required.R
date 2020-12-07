create_required_attr <- function(dat, n.new) {
  dat <- add_attr(dat, "active")
  dat <- add_attr(dat, "uid")
  dat <- append_required_attr(dat, n.new)

  return(dat)
}

append_required_attr <- function(dat, n.new) {
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- update_uids(dat, n)

  return(dat)
}

update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0 else last_uid
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + n.new
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}

