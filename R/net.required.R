create_required_attr <- function(dat, n.new) {
  required_attr <- list_required_attr()
  for (i in seq_along(required_attr)) {
    dat <- add_attr(dat, names(required_attr)[i])
  }

  append_required_attr(dat, n.new)

  return(dat)
}

append_required_attr <- function(dat, n.new) {
  required_attr <- list_required_attr()
  for (i in seq_along(required_attr)) {
    dat <- append_attr(
      dat,
      names(required_attr)[i], required_attr[i],
      n.new
    )
  }

  dat <- update_uids(dat, n)
}

update_uids <- function(dat) {
}

is_attr_required <- function(item) {
  return(item %in% names(list_required_attr))
}

list_required_attr <- function(dat) {
  required_attr <- list(
    active = 1,
    uid = NA
  )
}

get_next_uids <- function(dat, n) {
  last_uid <- dat[["_last_uid"]]
  last_uid <- if (is.null(last_uid)) 0 else last_uid
  next_uids <- seq_len(n) + last_uid

  return(next_uids)
}
