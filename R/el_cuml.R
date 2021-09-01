get_edgelist <- function(dat, network) {
  if (get_control(dat, "tergmLite")) {
    el <- dat[["el"]][[network]]
  } else {
    el <- network::as.edgelist(dat[["nw"]][[network]])
  }

  return(el)
}

get_cumulative_edgelist <- function(dat, network) {
  el_cuml <- dat[["el_cuml"]][[network]]

  if(is.null(el_cuml)) {
    el_cuml <- data.frame(
      head = numeric(0),
      tail = numeric(0),
      start = numeric(0),
      stop = numeric(0)
    )
  }

  return(el_cuml)
}

update_cumulative_edgelist <- function(dat, network) {
  truncate.el_cuml <- get_control(
    dat, "truncate.el_cuml", override.null.error = FALSE)
  truncate.el_cuml <- if (is.na(truncate.el_cuml)) Inf else truncate.el_cuml

  el <- get_edgelist(dat, network)
  el_cuml <- get_edgelist(dat, network)

  el <- data.frame(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = T
  )

  el_cuml <- merge(el_cuml, el, by = c("head", "tail"), all = T)
  el_cuml[is.na(el_cuml[["start"]]), ]$start <- at
  el_cuml[is.na(el_cuml[["current"]]), ]$stop <- at

  if (truncate.el_cuml != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate.el_cuml, ]
  }

  dat[["el_cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]

  return(dat)
}

update_cumulative_edgelist <- function(dat, at, network) {
  truncate.el_cuml <- get_control(
    dat, "truncate.el_cuml", override.null.error = FALSE)
  truncate.el_cuml <- if (is.na(truncate.el_cuml)) Inf else truncate.el_cuml

  el <- get_edgelist(dat, network)
  el_cuml <- get_edgelist(dat, network)

  el <- data.frame(
    head = get_unique_ids(dat, el[, 1]),
    tail = get_unique_ids(dat, el[, 2]),
    current = T
  )

  el_cuml <- merge(el_cuml, el, by = c("head", "tail"), all = TRUE)

  new_edges <- is.na(el_cuml[["start"]])
  if (any(new_edges)) {
    el_cuml[new_edges, ]$start <- at
  }

  terminated_edges <- is.na(el_cuml[["current"]])
  if (any(terminated_edges)) {
    el_cuml[terminated_edges, ]$stop <- at - 1
  }

  if (truncate.el_cuml != Inf) {
    rel.age <- at - el_cuml[["stop"]]
    rel.age <- ifelse(is.na(rel.age), 0, rel.age)
    el_cuml <- el_cuml[rel.age <= truncate.el_cuml, ]
  }

  dat[["el_cuml"]][[network]] <- el_cuml[, c("head", "tail", "start", "stop")]

  return(dat)
}
