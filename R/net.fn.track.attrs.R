init_track_attrs <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  nsteps <- get_control(dat, "nsteps")
  dat$run$tracked_attributes <- padded_vector(list(), nsteps)
  stored <- list()
  unique_ids <- get_unique_ids(dat)
  for (item in c("active", tracked_items)) {
    stored[[item]]$uid <- unique_ids
    stored[[item]]$value <- get_attr(dat, item)
  }
  dat$run$tracked_attributes[[get_current_timestep(dat)]] <- stored
  return(dat)
}

start_track_attrs_step <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  if (length(tracked_items) == 0)
    return(dat)

  dat$run$track_refs <- get_attr_list(dat, c("unique_id", tracked_items))

  return(dat)
}

stop_track_attrs_step <- function(dat) {
  if (is.null(dat$run$track_refs)) {
    return(dat)
  } else if (is.null(dat$run$tracked_attributes)) {
    dat <- init_track_attrs(dat)
    return(dat)
  } else if (length(dat$run$tracked_attributes) < get_control(dat, "nsteps")) {
    dat$run$tracked_attributes <-
      padded_vector(dat$run$tracked_attributes, get_control(dat, "nsteps"))
  }

  tracked_items <- get_control(dat, "tracked.attributes")
  refs_attrs <- dat$run$track_refs
  cur_attrs <- get_attr_list(dat, c("unique_id", tracked_items))
  stored <- list()

  # old nodes - store end time
  departed_uid <- setdiff(refs_attrs$unique_id, cur_attrs$unique_id)
  # new nodes - store all tracked
  new_uid <- setdiff(cur_attrs$unique_id, refs_attrs$unique_id)
  new_pos <- get_posit_ids(dat, new_uid)

  # Store arrivals and exits (`active`)
  stored$active <- list(
    uid = c(departed_uid, new_uid),
    value = c(rep(0L, length(departed_uid)), rep(1L, length(new_uid)))
  )

  # cur nodes - store diff
  cur_uid <- intersect(cur_attrs$unique_id, refs_attrs$unique_id)
  cur_ref_pos <- base::match(cur_uid, refs_attrs$unique_id)
  cur_pos <- get_posit_ids(dat, cur_uid)

  for (item in tracked_items) {
    changed_pos <- which(refs_attrs[[item]][cur_ref_pos] != cur_attrs[[item]][cur_pos])
    all_pos <- c(changed_pos, new_pos)
    stored[[item]]$uid <- get_unique_ids(dat, all_pos)
    stored[[item]]$value <- cur_attrs[[item]][all_pos]
  }

  dat$run$tracked_attributes[[get_current_timestep(dat)]] <- stored
  dat$run$track_refs <- NULL
  return(dat)
}
