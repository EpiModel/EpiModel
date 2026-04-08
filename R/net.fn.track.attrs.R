#' @title Snapshot Tracked Attributes at Step Start
#'
#' @description
#' Captures the current values of all attributes listed in
#' `control$tracked.attributes` along with their `unique_id`s. This snapshot
#' is stored in `dat$run$tracked_attrs_ref` and used by
#' `tracked_attrs_record` at the end of the step to compute deltas.
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @seealso [tracked_attrs_record()], [record_attr_history()]
#' @keywords internal
#' @noRd
tracked_attrs_set_ref <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  if (length(tracked_items) == 0)
    return(dat)
  dat$run$tracked_attrs_ref <- get_attr_list(dat, c("unique_id", tracked_items))
  return(dat)
}

#' @title Record Tracked Attribute Changes at Step End
#'
#' @description
#' Compares the current attribute values against the reference snapshot taken
#' by `tracked_attrs_set_ref` and records only the changes via
#' `record_attr_history`. On its first call (when `dat$run$tracking_attrs`
#' is `FALSE`), records a full initial snapshot instead of a delta.
#'
#' Three types of changes are detected:
#' - **Departures**: nodes present in the reference but absent now (recorded
#'   as `active = 0L` if `"active"` is tracked).
#' - **Arrivals**: nodes absent from the reference but present now (all
#'   tracked attribute values recorded).
#' - **Value changes**: persisting nodes whose attribute value differs from
#'   the reference, including `NA` transitions.
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @seealso [tracked_attrs_set_ref()], [record_attr_history()]
#' @keywords internal
#' @noRd
tracked_attrs_record <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  if (length(tracked_items) == 0) {
    return(dat)
  } else if (!dat$run$tracking_attrs) { # Initialize Tracking
    unique_ids <- get_unique_ids(dat)
    for (item in tracked_items) {
      value <- get_attr(dat, item)
      dat <- record_attr_history(dat, item, value, unique_ids = unique_ids)
    }
    dat$run$tracking_attrs <- TRUE
    return(dat)
  }

  ref_attrs <- dat$run$tracked_attrs_ref
  cur_attrs <- get_attr_list(dat, c("unique_id", tracked_items))

  # old nodes - store end time
  if ("active" %in% tracked_items) {
    departed_uid <- setdiff(ref_attrs$unique_id, cur_attrs$unique_id)
    dat <- record_attr_history(dat, "active", 0L, unique_ids = departed_uid)
  }

  # new nodes - store all tracked
  new_uid <- setdiff(cur_attrs$unique_id, ref_attrs$unique_id)
  new_pos <- base::match(new_uid, cur_attrs$unique_id)
  for (item in tracked_items) {
    value <- cur_attrs[[item]][new_pos]
    dat <- record_attr_history(dat, item, value, unique_ids = new_uid)
  }
  # persisting nodes - store diff
  pers_uid <- intersect(cur_attrs$unique_id, ref_attrs$unique_id)
  pers_ref_pos <- base::match(pers_uid, ref_attrs$unique_id)
  pers_cur_pos <- base::match(pers_uid, cur_attrs$unique_id)

  for (item in tracked_items) {
    ref_v <- ref_attrs[[item]][pers_ref_pos]
    cur_v <- cur_attrs[[item]][pers_cur_pos]
    # Must take into account NA -> non_NA and vice-versa
    chng_pos <- which(ref_v != cur_v | is.na(ref_v) != is.na(cur_v))
    chng_uid <- cur_attrs[["unique_id"]][pers_cur_pos][chng_pos]
    chng_value <- cur_v[chng_pos]
    dat <- record_attr_history(dat, item, chng_value, unique_ids = chng_uid)
  }
  dat$run$tracked_attrs_ref <- NULL
  return(dat)
}
