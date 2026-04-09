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
tracked_attrs_set_ref <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  if (length(tracked_items) == 0)
    return(dat)
  ref_items <- unique(c("unique_id", "active", tracked_items))
  dat$run$tracked_attrs_ref <- get_attr_list(dat, ref_items)
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
tracked_attrs_record <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  if (length(tracked_items) == 0) {
    return(dat)
  } else if (!dat$run$tracking_attrs) { # Initialize Tracking
    unique_ids <- get_unique_ids(dat)
    for (item in unique(c("active", tracked_items))) {
      value <- get_attr(dat, item)
      dat <- record_attr_history(dat, item, value, unique_ids = unique_ids)
    }
    dat$run$tracking_attrs <- TRUE
    return(dat)
  }

  ref_attrs <- dat$run$tracked_attrs_ref
  tracked_items <- unique(c("active", tracked_items))
  cur_attrs <- get_attr_list(dat, c("unique_id", tracked_items))

  # old nodes - store end time
  departed_uid <- setdiff(ref_attrs$unique_id, cur_attrs$unique_id)
  dat <- record_attr_history(dat, "active", 0L, unique_ids = departed_uid)

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

#' @title Get Tracked Attributes at a Given Time Step
#'
#' @description
#' Reconstructs the state of all tracked attributes at a specific time step
#' from the delta-compressed attribute history. Returns a `tibble` with one
#' row per active node and one column per tracked attribute.
#'
#' @param sim An `EpiModel` object of class `netsim`.
#' @param at The time step at which to reconstruct attributes.
#' @param sim_num The simulation number to use (default 1).
#'
#' @return A `tibble` with columns `unique_id`, `active`, and one column per
#'   tracked attribute, containing the most recent value at or before `at`
#'   for each active node.
#'
#' @seealso [record_attr_history()], [get_attr_history()]
#' @export
get_attr_at <- function(sim, at, sim_num = 1) {
  if (!inherits(sim, "netsim"))
    stop("`sim` must be of class netsim")
  if (length(sim$control$tracked.attributes) < 1)
    stop("At least the `active` attribute must be tracked. ",
         "Check control.net `tracked.attributes` settings.")
  if (sim_num > sim$control$nsims || sim_num < 1)
    stop("Specify a single sim_num between 1 and ", sim$control$nsims)
  if (at < 1 || at > sim$control$nsteps)
    stop("Specify at between 1 and ", sim$control$nsteps)

  attr_hist <- get_attr_history(sim)

  # Get active UIDs at that time
  spells <- get_nodes_spell(attr_hist$active[attr_hist$active$sim == sim_num, ])
  present_uids <- spells$uid[at >= spells$onset & at < spells$terminus]

  d_attrs_at <- dplyr::tibble(unique_id = present_uids)

  for (item in names(attr_hist)) {
    d <- attr_hist[[item]] |>
      dplyr::filter(
        sim == sim_num,
        time <= at,
        uids %in% present_uids
      ) |>
      dplyr::group_by(uids) |>
      dplyr::filter(time == max(time)) |>
      dplyr::ungroup() |>
      dplyr::select(uids, values)
    names(d) <- c("unique_id", item)
    d_attrs_at <- dplyr::left_join(d_attrs_at, d, by = "unique_id")
  }

  d_attrs_at
}

#' @title Compute Node Activity Spells from Active Attribute History
#'
#' @description
#' Derives onset and terminus times for each node from the `active` attribute
#' history records. Uses the earliest `active == 1` record as onset and the
#' earliest `active == 0` record as terminus. Nodes that never departed get
#' `terminus = Inf`.
#'
#' @param d_active A `data.frame` of active attribute history records for a
#'   single simulation, with columns `time`, `uids`, and `values`.
#'
#' @return A `data.frame`:
#'   \describe{
#'     \item{uid}{Integer unique IDs.}
#'     \item{onset}{Numeric onset times.}
#'     \item{terminus}{Numeric terminus times (`Inf` if never departed).}
#'   }
#'
#' @keywords internal
#' @noRd
get_nodes_spell <- function(d_active) {
  d_on <- d_active[d_active$values == 1, ]
  onset <- tapply(d_on$time, d_on$uids, min)

  # Terminus: earliest time each uid appears with active == 0
  d_off <- d_active[d_active$values == 0, ]
  if (nrow(d_off) > 0) {
    terminus <- tapply(d_off$time, d_off$uids, min)
  } else {
    terminus <- numeric(0)
  }

  # Match terminus to onset uids; nodes that never departed get Inf
  all_uids <- as.integer(names(onset))
  term_vals <- terminus[as.character(all_uids)]
  term_vals[is.na(term_vals)] <- Inf

  data.frame(
    uid = all_uids,
    onset = as.numeric(onset),
    terminus = as.numeric(term_vals)
  )
}
