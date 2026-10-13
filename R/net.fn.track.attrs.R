#' @title Snapshot Tracked Attributes at Step Start
#'
#' @description
#' Captures the current values of all attributes listed in
#' `control$tracked.attributes` and `control$tracked.attributes.once` along
#' with their `unique_id`s. This snapshot is stored in
#' `dat$run$tracked_attrs_ref` and used by `tracked_attrs_record` at the end of
#' the step to compute deltas.
#' On its first call (when `dat$run$tracking_attrs` is `FALSE`), records a full
#' initial snapshot instead for the step `at - 1L`. Attributes that do not
#' exist yet (e.g., created by a module on its first call) are left out of the
#' snapshot and recorded in full at the end of the step in which they appear.
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @seealso [tracked_attrs_record()], [record_attr_history()]
#' @keywords internal
tracked_attrs_set_ref <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  tracked_items_once <- get_control(dat, "tracked.attributes.once")
  if (length(tracked_items) == 0 && length(tracked_items_once) == 0) {
    return(dat)
  }

  recorded_items <- unique(c("active", tracked_items, tracked_items_once))
  present_items <- intersect(recorded_items, names(dat$run$attr))

  if (!dat$run$tracking_attrs) { # Initialize Tracking
    unique_ids <- get_unique_ids(dat)
    for (item in present_items) {
      value <- get_attr(dat, item)
      dat <- record_attr_history(
        dat, item, value,
        at = get_current_timestep(dat) - 1L,
        unique_ids = unique_ids
      )
    }
    dat$run$tracking_attrs <- TRUE
  }

  ref_items <- unique(c("unique_id", present_items))
  dat$run$tracked_attrs_ref <- get_attr_list(dat, ref_items)
  return(dat)
}

#' @title Record Tracked Attribute Changes at Step End
#'
#' @description
#' Compares the current attribute values against the reference snapshot taken
#' by `tracked_attrs_set_ref` and records only the changes via
#' `record_attr_history`.
#'
#' Four types of changes are detected:
#' - **New attributes**: attributes absent from the reference but present now
#'   (e.g., created by a module on its first call), recorded for all nodes.
#' - **Departures**: nodes present in the reference but absent now (recorded
#'   as `active = 0L`).
#' - **Arrivals**: nodes absent from the reference but present now (all
#'   tracked attribute values recorded).
#' - **Value changes**: persisting nodes whose attribute value differs from
#'   the reference, including `NA` transitions.
#'
#' A tracked attribute that does not exist at the end of a step is an error.
#'
#' @inheritParams recovery.net
#' @inherit recovery.net return
#'
#' @seealso [tracked_attrs_set_ref()], [record_attr_history()]
#' @keywords internal
tracked_attrs_record <- function(dat) {
  tracked_items <- get_control(dat, "tracked.attributes")
  tracked_items_once <- get_control(dat, "tracked.attributes.once")
  if (length(tracked_items) == 0 && length(tracked_items_once) == 0) {
    return(dat)
  }
  tracked_items <- unique(c("active", tracked_items))
  recorded_items <- unique(c(tracked_items, tracked_items_once))

  missing_items <- setdiff(recorded_items, names(dat$run$attr))
  if (length(missing_items) > 0) {
    stop(
      "`tracked.attributes` and `tracked.attributes.once` name attributes ",
      "that do not exist at the end of time step ", get_current_timestep(dat),
      ": ", paste.and(paste0("`", missing_items, "`")), ". A tracked ",
      "attribute must exist from the end of the first time step on.",
      call. = FALSE
    )
  }

  ref_attrs <- dat$run$tracked_attrs_ref
  cur_attrs <- get_attr_list(dat, unique(c("unique_id", recorded_items)))

  # new attributes - store all nodes
  new_items <- setdiff(recorded_items, names(ref_attrs))
  for (item in new_items) {
    dat <- record_attr_history(
      dat, item, cur_attrs[[item]],
      unique_ids = cur_attrs[["unique_id"]]
    )
  }
  recorded_items <- setdiff(recorded_items, new_items)
  tracked_items <- setdiff(tracked_items, new_items)

  # old nodes - store end time
  departed_uid <- setdiff(ref_attrs$unique_id, cur_attrs$unique_id)
  if (length(departed_uid) > 0L) {
    dat <- record_attr_history(dat, "active", 0L, unique_ids = departed_uid)
  }

  # new nodes - store all tracked
  new_uid <- setdiff(cur_attrs$unique_id, ref_attrs$unique_id)
  if (length(new_uid) > 0L) {
    new_pos <- base::match(new_uid, cur_attrs$unique_id)
    for (item in recorded_items) {
      value <- cur_attrs[[item]][new_pos]
      dat <- record_attr_history(dat, item, value, unique_ids = new_uid)
    }
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
    if (length(chng_pos) == 0L) {
      next
    }
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
#' @param compute_age If `TRUE` and `age` is tracked, compute the age of nodes
#'        at `at` from their last recorded age, assuming age in years and 364
#'        days in a year. This needs a `time.unit` parameter (the number of days
#'        in a time step) in `sim$param`, passed to [`param.net`] or set by the
#'        initialization module (as EpiModelHIV does from its `epistats`).
#'
#' @return A `tibble` with columns `unique_id`, `active`, and one column per
#'   tracked attribute, containing the most recent value at or before `at`
#'   for each active node.
#'
#' @seealso [record_attr_history()], [get_attr_history()]
#' @export
get_attr_at <- function(sim, at, sim_num = 1, compute_age = FALSE) {
  if (!inherits(sim, "netsim"))
    stop("`sim` must be of class netsim")
  if (length(sim$control$tracked.attributes) == 0 &&
        length(sim$control$tracked.attributes.once) == 0) {
    stop("At least the `active` attribute must be tracked. ",
         "Check control.net `tracked.attributes` settings.")
  }
  if (sim_num > sim$control$nsims || sim_num < 1)
    stop("Specify a single sim_num between 1 and ", sim$control$nsims)
  # the history holds the time steps labeled from `epi_start(sim)`, from the
  # first step simulated (`control$start` in a restarted simulation) on
  first <- max(sim$control$start, epi_start(sim))
  last <- epi_start(sim) + sim$control$nsteps - 1
  if (length(at) != 1 || !is.numeric(at) || at < first || at > last)
    stop("Specify a single at between ", first, " and ", last)
  if (compute_age && is.null(sim$param$time.unit))
    stop("`compute_age = TRUE` needs a `time.unit` parameter (the number of ",
         "days in a time step) in `param`")

  attr_hist <- get_attr_history(sim)
  check_attr_history(attr_hist, sim_num)

  # Get active UIDs at that time
  spells <- get_nodes_spell(attr_hist$active[attr_hist$active$sim == sim_num, ])
  present_uids <- spells$uid[at >= spells$onset & at < spells$terminus]

  d_attrs_at <- dplyr::tibble(unique_id = present_uids)

  for (item in names(attr_hist)) {
    if (item == "age" && compute_age) next
    d <- attr_hist[[item]] |>
      dplyr::filter(
        .data$sim == sim_num,
        .data$time <= at,
        .data$uids %in% present_uids
      ) |>
      dplyr::group_by(.data$uids) |>
      dplyr::slice_max(.data$time, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select("uids", "values")
    names(d) <- c("unique_id", item)
    d_attrs_at <- dplyr::left_join(d_attrs_at, d, by = "unique_id")
  }

  if (compute_age && "age" %in% names(attr_hist)) {
    year_steps <- 364 / sim$param$time.unit
    d <- attr_hist[["age"]] |>
      dplyr::filter(
        .data$sim == sim_num,
        .data$time <= at,
        .data$uids %in% present_uids
      ) |>
      dplyr::group_by(.data$uids) |>
      dplyr::slice_max(.data$time, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(age = .data$values + (at - .data$time) / year_steps) |>
      dplyr::select(unique_id = "uids", "age")
    d_attrs_at <- dplyr::left_join(d_attrs_at, d, by = "unique_id")
  }

  d_attrs_at
}

# Stop when a `netsim` object holds no attribute history for a simulation
# @noRd
check_attr_history <- function(attr_hist, sim_num) {
  if (is.null(attr_hist$active) || !any(attr_hist$active$sim == sim_num)) {
    stop("There is no attribute history for simulation ", sim_num, ". ",
         "`make_restart_point()` drops it, and so does `merge()` with ",
         "`keep.attr.history = FALSE`.", call. = FALSE)
  }
  invisible(attr_hist)
}

#' @title Compute Node Activity Spells from Active Attribute History
#'
#' @description
#' Derives the activity spells of each node from the `active` attribute
#' history records. A spell starts at a record where the node becomes active
#' (`active == 1`) and ends at the next record where it becomes inactive
#' (`active == 0`), or at `Inf` if there is none. When a node has several
#' records at the same time step, the last one is used. A node that becomes
#' active again has one spell per period of activity.
#'
#' @param d_active A `data.frame` of active attribute history records for a
#'   single simulation, with columns `time`, `uids`, and `values`.
#'
#' @return A `data.frame` with one row per spell, ordered by `uid` and `onset`:
#'   \describe{
#'     \item{uid}{Unique IDs.}
#'     \item{onset}{Numeric onset times.}
#'     \item{terminus}{Numeric terminus times (`Inf` if never departed).}
#'   }
#'
#' @keywords internal
#' @noRd
get_nodes_spell <- function(d_active) {
  # Unique IDs are matched by value and never by name, as a double ID such as
  # 1e+05 does not print as the integer 100000
  d <- d_active[order(d_active$uids, d_active$time),
                c("uids", "time", "values")]
  d <- d[!duplicated(d[c("uids", "time")], fromLast = TRUE), , drop = FALSE]
  n <- nrow(d)
  # keep the changes of state only
  chg <- c(TRUE, d$uids[-1] != d$uids[-n] | d$values[-1] != d$values[-n])
  d <- d[chg[seq_len(n)], , drop = FALSE]
  n <- nrow(d)
  on <- which(d$values == 1)
  nxt <- pmin(on + 1L, n)
  ends <- on < n & d$uids[nxt] == d$uids[on]

  dplyr::tibble(
    uid = d$uids[on],
    onset = as.numeric(d$time[on]),
    terminus = as.numeric(ifelse(ends, d$time[nxt], Inf))
  )
}

# The terminus of the activity spell of each node `uid` that holds time `at`,
# from the output of `get_nodes_spell()`. For a node with one spell, that
# spell's terminus; `NA` for a node without a spell.
# @noRd
get_spell_terminus <- function(spells, uid, at) {
  pos <- match(uid, spells$uid)
  terminus <- spells$terminus[pos]
  # nodes with more than one spell
  multi <- unique(spells$uid[duplicated(spells$uid)])
  for (u in multi) {
    idx <- which(uid == u)
    s <- spells[spells$uid == u, ]
    k <- findInterval(at[idx], s$onset)
    terminus[idx] <- s$terminus[ifelse(k > 0, k, NA)]
  }
  terminus
}
