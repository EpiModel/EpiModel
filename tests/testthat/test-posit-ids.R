context("Unique ID lookup")

test_that("get_posit_ids gives what match() gives", {
  dat <- create_dat_object(control = list(nsteps = 10))
  dat <- append_core_attr(dat, at = 1, n.new = 5000)
  set.seed(1)
  departed_ids <- sort(sample.int(5000L, 1000L))
  dat <- delete_attr(dat, departed_ids)
  all_unique_ids <- get_attr(dat, "unique_id")
  expect_type(all_unique_ids, "integer")

  id_sets <- list(
    few = sample(all_unique_ids, 20),
    many = sample(all_unique_ids, 3000),
    all = rev(all_unique_ids),
    departed = c(sample(all_unique_ids, 10), departed_ids[1:5]),
    out_of_range = c(0L, -3L, 6000L, all_unique_ids[c(1, 4000)]),
    repeated = rep(all_unique_ids[c(7, 7, 300)], 2),
    none = integer(0),
    missing = c(all_unique_ids[3], NA, all_unique_ids[9]),
    not_whole = c(all_unique_ids[5] + 0.5, all_unique_ids[5], Inf, -Inf)
  )
  for (set_name in names(id_sets)) {
    unique_ids <- id_sets[[set_name]]
    expect_identical(suppressWarnings(get_posit_ids(dat, unique_ids)),
                     match(unique_ids, all_unique_ids), info = set_name)
    expect_identical(
      suppressWarnings(get_posit_ids(dat, as.numeric(unique_ids))),
      match(unique_ids, all_unique_ids), info = set_name
    )
  }
})

test_that("get_posit_ids and is_active_unique_ids after departures", {
  dat <- create_dat_object(control = list(nsteps = 10))
  dat <- append_core_attr(dat, at = 1, n.new = 500)
  dat <- delete_attr(dat, c(21:30, 400:420))
  dat <- append_core_attr(dat, at = 2, n.new = 100)
  all_unique_ids <- get_attr(dat, "unique_id")

  expect_warning(
    posit_ids <- get_posit_ids(dat, c(19L, 25L, 410L, 421L)),
    "deactivated"
  )
  expect_identical(posit_ids, c(19L, NA, NA, match(421L, all_unique_ids)))
  expect_identical(is_active_unique_ids(dat, c(19L, 25L, 410L, 421L)),
                   c(TRUE, FALSE, FALSE, TRUE))
})
