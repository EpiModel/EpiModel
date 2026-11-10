context("Unique ID lookup")

# get_posit_ids() looks unique IDs up with match_unique_ids(), which has to
# return what match() returns for every input.

test_that("match_unique_ids gives what match() gives", {
  set.seed(1)
  uid <- sort(sample.int(5000L, 4000L))
  gone <- setdiff(seq_len(5000L), uid)
  cases <- list(
    few = sample(uid, 20),
    many = sample(uid, 3000),
    departed = c(sample(uid, 10), gone[1:5]),
    out.of.range = c(0L, -3L, 6000L, uid[1], uid[4000]),
    repeated = rep(uid[c(7, 7, 300)], 2),
    none = integer(0),
    missing = c(uid[3], NA, uid[9])
  )
  for (nm in names(cases)) {
    ids <- cases[[nm]]
    for (table in list(uid, as.numeric(uid))) {
      expect_identical(match_unique_ids(ids, table), match(ids, table),
                       info = nm)
      expect_identical(match_unique_ids(as.numeric(ids), table),
                       match(as.numeric(ids), table), info = nm)
    }
  }
  # values that are not whole numbers
  odd <- c(uid[5] + 0.5, uid[5], Inf, -Inf, NaN)
  expect_identical(match_unique_ids(odd[1:4], uid), match(odd[1:4], uid))
  expect_identical(match_unique_ids(odd, uid), match(odd, uid))
  # a table out of order or holding NA, an empty table, and IDs of another type
  shuffled <- sample(uid)
  expect_identical(match_unique_ids(cases$few, shuffled),
                   match(cases$few, shuffled))
  with.na <- c(uid[1:10], NA, uid[11:20])
  expect_identical(match_unique_ids(uid[c(3, 15)], with.na),
                   match(uid[c(3, 15)], with.na))
  expect_identical(match_unique_ids(1:3, integer(0)), match(1:3, integer(0)))
  expect_identical(match_unique_ids(as.character(uid[1:3]), uid),
                   match(as.character(uid[1:3]), uid))
})

test_that("get_posit_ids and is_active_unique_ids after departures", {
  dat <- create_dat_object(control = list(nsteps = 10))
  dat <- append_core_attr(dat, at = 1, n.new = 500)
  dat <- delete_attr(dat, c(21:30, 400:420))
  dat <- append_core_attr(dat, at = 2, n.new = 100)
  uid <- get_attr(dat, "unique_id")
  expect_type(uid, "integer")

  # a few IDs, many IDs, and all of them
  expect_identical(get_posit_ids(dat, uid[c(5, 480, 20, 21)]),
                   c(5L, 480L, 20L, 21L))
  expect_identical(get_posit_ids(dat, rev(uid)), rev(seq_along(uid)))
  expect_identical(get_posit_ids(dat, as.numeric(uid[1:300])), 1:300)

  # departed nodes
  expect_warning(pid <- get_posit_ids(dat, c(19L, 25L, 410L, 421L)),
                 "deactivated")
  expect_identical(pid, c(19L, NA, NA, match(421L, uid)))
  expect_identical(is_active_unique_ids(dat, c(19L, 25L, 410L, 421L)),
                   c(TRUE, FALSE, FALSE, TRUE))
})
