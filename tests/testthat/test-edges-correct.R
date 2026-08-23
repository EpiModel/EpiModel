# Helpers ---------------------------------------------------------------------
#
# `edges_correct()` reads `active` (and `group`), the `resimulate.network` and
# `edges.correct.attr` controls, the `groups` param, and the previous-step
# counts kept in `dat$run`. It writes the first formation coefficient of every
# network. These build the smallest object carrying all of that.

mk_ec_dat <- function(n_active, n_total = max(n_active, 10), n_eligible = NULL,
                      eligible = NULL, n_group2 = 0, groups = 1,
                      attr.name = NULL, run = list(), resim = TRUE) {
  control <- list(resimulate.network = resim)
  if (!is.null(attr.name)) {
    control$edges.correct.attr <- attr.name
  }
  dat <- create_dat_object(param = list(groups = groups),
                           control = control, run = run)
  # The attribute list starts empty, so the first write sets the length.
  dat <- set_attr(dat, "active",
                  c(rep(1, n_active), rep(0, n_total - n_active)),
                  override.length.check = TRUE)
  if (is.null(eligible) && !is.null(n_eligible)) {
    eligible <- c(rep(1, n_eligible), rep(0, n_total - n_eligible))
  }
  if (!is.null(eligible)) {
    dat <- set_attr(dat, "active.sex", eligible)
  }
  if (groups == 2) {
    dat <- set_attr(dat, "group",
                    c(rep(1, n_total - n_group2), rep(2, n_group2)))
  }
  dat$num.nw <- 2L
  dat$nwparam <- list(list(coef.form = c(-3.5, 1, 2)),
                      list(coef.form = c(-4.25, 5)))
  dat
}

# The first formation coefficient of each network, which is what the correction
# moves.
coef_form <- function(dat) {
  vapply(dat$nwparam, function(x) x$coef.form[1], numeric(1))
}


test_that("edges_correct is unchanged by default", {
  # The function exactly as it stood before `edges.correct.attr` was added.
  # Backwards compatibility is asserted against this rather than against
  # hand-computed numbers, so the test cannot drift away from the old behavior.
  edges_correct_legacy <- function(dat, at) {
    resimulate.network <- get_control(dat, "resimulate.network")
    groups <- get_param(dat, "groups")
    active <- get_attr(dat, "active")
    if (resimulate.network == TRUE) {
      if (groups == 1) {
        old.num <- dat$run$num
        new.num <- sum(active == 1)
        adjustment <- log(old.num) - log(new.num)
      }
      if (groups == 2) {
        old.num.g1 <- dat$run$num
        old.num.g2 <- dat$run$num.g2
        group <- get_attr(dat, "group")
        new.num.g1 <- sum(active == 1 & group == 1)
        new.num.g2 <- sum(active == 1 & group == 2)
        adjustment <-
          log(2 * old.num.g1 * old.num.g2 / (old.num.g1 + old.num.g2)) -
          log(2 * new.num.g1 * new.num.g2 / (new.num.g1 + new.num.g2))
      }
      for (network in seq_len(dat$num.nw)) {
        dat$nwparam[[network]]$coef.form[1] <-
          dat$nwparam[[network]]$coef.form[1] + adjustment
      }
    }
    return(dat)
  }

  # One group, across a grid of growing and shrinking populations.
  for (n_active in c(50, 100, 137, 400)) {
    for (old_num in c(50, 100, 137, 400)) {
      dat <- mk_ec_dat(n_active = n_active, run = list(num = old_num))
      # The whole network-parameter list, not just the coefficient that moves,
      # so an unintended write anywhere in it fails the test.
      expect_equal(edges_correct(dat, 2)$nwparam,
                   edges_correct_legacy(dat, 2)$nwparam)
    }
  }

  # Two groups.
  for (g1 in c(40, 90)) {
    for (g2 in c(60, 110)) {
      dat <- mk_ec_dat(n_active = 100, n_group2 = 80, groups = 2,
                       run = list(num = g1, num.g2 = g2))
      expect_equal(edges_correct(dat, 2)$nwparam,
                   edges_correct_legacy(dat, 2)$nwparam)
    }
  }

  # And it stays a no-op when the network is not resimulated.
  dat <- mk_ec_dat(n_active = 100, run = list(num = 50), resim = FALSE)
  expect_equal(coef_form(edges_correct(dat, 2)), coef_form(dat))
})

test_that("the default correction still preserves mean degree over a run", {
  # The property the correction exists to produce, with no eligibility
  # attribute in play: a population that grows by half moves the coefficient by
  # exactly log(1 / 1.5).
  dat <- mk_ec_dat(n_active = 150, run = list(num = 100))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]],
               log(100) - log(150))
})

test_that("edges.correct.attr counts only the tie-eligible nodes", {
  # 100 active, of whom 60 can form ties; 70 could on the previous step.
  dat <- mk_ec_dat(n_active = 100, n_eligible = 60,
                   attr.name = "active.sex",
                   run = list(num = 100, num.elig = 70))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]],
               log(70) - log(60))

  # The excluded band must not enter the count. Here the population grows from
  # 100 to 130 while the eligible count holds at 100, which is the case the
  # control exists for: the coefficient must not move.
  dat <- mk_ec_dat(n_active = 130, n_eligible = 100,
                   attr.name = "active.sex",
                   run = list(num = 100, num.elig = 100))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]], 0)

  # Without the control, the same dat dilutes by the excluded share.
  dat_default <- mk_ec_dat(n_active = 130, n_eligible = 100,
                           run = list(num = 100, num.elig = 100))
  expect_equal(coef_form(edges_correct(dat_default, 2))[[1]] -
                 coef_form(dat_default)[[1]],
               log(100) - log(130))
})

test_that("NA is treated as ineligible", {
  active <- rep(1, 10)
  eligible <- c(rep(1, 6), rep(NA, 2), rep(0, 2))
  dat <- mk_ec_dat(n_active = 10, attr.name = "active.sex",
                   eligible = eligible, run = list(num = 10, num.elig = 6))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]], 0)
})

test_that("edges.correct.attr works with two groups", {
  dat <- mk_ec_dat(n_active = 200, n_group2 = 100, groups = 2,
                   n_eligible = 150, attr.name = "active.sex",
                   run = list(num = 90, num.g2 = 90,
                              num.elig = 80, num.elig.g2 = 40))
  # 150 eligible over 200 nodes: nodes 1:100 are group 1, 101:200 group 2, and
  # the first 150 are eligible, so 100 eligible in group 1 and 50 in group 2.
  expected <- log(2 * 80 * 40 / (80 + 40)) - log(2 * 100 * 50 / (100 + 50))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]],
               expected)
})

test_that("a run with no stored eligible count treats the step as a no-op", {
  # A simulation restarted from output saved before the control existed.
  dat <- mk_ec_dat(n_active = 130, n_eligible = 100,
                   attr.name = "active.sex", run = list(num = 100))
  expect_equal(coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]], 0)
})

test_that("an empty count warns and leaves the coefficients alone", {
  # Previously this added -Inf or NaN to the edges coefficient, after which
  # every proposed tie is rejected for the rest of the run: a slow collapse
  # rather than an error.
  dat <- mk_ec_dat(n_active = 100, n_eligible = 0, attr.name = "active.sex",
                   run = list(num = 100, num.elig = 100))
  expect_warning(out <- edges_correct(dat, 2), "non-finite")
  expect_equal(coef_form(out), coef_form(dat))

  dat <- mk_ec_dat(n_active = 0, run = list(num = 100))
  expect_warning(out <- edges_correct(dat, 2), "non-finite")
  expect_equal(coef_form(out), coef_form(dat))
})

test_that("the adjustment is finite at a realistic population size", {
  # Population counts are integers, and a product of two of them overflows
  # 32-bit arithmetic above roughly 46,000 per factor. R returns NA rather
  # than promoting, and an NA reaching the edges coefficient collapses the
  # network silently, so the two-group harmonic mean is checked at scale.
  dat <- mk_ec_dat(n_active = 200000, n_group2 = 100000, groups = 2,
                   run = list(num = 99000, num.g2 = 99000))
  adj <- coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]]
  expect_true(is.finite(adj))

  dat <- mk_ec_dat(n_active = 122875, n_eligible = 99361,
                   attr.name = "active.sex",
                   run = list(num = 122000, num.elig = 99915))
  adj <- coef_form(edges_correct(dat, 2))[[1]] - coef_form(dat)[[1]]
  expect_true(is.finite(adj))
  expect_equal(adj, log(99915) - log(99361))
})

test_that("a mis-specified edges.correct.attr errors clearly", {
  dat <- mk_ec_dat(n_active = 100, n_eligible = 60, attr.name = "active.sex",
                   run = list(num = 100, num.elig = 70))
  dat <- set_attr(dat, "active.sex", rep(1, 5), override.length.check = TRUE)
  expect_error(edges_correct(dat, 2), "has length 5")

  dat <- mk_ec_dat(n_active = 100, n_eligible = 60, attr.name = "missing.attr",
                   run = list(num = 100, num.elig = 70))
  expect_error(edges_correct(dat, 2), "missing.attr")
})

test_that("update_sim_num maintains the eligible counts", {
  # Without the control it must not create them at all, so output saved by an
  # unrelated model is unchanged.
  dat <- mk_ec_dat(n_active = 130, n_eligible = 100)
  out <- update_sim_num(dat)
  expect_equal(out$run$num, 130)
  expect_null(out$run$num.elig)

  dat <- mk_ec_dat(n_active = 130, n_eligible = 100, attr.name = "active.sex")
  out <- update_sim_num(dat)
  expect_equal(out$run$num, 130)
  expect_equal(out$run$num.elig, 100)

  dat <- mk_ec_dat(n_active = 200, n_group2 = 100, groups = 2,
                   n_eligible = 150, attr.name = "active.sex")
  out <- update_sim_num(dat)
  expect_equal(out$run$num, 100)
  expect_equal(out$run$num.g2, 100)
  expect_equal(out$run$num.elig, 100)
  expect_equal(out$run$num.elig.g2, 50)
})

test_that("control.net carries edges.correct.attr and defaults to NULL", {
  ctrl <- control.net(type = "SI", nsteps = 5)
  expect_null(ctrl$edges.correct.attr)

  ctrl <- control.net(type = "SI", nsteps = 5,
                      edges.correct.attr = "active.sex")
  expect_equal(ctrl$edges.correct.attr, "active.sex")
})

test_that("the correction holds mean degree across a filling excluded band", {
  # An end-to-end property check: over 300 steps the total population grows by
  # a quarter while the eligible count holds, and the cumulative drift must be
  # zero with the control set and log(1 / 1.25) without it.
  drift_on <- 0
  drift_off <- 0
  total <- 100000
  prev_elig <- 100000
  for (i in seq_len(300)) {
    prev_total <- total
    total <- total + 100
    on <- mk_ec_dat(n_active = total, n_eligible = 100000,
                    attr.name = "active.sex",
                    run = list(num = prev_total, num.elig = prev_elig))
    off <- mk_ec_dat(n_active = total, n_eligible = 100000,
                     run = list(num = prev_total, num.elig = prev_elig))
    drift_on <- drift_on + coef_form(edges_correct(on, 2))[[1]] -
      coef_form(on)[[1]]
    drift_off <- drift_off + coef_form(edges_correct(off, 2))[[1]] -
      coef_form(off)[[1]]
    prev_elig <- update_sim_num(on)$run$num.elig
  }
  expect_equal(drift_on, 0)
  expect_equal(drift_off, log(100000) - log(total))
  expect_equal(total, 130000)
})
