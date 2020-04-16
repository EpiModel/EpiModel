context("`dat` object getters and setters")

test_that("`dat` getters and setter", {

  dat <- list(
    attr = list(
      active = rbinom(100, 1, 0.9)
    ),
    epi = list(),
    param = list(),
    init = list(),
    control = list(
      nsteps = 150
    )
  )

  ## Attr tests
  dat <- add_attr(dat, "age")
  expect_equal(dat$attr$age, rep(NA, length(dat$attr$active)))
  expect_error(dat <- add_attr(dat, "age"))

  expect_error(set_attr(dat, "age", 4))

  new_ages <- runif(length(dat$attr$active))
  dat <- set_attr(dat, "age", new_ages)
  expect_equal(dat$attr$age, new_ages)

  expect_equal(get_attr(dat, "age"), new_ages)
  expect_equal(get_attr(dat, "age", c(1, 5)), new_ages[c(1, 5)])
  expect_equal(get_attr(dat, "age", new_ages > 0.5), new_ages[new_ages > 0.5])

  expect_error(get_attr(dat, "age_absent"))
  expect_null(get_attr(dat, "age_absent", override.null.error = TRUE))
  expect_error(get_attr(dat, "age", c(1, 1000)))
  expect_error(get_attr(dat, "age", c(TRUE, FALSE)))

  expect_silent(dat <- set_attr(dat, "status",
                          rbinom(length(dat$attr$active), 1, 0.4)))

  expect_equal(get_attr_list(dat), dat$attr)
  expect_equal(get_attr_list(dat, c("age", "status")),
               dat$attr[c("age", "status")])

  expect_error(get_attr_list(dat, "sex"))

  expect_error(dat <- append_attr(dat, "status",
                          rbinom(length(dat$attr$active), 10)))

  expect_error(dat <- append_attr(dat, "status", 1, -1))

  dat <- append_attr(dat, "active", 1, 10)
  expect_length(get_attr(dat, "active"), 110)

  dat <- append_attr(dat, "status", sample(0:1, 10, TRUE), 10)
  expect_length(get_attr(dat, "status"), 110)

  ## Epi tests
  dat <- add_epi(dat, "i")
  expect_equal(dat$epi$i, rep(NA, dat$control$nsteps))

  expect_error(set_epi_at(dat, "i", c(1, 4), 4))

  new_i <- runif(dat$control$nsteps)

  dat <- set_epi(dat, "i", new_i)
  expect_equal(dat$epi$i, new_i)

  expect_equal(get_epi(dat, "i", c(1, 100)), dat$epi$i[c(1, 100)])
  expect_equal(get_epi(dat, "i", dat$epi$i > 0.5), dat$epi$i[dat$epi$i > 0.5])

  expect_error(get_epi(dat, "age_absent"))
  expect_null(get_epi(dat, "age_absent", override.null.error = TRUE))
  expect_error(get_epi(dat, "i", c(1, 300)))
  expect_error(get_epi(dat, "i", c(TRUE, FALSE)))

  dat$control$nsteps <- 200
  dat <- set_epi_at(dat, "i", 160, 8)
  expect_length(dat$epi$i, 200)

  expect_silent(dat <- set_epi(dat, "s", rbinom(dat$control$nsteps, 1, 0.4)))

  expect_equal(get_epi_list(dat), dat$epi)
  expect_equal(get_epi_list(dat, c("i", "s")), dat$epi[c("i", "s")])

  expect_error(get_epi_list(dat, "r"))

  # param, init, control tests
  dat <- add_param(dat, "x")

  dat <- add_init(dat, "x")
  dat <- add_control(dat, "x")

  expect_equal(dat$param$x, NA)
  expect_equal(dat$init$x, NA)
  expect_equal(dat$control$x, NA)

  expect_silent(dat <- set_param(dat, "y", 4))
  expect_silent(dat <- set_init(dat, "y", 4))
  expect_silent(dat <- set_control(dat, "y", 4))

  dat <- set_param(dat, "y", 5)
  dat <- set_init(dat, "y", 5)
  dat <- set_control(dat, "y", 5)

  expect_equal(dat$param$y, 5)
  expect_equal(dat$init$y, 5)
  expect_equal(dat$control$y, 5)

  expect_equal(get_param(dat, "y"), 5)
  expect_equal(get_init(dat, "y"), 5)
  expect_equal(get_control(dat, "y"), 5)

  expect_error(get_param(dat, "z"))
  expect_error(get_init(dat, "z"))
  expect_error(get_control(dat, "z"))

  expect_null(get_param(dat, "z", override.null.error = TRUE))
  expect_null(get_init(dat, "z", override.null.error = TRUE))
  expect_null(get_control(dat, "z", override.null.error = TRUE))

  expect_equal(get_param_list(dat), dat$param)
  expect_equal(get_init_list(dat), dat$init)
  expect_equal(get_control_list(dat), dat$control)

  expect_equal(get_param_list(dat, "x"), dat$param["x"])
  expect_equal(get_init_list(dat, "x"), dat$init["x"])
  expect_equal(get_control_list(dat, "x"), dat$control["x"])

  expect_error(get_param_list(dat, "z"))
  expect_error(get_init_list(dat, "z"))
  expect_error(get_control_list(dat, "z"))

})
