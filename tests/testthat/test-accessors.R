context("Accessors: Getters and Setters")

test_that("`dat` getters and setter", {

  n_nodes <- 100
  dat <- create_dat_object(control = list(nsteps = 150))
  dat <- append_core_attr(dat, at = 1,  n.new = 100)

  ## Attr tests
  dat <- add_attr(dat, "age")
  expect_equal(get_attr(dat, "age"), rep(NA, n_nodes))
  expect_error(dat <- add_attr(dat, "age"))

  expect_error(set_attr(dat, "age", 4))

  new_ages <- runif(n_nodes)
  dat <- set_attr(dat, "age", new_ages)
  expect_equal(get_attr(dat, "age"), new_ages)

  dat <- set_attr(dat, "age2", new_ages)
  expect_silent(dat <- set_attr(dat, "age2", rep(new_ages, 2),
                                override.length.check = TRUE))
  expect_length(get_attr(dat, "age2"), 2 * length(new_ages))


  expect_equal(get_attr(dat, "age"), new_ages)
  expect_equal(get_attr(dat, "age", c(1, 5)), new_ages[c(1, 5)])

  expect_error(get_attr(dat, "age_absent"))
  expect_null(get_attr(dat, "age_absent", override.null.error = TRUE))
  expect_error(get_attr(dat, "age", c(1, 1000)))
  expect_error(get_attr(dat, "age", c(TRUE, FALSE)))

  expect_silent(dat <- set_attr(dat, "status", rbinom(n_nodes, 1, 0.4)))

  expect_silent(dat <- set_attr(dat, "age", 2, posit_ids = 1:4))
  expect_equal(get_attr(dat, "age", 1:4), rep(2, 4))

  expect_error(dat <- set_attr(dat, "age", c(1, 2), posit_ids = 1:4))
  expect_error(dat <- set_attr(dat, "age", 1, posit_ids = c(1, 1000)))

  expect_error(dat <- set_attr(dat, "age", 1, posit_ids = TRUE))
  expect_error(dat <- set_attr(dat, "age", 1, posit_ids = "a"))


  expect_error(get_attr_list(dat, "sex"))

  expect_error(dat <- append_attr(dat, "status", rbinom(n_nodes, 10)))

  expect_error(dat <- append_attr(dat, "status", 1, -1))

  dat <- append_attr(dat, "active", 1, 10)
  expect_length(get_attr(dat, "active"), 110)

  dat <- append_attr(dat, "status", sample(0:1, 10, TRUE), 10)
  expect_length(get_attr(dat, "status"), 110)

  ## Epi tests
  dat <- add_epi(dat, "i")
  expect_equal(dat$epi$i, rep(NA_integer_, get_control(dat, "nsteps")))

  expect_error(set_epi(dat, "i", c(1, 4), 4))

  dat <- set_epi(dat, "i", 150, 10)
  expect_equal(get_epi(dat, "i")[150], 10)

  dat <- set_epi(dat, "s", 110, 10)
  expect_equal(get_epi(dat, "s")[110], 10)

  expect_equal(get_epi(dat, "i", c(1, 100)), dat$epi$i[c(1, 100)])

  expect_error(get_epi(dat, "age_absent"))
  expect_null(get_epi(dat, "age_absent", override.null.error = TRUE))
  expect_error(get_epi(dat, "i", c(1, 300)))
  expect_error(get_epi(dat, "i", c(TRUE, FALSE)))

  dat$control$nsteps <- 200
  dat <- set_epi(dat, "i", 160, 8)
  expect_length(dat$epi$i, 200)

  expect_equal(get_epi_list(dat, c("i", "s")), get_epi_list(dat)[c("i", "s")])

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

  expect_equal(get_param_list(dat, "x"), get_param_list(dat)["x"])
  expect_equal(get_init_list(dat, "x"), get_init_list(dat)["x"])
  expect_equal(get_control_list(dat, "x"), get_control_list(dat)["x"])

  expect_error(get_param_list(dat, "z"))
  expect_error(get_init_list(dat, "z"))
  expect_error(get_control_list(dat, "z"))

})

test_that("Net core attributes", {

  dat <- create_dat_object(control = list(nsteps = 150))

  # Append the first nodes (empty list before)
  dat <- append_core_attr(dat, at = 1,  n.new = 100)
  expect_equal(get_attr(dat, "active"), rep(1, 100))
  expect_equal(get_attr(dat, "unique_id"), 1:100)

  # Remove some nodes to check if unique_ids are unique
  dat <- delete_attr(dat, 21:30)
  dat <- append_core_attr(dat, at = 2, n.new = 100)
  expect_equal(get_attr(dat, "active"), rep(1, 190))
  expect_equal(get_attr(dat, "unique_id"), c(1:20, 31:200))
  expect_type(get_attr(dat, "unique_id"), "integer")

  # Test unique_ids posit_ids converters
  expect_equal(
    get_attr(dat, "unique_id"),
    get_unique_ids(dat, seq_along(get_attr(dat, "active")))
  )
  expect_equal(
    seq_along(get_attr(dat, "active")),
    get_posit_ids(dat, get_attr(dat, "unique_id"))
  )
  expect_warning(get_posit_ids(dat, 25:35))
})

context("Accessors: Attribute Copying")

################################################################################

test_that("Copying attributes from network to attribute list", {
  skip_on_cran()

  num1 <- num2 <- 500
  nw <- network_initialize(num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  nw <- set_vertex_attribute(nw, "race", sample(c("B", "W"), num1 + num2,
                                                replace = TRUE))
  nw <- set_vertex_attribute(nw, "region", sample(1:4, num1 + num2,
                                                  replace = TRUE))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                     act.rate = 5)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2, tergmLite = FALSE,
                         raw.output = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)

  for (simno in c(1, 2)) {
    dat <- sim[[simno]]
    # Character attribute
    dat.attr <- prop.table(table(get_attr(dat, "race")))
    nw.attr <- prop.table(table(get_vertex_attribute(dat$run$nw[[1]], "race")))

    expect_equal(dat.attr, nw.attr)

    # Numeric attribute
    dat.attr <- prop.table(table(get_attr(dat, "region")))
    nw.attr <- prop.table(table(get_vertex_attribute(dat$run$nw[[1]], "region")))

    expect_equal(dat.attr, nw.attr)
  }
})

test_that("set_attr and append_attr keep one value per node", {
  n_nodes <- 10
  dat <- create_dat_object(control = list(nsteps = 10))
  dat <- append_core_attr(dat, at = 1, n.new = n_nodes)

  ## Creating an attribute through posit_ids fills the other nodes with NA
  d1 <- set_attr(dat, "new", c(5, 6), posit_ids = c(1, 2))
  expect_length(get_attr(d1, "new"), n_nodes)
  expect_equal(get_attr(d1, "new"), c(5, 6, rep(NA, n_nodes - 2)))
  expect_true(check_attr_lengths(d1))

  ## The length does not depend on which nodes were selected
  d2 <- set_attr(dat, "new", c(5, 6), posit_ids = c(4, 7))
  expect_length(get_attr(d2, "new"), n_nodes)
  expect_equal(which(!is.na(get_attr(d2, "new"))), c(4, 7))

  ## Writing to an attribute that already exists is unchanged
  d3 <- set_attr(dat, "age", seq_len(n_nodes))
  d3 <- set_attr(d3, "age", 99, posit_ids = 3)
  expect_length(get_attr(d3, "age"), n_nodes)
  expect_equal(get_attr(d3, "age")[3], 99)

  ## Appending an attribute that does not exist yet pads the older nodes
  d4 <- append_core_attr(dat, at = 2, n.new = 2)
  d4 <- append_attr(d4, "brandnew", 99, 2)
  expect_length(get_attr(d4, "brandnew"), n_nodes + 2)
  expect_equal(get_attr(d4, "brandnew"), c(rep(NA, n_nodes), 99, 99))
  expect_true(check_attr_lengths(d4))

  ## The padding does not change the type of the appended values
  expect_type(get_attr(append_attr(d4, "int", 3L, 2), "int"), "integer")
  expect_type(get_attr(append_attr(d4, "chr", "a", 2), "chr"), "character")

  ## Core attributes are still created from an empty attribute list
  d5 <- create_dat_object(control = list(nsteps = 10))
  d5 <- append_core_attr(d5, at = 1, n.new = 4)
  expect_equal(get_attr(d5, "unique_id"), 1:4)
  expect_length(get_attr(d5, "active"), 4)
  expect_true(check_attr_lengths(d5))
})

test_that("check_attr_lengths flags attributes of the wrong length", {
  dat <- create_dat_object(control = list(nsteps = 10))
  dat <- append_core_attr(dat, at = 1, n.new = 10)
  expect_true(check_attr_lengths(dat))

  bad <- dat
  bad$run$attr$broken <- 1:3
  expect_error(check_attr_lengths(bad), "`broken`: 3")

  no.active <- dat
  no.active$run$attr$active <- NULL
  expect_error(check_attr_lengths(no.active), "`active` nodal attribute is missing")
})
