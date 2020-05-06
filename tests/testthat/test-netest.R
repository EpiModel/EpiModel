context("Network estimation functions")


# Test formation models ---------------------------------------------------

test_that("netest works for edges only model", {
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")
})


test_that("netest works for edges + nodematch model", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")
})


test_that("netest works with offset.coef terms", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "role", rep(c("I", "V", "R"), c(10, 80, 10)))
  est <- netest(nw, formation = ~edges + offset(nodematch('role', diff = TRUE, keep = 1:2)),
                coef.form = c(-Inf, -Inf), target.stats = c(40),
                coef.diss = dissolution_coefs(~offset(edges), 52 * 2, 0.0009),
                verbose = FALSE)
  expect_is(est, "netest")
})



# Test dissolution models -------------------------------------------------

test_that("netest works for heterogeneous dissolutions", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(50, 20),
                coef.diss = dissolution_coefs(~offset(edges) +
                                                offset(nodematch("race")), c(10, 20)),
                verbose = FALSE
  )
  expect_is(est, "netest")
})

test_that("netest diss_check flags bad models", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))

  formation <- ~edges + nodematch("race")
  dissolution <- ~offset(edges) + offset(nodefactor("race"))
  coef.diss <- dissolution_coefs(dissolution, c(10, 20))
  expect_error(netest(nw, formation, target.stats = c(50, 20),
               coef.diss, verbose = FALSE),
               "Dissolution model is not a subset of formation model.")

  formation <- ~edges + nodematch("race", diff = TRUE)
  dissolution <- ~offset(edges) + offset(nodematch("race"))
  coef.diss <- dissolution_coefs(dissolution, c(10, 20))
  expect_error(netest(nw, formation, target.stats = c(50, 20), coef.diss, verbose = FALSE),
               "Term options for one or more terms in dissolution model")
})


# Other tests -----------------------------------------------------------------------

test_that("Error if incorrect coef.diss parameter", {
  nw <- network_initialize(n = 50)
  coef.diss <- 1
  expect_error(netest(nw, formation = ~edges, target.stats = 25, coef.diss = coef.diss,
                verbose = FALSE))
})

test_that("update_dissolution tests", {
  nw <- network_initialize(n = 1000)

  diss.300 <- dissolution_coefs(~offset(edges), 300, 0.001)
  diss.200 <- dissolution_coefs(~offset(edges), 200, 0.001)

  est300 <- netest(nw = nw,
                  formation = ~edges,
                  target.stats = 500,
                  coef.diss = diss.300)

  est200 <- netest(nw = nw,
                  formation = ~edges,
                  target.stats = 500,
                  coef.diss = diss.200)

  est200.compare <- update_dissolution(est300, diss.200)
  expect_true(round(as.numeric(est200$coef.form, 3)) == round(as.numeric(est200.compare$coef.form, 3)))

  expect_error(update_dissolution(1, diss.200), "old.netest must be an object of class netest")
  expect_error(update_dissolution(est300, 1), "new.coef.diss must be an object of class disscoef")
  est300$edapprox <- FALSE
  expect_error(update_dissolution(est300, diss.200), "Edges dissolution approximation must be used")
})

# STERGM ----------------------------------------------------------------------------

test_that("Basic STERGM fit", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                edapprox = FALSE, verbose = FALSE)
  expect_is(est, "netest")
  expect_true(!est$edapprox)
})
