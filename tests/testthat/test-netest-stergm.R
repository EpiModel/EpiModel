context("Network STERGM estimation")

test_that("Basic STERGM fit", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                edapprox = FALSE, verbose = FALSE)
  expect_is(est, "netest")
  expect_true(!est$edapprox)
})
