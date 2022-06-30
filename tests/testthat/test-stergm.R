
context("Full STERGM Workflow")

test_that("Full STERGM", {
  skip_on_cran()
  skip_on_os("windows")
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                edapprox = FALSE, verbose = FALSE)
  
  for (trim in c(FALSE, TRUE)) {
    if (trim == TRUE) {
      est2 <- trim_netest(est)
    } else {
      est2 <- est
    }

    # one core test
    dx <- netdx(est2, nsims = 1, nsteps = 10, verbose = FALSE)
    expect_is(dx, "netdx")
    expect_true(!dx$edapprox)
    expect_true(colnames(dx$stats[[1]]) == "edges")
    
    # parallel test
    dx <- netdx(est2, nsims = 2, nsteps = 10, ncores = 2, verbose = FALSE)
    expect_is(dx, "netdx")
    expect_true(dx$nsims == 2)
    expect_is(dx$nw, "network")
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
    mod <- netsim(est2, param, init, control)
    expect_is(mod, "netsim")
  }
})
