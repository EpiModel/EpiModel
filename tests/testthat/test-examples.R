
test_that("init_tergmLite", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.1.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    str(dat, max.level = 1)
    
    # Element added is el (edgelist representation of network)...
    dat$el
    
    # ... and nw is now a networkLite
    dat$nw[[1]]
  }
})


test_that("networkLite", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.1.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)

    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)

    # Conversion to networkLite class format
    nwl <- networkLite(dat$el[[1]], dat$attr)
    nwl
  }
})


test_that("add_vertices", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.1.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)

    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)

    # Check current network size
    attributes(dat$el[[1]])$n

    # Add 10 vertices
    dat$el[[1]] <- add_vertices(dat$el[[1]], 10)

    # Check new network size
    attributes(dat$el[[1]])$n
  }
})


test_that("delete_vertices", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.1.0") {
    set.seed(12345)
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)

    # Set seed for reproducibility
    set.seed(123456)

    # networkLite representation structure after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)

    # Current edges
    head(dat$el[[1]], 20)

    # Remove nodes 1 and 2
    nodes.to.delete <- 1:2
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes.to.delete)

    # Newly permuted edges
    head(dat$el[[1]], 20)
  }
})
