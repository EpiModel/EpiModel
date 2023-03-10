context("Mutate")

test_that("mutate_epi.netsim", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
  init <- init.net(i.num = 1, i.num.g2 = 0)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         verbose = FALSE)
  mod1 <- netsim(est1, param, init, control)

  mod1 <- mutate_epi(mod1, i.prev = i.num / num,
                     i.prev.g2 = i.num.g2 / num.g2)
  expect_equal(names(mod1$epi), c("sim.num", "sim.num.g2", "s.num", "i.num",
                                  "num", "s.num.g2", "i.num.g2", "num.g2",
                                  "si.flow", "si.flow.g2", "i.prev",
                                  "i.prev.g2"))

  # Add incidence rate per 100 person years (assume time step = 1 week)
  mod1 <- mutate_epi(mod1, ir100 = 5200 * (si.flow + si.flow.g2) /
                       (s.num + s.num.g2))
  df <- as.data.frame(mod1)
  expect_true("ir100" %in% names(df))
})

test_that("mutate_epi.dcm", {

  param <- param.dcm(inf.prob = 0.3, inf.prob.g2 = 0.15, balance = "g1")
  init <- init.dcm(s.num = 1000, s.num.g2 = 1000, i.num = 1, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 100)
  mod1 <- dcm(param, init, control)

  mod1 <- mutate_epi(mod1, i.prev = i.num / num,
                     i.prev.g2 = i.num.g2 / num.g2)
  expect_equal(names(mod1$epi), c("s.num", "s.num.g2", "i.num", "i.num.g2",
                                  "si.flow", "si.flow.g2", "num", "num.g2",
                                  "i.prev", "i.prev.g2"))

  # Add incidence rate per 100 person years (assume time step = 1 week)
  mod1 <- mutate_epi(mod1, ir100 = 5200 * (si.flow + si.flow.g2) /
                       (s.num + s.num.g2))
  df <- as.data.frame(mod1)
  expect_true("ir100" %in% names(df))
})

test_that("mutate_epi.dcm, new mod", {

  syph <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      # 1. track the total population size
      num <- s.num + i.num + r.num

      # 2. define lambda, mu, gamma, and sigma
      ce <- R0 / dur.inf
      lambda <- ce * i.num / num
      mu <- 1 / life.expt
      gamma <- 1 / dur.inf
      sigma <- 1 / dur.imm

      # 3. Write out the four differential equations
      dS <- -lambda * s.num + mu * num - mu * s.num + sigma * r.num
      dI <- lambda * s.num - gamma * i.num - mu * i.num
      dR <- gamma * i.num - mu * r.num - sigma * r.num

      # 4. Outputs
      list(c(dS,
             dI,
             dR))
    })
  }

  param <- param.dcm(R0 = 1.5, life.expt = 365 * 30,
                     dur.inf = 60, dur.imm = 365 * c(8, 10, 12))
  init <- init.dcm(s.num = 1e5, i.num = 1, r.num = 0)
  control <- control.dcm(nsteps = 365, new.mod = syph)
  sim <- dcm(param, init, control)

  sim <- mutate_epi(sim, num = s.num + i.num + r.num)
  sim <- mutate_epi(sim, prev = i.num / num)

  df <- as.data.frame(sim)
  expect_true("num" %in% names(df))
  expect_true("prev" %in% names(df))
  expect_true("num" %in% names(sim$epi))
  expect_true("prev" %in% names(sim$epi))
})

test_that("mutate DCM with constant", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500)
  mod1 <- dcm(param, init, control)
  mod1 <- mutate_epi(mod1, prev = i.num / num,
                     cm = 3)
  expect_true(all(sapply(mod1$epi, class) == "data.frame"))
  expect_true(length(unique(sapply(mod1$epi, nrow))) == 1)

  # by itself
  mod1 <- dcm(param, init, control)
  mod1 <- mutate_epi(mod1, cm = 3)
  expect_true(all(sapply(mod1$epi, class) == "data.frame"))
  expect_true(length(unique(sapply(mod1$epi, nrow))) == 1)
})
