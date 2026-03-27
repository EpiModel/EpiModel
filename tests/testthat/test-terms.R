context("EpiModel Terms")

nw <- network_initialize(n = 50)
age <- runif(50)
sex <- rep(c(0,1), length.out = 50)
nw %v% "age" <- age
nw %v% "sex" <- sex

test_that("netest works for EpiModel terms", {
  skip_on_cran()
  est <- netest(nw, formation = ~edges + absdiffby("age", "sex", 1) + nodematch("sex"), target.stats = c(25, 25, 0),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")

  est2 <- netest(nw, formation = ~edges + absdiffnodemix("age", "sex"), target.stats = c(40, 5, 10, 5),
                 coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                 verbose = FALSE)
  expect_is(est2, "netest")
})

test_that("EpiModel terms produce correct summary statistics", {
  skip_on_cran()
  nw1 <- san(nw ~ edges + offset(nodematch("sex")), target.stats = c(30), offset.coef = c(-Inf))
    
  el1 <- as.edgelist(nw1)
  
  expect_equal(summary(nw1 ~ absdiffby("age", "sex", 2.3)),
               sum(abs(age[el1[,1]] - 2.3*sex[el1[,1]] - age[el1[,2]] + 2.3*sex[el1[,2]])),
               check.attributes = FALSE)
  
  nw2 <- san(nw ~ edges, target.stats = c(40))

  el2 <- as.edgelist(nw2)

  expect_equal(summary(nw2 ~ absdiffnodemix("age", "sex")),
               c(sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 0 & sex[el2[,2]] == 0)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(((sex[el2[,1]] == 0) + (sex[el2[,2]] == 0)) == 1)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 1 & sex[el2[,2]] == 1))),
               check.attributes = FALSE)
})

test_that("EpiModel terms produce correct change statistics", {
  skip_on_cran()
  nw1 <- simulate(nw ~ edges + offset(nodematch("sex")),
                  coef = c(-3, -Inf),
                  monitor = ~absdiffby("age", "sex", 2.3))
  
  stats1 <- attr(nw1, "stats")
  
  el1 <- as.edgelist(nw1)
  
  expect_equal(stats1[-c(1,2)],
               sum(abs(age[el1[,1]] - 2.3*sex[el1[,1]] - age[el1[,2]] + 2.3*sex[el1[,2]])),
               check.attributes = FALSE)
  
  nw2 <- simulate(nw ~ edges,
                  coef = c(-3),
                  monitor = ~absdiffnodemix("age", "sex"))
  
  stats2 <- attr(nw2, "stats")
  
  el2 <- as.edgelist(nw2)

  expect_equal(stats2[-c(1)],
               c(sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 0 & sex[el2[,2]] == 0)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(((sex[el2[,1]] == 0) + (sex[el2[,2]] == 0)) == 1)),
                 sum(abs(age[el2[,1]] - age[el2[,2]])*(sex[el2[,1]] == 1 & sex[el2[,2]] == 1))),
               check.attributes = FALSE)
})

context("fuzzynodematch Term")

test_that("fuzzynodematch works as intended", {
  skip_on_cran()
  n <- 1000L
  bip <- 400L
  nv <- 10L
  nvm <- 2000L
  prob <- 0.1

  for (directed in list(FALSE, TRUE)) {
    for (bipartite in list(FALSE, bip)) {
      for (split in list("|", ".")) {
        for (binary in list(FALSE, TRUE)) {
          if (directed == TRUE && bipartite == bip) {
            next
          }

          vids <- as.character(sample(seq_len(nvm), nv, FALSE))
          vids <- paste0(rep(c("a","b",""), length.out = length(vids)), vids)

          vcs <- matrix(as.logical(rbinom(n*nv, 1L, prob)), nrow = n)

          duplicate <- sample(c(FALSE, TRUE), n, TRUE)

          attr <- character(n)
          for (i in seq_along(attr)) {
            charvec <- vids[vcs[i,]]
            if (duplicate[i] == TRUE) {
              charvec <- c(charvec, sample(charvec, length(charvec), TRUE))
            }
            charvec <- sample(charvec)
            attr[i] <- paste(charvec, collapse = split)
          }

          nw <- network.initialize(n, directed = directed, bipartite = bipartite)
          nw %v% "attr" <- attr

          el <- as.edgelist(san(nw ~ edges, target = c(n)))

          toggles <- rbind(el, el)
          toggles <- toggles[sample(seq_len(NROW(toggles))), , drop = FALSE]
          toggles <- cbind(seq_len(NROW(toggles)), toggles)

          changes <- cbind(toggles, 1L)
          for (i in seq_len(NROW(changes))) {
            if(min(which(changes[,2L] == changes[i,2L] & changes[,3L] == changes[i,3L])) < i) {
              changes[i,4L] <- 0L
            }
          }

          gf_stats <- tergm.godfather(nw ~ fuzzynodematch(~attr, split = split, binary = binary), toggles = toggles, stats.start = TRUE)

          manual_stats <- integer(NROW(toggles))
          for (i in seq_along(manual_stats)) {
            manual_stats[i] <- sum(vcs[toggles[i,2L],]*vcs[toggles[i,3L],])
            if (binary == TRUE) {
              manual_stats[i] <- as.integer(manual_stats[i] > 0)
            }
            if (changes[i,4L] == 0L) {
              manual_stats[i] <- -manual_stats[i]
            }
          }
          manual_stats <- cumsum(c(0L, manual_stats))

          expect_identical(manual_stats, as.integer(gf_stats))
        }
      }
    }
  }
})
