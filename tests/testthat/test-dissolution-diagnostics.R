context("dissolution diagnostics")

test_that("simulation diagnostics work as expected", {

  ## direct implementation of dissolution statistics based on tedgelist
  simulate_diss_stats <- function(est, nsteps, dyad_indexer) {
    nws1 <- simulate(est$formula,
                     coef = est$coef.form.crude,
                     basis = est$newnetwork,
                     dynamic = FALSE)
  
    nws2 <- simulate(~Form(est$formation) + Persist(est$coef.diss$dissolution),
                     coef = c(est$coef.form, est$coef.diss$coef.crude),
                     basis = nws1,
                     time.slices = nsteps,
                     dynamic = TRUE)
  
    if (est$coef.diss$diss.model.type == "nodefactor") {
      dissolution <- ~offset(edges)
      durs <- mean(est$coef.diss$duration)
    } else {
      dissolution <- est$coef.diss$dissolution
      durs <- est$coef.diss$duration
    }
    
    sim.df <- as.data.frame(nws2)
  
    edgecounts <- matrix(0, nrow = nsteps, ncol = length(durs))
    edgeages <- matrix(0, nrow = nsteps, ncol = length(durs))
    edgediss <- matrix(0, nrow = nsteps, ncol = length(durs))
      
    dyad_types <- dyad_indexer(sim.df$tail, sim.df$head, nws1)
    for (i in seq_len(NROW(sim.df))) {
      onset <- max(sim.df$onset[i], 1)
      terminus <- sim.df$terminus[i]
      
      if (terminus > 1) {
        edgecounts[onset - 1 + seq_len(terminus - onset), dyad_types[i]] <- 
          edgecounts[onset - 1 + seq_len(terminus - onset), dyad_types[i]] + 1
        edgeages[onset - 1 + seq_len(terminus - onset), dyad_types[i]] <- 
          edgeages[onset - 1 + seq_len(terminus - onset), dyad_types[i]] + 
          seq_len(terminus - onset) + (sim.df$onset[i] == 0)
      }
      
      if (terminus <= nsteps) {
        edgediss[terminus, dyad_types[i]] <- edgediss[terminus, dyad_types[i]] + 1
      }
    }
  
    init.edgecounts <- summary(dissolution, basis = nws2, at = 0)
    if (length(durs) > 1) {
      init.edgecounts[1] <- init.edgecounts[1] - sum(init.edgecounts[-1])  
    }
    edgecounts <- rbind(init.edgecounts, edgecounts)
    
    edgeagesimputed <- edgeages
    
    wti <- which(sim.df$onset == 0 & sim.df$terminus > 1)
    for (index in wti) {
      edgeagesimputed[seq_len(sim.df$terminus[index] - 1), dyad_types[index]] <- 
        edgeagesimputed[seq_len(sim.df$terminus[index] - 1), dyad_types[index]] + 
        rgeom(1, 1/durs[dyad_types[index]])
    }
    
    pages <- array(edgeages/edgecounts[-1,,drop=FALSE],
                   dim = c(nsteps,length(durs),1))
    pages_imptd <- array(edgeagesimputed/edgecounts[-1,,drop=FALSE],
                         dim = c(nsteps,length(durs),1))
    prop.diss <- array(edgediss/edgecounts[-NROW(edgecounts),,drop=FALSE],
                       dim = c(nsteps,length(durs),1))
    
    pages[is.nan(pages)] <- 0
    pages_imptd[is.nan(pages_imptd)] <- 0
    prop.diss[is.nan(prop.diss)] <- 0
    
    list(pages = pages,
         pages_imptd = pages_imptd,
         prop.diss = prop.diss)
  }

  net_size <- 100
  nsteps <- 10
  nsims <- 2
  seed <- 0
  
  coefs.diss <- list(dissolution_coefs(~offset(edges), 10, 0),
                     dissolution_coefs(~offset(edges)+offset(nodemix("attr", levels = 2:3, levels2 = c(2,1))), c(8,5,7), 0),
                     dissolution_coefs(~offset(edges)+offset(nodematch("attr", diff = TRUE, levels = c(3,1,2))), c(31, 22, 3, 4), 0),
                     suppressWarnings(dissolution_coefs(~offset(edges)+offset(nodefactor("attr", levels = c(3,1))), c(31, 3, 4), 0)))
  
  dyad_indexers <- list(function(tails, heads, nw) rep(1, length(tails)),
                        function(tails, heads, nw) {
                          attr <- nw %v% "attr"
                          tailattr <- attr[tails]
                          headattr <- attr[heads]
                          indices <- rep(1, length(tails))
                          indices[(tailattr == 2 & headattr == 3) | (tailattr == 3 & headattr == 2)] <- 2
                          indices[tailattr == 2 & headattr == 2] <- 3
                          indices
                        },
                        function(tails, heads, nw) {
                          attr <- nw %v% "attr"
                          tailattr <- attr[tails]
                          headattr <- attr[heads]
                          indices <- rep(1, length(tails))
                          indices[tailattr == 3 & headattr == 3] <- 2
                          indices[tailattr == 1 & headattr == 1] <- 3
                          indices[tailattr == 2 & headattr == 2] <- 4
                          indices
                        },
                        function(tails, heads, nw) rep(1, length(tails)))
  
  formations <- list(~edges,
                     ~edges + nodemix("attr", levels = 2:3, levels2 = c(2,1)),
                     ~edges + nodematch("attr", diff = TRUE, levels = c(3,1,2)),
                     ~edges + nodefactor("attr", levels = c(3,1)))
  
  targets <- lapply(list(c(100),
                         c(100, 2*100/9, 2*100/9),
                         c(100, 100/9, 100/9, 100/9),
                         c(100, 200/3, 200/3)),
                    as.integer)
  
  for (nested.edapprox in list(FALSE, TRUE)) { 
    for (init.edges in list(FALSE, TRUE)) {
      for (index in seq_along(coefs.diss)) {
        coef.diss <- coefs.diss[[index]]
        dyad_indexer <- dyad_indexers[[index]]
        if (nested.edapprox == TRUE) {
          formation <- formations[[index]]
          target.stats <- targets[[index]]
        } else {
          formation <- ~edges
          target.stats <- c(100)
        }

        nw <- network.initialize(net_size, directed = FALSE)
        nw %v% "attr" <- rep(1:3, length.out = net_size)
                          
        est <- suppressWarnings(netest(nw, 
                                       formation = formation,
                                       coef.diss = coef.diss,
                                       target.stats = target.stats,
                                       nested.edapprox = nested.edapprox))
        
        if (init.edges == FALSE) {
          est$newnetwork[,] <- FALSE
        }
        
        set.seed(seed)
        dx <- netdx(est, nsims = nsims, nsteps = nsteps)
        
        set.seed(seed)
        ds <- list()
        for (i in seq_len(nsims)) {
          ds[[i]] <- simulate_diss_stats(est, nsteps, dyad_indexer)
        }
        
        if (coef.diss$diss.model.type == "nodefactor") {
          durs <- mean(coef.diss$duration)
        } else {
          durs <- coef.diss$duration
        }
        
        pages <- array(unlist(lapply(ds, `[[`, "pages")),
                       dim = c(nsteps,length(durs),nsims))
        pages_imptd <- array(unlist(lapply(ds, `[[`, "pages_imptd")),
                             dim = c(nsteps,length(durs),nsims))
        prop.diss <- array(unlist(lapply(ds, `[[`, "prop.diss")),
                           dim = c(nsteps,length(durs),nsims))
        
        expect_equal(pages, dx$pages)
        expect_equal(pages_imptd, dx$pages_imptd)
        expect_equal(prop.diss, dx$prop.diss)
      }
    }
  }
})

test_that("tedgelist_to_toggles functions as expected", {
  logit <- function(p) log(p/(1-p))
  density <- 1/50
  D <- 10

  for (init.edges in list(FALSE, TRUE)) {
    nw <- network.initialize(100, directed = FALSE)
    if (init.edges == TRUE) {
      nw <- san(nw ~ edges, target.stats = c(100))
    }

    set.seed(0)
    nwd <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(logit(density) - log(D), log(D - 1)),
                    time.slices = 10,
                    output = "networkDynamic",
                    dynamic = TRUE)
    toggles <- tedgelist_to_toggles(as.data.frame(nwd))
  
    set.seed(0)
    changes <- simulate(nw ~ Form(~edges) + Persist(~edges),
                        coef = c(logit(density) - log(D), log(D - 1)),
                        time.slices = 10,
                        output = "changes",
                        dynamic = TRUE)
    
    if (init.edges == TRUE) {
      changes <- rbind(cbind(0L, as.edgelist(nw), 1L),
                       changes)
    }
    
    toggles2 <- changes[,-4L,drop=FALSE]
    
    toggles <- toggles[order(toggles[,1], toggles[,2], toggles[,3]),,drop=FALSE]
    toggles2 <- toggles2[order(toggles2[,1], toggles2[,2], toggles2[,3]),,drop=FALSE]
    expect_identical(unname(toggles), unname(toggles2))
  }
})
