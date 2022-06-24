context("fuzzynodematch")

test_that("fuzzynodematch works as intended", {
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
