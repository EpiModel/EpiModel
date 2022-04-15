context("fuzzynodematch")

test_that("fuzzynodematch works as intended", {
  n <- 10000L
  bip <- 4000L
  nv <- 10L
  nvm <- 1000L
  prob <- 0.05
  
  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip)) {
      if(directed && bipartite) {
        next
      }
        
      vids <- as.character(sample(seq_len(nvm), nv, FALSE))
      for(i in seq_along(vids)) {
        vids[i] <- paste0(c("v", rep("0", 4L - nchar(vids[i])), vids[i]), collapse = "")
      }
      
      vcs <- matrix(as.logical(rbinom(n*nv, 1L, prob)), nrow = n)
      
      attr <- character(n)
      for(i in seq_along(attr)) {
        attr[i] <- paste(vids[vcs[i,]], collapse = "|")
        if(nchar(attr[i]) == 0L) {
          attr[i] <- paste0(c("ego", rep("0", 5L - nchar(i)), i), collapse = "")
        }
      }
      
      nw <- network.initialize(n, directed = directed, bipartite = bipartite)
      nw %v% "attr" <- attr
      
      el <- as.edgelist(san(nw ~ edges, target = c(n)))
      
      toggles <- rbind(el, el)
      toggles <- toggles[sample(seq_len(NROW(toggles))), , drop = FALSE]
      toggles <- cbind(seq_len(NROW(toggles)), toggles)  
      
      changes <- cbind(toggles, 1L)
      for(i in seq_len(NROW(changes))) {
        if(min(which(changes[,2L] == changes[i,2L] & changes[,3L] == changes[i,3L])) < i) {
          changes[i,4L] <- 0L
        }
      }
      
      for(binary in list(FALSE, TRUE)) {
        gf_stats <- tergm.godfather(nw ~ fuzzynodematch(~attr, binary), toggles = toggles, stats.start = TRUE)
      
        manual_stats <- integer(NROW(toggles))
        for(i in seq_along(manual_stats)) {
          manual_stats[i] <- sum(vcs[toggles[i,2L],]*vcs[toggles[i,3L],])
          if(binary) {
            manual_stats[i] <- as.integer(manual_stats[i] > 0)
          }
          if(changes[i,4L] == 0L) {
            manual_stats[i] <- -manual_stats[i]
          }
        }
        manual_stats <- cumsum(c(0L, manual_stats))
        
        expect_identical(manual_stats, as.integer(gf_stats))
      }
    }
  }
})
