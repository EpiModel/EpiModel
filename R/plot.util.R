
## Non-exported plot utilities

draw_qnts <- function(x, y, qnts, qnts.pal, loc = "epi") {

  lcomp <- length(y)
  for (j in seq_len(lcomp)) {
    quants <- c((1-qnts)/2, 1-((1-qnts)/2))
    qnt.prev <- apply(x[[loc]][[y[j]]], 1,
                      function(x) quantile(x, c(quants[1], quants[2])))
    xx <- c(1:(ncol(qnt.prev)), (ncol(qnt.prev)):1)
    yy <- c(qnt.prev[1,], rev(qnt.prev[2,]))
    polygon(xx, yy, col = qnts.pal[j], border = NA)
  }

}
