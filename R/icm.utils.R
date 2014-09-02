
#' @title Get Epidemic Output from icm Model
#'
#' @description This function provides all active model state sizes from
#'              the network at the specified time step, output to a list of
#'              vectors.
#'
#' @param all master data list object.
#' @param at current time step.
#'
#' @export
#' @keywords internal
#'
get_prev.icm <- function(all, at) {

  if (at == 1) {
    all$out <- list()
    all$out$s.num <- sum(all$attr$active == 1 &
                         all$attr$status == "s" &
                         all$attr$group == 1)
    all$out$i.num <- sum(all$attr$active == 1 &
                         all$attr$status == "i" &
                         all$attr$group == 1)
    all$out$num <- all$out$s.num + all$out$i.num
    if (all$control$type == "SIR") {
      all$out$r.num <- sum(all$attr$active == 1 &
                           all$attr$status == "r" &
                           all$attr$group ==1)
      all$out$num <- all$out$s.num +
                     all$out$i.num +
                     all$out$r.num
    }
    if (all$param$groups == 2) {
      all$out$s.num.g2 <- sum(all$attr$active == 1 &
                              all$attr$status == "s" &
                              all$attr$group == 2)
      all$out$i.num.g2 <- sum(all$attr$active == 1 &
                              all$attr$status == "i" &
                              all$attr$group == 2)
      all$out$num.g2 <- all$out$s.num.g2 + all$out$i.num.g2
      if (all$control$type == "SIR") {
        all$out$r.num.g2 <- sum(all$attr$active == 1 &
                                all$attr$status == "r" &
                                all$attr$group == 2)
        all$out$num.g2 <- all$out$s.num.g2 +
                          all$out$i.num.g2 +
                          all$out$r.num.g2
      }
    }
  } else {
    all$out$s.num[at] <- sum(all$attr$active == 1 &
                             all$attr$status == "s" &
                             all$attr$group == 1)
    all$out$i.num[at] <- sum(all$attr$active == 1 &
                             all$attr$status == "i" &
                             all$attr$group == 1)
    all$out$num[at] <- all$out$s.num[at] + all$out$i.num[at]
    if (all$control$type == "SIR") {
      all$out$r.num[at] <- sum(all$attr$active == 1 &
                               all$attr$status == "r" &
                               all$attr$group ==1)
      all$out$num[at] <- all$out$s.num[at] +
                         all$out$i.num[at] +
                         all$out$r.num[at]
    }
    if (all$param$groups == 2) {
      all$out$s.num.g2[at] <- sum(all$attr$active == 1 &
                                  all$attr$status == "s" &
                                  all$attr$group == 2)
      all$out$i.num.g2[at] <- sum(all$attr$active == 1 &
                                  all$attr$status == "i" &
                                  all$attr$group == 2)
      all$out$num.g2[at] <- all$out$s.num.g2[at] + all$out$i.num.g2[at]
      if (all$control$type == "SIR") {
        all$out$r.num.g2[at] <- sum(all$attr$active == 1 &
                                    all$attr$status == "r" &
                                    all$attr$group == 2)
        all$out$num.g2[at] <- all$out$s.num.g2[at] +
                              all$out$i.num.g2[at] +
                              all$out$r.num.g2[at]
      }
    }
  }

  return(all)
}
