
#' @title Get Epidemic Output from icm Model
#'
#' @description This function provides all active model state sizes from
#'              the network at the specified time step, output to a list of
#'              vectors.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
prevalence.icm <- function(dat, at) {

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "s" &
                           dat$attr$group == 1)
    dat$epi$i.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "i" &
                           dat$attr$group == 1)
    dat$epi$num <- dat$epi$s.num + dat$epi$i.num
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(dat$attr$active == 1 &
                             dat$attr$status == "r" &
                             dat$attr$group == 1)
      dat$epi$num <- dat$epi$s.num +
        dat$epi$i.num +
        dat$epi$r.num
    }
  } else {
    dat$epi$s.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "s" &
                               dat$attr$group == 1)
    dat$epi$i.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "i" &
                               dat$attr$group == 1)
    dat$epi$num[at] <- dat$epi$s.num[at] + dat$epi$i.num[at]
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$status == "r" &
                                 dat$attr$group == 1)
      dat$epi$num[at] <- dat$epi$s.num[at] +
        dat$epi$i.num[at] +
        dat$epi$r.num[at]
    }
  }

  return(dat)
}


#' @title Get Epidemic Output from icm Model
#'
#' @description This function provides all active model state sizes from
#'              the network at the specified time step, output to a list of
#'              vectors.
#'
#' @param dat Master data list object.
#' @param at Current time step.
#'
#' @export
#' @keywords internal
#'
prevalence.icm.bip <- function(dat, at) {

  if (at == 1) {
    dat$epi <- list()
    dat$epi$s.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "s" &
                           dat$attr$group == 1)
    dat$epi$i.num <- sum(dat$attr$active == 1 &
                           dat$attr$status == "i" &
                           dat$attr$group == 1)
    dat$epi$num <- dat$epi$s.num + dat$epi$i.num
    dat$epi$s.num.g2 <- sum(dat$attr$active == 1 &
                              dat$attr$status == "s" &
                              dat$attr$group == 2)
    dat$epi$i.num.g2 <- sum(dat$attr$active == 1 &
                              dat$attr$status == "i" &
                              dat$attr$group == 2)
    dat$epi$num.g2 <- dat$epi$s.num.g2 + dat$epi$i.num.g2
    if (dat$control$type == "SIR") {
      dat$epi$r.num <- sum(dat$attr$active == 1 &
                             dat$attr$status == "r" &
                             dat$attr$group == 1)
      dat$epi$num <- dat$epi$s.num +
        dat$epi$i.num +
        dat$epi$r.num
      dat$epi$r.num.g2 <- sum(dat$attr$active == 1 &
                                dat$attr$status == "r" &
                                dat$attr$group == 2)
      dat$epi$num.g2 <- dat$epi$s.num.g2 +
        dat$epi$i.num.g2 +
        dat$epi$r.num.g2
    }
  } else {
    dat$epi$s.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "s" &
                               dat$attr$group == 1)
    dat$epi$i.num[at] <- sum(dat$attr$active == 1 &
                               dat$attr$status == "i" &
                               dat$attr$group == 1)
    dat$epi$num[at] <- dat$epi$s.num[at] + dat$epi$i.num[at]
    dat$epi$s.num.g2[at] <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "s" &
                                  dat$attr$group == 2)
    dat$epi$i.num.g2[at] <- sum(dat$attr$active == 1 &
                                  dat$attr$status == "i" &
                                  dat$attr$group == 2)
    dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] + dat$epi$i.num.g2[at]
    if (dat$control$type == "SIR") {
      dat$epi$r.num[at] <- sum(dat$attr$active == 1 &
                                 dat$attr$status == "r" &
                                 dat$attr$group == 1)
      dat$epi$num[at] <- dat$epi$s.num[at] +
        dat$epi$i.num[at] +
        dat$epi$r.num[at]
      dat$epi$r.num.g2[at] <- sum(dat$attr$active == 1 &
                                    dat$attr$status == "r" &
                                    dat$attr$group == 2)
      dat$epi$num.g2[at] <- dat$epi$s.num.g2[at] +
        dat$epi$i.num.g2[at] +
        dat$epi$r.num.g2[at]
    }
  }

  return(dat)
}

