#'
#' @title Concurrency Microsimulation Model for HIV-1 Transmission Dynamics
#'
#' @description Simulates an HIV-1 epidemic in a population of men and women
#'              with purely heterosexual mixing under varying scenarios
#'              of sexual partnership concurrency.
#'
#' @param s.num.f number of initial susceptible females in the population.
#' @param i.num.f number of initial infected females in the population.
#' @param s.num.m number of initial susceptible males in the population.
#' @param i.num.m number of initial infected females in the population.
#' @param monog.f if \code{TRUE}, enforce a momentary degree constraint of
#'        monogamy for females (females not allowed concurrent partnerships).
#' @param monog.m if \code{TRUE}, enforce a momentary degree constraint of
#'        monogamy for males (males not allowed concurrent partnerships).
#' @param meandeg average momentary mean degree (number of current partnerships)
#'        in the population.
#' @param part.duration average length of partnerships in months.
#' @param nsteps number of time steps to simulate the model over. This must be a
#'        positive integer.
#' @param nsims number of simulations to run.
#' @param verbose if \code{TRUE}, print model progress to the console.
#'
#' @details
#' This function runs a microsimulation model of HIV-1 transmission in a purely
#' heterosexual context to investigate how the presence or absence of relational
#' concurrency in sexual partnerships affects the prevalence of disease
#' prevalence at the population level.
#'
#' The parameters of the model include the initial number of susceptible
#' females and males, the initial number of infected females and males, whether
#' males and females are allowed concurrency, the mean degree of all persons in
#' the population, and the average duration of partnerships (in momths). As the
#' four examples below show, all parameters except concurrency may be held
#' constant to test the effects of concurrency.
#'
#' @section Epidemiology:
#' HIV infection is simulated based on a four-stage disease progression model in
#' which persons transition from acute to latent to pre-AIDS to AIDS stages. These
#' transitions occur at deterministic intervals based on estimates of the average
#' time per stage.
#'
#' The transmission probability to uninfected partners varies by stage of the
#' infected partner: it is highest in the acute stage and lowest in the AIDS
#' stage when no sexual acts occur. See the Hollingsworth reference for further
#' details.
#'
#' @section Model Assumptions:
#' This model makes several simplifying assumptions about partnership formation
#' and dissolution, such as the phenomenon of all dissolving partnerships being
#' immediately replaced by a new partnership in the network. Additionally, the
#' user may specify whether concurrency is allowed, but not the level of
#' concurrency (it is calculated based on a binomial distribution with the
#' probability equal to the mean degree parameter). Therefore, this model serves
#' as an introduction to network modeling featured in the network class of
#' functions in \code{EpiModel}; there the user has much more control over the
#' network parameterization and simulation.
#'
#' @references
#' A web-based implementation of this model is available at
#' \url{http://statnet.org/apps/Conc}. The background and details of this model
#' are explained in a full tutorial on concurrency at
#' \url{http://statnet.org/trac/wiki/ConcurrencyIndex}.
#'
#' Hollingsworth TD, Anderson RM, Fraser C. HIV-1 transmission, by stage of
#' infection. Journal of Infectious Diseases. 2008; 198(5): 687-693.
#'
#' @keywords model internal
#' @export
#'
#' @examples
#' \dontrun{
#' # No concurrency model
#' no.conc <- conc_microsim(
#'    s.num.f = 1000,
#'    i.num.f = 50,
#'    s.num.m = 1000,
#'    i.num.m = 50,
#'    monog.f = TRUE,
#'    monog.m = TRUE,
#'    meandeg = 0.8,
#'    part.duration = 10,
#'    nsteps = 2000,
#'    nsims = 10,
#'    verbose = TRUE)
#'
#' # Male concurrency only model
#' male.conc <- conc_microsim(
#'    s.num.f = 1000,
#'    i.num.f = 50,
#'    s.num.m = 1000,
#'    i.num.m = 50,
#'    monog.f = TRUE,
#'    monog.m = FALSE,
#'    meandeg = 0.8,
#'    part.duration = 10,
#'    nsteps = 2000,
#'    nsims = 10,
#'    verbose = TRUE)
#'
#' # Female concurrency only model
#' feml.conc <- conc_microsim(
#'    s.num.f = 1000,
#'    i.num.f = 50,
#'    s.num.m = 1000,
#'    i.num.m = 50,
#'    monog.f = FALSE,
#'    monog.m = TRUE,
#'    meandeg = 0.8,
#'    part.duration = 10,
#'    nsteps = 2000,
#'    nsims = 10,
#'    verbose = TRUE)
#'
#' # Both sexes concurrency model
#' both.conc <- conc_microsim(
#'    s.num.f = 1000,
#'    i.num.f = 50,
#'    s.num.m = 1000,
#'    i.num.m = 50,
#'    monog.f = FALSE,
#'    monog.m = FALSE,
#'    meandeg = 0.8,
#'    part.duration = 10,
#'    nsteps = 2000,
#'    nsims = 10,
#'    verbose = TRUE)
#'
#' # Plot the results
#' par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2,1,0))
#' plot(no.conc, alpha=0.5, ylim=c(0, 0.5), main="No Concurrency")
#' plot(male.conc, alpha=0.5, ylim=c(0, 0.5), main="Male Concurrency")
#' plot(feml.conc, alpha=0.5, ylim=c(0, 0.5), main="Female Concurrency")
#' plot(both.conc, alpha=0.5, ylim=c(0, 0.5), main="Both Concurrency")
#'}
#'
conc_microsim <- function(s.num.f,
                          i.num.f,
                          s.num.m,
                          i.num.m,
                          monog.f = TRUE,
                          monog.m = TRUE,
                          meandeg,
                          part.duration,
                          nsteps,
                          nsims = 1,
                          verbose = TRUE) {


  # Probability of transmission per month for active relationship
  # From Hollingsworth et al. 2008
  betaByTime <- c(rep(0.2055, 3), rep(0.0088, 100),
                  rep(0.0614, 9), rep(0, 10))


  out <- list()
  for (sim in 1:nsims){
    ptm <- proc.time()

    # Basic calculations ------------------------------------------------------

    # Total pop size
    nFeml <- s.num.f + i.num.f
    nMale <- s.num.m + i.num.m
    n <- nFeml + nMale

    # Expected # of relationships ("edges") in the population at any time
    expectedEdges <- round(meandeg * n / 2)

    # Daily prob of dissolution for an existing relationship
    probDissolution <- 1 / part.duration

    # Length of time from HIV infection until death (in months)
    timeAidsDeath <- length(betaByTime)


    # Female Data Vectors -----------------------------------------------------

    # ID numbers for the females
    idsFeml <- 1:nFeml

    # Create a variable called "hivStatusFeml", with all set to 0 (for now)
    hivStatFeml <- rep(0, nFeml)

    # Create a variable called "infTimeFeml" (infection time), set to missing (for now)
    infTimeFeml <- rep(NA, nFeml)

    # Randomly set infection status and infection time for females
    if (i.num.f > 0) {
      # Sample from vector of female ids with size of initially infected females
      idsInfFeml <- sample(idsFeml, i.num.f)

      # Set the status for those to 1
      hivStatFeml[idsInfFeml] <- 1

      # Sample times backwards from present to length of infection
      infTimeFeml[idsInfFeml] <- sample(0:(-timeAidsDeath + 2),
                                        i.num.f, replace = TRUE)
    }


    # Male Data Vectors -------------------------------------------------------

    # Parallel to female data frame above
    idsMale <- 1:nMale
    hivStatMale <- rep(0, nMale)
    infTimeMale <- rep(NA, nMale)
    if (i.num.m > 0) {
      idsInfMale <- sample(idsMale, i.num.m)
      hivStatMale[idsInfMale] <- 1
      infTimeMale[idsInfMale] <- sample(0:(-timeAidsDeath + 2),
                                        i.num.m, replace = TRUE)
    }


    # Initial contact network -------------------------------------------------

    # If we are *not* enforcing female monogamy, sample from the female IDs with
    #   replacement, and assign them to a female edgelist (list of persons in
    #   active partnerships)
    if (monog.f == FALSE) {
      edgelistFeml <- sample(idsFeml, expectedEdges, replace = TRUE)
      #  otherwise, sample from the female IDs without replacement
    } else {
      edgelistFeml <- sample(idsFeml, expectedEdges, replace = FALSE)
    }

    # Same for males
    if (monog.m == FALSE) {
      edgelistMale <- sample(idsMale, expectedEdges, replace = TRUE)
    } else {
      edgelistMale <- sample(idsMale, expectedEdges, replace = FALSE)
    }

    # These lines provide a check to see a table of relationships
    #   per person in the population. When monogamy is enforced for
    #   either sex, the corresponding table should be limited to 0's and 1's;
    #   when it os not enforced, the table is not limited.
    # table(tabulate(edgelistFeml))
    # table(tabulate(edgelistMale))


    # Time loop ---------------------------------------------------------------

    for (tstep in 1:nsteps) {

      ## Transmissions ##
      # Vector with HIV status for the female partner in each relationship
      hivStatFemlPart <- hivStatFeml[edgelistFeml]

      # Same for males
      hivStatMalePart <- hivStatMale[edgelistMale]

      # Vector of IDs for *relationships* with serodiscordant positive female
      sdpFeml <- which(hivStatMalePart == 0 & hivStatFemlPart == 1)

      # IDs for *females* in SDPF relationships
      idsSdpFeml <- edgelistFeml[sdpFeml]

      # Vector that identifies how long these females were infected
      infTimeSdpFeml <- tstep - infTimeFeml[idsSdpFeml]

      # Probability of transmission, based on time since the female was infected
      probTransSdpFeml <- betaByTime[infTimeSdpFeml]

      # Random transmissions given transmission probability
      transSdpFeml <- rbinom(length(sdpFeml), 1, probTransSdpFeml)

      # Double-indexing (an R specialty!):
      # Vector of males in SDP female relationships with a transmission event
      newInfMale <- edgelistMale[sdpFeml[transSdpFeml == 1]]

      # Assign the newly infected males an hivStat of 1
      hivStatMale[newInfMale] <- 1

      # Assign the newly infected males an infection time of tstep
      infTimeMale[newInfMale] <- tstep

      # All parallel for serodiscordant with positive male (SDPM)
      sdpMale <- which(hivStatFemlPart == 0 & hivStatMalePart == 1)
      idsSdpMale <- edgelistMale[sdpMale]
      infTimeSdpMale <- tstep - infTimeMale[idsSdpMale]
      probTransSdpMale <- betaByTime[infTimeSdpMale]
      transSdpMale <- rbinom(length(sdpMale), 1, probTransSdpMale)
      newInfFeml <- edgelistFeml[sdpMale[transSdpMale == 1]]
      hivStatFeml[newInfFeml] <- 1
      infTimeFeml[newInfFeml] <- tstep


      ## Deaths to AIDS ##
      # Which females have been HIV+ long enough to die of AIDS?
      idsDeathAidsFeml <- which((tstep - infTimeFeml) == timeAidsDeath)
      # Which females have been HIV+ long enough to die of AIDS?
      idsDeathAidsMale <- which((tstep - infTimeMale) == timeAidsDeath)

      ## End of ties because of death ##
      # Determine the IDs of those relationships involving dying women
      edgesFemlDeath <- which(edgelistFeml %in% idsDeathAidsFeml)
      # Determine the IDs of those relationships involving dying women
      edgesMaleDeath <- which(edgelistMale %in% idsDeathAidsMale)
      # Combine and remove from edgelists
      edgesBothDeath <- c(edgesFemlDeath, edgesMaleDeath)
      if (length(edgesBothDeath > 0)){
        edgelistFeml <- edgelistFeml[-edgesBothDeath]
        edgelistMale <- edgelistMale[-edgesBothDeath]
      }

      ## End of other ties randomly ##
      # Flip a weighted coin for each remaning relationship
      edgesElig <- rbinom(length(edgelistMale), 1, probDissolution)
      # Make vector of those edges to break
      edgesBroken <- which(edgesElig == 1)
      # If there are any edges to break, then break them.
      if (length(edgesBroken > 0)) {
        edgelistMale <- edgelistMale[-edgesBroken]
        edgelistFeml <- edgelistFeml[-edgesBroken]
      }


      ## Add New Edges ##
      # Determine how many edges to add.
      # We assume same as # broken in this simple model.
      nNewTies <- expectedEdges - length(edgelistMale)

      # If there are edges to add,
      if (nNewTies > 0) {
        #  and if we are *not* enforcing female monogamy,
        if (monog.f == FALSE) {
          # sample from the female IDs with replacement,
          f <- sample(1:nFeml, nNewTies, replace = TRUE)
          #	else if we are enforcing female monogamy,
        } else {
          # determine which women do not currently have a relationship, and
          femlNoTies <- setdiff(1:nFeml, edgelistFeml)
          # sample from them without replacement
          f <- sample(femlNoTies, nNewTies, replace = FALSE)
        }
        # and if we are *not* enforcing male monogamy,
        if (monog.m == FALSE) {
          # sample from the male IDs with replacement,
          m <- sample(1:nMale, nNewTies, replace = TRUE)
          # else if we are enforcing male monogamy,
        } else {
          # determine which men do not currently have a relationship, and
          maleNoTies <- setdiff(1:nMale, edgelistMale)
          # sample from them without replacement
          m <- sample(maleNoTies, nNewTies, replace = FALSE)
        }
        # Concatenate that into the existing edgelists
        edgelistFeml <- c(edgelistFeml, f)
        edgelistMale <- c(edgelistMale, m)
      }


      ## Insert new births ##
      # Vector of the IDs of females who just died, replaced them with a new birth
      idsBirthFeml <- idsDeathAidsFeml
      # Reset their attributes
      hivStatFeml[idsBirthFeml] <- 0
      infTimeFeml[idsBirthFeml] <- NA

      # All parallel with males
      idsBirthMale <- idsDeathAidsMale
      hivStatMale[idsBirthMale] <- 0
      infTimeMale[idsBirthMale] <- NA


      ## Track prevalence ##
      # Calculate current prevalence by sex
      if (tstep == 1) {
        femlPrev <- mean(hivStatFeml)
        malePrev <- mean(hivStatMale)
      } else {
        femlPrev[tstep] <- mean(hivStatFeml)
        malePrev[tstep] <- mean(hivStatMale)
      }

    }
    # Progress tracker
    stepTime <- (proc.time()[3] - ptm[3])
    doneTime <- (proc.time()[3] - ptm[3]) * (nsims-sim)
    if (verbose == TRUE) {
      if (doneTime/60 < 60) {
        cat("SIM = ", sim, "/", nsims,
            " | SIM Time = ", round(stepTime/60, 1), " Min",
            " | SIM Done = ", round(doneTime/60), " Min", "\n", sep="")
      } else {
        cat("SIM = ", sim, "/", nsims,
            " | SIM Time = ", round(stepTime/60, 1), " Min",
            " | SIM Done = ", round(doneTime/60/60, 1), " Hrs", "\n", sep="")
      }
    }

    prevBoth <- matrix(c(femlPrev, malePrev), ncol = 2)
    out[[sim]] <- prevBoth

  } # end sim loop

  class(out) <- "conc_microsim"
  return(out)
}



#' @title Plot Values from a Concurrency Microsimulation Model
#'
#' @description Plots values from an concurrency microsimulation
#'              epidemic model simulated with \code{conc_microsim}.
#'
#' @param x an object of class \code{conc_microsim}.
#' @param xlim x-axis scale limits for plot, with default based on model time steps.
#' @param ylim y-axis scale limits for plot, with default based on range of
#'        simulated data.
#' @param alpha transparency level for simulation lines, where 0 = transparent
#'        and 1 = opaque (see \code{\link{transco}}).
#' @param lwd line width for simulation lines.
#' @param ... additional arguments to pass to main plot window (see
#'        \code{\link{plot.default}}).
#'
#' @details
#' This function plots the disease prevalence from an \code{conc_microsim}
#' model. The function is currently limited to individual simulation lines of
#' disease prevalence. Future releases will standardize the plotting options to
#' those available with stochastic epidemic models, similar to those in
#' \code{\link{plot.icm}}.
#'
#' @method plot conc_microsim
#' @keywords internal
#' @export
#'
plot.conc_microsim <- function(x,
                               xlim,
                               ylim,
                               alpha = 1,
                               lwd = 1,
                               ...) {

  nsims <- length(x)
  nsteps <- nrow(x[[1]])
  pal <- transco(brewer.pal(3, "Set1"), alpha)

  if (missing(ylim)) ylim <- 0:1
  if (missing(xlim)) xlim <- c(0, nsteps)

  plot(1, 1, type="n", xlim = xlim, ylim = ylim, bty="n",
       xlab="Time (months)", ylab="Infected", ...)
  for (i in 1:nsims) {
    lines(1:nsteps, x[[i]][,1], col = pal[1], lwd = lwd)
    lines(1:nsteps, x[[i]][,2], col = pal[2], lwd = lwd)
  }
  legend("topleft", c("Female", "Male"), lwd=3, col=pal[1:2], bty="n", cex=0.9)

}
