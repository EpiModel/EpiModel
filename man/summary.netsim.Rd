% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.netsim}
\alias{summary.netsim}
\title{Summary Model Statistics}
\usage{
\method{summary}{netsim}(object, at, digits = 3, ...)
}
\arguments{
\item{object}{An \code{EpiModel} object of class \code{netsim}.}

\item{at}{Time step for model statistics.}

\item{digits}{Number of significant digits to print.}

\item{...}{Additional summary function arguments.}
}
\description{
Extracts and prints model statistics simulated with
\code{netsim}.
}
\details{
This function provides summary statistics for the main epidemiological
outcomes (state and transition size and prevalence) from a \code{netsim}
model. Time-specific summary measures are provided, so it is necessary to
input a time of interest.
}
\examples{
\dontrun{
## SI Model without Network Feedback
# Initialize network and set network model parameters
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Estimate the ERGM models (see help for netest)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Parameters, initial conditions, and controls for model
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)

# Run the model simulation
mod <- netsim(est1, param, init, control)

summary(mod, at = 1)
summary(mod, at = 50)
summary(mod, at = 100)
}

}
\seealso{
\code{\link{netsim}}
}
\keyword{extract}
