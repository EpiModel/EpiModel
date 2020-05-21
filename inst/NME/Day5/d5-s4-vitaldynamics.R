
##
## Simple SI Model with Vital Dynamics
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Author: Samuel M. Jenness (Emory University)
## Date: August 2018
##

# Load EpiModel
library(EpiModel)
rm(list = ls())


# Vital Dynamics Setup ----------------------------------------------------

ages <- 0:85

# Rates per 100,000 for age groups: <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29,
#                                   30-34, 35-39, 40-44, 45-49, 50-54, 55-59,
#                                   60-64, 65-69, 70-74, 75-79, 80-84, 85+
# source: https://www.statista.com/statistics/241572/death-rate-by-age-and-sex-in-the-us/
mortality_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)
# rate per person, per week
mr_pp_pw <- mortality_rate / 1e5 / 52

# Build out a mortality rate vector
age_spans <- c(1, 4, rep(5, 16), 1)
mr_vec <- rep(mr_pp_pw, times = age_spans)
data.frame(ages, mr_vec)

par(mar = c(3,3,2,1), mgp = c(2,1,0))
barplot(mr_vec, col = "steelblue1", xlab = "age", ylab = "Death Rate")


# Network Model Estimation ------------------------------------------------

# Initialize the network
n <- 500
nw <- network.initialize(n, directed = FALSE)

# Set up ages
ageVec <- sample(ages, n, replace = TRUE)
nw <- set.vertex.attribute(nw, "age", ageVec)

# Define the formation model: edges
formation <- ~edges + absdiff("age")

# Input the appropriate target statistics for each term
mean_degree <- 0.8
edges <- mean_degree * (n/2)
avg.abs.age.diff <- 1.5
absdiff <- edges * avg.abs.age.diff

target.stats <- c(edges, absdiff)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(~offset(edges), 60, mean(mr_vec))
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 8, ncores = 8, nsteps = 500,
            nwstats.formula = ~edges + absdiff("age") + isolates + degree(0:5))
print(dx)
plot(dx, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Epidemic model parameters
param <- param.net(inf.prob = 0.15,
                   mortality.rates = mr_vec,
                   mortality.disease.mult = 2,
                   birth.rate = mean(mr_vec))

# Initial conditions
init <- init.net(i.num = 50)

# Read in the module functions
source("module-fx.R", echo = TRUE)

# Control settings
control <- control.net(type = "SI",
                       nsims = 1,
                       ncores = 1,
                       nsteps = 250,
                       aging.FUN = aging,
                       deaths.FUN = dfunc,
                       births.FUN = bfunc,
                       delete.nodes = TRUE,
                       depend = TRUE)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)
print(sim)

# Plot outcomes
par(mfrow = c(1,2))
plot(sim, main = "State Prevalences")
plot(sim, main = "State Sizes", sim.lines = TRUE,
     qnts = FALSE, mean.smooth = FALSE)

par(mfrow = c(1, 2))
plot(sim, y = "num", main = "Population Size", qnts = 1, ylim = c(450, 550))
plot(sim, y = "meanAge", main = "Mean Age", qnts = 1, ylim = c(35, 50))

par(mfrow = c(1, 2))
plot(sim, y = "d.flow", mean.smooth = TRUE, qnts = 1, main = "Deaths")
plot(sim, y = "b.flow", mean.smooth = TRUE, qnts = 1, main = "Births")

# Examine the data
df <- as.data.frame(sim)
head(df, 25)
