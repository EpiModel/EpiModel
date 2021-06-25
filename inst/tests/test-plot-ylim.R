
# Initialize network
nw <- network_initialize(n = 500)

# Define the formation model
formation <- ~edges + concurrent + degrange(from = 4)

# Input the appropriate target statistics for each term
target.stats <- c(175, 110, 0)

# Parameterize the dissolulation model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Review the arguments for the netest function
args(netest)

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 10, nsteps = 1000,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)
dx

# Use multiple cores
parallel::detectCores()

dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 5,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)
dx

# Plot the formation diagnostics
plot(dx)

nwstats1 <- get_nwstats(dx, sim = 1)
head(nwstats1, 20)

# Plot the dissolution diagnostics
par(mfrow = c(1, 2), mgp = c(2,1,0), mar = c(3,3,1,1))
plot(dx, type = "duration", duration.imputed = TRUE)
plot(dx, type = "duration", duration.imputed = FALSE)

plot(dx, type = "duration", sim.lines = TRUE, duration.imputed = TRUE, qnts = FALSE)
plot(dx, type = "duration", sim.lines = TRUE, duration.imputed = FALSE, qnts = FALSE)

dx <- netdx(est, nsims = 1, nsteps = 1000,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)
dx

plot(dx, type = "duration", duration.imputed = TRUE)
plot(dx, type = "duration", duration.imputed = FALSE)

plot(dx, type = "duration", sim.lines = TRUE, duration.imputed = TRUE, qnts = FALSE)
plot(dx, type = "duration", sim.lines = TRUE, duration.imputed = FALSE, qnts = FALSE)
