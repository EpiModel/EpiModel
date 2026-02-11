# Truncate Simulation Time Series

Left-truncates simulation epidemiological summary statistics and network
statistics at a specified time step.

## Usage

``` r
truncate_sim(x, at, reset.time)

# S3 method for class 'dcm'
truncate_sim(x, at, reset.time = TRUE)

# S3 method for class 'icm'
truncate_sim(x, at, reset.time = TRUE)

# S3 method for class 'netsim'
truncate_sim(x, at, reset.time = TRUE)
```

## Arguments

- x:

  Object of class `dcm`, `netsim`, or `icm`.

- at:

  Time step at which to left-truncate the time series.

- reset.time:

  If `TRUE`, the time step sequence in the truncated model will be reset
  to start at 1. If `FALSE`, the original time step values will be
  preserved. Default is `TRUE`.

## Value

The updated object of class `dcm`, `netsim`, or `icm`.

## Details

This function would be used when running a follow-up simulation from
time steps `b` to `c` after a burn-in period from time `a` to `b`, where
the final time window of interest for data analysis is `b` to `c` only.

## Examples

``` r
# DCM examples
param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
init <- init.dcm(s.num = 500, i.num = 1)
control <- control.dcm(type = "SI", nsteps = 200)
mod1 <- dcm(param, init, control)
plot(mod1)


# Reset time
mod2a <- truncate_sim(mod1, at = 150)
plot(mod2a)

head(as.data.frame(mod2a))
#>   time     s.num    i.num  si.flow num
#> 1    1 112.84481 388.1552 4.311261 501
#> 2    2 108.53355 392.4664 4.190826 501
#> 3    3 104.34273 396.6573 4.070349 501
#> 4    4 100.27238 400.7276 3.950125 501
#> 5    5  96.32225 404.6777 3.830430 501
#> 6    6  92.49182 408.5082 3.711524 501

# Do not reset time
mod2b <- truncate_sim(mod1, at = 150, reset.time = FALSE)
plot(mod2b)

head(as.data.frame(mod2b))
#>   time     s.num    i.num  si.flow num
#> 1  150 112.84481 388.1552 4.311261 501
#> 2  151 108.53355 392.4664 4.190826 501
#> 3  152 104.34273 396.6573 4.070349 501
#> 4  153 100.27238 400.7276 3.950125 501
#> 5  154  96.32225 404.6777 3.830430 501
#> 6  155  92.49182 408.5082 3.711524 501

# ICM example
param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SI", nsteps = 200, nsims = 1)
mod1 <- icm(param, init, control)

# Reset time
mod2a <- truncate_sim(mod1, at = 150)
plot(mod2a)

head(as.data.frame(mod2a))
#>   sim time s.num i.num num si.flow
#> 1   1    1    49   452 501       3
#> 2   1    2    48   453 501       1
#> 3   1    3    42   459 501       6
#> 4   1    4    41   460 501       1
#> 5   1    5    41   460 501       0
#> 6   1    6    41   460 501       0

# Do not reset time
mod2b <- truncate_sim(mod1, at = 150, reset.time = FALSE)
plot(mod2b)

head(as.data.frame(mod2b))
#>   sim time s.num i.num num si.flow
#> 1   1  150    49   452 501       3
#> 2   1  151    48   453 501       1
#> 3   1  152    42   459 501       6
#> 4   1  153    41   460 501       1
#> 5   1  154    41   460 501       0
#> 6   1  155    41   460 501       0
```
