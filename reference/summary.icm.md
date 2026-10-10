# Summary Model Statistics

Extracts and prints model statistics simulated with `icm`.

## Usage

``` r
# S3 method for class 'icm'
summary(object, at, digits = 3, ...)
```

## Arguments

- object:

  An `EpiModel` object of class `icm`.

- at:

  Time step for model statistics.

- digits:

  Number of significant digits to print.

- ...:

  Additional summary function arguments.

## Details

This function provides summary statistics for the main epidemiological
outcomes (state and transition size and prevalence) from an `icm` model.
Time-specific summary measures are provided, so it is necessary to input
a time of interest.

## See also

[`icm()`](https://epimodel.github.io/EpiModel/reference/icm.md)

## Examples

``` r
# \donttest{
## Stochastic ICM SI model with 3 simulations
param <- param.icm(inf.prob = 0.2, act.rate = 1)
init <- init.icm(s.num = 500, i.num = 1)
control <- control.icm(type = "SI", nsteps = 50,
                       nsims = 5, verbose = FALSE)
mod <- icm(param, init, control)
summary(mod, at = 25)
#> EpiModel Summary
#> =======================
#> Model class: icm
#> 
#> Simulation Details
#> -----------------------
#> Model type: SI
#> No. simulations: 5
#> No. time steps: 50
#> No. groups: 1
#> 
#> Model Statistics
#> ------------------------------
#> Time: 25 
#> ------------------------------ 
#>            mean      sd   pct
#> Suscept.  441.0  38.523  0.88
#> Infect.    60.0  38.523  0.12
#> Total     501.0   0.000  1.00
#> S -> I      7.6   5.983    NA
#> ------------------------------ 
summary(mod, at = 50)
#> EpiModel Summary
#> =======================
#> Model class: icm
#> 
#> Simulation Details
#> -----------------------
#> Model type: SI
#> No. simulations: 5
#> No. time steps: 50
#> No. groups: 1
#> 
#> Model Statistics
#> ------------------------------
#> Time: 50 
#> ------------------------------ 
#>            mean     sd    pct
#> Suscept.   37.2  20.56  0.074
#> Infect.   463.8  20.56  0.926
#> Total     501.0   0.00  1.000
#> S -> I      8.6   4.93     NA
#> ------------------------------ 
# }
```
