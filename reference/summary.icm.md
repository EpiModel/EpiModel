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

[`icm()`](http://epimodel.github.io/EpiModel/reference/icm.md)

## Examples

``` r
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
#> Suscept.  480.8  17.283  0.96
#> Infect.    20.2  17.283  0.04
#> Total     501.0   0.000  1.00
#> S -> I      3.0   4.243    NA
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
#>            mean       sd    pct
#> Suscept.  151.6  129.855  0.303
#> Infect.   349.4  129.855  0.697
#> Total     501.0    0.000  1.000
#> S -> I     16.2    4.266     NA
#> ------------------------------ 
```
