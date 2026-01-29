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
#>            mean      sd    pct
#> Suscept.  448.8  60.648  0.896
#> Infect.    52.2  60.648  0.104
#> Total     501.0   0.000  1.000
#> S -> I      9.6   8.792     NA
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
#>            mean      sd    pct
#> Suscept.   84.0  66.667  0.168
#> Infect.   417.0  66.667  0.832
#> Total     501.0   0.000  1.000
#> S -> I     13.6  10.877     NA
#> ------------------------------ 
```
