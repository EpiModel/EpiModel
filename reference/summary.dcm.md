# Summary Model Statistics

Extracts and prints model statistics solved with `dcm`.

## Usage

``` r
# S3 method for class 'dcm'
summary(object, at, run = 1, digits = 3, ...)
```

## Arguments

- object:

  An `EpiModel` object of class `dcm`.

- at:

  Time step for model statistics.

- run:

  Model run number, for `dcm` class models with multiple runs
  (sensitivity analyses).

- digits:

  Number of significant digits to print.

- ...:

  Additional summary function arguments (not used).

## Details

This function provides summary statistics for the main epidemiological
outcomes (state and transition size and prevalence) from a `dcm` model.
Time-specific summary measures are provided, so it is necessary to input
a time of interest. For multiple-run models (sensitivity analyses),
input a model run number. See examples below.

## See also

[`dcm()`](http://epimodel.github.io/EpiModel/reference/dcm.md)

## Examples

``` r
## Deterministic SIR model with varying act.rate
param <- param.dcm(inf.prob = 0.2, act.rate = 2:4, rec.rate = 1/3,
                   a.rate = 0.011, ds.rate = 0.01,
                   di.rate = 0.03, dr.rate = 0.01)
init <- init.dcm(s.num = 1000, i.num = 1, r.num = 0)
control <- control.dcm(type = "SIR", nsteps = 50)
mod <- dcm(param, init, control)
summary(mod, at = 25, run = 1)
#> EpiModel Summary
#> =======================
#> Model class: dcm
#> 
#> Simulation Summary
#> -----------------------
#> Model type: SIR
#> No. runs: 3
#> No. time steps:
#> No. groups: 1
#> 
#> Model Statistics
#> ------------------------------
#> Time: 25  Run: 1 
#> ------------------------------ 
#>                        n    pct
#> Suscept.        1010.999  0.987
#> Infect.            2.265  0.002
#> Recov.            11.293  0.011
#> Total           1024.557  1.000
#> S -> I             0.908     NA
#> I -> R             0.767     NA
#> Arrival ->        11.276     NA
#> S Departure ->    10.111     NA
#> I Departure ->     0.069     NA
#> R Departure ->     0.116     NA
#> ------------------------------ 
summary(mod, at = 25, run = 3)
#> EpiModel Summary
#> =======================
#> Model class: dcm
#> 
#> Simulation Summary
#> -----------------------
#> Model type: SIR
#> No. runs: 3
#> No. time steps:
#> No. groups: 1
#> 
#> Model Statistics
#> ------------------------------
#> Time: 25  Run: 3 
#> ------------------------------ 
#>                       n    pct
#> Suscept.        224.736  0.229
#> Infect.          79.014  0.081
#> Recov.          677.357  0.690
#> Total           981.106  1.000
#> S -> I           13.102     NA
#> I -> R           24.087     NA
#> Arrival ->       10.789     NA
#> S Departure ->    2.222     NA
#> I Departure ->    2.168     NA
#> R Departure ->    6.863     NA
#> ------------------------------ 
summary(mod, at = 26, run = 3)
#> EpiModel Summary
#> =======================
#> Model class: dcm
#> 
#> Simulation Summary
#> -----------------------
#> Model type: SIR
#> No. runs: 3
#> No. time steps:
#> No. groups: 1
#> 
#> Model Statistics
#> ------------------------------
#> Time: 26  Run: 3 
#> ------------------------------ 
#>                       n    pct
#> Suscept.        220.201  0.225
#> Infect.          65.862  0.067
#> Recov.          694.580  0.708
#> Total           980.642  1.000
#> S -> I           10.746     NA
#> I -> R           20.049     NA
#> Arrival ->       10.786     NA
#> S Departure ->    2.189     NA
#> I Departure ->    1.804     NA
#> R Departure ->    7.014     NA
#> ------------------------------ 
```
