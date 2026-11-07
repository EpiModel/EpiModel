# Check Mean Edge Durations Against the Departure Ceiling

Stops when a mean edge duration is at or above the departure ceiling and
warns when one is at or above 90% of it, for
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).

## Usage

``` r
check_departure_ceiling(duration, d.rate)
```

## Arguments

- duration:

  Vector of mean edge durations.

- d.rate:

  Departure rate, a single number.

## Value

`NULL`, invisibly.

## Details

The ceiling is `1 / (1 - (1 - d.rate)^2)`, the mean duration at which
every partnership ending is a departure; see the Details of
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).
The error and the warning carry the caller's call, so they name
[`dissolution_coefs()`](https://epimodel.github.io/EpiModel/reference/dissolution_coefs.md).
