# Apportion Using the Largest Remainder Method

Apportions a vector of values given a specified frequency distribution
of those values such that the length of the output is robust to rounding
and other instabilities.

## Usage

``` r
apportion_lr(vector.length, values, proportions, shuffled = FALSE)
```

## Arguments

- vector.length:

  Length for the output vector.

- values:

  Values for the output vector.

- proportions:

  Proportion distribution with one number for each value. This must sum
  to 1.

- shuffled:

  If `TRUE`, randomly shuffle the order of the vector.

## Value

A vector of length `vector.length` containing the apportioned values
from `values`.

## Examples

``` r
if (FALSE) { # \dontrun{
## Example 1: Without rounding
apportioned_vec_1 <- apportion_lr(4, c(1, 2, 3, 4, 5),
                                     c(0.25, 0, 0.25, 0.25, 0.25))

## Example 2: With rounding
apportioned_vec_2 <- apportion_lr(5, c(1, 2, 3, 4, 5),
                                     c(0.21, 0, 0.29, 0.25, 0.25))
} # }
```
