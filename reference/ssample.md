# Stable Sampling Function

Provides a sampling function useful for dynamic simulations, in which
the length of the input vector may be multiple lengths and the size of
the sample may be 0.

## Usage

``` r
ssample(x, size, replace = FALSE, prob = NULL)
```

## Arguments

- x:

  Either a vector of one or more elements from which to choose, or a
  positive integer.

- size:

  Non-negative integer giving the number of items to choose.

- replace:

  Should sampling be with replacement?

- prob:

  Vector of probability weights for obtaining the elements of the vector
  being sampled.

## Value

A vector containing the sampled value(s).
