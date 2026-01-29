# Create a Value Sampler for Random Parameters

This function returns a 0 argument function that can be used as a
generator function in the `random.params` argument of the
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
function.

## Usage

``` r
param_random(values, prob = NULL)
```

## Arguments

- values:

  A vector of values to sample from.

- prob:

  A vector of weights to use during sampling. If `NULL`, all values have
  the same probability of being picked (default = `NULL`).

## Value

A 0 argument generator function to sample one of the values from the
`values` vector.

## See also

[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
and
[`generate_random_params()`](http://epimodel.github.io/EpiModel/reference/generate_random_params.md)

## Examples

``` r
# Define function with equal sampling probability
a <- param_random(1:5)
a()
#> [1] 5

# Define function with unequal sampling probability
b <- param_random(1:5, prob = c(0.1, 0.1, 0.1, 0.1, 0.6))
b()
#> [1] 1
```
