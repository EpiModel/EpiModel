# Generate Values for Random Parameters

This function uses the generative functions in the `random.params` list
to create values for the parameters.

## Usage

``` r
generate_random_params(param, verbose = FALSE)
```

## Arguments

- param:

  The `param` argument received by the `netsim` functions.

- verbose:

  Should the function output the generated values (default = FALSE)?

## Value

A fully instantiated `param` list.

## `random.params`

The `random.params` argument to the
[`param.net()`](http://epimodel.github.io/EpiModel/reference/param.net.md)
function must be a named list of functions that each return a value that
can be used as the argument with the same name. In the example below,
`param_random` is a function factory provided by EpiModel for `act.rate`
and for `tx.halt.part.prob` we provide bespoke functions. A function
factory is a function that returns a new function (see
https://adv-r.hadley.nz/function-factories.html).

## Generator Functions

The functions used inside `random_params` must be 0 argument functions
returning a valid value for the parameter with the same name.

## `param_random_set`

The `random_params` list can optionally contain a `param_random_set`
element. It must be a `data.frame` of possible values to be used as
parameters.

The column names must correspond either to: the name of one parameter,
if this parameter is of size 1; or the name of one parameter with "\_1",
"*2", etc. appended, with the number representing the position of the
value, if this parameter is of size \> 1. This means that the parameter
names cannot contain any underscores "*" if you intend to use
`param_random_set`.

The point of the `param.random.set` `data.frame` is to allow the random
parameters to be correlated. To achieve this, a whole row of the
`data.frame` is selected for each simulation.

## Examples

``` r
if (FALSE) { # \dontrun{

## Example with only the generator function

# Define random parameter list
my_randoms <- list(
  act.rate = param_random(c(0.25, 0.5, 0.75)),
  tx.prob = function() rbeta(1, 1, 2),
  stratified.test.rate = function() c(
    rnorm(1, 0.05, 0.01),
    rnorm(1, 0.15, 0.03),
    rnorm(1, 0.25, 0.05)
  )
)

# Parameter model with fixed and random parameters
param <- param.net(inf.prob = 0.3, random.params = my_randoms)

# Below, `tx.prob` is set first to 0.3 then assigned a random value using
# the function from `my_randoms`. A warning notifying of this overwrite is
# therefore produced.
param <- param.net(tx.prob = 0.3, random.params = my_randoms)


# Parameters are drawn automatically in netsim by calling the function
# within netsim_loop. Demonstrating draws here but this is not used by
# end user.
paramDraw <- generate_random_params(param, verbose = TRUE)
paramDraw


## Addition of the `param.random.set` `data.frame`

# This function will generate sets of correlated parameters
 generate_correlated_params <- function() {
   param.unique <- runif(1)
   param.set.1 <- param.unique + runif(2)
   param.set.2 <- param.unique * rnorm(3)

   return(list(param.unique, param.set.1, param.set.2))
 }

 # Data.frame set of random parameters :
 correlated_params <- t(replicate(10, unlist(generate_correlated_params())))
 correlated_params <- as.data.frame(correlated_params)
 colnames(correlated_params) <- c(
   "param.unique",
   "param.set.1_1", "param.set.1_2",
   "param.set.2_1", "param.set.2_2", "param.set.2_3"
 )

# Define random parameter list with the `param.random.set` element
my_randoms <- list(
  act.rate = param_random(c(0.25, 0.5, 0.75)),
  param.random.set = correlated_params
)

# Parameter model with fixed and random parameters
param <- param.net(inf.prob = 0.3, random.params = my_randoms)

# Parameters are drawn automatically in netsim by calling the function
# within netsim_loop. Demonstrating draws here but this is not used by
# end user.
paramDraw <- generate_random_params(param, verbose = TRUE)
paramDraw

} # }
```
